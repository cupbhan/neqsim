package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryGlobalStabilityGate.RetainedPhaseLocalStability;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Classifies a stopped hydrocarbon-water two-to-three-phase continuation without promoting a numerical failure to a
 * physical endpoint.
 *
 * <p>
 * A retained-phase spinodal is accepted only when an independently stable boundary root and the first rejected boundary
 * root bracket a zero homogeneous-phase stationarity eigenvalue and {@link HydrocarbonWaterRetainedPhaseSpinodalSolver}
 * converges inside that bracket. A lower TPD branch, a minimum continuation step, or an unrefined negative eigenvalue
 * remains unresolved.
 * </p>
 */
public final class HydrocarbonWaterBoundaryTerminationClassifier {
  private final SystemInterface template;
  private int maximumIterations = 40;
  private int boundaryMaximumIterations = 80;
  private double boundaryResidualTolerance = 1.0e-8;
  private double boundaryFiniteDifferenceStep = 2.0e-5;
  private double curvatureFiniteDifferenceStep = 2.0e-4;
  private double pressureBracketToleranceBara = 1.0e-4;
  private double eigenvalueTolerance = 1.0e-6;

  /**
   * Creates a non-destructive termination classifier.
   *
   * @param template configured thermodynamic system
   */
  public HydrocarbonWaterBoundaryTerminationClassifier(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /**
   * Sets the retained-phase spinodal refinement controls.
   *
   * @param maximumIterations maximum spinodal bisections
   * @param boundaryMaximumIterations maximum fixed-pressure boundary iterations
   * @param boundaryResidualTolerance boundary residual tolerance
   * @param boundaryFiniteDifferenceStep boundary finite-difference step
   * @param curvatureFiniteDifferenceStep stationarity-Jacobian finite-difference step
   * @param pressureBracketToleranceBara terminal pressure-bracket width in bara
   * @param eigenvalueTolerance terminal absolute eigenvalue tolerance
   * @return this classifier
   */
  public HydrocarbonWaterBoundaryTerminationClassifier setNumericalControls(int maximumIterations,
      int boundaryMaximumIterations, double boundaryResidualTolerance, double boundaryFiniteDifferenceStep,
      double curvatureFiniteDifferenceStep, double pressureBracketToleranceBara, double eigenvalueTolerance) {
    if (maximumIterations < 1 || boundaryMaximumIterations < 1 || !positive(boundaryResidualTolerance)
        || !positive(boundaryFiniteDifferenceStep) || !positive(curvatureFiniteDifferenceStep)
        || !positive(pressureBracketToleranceBara) || !positive(eigenvalueTolerance)) {
      throw new IllegalArgumentException("invalid hydrocarbon-water termination controls");
    }
    this.maximumIterations = maximumIterations;
    this.boundaryMaximumIterations = boundaryMaximumIterations;
    this.boundaryResidualTolerance = boundaryResidualTolerance;
    this.boundaryFiniteDifferenceStep = boundaryFiniteDifferenceStep;
    this.curvatureFiniteDifferenceStep = curvatureFiniteDifferenceStep;
    this.pressureBracketToleranceBara = pressureBracketToleranceBara;
    this.eigenvalueTolerance = eigenvalueTolerance;
    return this;
  }

  /**
   * Classifies one stopped trace using its last accepted seed as fallback bracket evidence.
   *
   * @param trace pseudo-arclength trace with global-stability trials
   * @param acceptedSeed independently accepted boundary root adjacent to the traced end
   * @return immutable termination classification
   */
  public Result classify(TwoToThreePhasePseudoArcLengthTracer.Result trace,
      TwoToThreePhaseBoundaryPointSolver.Result acceptedSeed) {
    if (trace == null || acceptedSeed == null || !acceptedSeed.isConverged()) {
      throw new IllegalArgumentException("a trace and converged accepted seed are required");
    }
    if (trace.hasCompletedRequestedPoints()) {
      TwoToThreePhaseArcLengthCorrector.State last = lastState(trace);
      return Result.noTermination(last.getTemperatureK(), last.getPressureBara());
    }

    HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejected = null;
    TwoToThreePhaseBoundaryPointSolver.Result nearestAccepted = acceptedSeed;
    if (!trace.getGlobalStabilityEvidence().isEmpty()) {
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result terminalEvidence = trace.getGlobalStabilityEvidence()
          .get(trace.getGlobalStabilityEvidence().size() - 1);
      if (!terminalEvidence.isAccepted()) {
        rejected = terminalEvidence;
        RetainedPhaseLocalStability terminalUnstable = mostUnstableRetainedPhase(terminalEvidence);
        CandidatePhase destabilizingPhase = terminalUnstable == null ? null : terminalUnstable.getPhase();
        for (int index = trace.getGlobalStabilityEvidence().size() - 2; index >= 0; index--) {
          HydrocarbonWaterBoundaryGlobalStabilityGate.Result evidence = trace.getGlobalStabilityEvidence().get(index);
          if (evidence.isAccepted()
              && (destabilizingPhase == null || hasNonnegativeRetainedMode(evidence, destabilizingPhase))) {
            nearestAccepted = evidence.getBoundaryRoot();
            break;
          }
        }
      }
    }
    if (rejected == null) {
      TwoToThreePhaseArcLengthCorrector.State last = lastState(trace);
      return Result.numerical(last.getTemperatureK(), last.getPressureBara(),
          "continuation stopped without a rejected global-stability boundary root: " + trace.getTerminationReason()
              + (trace.getFailureMessage() == null ? "" : "; " + trace.getFailureMessage()));
    }
    return classify(nearestAccepted, rejected);
  }

  /**
   * Refines and classifies a stable-to-rejected boundary bracket.
   *
   * @param acceptedRoot last independently accepted boundary root
   * @param rejectedEvidence first rejected global-stability evaluation
   * @return immutable termination classification
   */
  public Result classify(TwoToThreePhaseBoundaryPointSolver.Result acceptedRoot,
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejectedEvidence) {
    if (acceptedRoot == null || !acceptedRoot.isConverged() || rejectedEvidence == null
        || rejectedEvidence.getBoundaryRoot() == null || rejectedEvidence.isAccepted()) {
      throw new IllegalArgumentException(
          "a converged accepted root and rejected global-stability evidence are required");
    }
    TwoToThreePhaseBoundaryPointSolver.Result rejectedRoot = rejectedEvidence.getBoundaryRoot();
    validateSameBoundaryFamily(acceptedRoot, rejectedRoot);
    RetainedPhaseLocalStability unstable = mostUnstableRetainedPhase(rejectedEvidence);
    if (unstable == null) {
      return Result.unresolved(rejectedRoot.getTemperatureK(), rejectedRoot.getPressureBara(), null, Double.NaN,
          "global stability failed without a negative retained-phase eigenvalue" + diagnosticSuffix(rejectedEvidence));
    }

    HydrocarbonWaterRetainedPhaseSpinodalSolver solver = new HydrocarbonWaterRetainedPhaseSpinodalSolver(template,
        acceptedRoot.getRetainedPhaseZero(), acceptedRoot.getRetainedPhaseOne(), acceptedRoot.getIncipientPhase(),
        unstable.getPhase()).setNumericalControls(maximumIterations, boundaryMaximumIterations,
            boundaryResidualTolerance, boundaryFiniteDifferenceStep, curvatureFiniteDifferenceStep,
            pressureBracketToleranceBara, eigenvalueTolerance);
    HydrocarbonWaterRetainedPhaseSpinodalSolver.Result refined = solver.solve(acceptedRoot, rejectedRoot);
    if (!refined.isConverged()) {
      return Result.unresolved(rejectedRoot.getTemperatureK(), rejectedRoot.getPressureBara(), unstable.getPhase(),
          Math.abs(unstable.getMinimumEigenvalue()),
          "retained-phase stability loss was detected but the zero-eigenvalue bracket did not " + "converge: "
              + refined.getFailureMessage() + diagnosticSuffix(rejectedEvidence));
    }

    HydrocarbonWaterRetainedPhaseSpinodalSolver.CurvatureState spinodal = refined.getSpinodal();
    String diagnostic = "retained " + unstable.getPhase()
        + " homogeneous stationarity eigenvalue was independently bracketed and refined to zero; "
        + "pressureBracketWidthBara=" + refined.getPressureBracketWidthBara() + "; iterations="
        + refined.getIterations();
    return Result.spinodal(spinodal.getState().getTemperatureK(), spinodal.getState().getPressureBara(),
        unstable.getPhase(), Math.abs(spinodal.getMinimumEigenvalue()), diagnostic, refined);
  }

  private static TwoToThreePhaseArcLengthCorrector.State lastState(TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    if (trace.getPoints().isEmpty()) {
      throw new IllegalArgumentException("termination classification requires at least one trace point");
    }
    return trace.getPoints().get(trace.getPoints().size() - 1);
  }

  private static RetainedPhaseLocalStability mostUnstableRetainedPhase(
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejected) {
    RetainedPhaseLocalStability selected = null;
    for (RetainedPhaseLocalStability candidate : rejected.getRetainedPhaseStability()) {
      if (!candidate.isAccepted() && Double.isFinite(candidate.getMinimumEigenvalue())
          && candidate.getMinimumEigenvalue() < 0.0
          && (selected == null || candidate.getMinimumEigenvalue() < selected.getMinimumEigenvalue())) {
        selected = candidate;
      }
    }
    return selected;
  }

  private static boolean hasNonnegativeRetainedMode(HydrocarbonWaterBoundaryGlobalStabilityGate.Result evidence,
      CandidatePhase phase) {
    for (RetainedPhaseLocalStability retained : evidence.getRetainedPhaseStability()) {
      if (retained.getPhase() == phase && Double.isFinite(retained.getMinimumEigenvalue())
          && retained.getMinimumEigenvalue() >= 0.0) {
        return true;
      }
    }
    return false;
  }

  private static void validateSameBoundaryFamily(TwoToThreePhaseBoundaryPointSolver.Result accepted,
      TwoToThreePhaseBoundaryPointSolver.Result rejected) {
    if (accepted.getRetainedPhaseZero() != rejected.getRetainedPhaseZero()
        || accepted.getRetainedPhaseOne() != rejected.getRetainedPhaseOne()
        || accepted.getIncipientPhase() != rejected.getIncipientPhase()) {
      throw new IllegalArgumentException("accepted and rejected roots belong to different boundaries");
    }
  }

  private static String diagnosticSuffix(HydrocarbonWaterBoundaryGlobalStabilityGate.Result evidence) {
    return evidence.getFailureMessage() == null ? "" : "; " + evidence.getFailureMessage();
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Mutually exclusive continuation termination classifications. */
  public enum Type {
    /** The requested continuation window ended before a branch endpoint was observed. */
    NO_TERMINATION_OBSERVED,
    /** A retained phase lost local composition stability at a refined zero eigenvalue. */
    RETAINED_PHASE_SPINODAL,
    /** Global stability failed, but the physical endpoint equations did not converge. */
    UNRESOLVED_GLOBAL_STABILITY_LIMIT,
    /** Continuation stopped for a numerical reason without a physical endpoint bracket. */
    NUMERICAL_TERMINATION
  }

  /** Immutable endpoint classification and refinement audit. */
  public static final class Result {
    private final Type type;
    private final double temperatureK;
    private final double pressureBara;
    private final CandidatePhase destabilizingPhase;
    private final double qualityMeasure;
    private final String diagnostic;
    private final HydrocarbonWaterRetainedPhaseSpinodalSolver.Result spinodalResult;

    private Result(Type type, double temperatureK, double pressureBara, CandidatePhase destabilizingPhase,
        double qualityMeasure, String diagnostic, HydrocarbonWaterRetainedPhaseSpinodalSolver.Result spinodalResult) {
      this.type = type;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.destabilizingPhase = destabilizingPhase;
      this.qualityMeasure = qualityMeasure;
      this.diagnostic = diagnostic;
      this.spinodalResult = spinodalResult;
    }

    private static Result noTermination(double temperatureK, double pressureBara) {
      return new Result(Type.NO_TERMINATION_OBSERVED, temperatureK, pressureBara, null, Double.NaN,
          "requested continuation points completed; no endpoint was observed", null);
    }

    private static Result spinodal(double temperatureK, double pressureBara, CandidatePhase destabilizingPhase,
        double qualityMeasure, String diagnostic, HydrocarbonWaterRetainedPhaseSpinodalSolver.Result spinodalResult) {
      return new Result(Type.RETAINED_PHASE_SPINODAL, temperatureK, pressureBara, destabilizingPhase, qualityMeasure,
          diagnostic, spinodalResult);
    }

    private static Result unresolved(double temperatureK, double pressureBara, CandidatePhase destabilizingPhase,
        double qualityMeasure, String diagnostic) {
      return new Result(Type.UNRESOLVED_GLOBAL_STABILITY_LIMIT, temperatureK, pressureBara, destabilizingPhase,
          qualityMeasure, diagnostic, null);
    }

    private static Result numerical(double temperatureK, double pressureBara, String diagnostic) {
      return new Result(Type.NUMERICAL_TERMINATION, temperatureK, pressureBara, null, Double.NaN, diagnostic, null);
    }

    /** @return termination classification */
    public Type getType() {
      return type;
    }

    /** @return true only for an independently refined physical boundary endpoint */
    public boolean isPhysicalEndpoint() {
      return type == Type.RETAINED_PHASE_SPINODAL && spinodalResult != null && spinodalResult.isConverged();
    }

    /** @return endpoint or last traced temperature in kelvin */
    public double getTemperatureK() {
      return temperatureK;
    }

    /** @return endpoint or last traced pressure in bara */
    public double getPressureBara() {
      return pressureBara;
    }

    /** @return destabilizing retained phase, or {@code null} when not identified */
    public CandidatePhase getDestabilizingPhase() {
      return destabilizingPhase;
    }

    /** @return absolute terminal eigenvalue or other type-specific quality measure */
    public double getQualityMeasure() {
      return qualityMeasure;
    }

    /** @return auditable classification diagnostic */
    public String getDiagnostic() {
      return diagnostic;
    }

    /** @return converged retained-phase spinodal evidence, or {@code null} */
    public HydrocarbonWaterRetainedPhaseSpinodalSolver.Result getSpinodalResult() {
      return spinodalResult;
    }

    /**
     * Reconstructs the rigorously corrected endpoint boundary root.
     *
     * @return corrected endpoint root, or {@code null} for an unresolved termination
     */
    public TwoToThreePhaseBoundaryPointSolver.Result getBoundaryRoot() {
      if (!isPhysicalEndpoint()) {
        return null;
      }
      HydrocarbonWaterRetainedPhaseSpinodalSolver.CurvatureState spinodal = spinodalResult.getSpinodal();
      return TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(spinodal.getState(),
          spinodal.getBoundaryResidual());
    }
  }
}
