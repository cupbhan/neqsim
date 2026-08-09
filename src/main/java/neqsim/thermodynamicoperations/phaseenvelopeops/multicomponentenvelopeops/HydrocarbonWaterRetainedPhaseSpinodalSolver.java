package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Locates a retained-phase spinodal where a secondary stationary TPD branch bifurcates from it. */
public final class HydrocarbonWaterRetainedPhaseSpinodalSolver {
  private final SystemInterface template;
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase incipientPhase;
  private final CandidatePhase bifurcatingPhase;
  private int maximumIterations = 40;
  private int boundaryMaximumIterations = 80;
  private double boundaryResidualTolerance = 1.0e-8;
  private double boundaryFiniteDifferenceStep = 2.0e-5;
  private double curvatureFiniteDifferenceStep = 2.0e-4;
  private double pressureBracketToleranceBara = 1.0e-4;
  private double eigenvalueTolerance = 1.0e-6;

  /** Creates a spinodal solver on one specified two-to-three-phase boundary family. */
  public HydrocarbonWaterRetainedPhaseSpinodalSolver(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase, CandidatePhase bifurcatingPhase) {
    if (template == null || retainedPhaseZero == null || retainedPhaseOne == null || incipientPhase == null
        || bifurcatingPhase == null) {
      throw new IllegalArgumentException("template and all phase identities are required");
    }
    if (bifurcatingPhase != retainedPhaseZero && bifurcatingPhase != retainedPhaseOne) {
      throw new IllegalArgumentException("the bifurcating phase must be one of the retained phases");
    }
    this.template = template.clone();
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.incipientPhase = incipientPhase;
    this.bifurcatingPhase = bifurcatingPhase;
  }

  /** Sets boundary-correction, curvature, pressure-bracket, and eigenvalue controls. */
  public HydrocarbonWaterRetainedPhaseSpinodalSolver setNumericalControls(int maximumIterations,
      int boundaryMaximumIterations, double boundaryResidualTolerance, double boundaryFiniteDifferenceStep,
      double curvatureFiniteDifferenceStep, double pressureBracketToleranceBara, double eigenvalueTolerance) {
    if (maximumIterations < 1 || boundaryMaximumIterations < 1 || !positive(boundaryResidualTolerance)
        || !positive(boundaryFiniteDifferenceStep) || !positive(curvatureFiniteDifferenceStep)
        || !positive(pressureBracketToleranceBara) || !positive(eigenvalueTolerance)) {
      throw new IllegalArgumentException("invalid retained-phase spinodal numerical controls");
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

  /** Refines two boundary roots whose retained-phase minimum TPD eigenvalues have opposite signs. */
  public Result solve(TwoToThreePhaseBoundaryPointSolver.Result first,
      TwoToThreePhaseBoundaryPointSolver.Result second) {
    validateRoot(first);
    validateRoot(second);
    TwoToThreePhaseArcLengthCorrector.State firstState = TwoToThreePhaseArcLengthCorrector.State.from(first);
    TwoToThreePhaseArcLengthCorrector.State secondState = TwoToThreePhaseArcLengthCorrector.State.from(second);
    CurvatureState firstCurvature = curvature(firstState, maximumResidual(first));
    CurvatureState secondCurvature = curvature(secondState, maximumResidual(second));
    if (!oppositeSigns(firstCurvature.minimumEigenvalue, secondCurvature.minimumEigenvalue)) {
      return Result.failure(firstCurvature, secondCurvature, "retained-phase minimum eigenvalues do not bracket zero");
    }

    CurvatureState lowerPressure = firstCurvature.state.getPressureBara() <= secondCurvature.state.getPressureBara()
        ? firstCurvature
        : secondCurvature;
    CurvatureState upperPressure = lowerPressure == firstCurvature ? secondCurvature : firstCurvature;
    TwoToThreePhaseArcLengthCorrector corrector = new TwoToThreePhaseArcLengthCorrector(template, retainedPhaseZero,
        retainedPhaseOne, incipientPhase)
        .setNumericalControls(boundaryMaximumIterations, boundaryResidualTolerance, boundaryFiniteDifferenceStep);
    CurvatureState best = Math.abs(lowerPressure.minimumEigenvalue) <= Math.abs(upperPressure.minimumEigenvalue)
        ? lowerPressure
        : upperPressure;
    int iterations = 0;
    String failureMessage = null;
    while (iterations < maximumIterations
        && upperPressure.state.getPressureBara() - lowerPressure.state.getPressureBara() > pressureBracketToleranceBara
        && Math.abs(best.minimumEigenvalue) > eigenvalueTolerance) {
      iterations++;
      double midpointPressure = 0.5 * (lowerPressure.state.getPressureBara() + upperPressure.state.getPressureBara());
      CurvatureState nearest = midpointPressure
          - lowerPressure.state.getPressureBara() <= upperPressure.state.getPressureBara() - midpointPressure
              ? lowerPressure
              : upperPressure;
      TwoToThreePhaseArcLengthCorrector.Result correction = corrector.correctAtPressure(nearest.state,
          midpointPressure);
      if (!correction.isConverged()) {
        CurvatureState alternate = nearest == lowerPressure ? upperPressure : lowerPressure;
        correction = corrector.correctAtPressure(alternate.state, midpointPressure);
      }
      if (!correction.isConverged()) {
        failureMessage = "fixed-pressure boundary correction failed during spinodal bisection: "
            + correction.getFailureMessage();
        break;
      }
      CurvatureState midpoint = curvature(correction.getState(), correction.getThermodynamicMaximumResidual());
      if (Math.abs(midpoint.minimumEigenvalue) < Math.abs(best.minimumEigenvalue)) {
        best = midpoint;
      }
      if (oppositeSigns(lowerPressure.minimumEigenvalue, midpoint.minimumEigenvalue)) {
        upperPressure = midpoint;
      } else if (oppositeSigns(midpoint.minimumEigenvalue, upperPressure.minimumEigenvalue)) {
        lowerPressure = midpoint;
      } else {
        failureMessage = "retained-phase eigenvalue lost its zero bracket during spinodal bisection";
        break;
      }
    }
    double bracketWidth = upperPressure.state.getPressureBara() - lowerPressure.state.getPressureBara();
    boolean converged = failureMessage == null
        && (bracketWidth <= pressureBracketToleranceBara || Math.abs(best.minimumEigenvalue) <= eigenvalueTolerance);
    if (!converged && failureMessage == null) {
      failureMessage = "retained-phase spinodal iteration limit reached";
    }
    return new Result(lowerPressure, upperPressure, best, iterations, bracketWidth, converged, failureMessage);
  }

  /** Evaluates the homogeneous retained-phase TPD curvature at one boundary root without refinement. */
  public CurvatureState analyze(TwoToThreePhaseBoundaryPointSolver.Result root) {
    validateRoot(root);
    return curvature(TwoToThreePhaseArcLengthCorrector.State.from(root), maximumResidual(root));
  }

  private CurvatureState curvature(TwoToThreePhaseArcLengthCorrector.State boundaryState, double boundaryResidual) {
    double[] bifurcatingComposition = bifurcatingPhase == retainedPhaseZero ? boundaryState.getPhaseZeroComposition()
        : boundaryState.getPhaseOneComposition();
    TwoToThreePhaseArcLengthCorrector.State homogeneous = TwoToThreePhaseArcLengthCorrector.State.create(
        retainedPhaseZero, retainedPhaseOne, bifurcatingPhase, boundaryState.getTemperatureK(),
        boundaryState.getPressureBara(), boundaryState.getBeta(), boundaryState.getPhaseZeroComposition(),
        boundaryState.getPhaseOneComposition(), bifurcatingComposition);
    IncipientPhaseStationarityJacobianAnalyzer.Result curvature = new IncipientPhaseStationarityJacobianAnalyzer(
        template, bifurcatingPhase, bifurcatingPhase).setFiniteDifferenceStep(curvatureFiniteDifferenceStep)
        .analyze(homogeneous);
    return new CurvatureState(boundaryState, curvature, boundaryResidual);
  }

  private static double maximumResidual(TwoToThreePhaseBoundaryPointSolver.Result root) {
    return Math.max(Math.abs(root.getRetainedFlashResidual()),
        Math.max(Math.abs(root.getTangentPlaneDistance()), Math.abs(root.getStationarityResidual())));
  }

  private void validateRoot(TwoToThreePhaseBoundaryPointSolver.Result root) {
    if (root == null || !root.isConverged() || root.getRetainedPhaseZero() != retainedPhaseZero
        || root.getRetainedPhaseOne() != retainedPhaseOne || root.getIncipientPhase() != incipientPhase) {
      throw new IllegalArgumentException("root does not match the configured two-to-three-phase boundary");
    }
  }

  private static boolean oppositeSigns(double first, double second) {
    return Double.isFinite(first) && Double.isFinite(second)
        && (first == 0.0 || second == 0.0 || Math.copySign(1.0, first) != Math.copySign(1.0, second));
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Boundary state with homogeneous retained-phase TPD curvature evidence. */
  public static final class CurvatureState {
    private final TwoToThreePhaseArcLengthCorrector.State state;
    private final IncipientPhaseStationarityJacobianAnalyzer.Result curvature;
    private final double minimumEigenvalue;
    private final double boundaryResidual;

    private CurvatureState(TwoToThreePhaseArcLengthCorrector.State state,
        IncipientPhaseStationarityJacobianAnalyzer.Result curvature, double boundaryResidual) {
      this.state = state;
      this.curvature = curvature;
      this.minimumEigenvalue = curvature.getBifurcationEigenvalue();
      this.boundaryResidual = boundaryResidual;
    }

    public TwoToThreePhaseArcLengthCorrector.State getState() {
      return state;
    }

    public IncipientPhaseStationarityJacobianAnalyzer.Result getCurvature() {
      return curvature;
    }

    public double getMinimumEigenvalue() {
      return minimumEigenvalue;
    }

    /** @return maximum thermodynamic residual of the corrected boundary state */
    public double getBoundaryResidual() {
      return boundaryResidual;
    }
  }

  /** Immutable retained-phase spinodal refinement result. */
  public static final class Result {
    private final CurvatureState lowerPressure;
    private final CurvatureState upperPressure;
    private final CurvatureState spinodal;
    private final int iterations;
    private final double pressureBracketWidthBara;
    private final boolean converged;
    private final String failureMessage;

    private Result(CurvatureState lowerPressure, CurvatureState upperPressure, CurvatureState spinodal, int iterations,
        double pressureBracketWidthBara, boolean converged, String failureMessage) {
      this.lowerPressure = lowerPressure;
      this.upperPressure = upperPressure;
      this.spinodal = spinodal;
      this.iterations = iterations;
      this.pressureBracketWidthBara = pressureBracketWidthBara;
      this.converged = converged;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CurvatureState first, CurvatureState second, String failureMessage) {
      CurvatureState best = Math.abs(first.minimumEigenvalue) <= Math.abs(second.minimumEigenvalue) ? first : second;
      return new Result(first, second, best, 0,
          Math.abs(second.state.getPressureBara() - first.state.getPressureBara()), false, failureMessage);
    }

    public boolean isConverged() {
      return converged && failureMessage == null;
    }

    public CurvatureState getLowerPressure() {
      return lowerPressure;
    }

    public CurvatureState getUpperPressure() {
      return upperPressure;
    }

    public CurvatureState getSpinodal() {
      return spinodal;
    }

    public int getIterations() {
      return iterations;
    }

    public double getPressureBracketWidthBara() {
      return pressureBracketWidthBara;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
