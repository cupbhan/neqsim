package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Separates ordinary two-to-three-phase continuation anchors from physical or unresolved endpoints.
 *
 * <p>
 * A zero third-phase TPD is not sufficient to define an ordinary boundary point. At a vanishing retained-phase fraction
 * the same equations also admit double-incidence and phase-coalescence solutions. Those states must be corrected with
 * the coupled three-phase equations and, when two phases coalesce, with independent mixture-criticality equations
 * before they receive a physical special-point label.
 * </p>
 */
public final class HydrocarbonWaterBoundaryEndpointClassifier {
  private final SystemInterface template;
  private double betaEndpointTolerance = 1.0e-6;
  private double distinctPhaseTolerance = 1.0e-5;
  private double residualTolerance = 1.0e-7;

  /** Creates a non-destructive endpoint classifier for one configured EOS and mixing rule. */
  public HydrocarbonWaterBoundaryEndpointClassifier(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets the phase-fraction, phase-distinctness, and coupled residual tolerances. */
  public HydrocarbonWaterBoundaryEndpointClassifier setTolerances(double betaEndpointTolerance,
      double distinctPhaseTolerance, double residualTolerance) {
    if (!positive(betaEndpointTolerance) || betaEndpointTolerance >= 0.5 || !positive(distinctPhaseTolerance)
        || !positive(residualTolerance)) {
      throw new IllegalArgumentException("invalid hydrocarbon-water endpoint tolerances");
    }
    this.betaEndpointTolerance = betaEndpointTolerance;
    this.distinctPhaseTolerance = distinctPhaseTolerance;
    this.residualTolerance = residualTolerance;
    return this;
  }

  /**
   * Classifies one strictly corrected zero-TPD root without ever promoting unresolved endpoint evidence.
   *
   * @param root specified retained-two-phase plus zero-third-phase-TPD solution
   * @return immutable classification and all independent correction evidence
   */
  public Result classify(TwoToThreePhaseBoundaryPointSolver.Result root) {
    if (root == null || !root.isConverged()) {
      throw new IllegalArgumentException("a converged strict boundary root is required");
    }
    double phaseZeroOneDistance = compositionDistance(root.getPhaseZeroComposition(), root.getPhaseOneComposition());
    double phaseZeroIncipientDistance = compositionDistance(root.getPhaseZeroComposition(),
        root.getIncipientComposition());
    double phaseOneIncipientDistance = compositionDistance(root.getPhaseOneComposition(),
        root.getIncipientComposition());
    double minimumDistance = Math.min(phaseZeroOneDistance,
        Math.min(phaseZeroIncipientDistance, phaseOneIncipientDistance));
    boolean endpointBeta = root.getBeta() <= betaEndpointTolerance || root.getBeta() >= 1.0 - betaEndpointTolerance;

    if (!endpointBeta && minimumDistance > distinctPhaseTolerance) {
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStability = new HydrocarbonWaterBoundaryGlobalStabilityGate(
          template).evaluate(root);
      if (!globalStability.isAccepted()) {
        return Result.metastable(root, phaseZeroOneDistance, phaseZeroIncipientDistance, phaseOneIncipientDistance,
            globalStability);
      }
      return Result.regular(root, phaseZeroOneDistance, phaseZeroIncipientDistance, phaseOneIncipientDistance,
          globalStability);
    }

    ThreePhasePointSolver.Result coupled = endpointBeta || minimumDistance <= distinctPhaseTolerance
        ? solveCoupledEndpoint(root)
        : null;
    if (coupled != null && coupled.isConverged()) {
      return Result.special(Classification.THREE_PHASE_POINT, SpecialPointType.THREE_PHASE_POINT, root,
          phaseZeroOneDistance, phaseZeroIncipientDistance, phaseOneIncipientDistance, coupled, null, null, null);
    }

    HydrocarbonWaterCriticalEndpointSolver.Result criticalZero = null;
    HydrocarbonWaterCriticalEndpointSolver.Result criticalOne = null;
    if (minimumDistance <= distinctPhaseTolerance
        || coupled != null && coupled.isResidualConverged() && !coupled.hasDistinctPhases()) {
      criticalZero = solveCriticalEndpoint(root, root.getRetainedPhaseZero());
      if (criticalZero != null && criticalZero.isPhysicalEndpoint()) {
        return Result.special(Classification.CRITICAL_END_POINT, SpecialPointType.CRITICAL_END_POINT, root,
            phaseZeroOneDistance, phaseZeroIncipientDistance, phaseOneIncipientDistance, coupled, criticalZero, null,
            null);
      }
      criticalOne = solveCriticalEndpoint(root, root.getRetainedPhaseOne());
      if (criticalOne != null && criticalOne.isPhysicalEndpoint()) {
        return Result.special(Classification.CRITICAL_END_POINT, SpecialPointType.CRITICAL_END_POINT, root,
            phaseZeroOneDistance, phaseZeroIncipientDistance, phaseOneIncipientDistance, coupled, criticalZero,
            criticalOne, null);
      }
      return Result.special(Classification.PHASE_COALESCENCE_CANDIDATE, SpecialPointType.PHASE_COALESCENCE_CANDIDATE,
          root, phaseZeroOneDistance, phaseZeroIncipientDistance, phaseOneIncipientDistance, coupled, criticalZero,
          criticalOne,
          "zero-TPD root contains duplicate phases but independent critical-endpoint conditions did not pass");
    }

    String failure = coupled == null ? "endpoint phase fraction has no coupled correction evidence"
        : coupled.getFailureMessage();
    return Result.special(Classification.UNRESOLVED_ENDPOINT, null, root, phaseZeroOneDistance,
        phaseZeroIncipientDistance, phaseOneIncipientDistance, coupled, null, null, failure);
  }

  private ThreePhasePointSolver.Result solveCoupledEndpoint(TwoToThreePhaseBoundaryPointSolver.Result root) {
    boolean phaseZeroIsMother = root.getBeta() >= 0.5;
    CandidatePhase mother = phaseZeroIsMother ? root.getRetainedPhaseZero() : root.getRetainedPhaseOne();
    CandidatePhase firstIncipient = phaseZeroIsMother ? root.getRetainedPhaseOne() : root.getRetainedPhaseZero();
    double[] firstComposition = phaseZeroIsMother ? root.getPhaseOneComposition() : root.getPhaseZeroComposition();
    return new ThreePhasePointSolver(template, mother, firstIncipient, root.getIncipientPhase())
        .setNumericalControls(80, residualTolerance, 2.0e-5)
        .solve(root.getTemperatureK(), root.getPressureBara(), firstComposition, root.getIncipientComposition());
  }

  private HydrocarbonWaterCriticalEndpointSolver.Result solveCriticalEndpoint(
      TwoToThreePhaseBoundaryPointSolver.Result root, CandidatePhase criticalPhase) {
    try {
      TwoToThreePhaseArcLengthCorrector.State state = TwoToThreePhaseArcLengthCorrector.State.create(
          root.getRetainedPhaseZero(), root.getRetainedPhaseOne(), root.getIncipientPhase(), root.getTemperatureK(),
          root.getPressureBara(), root.getBeta(), root.getPhaseZeroComposition(), root.getPhaseOneComposition(),
          root.getIncipientComposition());
      return new HydrocarbonWaterCriticalEndpointSolver(template, root.getRetainedPhaseZero(),
          root.getRetainedPhaseOne(), root.getIncipientPhase(), criticalPhase)
          .setNumericalControls(24, residualTolerance, Math.min(1.0e-9, residualTolerance), 1.0e-3, 2.0e-4,
              residualTolerance, 10.0 * distinctPhaseTolerance)
          .solve(state);
    } catch (RuntimeException error) {
      return null;
    }
  }

  private static double compositionDistance(double[] first, double[] second) {
    if (first.length != second.length) {
      return Double.POSITIVE_INFINITY;
    }
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Mutually exclusive physical/numerical classifications for a corrected zero-TPD root. */
  public enum Classification {
    REGULAR_BOUNDARY_POINT, THREE_PHASE_POINT, CRITICAL_END_POINT, PHASE_COALESCENCE_CANDIDATE,
    METASTABLE_BOUNDARY_ROOT, UNRESOLVED_ENDPOINT
  }

  /** Immutable classification with coupled and critical-equation audit evidence. */
  public static final class Result {
    private final Classification classification;
    private final SpecialPointType specialPointType;
    private final TwoToThreePhaseBoundaryPointSolver.Result boundaryRoot;
    private final double phaseZeroOneDistance;
    private final double phaseZeroIncipientDistance;
    private final double phaseOneIncipientDistance;
    private final ThreePhasePointSolver.Result coupledThreePhasePoint;
    private final HydrocarbonWaterCriticalEndpointSolver.Result criticalPhaseZeroResult;
    private final HydrocarbonWaterCriticalEndpointSolver.Result criticalPhaseOneResult;
    private final HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStabilityResult;
    private final String failureMessage;

    private Result(Classification classification, SpecialPointType specialPointType,
        TwoToThreePhaseBoundaryPointSolver.Result boundaryRoot, double phaseZeroOneDistance,
        double phaseZeroIncipientDistance, double phaseOneIncipientDistance,
        ThreePhasePointSolver.Result coupledThreePhasePoint,
        HydrocarbonWaterCriticalEndpointSolver.Result criticalPhaseZeroResult,
        HydrocarbonWaterCriticalEndpointSolver.Result criticalPhaseOneResult,
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStabilityResult, String failureMessage) {
      this.classification = classification;
      this.specialPointType = specialPointType;
      this.boundaryRoot = boundaryRoot;
      this.phaseZeroOneDistance = phaseZeroOneDistance;
      this.phaseZeroIncipientDistance = phaseZeroIncipientDistance;
      this.phaseOneIncipientDistance = phaseOneIncipientDistance;
      this.coupledThreePhasePoint = coupledThreePhasePoint;
      this.criticalPhaseZeroResult = criticalPhaseZeroResult;
      this.criticalPhaseOneResult = criticalPhaseOneResult;
      this.globalStabilityResult = globalStabilityResult;
      this.failureMessage = failureMessage;
    }

    private static Result regular(TwoToThreePhaseBoundaryPointSolver.Result root, double phaseZeroOneDistance,
        double phaseZeroIncipientDistance, double phaseOneIncipientDistance,
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStabilityResult) {
      return new Result(Classification.REGULAR_BOUNDARY_POINT, null, root, phaseZeroOneDistance,
          phaseZeroIncipientDistance, phaseOneIncipientDistance, null, null, null, globalStabilityResult, null);
    }

    private static Result metastable(TwoToThreePhaseBoundaryPointSolver.Result root, double phaseZeroOneDistance,
        double phaseZeroIncipientDistance, double phaseOneIncipientDistance,
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStabilityResult) {
      return new Result(Classification.METASTABLE_BOUNDARY_ROOT, null, root, phaseZeroOneDistance,
          phaseZeroIncipientDistance, phaseOneIncipientDistance, null, null, null, globalStabilityResult,
          globalStabilityResult.getFailureMessage());
    }

    private static Result special(Classification classification, SpecialPointType specialPointType,
        TwoToThreePhaseBoundaryPointSolver.Result root, double phaseZeroOneDistance, double phaseZeroIncipientDistance,
        double phaseOneIncipientDistance, ThreePhasePointSolver.Result coupledThreePhasePoint,
        HydrocarbonWaterCriticalEndpointSolver.Result criticalPhaseZeroResult,
        HydrocarbonWaterCriticalEndpointSolver.Result criticalPhaseOneResult, String failureMessage) {
      return new Result(classification, specialPointType, root, phaseZeroOneDistance, phaseZeroIncipientDistance,
          phaseOneIncipientDistance, coupledThreePhasePoint, criticalPhaseZeroResult, criticalPhaseOneResult, null,
          failureMessage);
    }

    public Classification getClassification() {
      return classification;
    }

    /** @return true only for a finite-fraction, pairwise-distinct ordinary continuation anchor */
    public boolean isOrdinaryBoundaryPoint() {
      return classification == Classification.REGULAR_BOUNDARY_POINT;
    }

    public SpecialPointType getSpecialPointType() {
      return specialPointType;
    }

    public TwoToThreePhaseBoundaryPointSolver.Result getBoundaryRoot() {
      return boundaryRoot;
    }

    public double getPhaseZeroOneDistance() {
      return phaseZeroOneDistance;
    }

    public double getPhaseZeroIncipientDistance() {
      return phaseZeroIncipientDistance;
    }

    public double getPhaseOneIncipientDistance() {
      return phaseOneIncipientDistance;
    }

    public double getMinimumPhaseDistance() {
      return Math.min(phaseZeroOneDistance, Math.min(phaseZeroIncipientDistance, phaseOneIncipientDistance));
    }

    public ThreePhasePointSolver.Result getCoupledThreePhasePoint() {
      return coupledThreePhasePoint;
    }

    public HydrocarbonWaterCriticalEndpointSolver.Result getCriticalPhaseZeroResult() {
      return criticalPhaseZeroResult;
    }

    public HydrocarbonWaterCriticalEndpointSolver.Result getCriticalPhaseOneResult() {
      return criticalPhaseOneResult;
    }

    /** @return global stability evidence for regular or metastable finite-fraction roots, otherwise null */
    public HydrocarbonWaterBoundaryGlobalStabilityGate.Result getGlobalStabilityResult() {
      return globalStabilityResult;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
