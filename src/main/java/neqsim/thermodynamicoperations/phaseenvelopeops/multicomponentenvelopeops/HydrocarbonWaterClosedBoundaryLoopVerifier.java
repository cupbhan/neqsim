package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.Branch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EnvelopeReport;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EvidencePoint;

/** Independently re-flashes and stability-verifies every unique state of one closed hydrocarbon-water boundary. */
public final class HydrocarbonWaterClosedBoundaryLoopVerifier {
  private final SystemInterface template;

  /** Creates a verifier from a fully configured EOS and mixing-rule template. */
  public HydrocarbonWaterClosedBoundaryLoopVerifier(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /**
   * Re-solves the retained flash and all physical stationary-phase trials at every unique loop state.
   *
   * <p>
   * The final state of a closed loop is an exact duplicate of the first and reuses its independent evidence. This
   * preserves an exact topological closure without performing the same thermodynamic calculation twice.
   * </p>
   */
  public Result verify(List<State> closedLoop) {
    if (closedLoop == null || closedLoop.size() < 4
        || !sameState(closedLoop.get(0), closedLoop.get(closedLoop.size() - 1))) {
      throw new IllegalArgumentException("verification requires a state- and identity-closed loop");
    }
    BoundaryFamily family = boundaryFamily(closedLoop.get(0));
    List<PointResult> pointResults = new ArrayList<PointResult>();
    List<EvidencePoint> evidence = new ArrayList<EvidencePoint>();
    int acceptedPointCount = 0;
    double maximumRetainedFlashResidual = 0.0;
    double maximumTargetTangentPlaneDistance = 0.0;
    double maximumTargetStationarityResidual = 0.0;
    double minimumNonTrivialTangentPlaneDistance = Double.POSITIVE_INFINITY;
    HydrocarbonWaterBoundaryGlobalStabilityGate globalGate = new HydrocarbonWaterBoundaryGlobalStabilityGate(template);

    for (int index = 0; index < closedLoop.size() - 1; index++) {
      State state = closedLoop.get(index);
      if (!sameTopology(closedLoop.get(0), state)) {
        throw new IllegalArgumentException("every loop state must share one hydrocarbon-water topology");
      }
      TwoToThreePhaseBoundaryPointSolver.Result root = TwoToThreePhaseBoundaryPointSolver.Result
          .fromContinuationState(state, 0.0);
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result global = globalGate.evaluate(root);
      SpecifiedTwoPhaseFlashSolver.Result retained = global.getRetainedFlash();
      Candidate target = global.getTargetCandidate();
      EvidencePoint point = evidence(state, retained, target);
      PointResult result = new PointResult(index, state, global, point);
      pointResults.add(result);
      evidence.add(point);
      if (result.isAccepted()) {
        acceptedPointCount++;
      }
      maximumRetainedFlashResidual = maximumFinite(maximumRetainedFlashResidual, point.getRetainedFlashResidual());
      maximumTargetTangentPlaneDistance = maximumFinite(maximumTargetTangentPlaneDistance,
          Math.abs(point.getTangentPlaneDistance()));
      maximumTargetStationarityResidual = maximumFinite(maximumTargetStationarityResidual,
          point.getStationarityResidual());
      minimumNonTrivialTangentPlaneDistance = Math.min(minimumNonTrivialTangentPlaneDistance,
          global.getMinimumNonTrivialTangentPlaneDistance());
    }
    evidence.add(evidence.get(0));
    Branch branch = new Branch("independently-verified-closed-hydrocarbon-water-loop", family.getDefinition(),
        family.getRetainedPhaseZero(), family.getRetainedPhaseOne(), family.getIncipientPhase(), evidence, null, null);
    EnvelopeReport quality = new TwoToThreePhaseBoundaryQualityGate().validate(Collections.singletonList(branch));
    boolean allUniquePointsAccepted = acceptedPointCount == closedLoop.size() - 1;
    return new Result(pointResults, evidence, quality, acceptedPointCount, allUniquePointsAccepted,
        maximumRetainedFlashResidual, maximumTargetTangentPlaneDistance, maximumTargetStationarityResidual,
        minimumNonTrivialTangentPlaneDistance, allUniquePointsAccepted && quality.isAccepted());
  }

  private static EvidencePoint evidence(State state, SpecifiedTwoPhaseFlashSolver.Result retained, Candidate target) {
    if (retained == null || !retained.isConverged() || target == null || !target.isConverged()) {
      return new EvidencePoint(state.getTemperatureK(), state.getPressureBara(), state.getBeta(),
          state.getPhaseZeroComposition(), state.getPhaseOneComposition(), state.getIncipientComposition(),
          Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
    }
    return new EvidencePoint(state.getTemperatureK(), state.getPressureBara(), retained.getBeta(),
        retained.getPhaseZeroComposition(), retained.getPhaseOneComposition(), target.getComposition(),
        retained.getMaximumResidual(), target.getTangentPlaneDistance(), target.getStationarityResidual());
  }

  private static BoundaryFamily boundaryFamily(State state) {
    for (BoundaryFamily family : BoundaryFamily.values()) {
      if (family.getRetainedPhaseZero() == state.getRetainedPhaseZero()
          && family.getRetainedPhaseOne() == state.getRetainedPhaseOne()
          && family.getIncipientPhase() == state.getIncipientPhase()) {
        return family;
      }
    }
    throw new IllegalArgumentException("unsupported hydrocarbon-water continuation topology");
  }

  private static boolean sameState(State first, State second) {
    return sameTopology(first, second) && Math.abs(first.getTemperatureK() - second.getTemperatureK()) <= 1.0e-7
        && Math.abs(first.getPressureBara() - second.getPressureBara()) <= 1.0e-7
        && compositionDistance(first.getPhaseZeroComposition(), second.getPhaseZeroComposition()) <= 1.0e-7
        && compositionDistance(first.getPhaseOneComposition(), second.getPhaseOneComposition()) <= 1.0e-7
        && compositionDistance(first.getIncipientComposition(), second.getIncipientComposition()) <= 1.0e-7;
  }

  private static boolean sameTopology(State first, State second) {
    return first != null && second != null && first.getRetainedPhaseZero() == second.getRetainedPhaseZero()
        && first.getRetainedPhaseOne() == second.getRetainedPhaseOne()
        && first.getIncipientPhase() == second.getIncipientPhase();
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

  private static double maximumFinite(double current, double candidate) {
    return Double.isFinite(current) && Double.isFinite(candidate) ? Math.max(current, candidate)
        : Double.POSITIVE_INFINITY;
  }

  /** Independent evidence and verdict for one unique loop state. */
  public static final class PointResult {
    private final int index;
    private final State inputState;
    private final HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStability;
    private final EvidencePoint evidence;

    private PointResult(int index, State inputState, HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStability,
        EvidencePoint evidence) {
      this.index = index;
      this.inputState = inputState;
      this.globalStability = globalStability;
      this.evidence = evidence;
    }

    public int getIndex() {
      return index;
    }

    public State getInputState() {
      return inputState;
    }

    public HydrocarbonWaterBoundaryGlobalStabilityGate.Result getGlobalStability() {
      return globalStability;
    }

    public EvidencePoint getEvidence() {
      return evidence;
    }

    public boolean isAccepted() {
      return globalStability != null && globalStability.isAccepted()
          && Double.isFinite(evidence.getRetainedFlashResidual()) && Double.isFinite(evidence.getTangentPlaneDistance())
          && Double.isFinite(evidence.getStationarityResidual());
    }
  }

  /** Full-loop independent evidence inventory plus the existing topology and continuity quality-gate report. */
  public static final class Result {
    private final List<PointResult> pointResults;
    private final List<EvidencePoint> evidence;
    private final EnvelopeReport qualityReport;
    private final int acceptedUniquePointCount;
    private final boolean allUniquePointsAccepted;
    private final double maximumRetainedFlashResidual;
    private final double maximumTargetTangentPlaneDistance;
    private final double maximumTargetStationarityResidual;
    private final double minimumNonTrivialTangentPlaneDistance;
    private final boolean internalQualityEligible;

    private Result(List<PointResult> pointResults, List<EvidencePoint> evidence, EnvelopeReport qualityReport,
        int acceptedUniquePointCount, boolean allUniquePointsAccepted, double maximumRetainedFlashResidual,
        double maximumTargetTangentPlaneDistance, double maximumTargetStationarityResidual,
        double minimumNonTrivialTangentPlaneDistance, boolean internalQualityEligible) {
      this.pointResults = Collections.unmodifiableList(new ArrayList<PointResult>(pointResults));
      this.evidence = Collections.unmodifiableList(new ArrayList<EvidencePoint>(evidence));
      this.qualityReport = qualityReport;
      this.acceptedUniquePointCount = acceptedUniquePointCount;
      this.allUniquePointsAccepted = allUniquePointsAccepted;
      this.maximumRetainedFlashResidual = maximumRetainedFlashResidual;
      this.maximumTargetTangentPlaneDistance = maximumTargetTangentPlaneDistance;
      this.maximumTargetStationarityResidual = maximumTargetStationarityResidual;
      this.minimumNonTrivialTangentPlaneDistance = minimumNonTrivialTangentPlaneDistance;
      this.internalQualityEligible = internalQualityEligible;
    }

    public List<PointResult> getPointResults() {
      return pointResults;
    }

    public List<EvidencePoint> getEvidence() {
      return evidence;
    }

    public EnvelopeReport getQualityReport() {
      return qualityReport;
    }

    public int getAcceptedUniquePointCount() {
      return acceptedUniquePointCount;
    }

    public boolean areAllUniquePointsAccepted() {
      return allUniquePointsAccepted;
    }

    public double getMaximumRetainedFlashResidual() {
      return maximumRetainedFlashResidual;
    }

    public double getMaximumTargetTangentPlaneDistance() {
      return maximumTargetTangentPlaneDistance;
    }

    public double getMaximumTargetStationarityResidual() {
      return maximumTargetStationarityResidual;
    }

    public double getMinimumNonTrivialTangentPlaneDistance() {
      return minimumNonTrivialTangentPlaneDistance;
    }

    public boolean isInternalQualityEligible() {
      return internalQualityEligible;
    }
  }
}
