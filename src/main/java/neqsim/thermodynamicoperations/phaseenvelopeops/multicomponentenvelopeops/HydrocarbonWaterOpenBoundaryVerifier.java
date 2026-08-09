package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.Branch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EnvelopeReport;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EvidencePoint;

/** Independently verifies an open hydrocarbon-water boundary from a declared domain exit to a physical endpoint. */
public final class HydrocarbonWaterOpenBoundaryVerifier {
  private final SystemInterface template;

  /** Creates a verifier from the exact EOS, mixing rule, composition, and parameter template under test. */
  public HydrocarbonWaterOpenBoundaryVerifier(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /**
   * Re-solves every ordinary state and independently refines the high-pressure retained-phase spinodal.
   *
   * @param orderedStates states ordered from the low-pressure domain exit toward the physical endpoint
   * @param minimumPressureBara declared lower pressure of the verified computational domain
   * @param firstRejectedHighState first corrected boundary state rejected beyond the high-pressure endpoint
   * @return immutable point, endpoint, and branch-quality evidence
   */
  public Result verify(List<State> orderedStates, double minimumPressureBara, State firstRejectedHighState) {
    if (orderedStates == null || orderedStates.size() < 3 || !Double.isFinite(minimumPressureBara)
        || minimumPressureBara <= 0.0 || firstRejectedHighState == null) {
      throw new IllegalArgumentException(
          "an ordered open branch, positive domain pressure, and rejected end state are required");
    }
    List<State> branchStates = new ArrayList<State>(orderedStates);
    boolean startsAtDomainExit = branchStates.get(0).getPressureBara() <= minimumPressureBara
        && branchStates.get(1).getPressureBara() > minimumPressureBara;
    int lastIndex = branchStates.size() - 1;
    boolean endsAtDomainExit = branchStates.get(lastIndex).getPressureBara() <= minimumPressureBara
        && branchStates.get(lastIndex - 1).getPressureBara() > minimumPressureBara;
    if (!startsAtDomainExit && endsAtDomainExit) {
      Collections.reverse(branchStates);
      startsAtDomainExit = true;
    }
    State reference = branchStates.get(0);
    boolean domainExitVerified = startsAtDomainExit;
    if (!domainExitVerified) {
      throw new IllegalArgumentException("the first two states do not bracket the declared low-pressure domain exit");
    }
    BoundaryFamily family = boundaryFamily(reference);
    HydrocarbonWaterBoundaryGlobalStabilityGate globalGate = new HydrocarbonWaterBoundaryGlobalStabilityGate(template);
    List<PointResult> pointResults = new ArrayList<PointResult>();
    List<EvidencePoint> evidence = new ArrayList<EvidencePoint>();
    int acceptedPointCount = 0;
    double maximumRetainedFlashResidual = 0.0;
    double maximumTargetTangentPlaneDistance = 0.0;
    double maximumTargetStationarityResidual = 0.0;
    double minimumNonTrivialTangentPlaneDistance = Double.POSITIVE_INFINITY;

    for (int index = 0; index < branchStates.size(); index++) {
      State state = branchStates.get(index);
      if (!sameTopology(reference, state)) {
        throw new IllegalArgumentException("every open-branch state must share one hydrocarbon-water topology");
      }
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result global = globalGate
          .evaluate(TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(state, 0.0));
      EvidencePoint point = evidence(state, global.getRetainedFlash(), global.getTargetCandidate());
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

    if (!sameTopology(reference, firstRejectedHighState)) {
      throw new IllegalArgumentException("the rejected endpoint state belongs to a different boundary topology");
    }
    TwoToThreePhaseBoundaryPointSolver.Result acceptedRoot = TwoToThreePhaseBoundaryPointSolver.Result
        .fromContinuationState(branchStates.get(branchStates.size() - 1), 0.0);
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejectedEvidence = globalGate
        .evaluate(TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(firstRejectedHighState, 0.0));
    if (rejectedEvidence.isAccepted()) {
      throw new IllegalArgumentException("the supplied high-pressure endpoint state is not independently rejected");
    }
    HydrocarbonWaterBoundaryTerminationClassifier.Result termination = new HydrocarbonWaterBoundaryTerminationClassifier(
        template).classify(acceptedRoot, rejectedEvidence);
    if (!termination.isPhysicalEndpoint()
        || termination.getType() != HydrocarbonWaterBoundaryTerminationClassifier.Type.RETAINED_PHASE_SPINODAL) {
      throw new IllegalArgumentException(
          "the high-pressure endpoint did not refine to a retained-phase spinodal: " + termination.getDiagnostic());
    }
    TwoToThreePhaseBoundaryPointSolver.Result endpointRoot = termination.getBoundaryRoot();
    EvidencePoint endpointEvidence = evidence(endpointRoot);
    evidence.add(endpointEvidence);
    maximumRetainedFlashResidual = maximumFinite(maximumRetainedFlashResidual,
        endpointEvidence.getRetainedFlashResidual());
    maximumTargetTangentPlaneDistance = maximumFinite(maximumTargetTangentPlaneDistance,
        Math.abs(endpointEvidence.getTangentPlaneDistance()));
    maximumTargetStationarityResidual = maximumFinite(maximumTargetStationarityResidual,
        endpointEvidence.getStationarityResidual());

    Branch branch = new Branch("independently-verified-open-hydrocarbon-water-branch", family.getDefinition(),
        family.getRetainedPhaseZero(), family.getRetainedPhaseOne(), family.getIncipientPhase(), evidence,
        SpecialPointType.DOMAIN_EXIT, SpecialPointType.RETAINED_PHASE_SPINODAL, null, termination);
    EnvelopeReport quality = new TwoToThreePhaseBoundaryQualityGate().validate(Collections.singletonList(branch));
    boolean allOrdinaryPointsAccepted = acceptedPointCount == branchStates.size();
    boolean internalQualityEligible = domainExitVerified && allOrdinaryPointsAccepted && quality.isAccepted();
    return new Result(pointResults, evidence, quality, acceptedPointCount, allOrdinaryPointsAccepted,
        domainExitVerified, minimumPressureBara, termination, rejectedEvidence, maximumRetainedFlashResidual,
        maximumTargetTangentPlaneDistance, maximumTargetStationarityResidual, minimumNonTrivialTangentPlaneDistance,
        internalQualityEligible);
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

  private static EvidencePoint evidence(TwoToThreePhaseBoundaryPointSolver.Result root) {
    return new EvidencePoint(root.getTemperatureK(), root.getPressureBara(), root.getBeta(),
        root.getPhaseZeroComposition(), root.getPhaseOneComposition(), root.getIncipientComposition(),
        root.getRetainedFlashResidual(), root.getTangentPlaneDistance(), root.getStationarityResidual());
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

  private static boolean sameTopology(State first, State second) {
    return first != null && second != null && first.getRetainedPhaseZero() == second.getRetainedPhaseZero()
        && first.getRetainedPhaseOne() == second.getRetainedPhaseOne()
        && first.getIncipientPhase() == second.getIncipientPhase();
  }

  private static double maximumFinite(double current, double candidate) {
    return Double.isFinite(current) && Double.isFinite(candidate) ? Math.max(current, candidate)
        : Double.POSITIVE_INFINITY;
  }

  /** Independent evidence for one ordinary input state. */
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

  /** Full open-branch point, endpoint, and topology-quality inventory. */
  public static final class Result {
    private final List<PointResult> pointResults;
    private final List<EvidencePoint> evidence;
    private final EnvelopeReport qualityReport;
    private final int acceptedPointCount;
    private final boolean allOrdinaryPointsAccepted;
    private final boolean domainExitVerified;
    private final double minimumPressureBara;
    private final HydrocarbonWaterBoundaryTerminationClassifier.Result highTermination;
    private final HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejectedHighEvidence;
    private final double maximumRetainedFlashResidual;
    private final double maximumTargetTangentPlaneDistance;
    private final double maximumTargetStationarityResidual;
    private final double minimumNonTrivialTangentPlaneDistance;
    private final boolean internalQualityEligible;

    private Result(List<PointResult> pointResults, List<EvidencePoint> evidence, EnvelopeReport qualityReport,
        int acceptedPointCount, boolean allOrdinaryPointsAccepted, boolean domainExitVerified,
        double minimumPressureBara, HydrocarbonWaterBoundaryTerminationClassifier.Result highTermination,
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejectedHighEvidence, double maximumRetainedFlashResidual,
        double maximumTargetTangentPlaneDistance, double maximumTargetStationarityResidual,
        double minimumNonTrivialTangentPlaneDistance, boolean internalQualityEligible) {
      this.pointResults = Collections.unmodifiableList(new ArrayList<PointResult>(pointResults));
      this.evidence = Collections.unmodifiableList(new ArrayList<EvidencePoint>(evidence));
      this.qualityReport = qualityReport;
      this.acceptedPointCount = acceptedPointCount;
      this.allOrdinaryPointsAccepted = allOrdinaryPointsAccepted;
      this.domainExitVerified = domainExitVerified;
      this.minimumPressureBara = minimumPressureBara;
      this.highTermination = highTermination;
      this.rejectedHighEvidence = rejectedHighEvidence;
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

    public int getAcceptedPointCount() {
      return acceptedPointCount;
    }

    public boolean areAllOrdinaryPointsAccepted() {
      return allOrdinaryPointsAccepted;
    }

    public boolean isDomainExitVerified() {
      return domainExitVerified;
    }

    public double getMinimumPressureBara() {
      return minimumPressureBara;
    }

    public HydrocarbonWaterBoundaryTerminationClassifier.Result getHighTermination() {
      return highTermination;
    }

    public HydrocarbonWaterBoundaryGlobalStabilityGate.Result getRejectedHighEvidence() {
      return rejectedHighEvidence;
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
