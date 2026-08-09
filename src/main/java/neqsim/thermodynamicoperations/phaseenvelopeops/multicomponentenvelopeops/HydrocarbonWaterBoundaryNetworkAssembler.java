package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.EndpointCandidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryEndpointClassifier.Classification;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EvidencePoint;

/**
 * Joins corrected hydrocarbon-water boundary branches to independently proven physical endpoints.
 *
 * <p>
 * The assembler is deliberately stricter than a geometric line join. A branch end can be attached only when the
 * endpoint supports the same {@link BoundaryFamily}, its temperature and pressure are within configured transformed
 * distances, and the endpoint carries converged solver evidence. Retained-phase spinodals are taken only from the
 * branch termination classifier. Isolated branches, unresolved endpoint candidates, duplicate endpoint evidence, and
 * incomplete endpoint degrees are preserved as network violations instead of being discarded.
 * </p>
 */
public final class HydrocarbonWaterBoundaryNetworkAssembler {
  private double maximumRelativeTemperatureDistance = 0.05;
  private double maximumLogPressureDistance = 0.25;
  private double duplicateRelativeTemperatureTolerance = 1.0e-6;
  private double duplicateLogPressureTolerance = 1.0e-6;

  /** Sets endpoint attachment and duplicate-detection tolerances in transformed PT space. */
  public HydrocarbonWaterBoundaryNetworkAssembler setAttachmentTolerances(double maximumRelativeTemperatureDistance,
      double maximumLogPressureDistance, double duplicateRelativeTemperatureTolerance,
      double duplicateLogPressureTolerance) {
    if (!positive(maximumRelativeTemperatureDistance) || !positive(maximumLogPressureDistance)
        || !positive(duplicateRelativeTemperatureTolerance) || !positive(duplicateLogPressureTolerance)) {
      throw new IllegalArgumentException("invalid hydrocarbon-water network attachment tolerances");
    }
    this.maximumRelativeTemperatureDistance = maximumRelativeTemperatureDistance;
    this.maximumLogPressureDistance = maximumLogPressureDistance;
    this.duplicateRelativeTemperatureTolerance = duplicateRelativeTemperatureTolerance;
    this.duplicateLogPressureTolerance = duplicateLogPressureTolerance;
    return this;
  }

  /**
   * Builds and validates a complete branch network.
   *
   * @param assembly preserved branch assembly and withheld endpoint candidates
   * @param independentEndpoints independently solved CEP or simultaneous-incidence evidence
   * @return immutable branch network, endpoint attachments, quality report, and every violation
   */
  public Result assemble(HydrocarbonWaterBoundaryBranchAssembler.Result assembly,
      List<EndpointEvidence> independentEndpoints) {
    if (assembly == null || independentEndpoints == null) {
      throw new IllegalArgumentException("branch assembly and endpoint evidence list are required");
    }
    List<String> violations = new ArrayList<String>();
    List<EndpointEvidence> endpoints = new ArrayList<EndpointEvidence>();
    for (EndpointEvidence endpoint : independentEndpoints) {
      if (endpoint == null) {
        violations.add("NULL_ENDPOINT_EVIDENCE");
      } else {
        endpoints.add(endpoint.copyWithoutAttachments());
      }
    }
    preserveEndpointCandidates(assembly.getEndpointCandidates(), endpoints, violations);
    detectDuplicateEndpointEvidence(endpoints, violations);

    List<MutableBranch> mutableBranches = new ArrayList<MutableBranch>();
    Set<String> branchIdentifiers = new HashSet<String>();
    for (AssembledBranch branch : assembly.getBranches()) {
      if (branch == null) {
        violations.add("NULL_ASSEMBLED_BRANCH");
        continue;
      }
      if (!branchIdentifiers.add(branch.getIdentifier())) {
        violations.add("DUPLICATE_BRANCH_ID:" + branch.getIdentifier());
      }
      MutableBranch mutable = new MutableBranch(branch);
      mutableBranches.add(mutable);
      attachClassifiedTermination(mutable, End.START, branch.getBackwardTermination(), endpoints, violations);
      attachClassifiedTermination(mutable, End.END, branch.getForwardTermination(), endpoints, violations);
    }

    attachTargetBranchMerges(mutableBranches, endpoints, violations);

    List<AttachmentTrial> trials = discoverAttachmentTrials(mutableBranches, endpoints);
    Collections.sort(trials, Comparator.comparingDouble(AttachmentTrial::getDistance));
    for (AttachmentTrial trial : trials) {
      if (trial.branch.isAttached(trial.end)
          || trial.endpoint.getAttachmentCount() >= trial.endpoint.getRequiredDegree()
          || trial.endpoint.isAttachedTo(trial.branch.branch.getIdentifier(), trial.end)) {
        continue;
      }
      trial.branch.attach(trial.end, trial.endpoint, trial.distance);
      trial.endpoint.recordAttachment(trial.branch.branch.getIdentifier(), trial.branch.branch.getFamily(), trial.end);
    }

    validateEndpointDegrees(endpoints, mutableBranches, violations);
    List<NetworkBranch> branches = new ArrayList<NetworkBranch>();
    List<TwoToThreePhaseBoundaryQualityGate.Branch> gateBranches = new ArrayList<TwoToThreePhaseBoundaryQualityGate.Branch>();
    for (MutableBranch mutable : mutableBranches) {
      NetworkBranch branch = mutable.freeze();
      branches.add(branch);
      if (branch.isIsolated()) {
        violations.add("ISOLATED_BRANCH:" + branch.getIdentifier());
        continue;
      }
      if (branch.getStartAttachment() == null) {
        violations.add("UNATTACHED_BRANCH_ENDPOINT:" + branch.getIdentifier() + ":START");
      }
      if (branch.getEndAttachment() == null) {
        violations.add("UNATTACHED_BRANCH_ENDPOINT:" + branch.getIdentifier() + ":END");
      }
      if (branch.getStartAttachment() != null && branch.getEndAttachment() != null && branch.getStartAttachment()
          .getEndpointIdentifier().equals(branch.getEndAttachment().getEndpointIdentifier())) {
        violations.add("SAME_ENDPOINT_ON_BOTH_BRANCH_ENDS:" + branch.getIdentifier());
      }
      gateBranches.add(branch.toQualityGateBranch());
    }

    TwoToThreePhaseBoundaryQualityGate.EnvelopeReport qualityReport = null;
    if (gateBranches.isEmpty()) {
      violations.add("NO_CONTINUATION_BRANCH");
    } else {
      qualityReport = new TwoToThreePhaseBoundaryQualityGate().validate(gateBranches);
      if (!qualityReport.isAccepted()) {
        violations.add("BOUNDARY_QUALITY_GATE_REJECTED");
      }
    }
    boolean engineeringEligible = violations.isEmpty() && qualityReport != null && qualityReport.isAccepted();
    return new Result(branches, endpoints, assembly.getEndpointCandidates().size(), qualityReport, violations,
        engineeringEligible);
  }

  private void preserveEndpointCandidates(List<EndpointCandidate> candidates, List<EndpointEvidence> endpoints,
      List<String> violations) {
    for (int candidateIndex = 0; candidateIndex < candidates.size(); candidateIndex++) {
      EndpointCandidate candidate = candidates.get(candidateIndex);
      if (candidate == null) {
        violations.add("NULL_WITHHELD_ENDPOINT_CANDIDATE");
        continue;
      }
      HydrocarbonWaterBoundaryEndpointClassifier.Result classification = candidate.getClassification();
      EndpointEvidence evidence = endpointEvidence(candidate, candidateIndex);
      if (evidence == null) {
        violations.add("UNRESOLVED_ENDPOINT_CANDIDATE:" + candidate.getFamily().name() + ":"
            + classification.getClassification().name());
      } else {
        endpoints.add(evidence);
      }
    }
  }

  private static EndpointEvidence endpointEvidence(EndpointCandidate candidate, int candidateIndex) {
    HydrocarbonWaterBoundaryEndpointClassifier.Result classification = candidate.getClassification();
    if (classification.getClassification() == Classification.CRITICAL_END_POINT) {
      HydrocarbonWaterCriticalEndpointSolver.Result critical = physicalCriticalResult(classification);
      if (critical != null) {
        return EndpointEvidence.criticalEndpoint("candidate-cep-" + candidate.getFamily().name() + "-" + candidateIndex,
            candidate.getFamily(), critical, EvidenceSource.ENDPOINT_CLASSIFIER);
      }
    }
    if (classification.getClassification() == Classification.THREE_PHASE_POINT
        && classification.getCoupledThreePhasePoint() != null
        && classification.getCoupledThreePhasePoint().isConverged()) {
      return EndpointEvidence.threePhasePoint(
          "candidate-three-phase-" + candidate.getFamily().name() + "-" + candidateIndex,
          classification.getCoupledThreePhasePoint(), EvidencePoint.from(classification.getBoundaryRoot()),
          EvidenceSource.ENDPOINT_CLASSIFIER);
    }
    return null;
  }

  private static HydrocarbonWaterCriticalEndpointSolver.Result physicalCriticalResult(
      HydrocarbonWaterBoundaryEndpointClassifier.Result classification) {
    HydrocarbonWaterCriticalEndpointSolver.Result zero = classification.getCriticalPhaseZeroResult();
    if (zero != null && zero.isPhysicalEndpoint()) {
      return zero;
    }
    HydrocarbonWaterCriticalEndpointSolver.Result one = classification.getCriticalPhaseOneResult();
    return one != null && one.isPhysicalEndpoint() ? one : null;
  }

  private void detectDuplicateEndpointEvidence(List<EndpointEvidence> endpoints, List<String> violations) {
    for (int firstIndex = 0; firstIndex < endpoints.size(); firstIndex++) {
      EndpointEvidence first = endpoints.get(firstIndex);
      for (int secondIndex = firstIndex + 1; secondIndex < endpoints.size(); secondIndex++) {
        EndpointEvidence second = endpoints.get(secondIndex);
        if (first.type == second.type && first.compatibleFamilies.equals(second.compatibleFamilies)
            && relativeTemperatureDistance(first.temperatureK,
                second.temperatureK) <= duplicateRelativeTemperatureTolerance
            && logPressureDistance(first.pressureBara, second.pressureBara) <= duplicateLogPressureTolerance) {
          violations.add("DUPLICATE_ENDPOINT_EVIDENCE:" + first.identifier + ":" + second.identifier);
        }
      }
    }
  }

  private void attachClassifiedTermination(MutableBranch branch, End end,
      HydrocarbonWaterBoundaryTerminationClassifier.Result termination, List<EndpointEvidence> endpoints,
      List<String> violations) {
    if (termination == null || !termination.isPhysicalEndpoint()) {
      return;
    }
    if (termination.getType() != HydrocarbonWaterBoundaryTerminationClassifier.Type.RETAINED_PHASE_SPINODAL) {
      violations.add("UNSUPPORTED_PHYSICAL_TERMINATION:" + branch.branch.getIdentifier() + ":" + end.name());
      return;
    }
    EndpointEvidence endpoint = EndpointEvidence.spinodal(
        branch.branch.getIdentifier() + "-" + end.name().toLowerCase() + "-spinodal", branch.branch.getFamily(),
        termination);
    endpoints.add(endpoint);
    branch.attach(end, endpoint, 0.0);
    endpoint.recordAttachment(branch.branch.getIdentifier(), branch.branch.getFamily(), end);
  }

  private List<AttachmentTrial> discoverAttachmentTrials(List<MutableBranch> branches,
      List<EndpointEvidence> endpoints) {
    List<AttachmentTrial> trials = new ArrayList<AttachmentTrial>();
    for (EndpointEvidence endpoint : endpoints) {
      if (endpoint.source == EvidenceSource.TERMINATION_CLASSIFIER || endpoint.requestedBranchEnd != null) {
        continue;
      }
      for (MutableBranch branch : branches) {
        if (branch.branch.isIsolated() || !endpoint.compatibleFamilies.contains(branch.branch.getFamily())) {
          continue;
        }
        addTrial(trials, branch, End.START, endpoint);
        addTrial(trials, branch, End.END, endpoint);
      }
    }
    return trials;
  }

  private void attachTargetBranchMerges(List<MutableBranch> branches, List<EndpointEvidence> endpoints,
      List<String> violations) {
    List<TargetBranchMergeTrial> trials = new ArrayList<TargetBranchMergeTrial>();
    for (EndpointEvidence endpoint : endpoints) {
      if (endpoint.requestedBranchEnd == null) {
        continue;
      }
      for (MutableBranch branch : branches) {
        TargetBranchMergeTrial trial = targetBranchMergeTrial(branch, endpoint);
        if (trial != null) {
          trials.add(trial);
        }
      }
    }
    Collections.sort(trials, Comparator.comparingDouble(TargetBranchMergeTrial::getDistance));
    for (TargetBranchMergeTrial trial : trials) {
      if (trial.endpoint.getAttachmentCount() >= trial.endpoint.getRequiredDegree()
          || trial.branch.isAttached(trial.endpoint.requestedBranchEnd)) {
        continue;
      }
      trial.branch.applyTargetBranchMergeTrim(trial.trim, trial.endpoint, trial.distance);
      trial.endpoint.recordAttachment(trial.branch.branch.getIdentifier(), trial.branch.branch.getFamily(),
          trial.endpoint.requestedBranchEnd);
    }
  }

  private TargetBranchMergeTrial targetBranchMergeTrial(MutableBranch branch, EndpointEvidence endpoint) {
    if (branch.branch.isIsolated() || branch.isAttached(endpoint.requestedBranchEnd)
        || !endpoint.compatibleFamilies.contains(branch.branch.getFamily()) || endpoint.excludedSidePoint == null
        || endpoint.exactBoundaryPoint == null || branch.points.size() < 3) {
      return null;
    }
    TargetBranchMergeTrimResult trim = trimTargetBranchMerge(branch.points, endpoint);
    return trim.isAccepted() ? new TargetBranchMergeTrial(branch, endpoint, trim, trim.getTransformedDistance()) : null;
  }

  /**
   * Removes the stationary-identity side beyond one independently proven target-branch merge.
   *
   * <p>
   * The excluded-side sample identifies which of two locally adjacent roots lies on the metastable continuation. This
   * avoids selecting the wrong side of a tight PT fold merely because both geometric segments pass close to the merge
   * coordinate. The exact merged root replaces the retained endpoint after the excluded side is removed.
   * </p>
   *
   * @param orderedPoints continuation evidence in branch order
   * @param endpoint independently proven target-branch merge evidence
   * @return immutable retained evidence and trimming diagnostics
   */
  public TargetBranchMergeTrimResult trimTargetBranchMerge(List<EvidencePoint> orderedPoints,
      EndpointEvidence endpoint) {
    if (orderedPoints == null || endpoint == null || endpoint.type != SpecialPointType.TARGET_BRANCH_MERGE
        || endpoint.source != EvidenceSource.SECONDARY_STATIONARY_BRANCH_TRACKER || endpoint.requestedBranchEnd == null
        || endpoint.excludedSidePoint == null || endpoint.exactBoundaryPoint == null || orderedPoints.size() < 3) {
      return TargetBranchMergeTrimResult.rejected("invalid target-branch merge trimming evidence");
    }
    int excludedIndex = closestEvidencePointIndex(orderedPoints, endpoint.excludedSidePoint);
    int retainedIndex = endpoint.requestedBranchEnd == End.END ? excludedIndex - 1 : excludedIndex + 1;
    if (retainedIndex < 0 || retainedIndex >= orderedPoints.size()) {
      return TargetBranchMergeTrimResult.rejected("excluded-side sample does not leave a retained neighbor");
    }
    EvidencePoint retained = orderedPoints.get(retainedIndex);
    EvidencePoint excluded = orderedPoints.get(excludedIndex);
    double transformedPtDistance = pointToSegmentDistance(Math.log(endpoint.temperatureK),
        Math.log(retained.getTemperatureK()), Math.log(excluded.getTemperatureK()), Math.log(endpoint.pressureBara),
        Math.log(retained.getPressureBara()), Math.log(excluded.getPressureBara()));
    double endpointCompositionDistance = Math.min(compositionDistance(endpoint.exactBoundaryPoint, retained),
        compositionDistance(endpoint.exactBoundaryPoint, excluded));
    double excludedCompositionDistance = compositionDistance(endpoint.excludedSidePoint, excluded);
    if (transformedPtDistance > Math.max(maximumRelativeTemperatureDistance, maximumLogPressureDistance)
        || endpointCompositionDistance > 0.35 || excludedCompositionDistance > 0.35) {
      return TargetBranchMergeTrimResult.rejected("target-branch merge does not match this continuation branch");
    }
    List<EvidencePoint> retainedPoints = new ArrayList<EvidencePoint>();
    if (endpoint.requestedBranchEnd == End.END) {
      if (excludedIndex < 2) {
        return TargetBranchMergeTrimResult.rejected("target-branch merge would leave fewer than two points");
      }
      retainedPoints.addAll(orderedPoints.subList(0, excludedIndex));
    } else {
      if (excludedIndex + 2 >= orderedPoints.size()) {
        return TargetBranchMergeTrimResult.rejected("target-branch merge would leave fewer than two points");
      }
      retainedPoints.addAll(orderedPoints.subList(excludedIndex + 1, orderedPoints.size()));
    }
    replaceOrAddEndpoint(retainedPoints, endpoint.requestedBranchEnd, endpoint.exactBoundaryPoint);
    int trimmedPointCount = endpoint.requestedBranchEnd == End.END ? orderedPoints.size() - excludedIndex
        : excludedIndex + 1;
    double distance = Math.max(
        transformedPtDistance / Math.max(maximumRelativeTemperatureDistance, maximumLogPressureDistance),
        Math.max(endpointCompositionDistance, excludedCompositionDistance) / 0.35);
    return TargetBranchMergeTrimResult.accepted(retainedPoints, trimmedPointCount, excludedIndex, distance);
  }

  private static void replaceOrAddEndpoint(List<EvidencePoint> points, End end, EvidencePoint exact) {
    int endpointIndex = end == End.START ? 0 : points.size() - 1;
    EvidencePoint current = points.get(endpointIndex);
    if (relativeTemperatureDistance(current.getTemperatureK(), exact.getTemperatureK()) <= 1.0e-10
        && logPressureDistance(current.getPressureBara(), exact.getPressureBara()) <= 1.0e-10) {
      points.set(endpointIndex, exact);
    } else if (end == End.START) {
      points.add(0, exact);
    } else {
      points.add(exact);
    }
  }

  private static int closestEvidencePointIndex(List<EvidencePoint> points, EvidencePoint target) {
    int closestIndex = -1;
    double closestDistance = Double.POSITIVE_INFINITY;
    for (int index = 0; index < points.size(); index++) {
      EvidencePoint point = points.get(index);
      double ptDistance = Math.abs(Math.log(point.getTemperatureK() / target.getTemperatureK()))
          + Math.abs(Math.log(point.getPressureBara() / target.getPressureBara()));
      double distance = ptDistance + compositionDistance(point, target);
      if (distance < closestDistance) {
        closestDistance = distance;
        closestIndex = index;
      }
    }
    return closestIndex;
  }

  private static double pointToSegmentDistance(double targetX, double firstX, double secondX, double targetY,
      double firstY, double secondY) {
    double deltaX = secondX - firstX;
    double deltaY = secondY - firstY;
    double squaredLength = deltaX * deltaX + deltaY * deltaY;
    if (!(squaredLength > 0.0)) {
      return Math.hypot(targetX - firstX, targetY - firstY);
    }
    double parameter = ((targetX - firstX) * deltaX + (targetY - firstY) * deltaY) / squaredLength;
    parameter = Math.max(0.0, Math.min(1.0, parameter));
    double projectedX = firstX + parameter * deltaX;
    double projectedY = firstY + parameter * deltaY;
    return Math.hypot(targetX - projectedX, targetY - projectedY);
  }

  private static double compositionDistance(EvidencePoint first, EvidencePoint second) {
    return Math.max(compositionDistance(first.getPhaseZeroComposition(), second.getPhaseZeroComposition()),
        Math.max(compositionDistance(first.getPhaseOneComposition(), second.getPhaseOneComposition()),
            compositionDistance(first.getIncipientComposition(), second.getIncipientComposition())));
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

  private void addTrial(List<AttachmentTrial> trials, MutableBranch branch, End end, EndpointEvidence endpoint) {
    if (branch.isAttached(end)) {
      return;
    }
    EvidencePoint point = branch.endpointPoint(end);
    double relativeTemperature = relativeTemperatureDistance(point.getTemperatureK(), endpoint.temperatureK);
    double logPressure = logPressureDistance(point.getPressureBara(), endpoint.pressureBara);
    if (relativeTemperature <= maximumRelativeTemperatureDistance && logPressure <= maximumLogPressureDistance) {
      double distance = Math.max(relativeTemperature / maximumRelativeTemperatureDistance,
          logPressure / maximumLogPressureDistance);
      trials.add(new AttachmentTrial(branch, end, endpoint, distance));
    }
  }

  private static void validateEndpointDegrees(List<EndpointEvidence> endpoints, List<MutableBranch> branches,
      List<String> violations) {
    for (EndpointEvidence endpoint : endpoints) {
      boolean compatibleBranchExists = false;
      for (MutableBranch branch : branches) {
        compatibleBranchExists |= !branch.branch.isIsolated()
            && endpoint.compatibleFamilies.contains(branch.branch.getFamily());
      }
      if (!compatibleBranchExists) {
        violations.add("NO_COMPATIBLE_BRANCH_FOR_ENDPOINT:" + endpoint.identifier);
      }
      if (endpoint.getAttachmentCount() != endpoint.requiredDegree) {
        violations.add("ENDPOINT_DEGREE_MISMATCH:" + endpoint.identifier + ":" + endpoint.getAttachmentCount() + "/"
            + endpoint.requiredDegree);
      }
      Set<BoundaryFamily> attachedFamilies = new HashSet<BoundaryFamily>();
      for (EndpointUse attachment : endpoint.attachments) {
        if (!attachedFamilies.add(attachment.family)) {
          violations
              .add("DUPLICATE_ENDPOINT_FAMILY_ATTACHMENT:" + endpoint.identifier + ":" + attachment.family.name());
        }
      }
    }
  }

  private static double relativeTemperatureDistance(double first, double second) {
    return Math.abs(first - second) / Math.max(1.0, Math.min(first, second));
  }

  private static double logPressureDistance(double first, double second) {
    return Math.abs(Math.log(first / second));
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  private static BoundaryFamily family(CandidatePhase first, CandidatePhase second) {
    if (contains(first, second, CandidatePhase.GAS, CandidatePhase.OIL)) {
      return BoundaryFamily.GO_TO_GOW;
    }
    if (contains(first, second, CandidatePhase.GAS, CandidatePhase.AQUEOUS)) {
      return BoundaryFamily.GW_TO_GOW;
    }
    if (contains(first, second, CandidatePhase.OIL, CandidatePhase.AQUEOUS)) {
      return BoundaryFamily.OW_TO_GOW;
    }
    throw new IllegalArgumentException("unsupported retained phase pair " + first + "/" + second);
  }

  private static boolean contains(CandidatePhase first, CandidatePhase second, CandidatePhase expectedFirst,
      CandidatePhase expectedSecond) {
    return first == expectedFirst && second == expectedSecond || first == expectedSecond && second == expectedFirst;
  }

  /** Branch endpoint orientation after evidence points are ordered from backward to forward continuation. */
  public enum End {
    START, END
  }

  /** Independent numerical source required before a physical endpoint may enter the network. */
  public enum EvidenceSource {
    CRITICAL_ENDPOINT_SOLVER, RETAINED_PAIR_CRITICAL_ENDPOINT_SOLVER, THREE_PHASE_POINT_SOLVER, ENDPOINT_CLASSIFIER,
    TERMINATION_CLASSIFIER, DOMAIN_BOUND_CLASSIFIER, SECONDARY_STATIONARY_BRANCH_TRACKER
  }

  /** Immutable physical endpoint evidence with expected network degree and actual attachments. */
  public static final class EndpointEvidence {
    private final String identifier;
    private final SpecialPointType type;
    private final double temperatureK;
    private final double pressureBara;
    private final Set<BoundaryFamily> compatibleFamilies;
    private final int requiredDegree;
    private final EvidenceSource source;
    private final EvidencePoint exactBoundaryPoint;
    private final HydrocarbonWaterBoundaryTerminationClassifier.Result termination;
    private final EvidencePoint excludedSidePoint;
    private final End requestedBranchEnd;
    private final List<EndpointUse> attachments = new ArrayList<EndpointUse>();

    private EndpointEvidence(String identifier, SpecialPointType type, double temperatureK, double pressureBara,
        Set<BoundaryFamily> compatibleFamilies, int requiredDegree, EvidenceSource source,
        EvidencePoint exactBoundaryPoint, HydrocarbonWaterBoundaryTerminationClassifier.Result termination) {
      this(identifier, type, temperatureK, pressureBara, compatibleFamilies, requiredDegree, source, exactBoundaryPoint,
          termination, null, null);
    }

    private EndpointEvidence(String identifier, SpecialPointType type, double temperatureK, double pressureBara,
        Set<BoundaryFamily> compatibleFamilies, int requiredDegree, EvidenceSource source,
        EvidencePoint exactBoundaryPoint, HydrocarbonWaterBoundaryTerminationClassifier.Result termination,
        EvidencePoint excludedSidePoint, End requestedBranchEnd) {
      if (identifier == null || identifier.trim().isEmpty() || type == null || !positive(temperatureK)
          || !positive(pressureBara) || compatibleFamilies == null || compatibleFamilies.isEmpty() || requiredDegree < 1
          || source == null) {
        throw new IllegalArgumentException("invalid physical endpoint evidence");
      }
      this.identifier = identifier;
      this.type = type;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.compatibleFamilies = Collections.unmodifiableSet(EnumSet.copyOf(compatibleFamilies));
      this.requiredDegree = requiredDegree;
      this.source = source;
      this.exactBoundaryPoint = exactBoundaryPoint;
      this.termination = termination;
      this.excludedSidePoint = excludedSidePoint;
      this.requestedBranchEnd = requestedBranchEnd;
    }

    /** Creates endpoint evidence from the retained-two-phase critical endpoint equations. */
    public static EndpointEvidence criticalEndpoint(String identifier, BoundaryFamily family,
        HydrocarbonWaterCriticalEndpointSolver.Result endpoint) {
      return criticalEndpoint(identifier, family, endpoint, EvidenceSource.CRITICAL_ENDPOINT_SOLVER);
    }

    private static EndpointEvidence criticalEndpoint(String identifier, BoundaryFamily family,
        HydrocarbonWaterCriticalEndpointSolver.Result endpoint, EvidenceSource source) {
      if (family == null || endpoint == null || !endpoint.isPhysicalEndpoint()
          || endpoint.getRetainedPhaseZero() != family.getRetainedPhaseZero()
          || endpoint.getRetainedPhaseOne() != family.getRetainedPhaseOne()
          || endpoint.getIncipientPhase() != family.getIncipientPhase()) {
        throw new IllegalArgumentException("critical endpoint does not prove the requested boundary family");
      }
      return new EndpointEvidence(identifier, SpecialPointType.CRITICAL_END_POINT, endpoint.getTemperatureK(),
          endpoint.getPressureBara(), EnumSet.of(family), 1, source, endpoint.toEvidencePoint(), null);
    }

    /** Creates endpoint evidence where a retained pair coalesces at a distinct zero-TPD third phase. */
    public static EndpointEvidence retainedPairCriticalEndpoint(String identifier, BoundaryFamily family,
        HydrocarbonWaterRetainedPairCriticalEndpointSolver.Result endpoint) {
      if (family == null || endpoint == null || !endpoint.isPhysicalEndpoint()
          || endpoint.getRetainedPhaseZero() != family.getRetainedPhaseZero()
          || endpoint.getRetainedPhaseOne() != family.getRetainedPhaseOne()
          || endpoint.getIncipientPhase() != family.getIncipientPhase()) {
        throw new IllegalArgumentException("retained-pair endpoint does not prove the requested boundary family");
      }
      return new EndpointEvidence(identifier, SpecialPointType.CRITICAL_END_POINT, endpoint.getTemperatureK(),
          endpoint.getPressureBara(), EnumSet.of(family), 1, EvidenceSource.RETAINED_PAIR_CRITICAL_ENDPOINT_SOLVER,
          null, null);
    }

    /** Creates a two-branch junction from simultaneous incidence of two phases from one mother phase. */
    public static EndpointEvidence threePhasePoint(String identifier, ThreePhasePointSolver.Result endpoint) {
      return threePhasePoint(identifier, endpoint, null, EvidenceSource.THREE_PHASE_POINT_SOLVER);
    }

    /**
     * Creates a declared-domain endpoint only when the trace actually reaches one configured PT bound.
     *
     * @param identifier endpoint identifier
     * @param family boundary family traced to the domain edge
     * @param trace terminated pseudo-arclength trace
     * @param minimumTemperatureK lower solver temperature bound
     * @param maximumTemperatureK upper solver temperature bound
     * @param minimumPressureBara lower solver pressure bound
     * @param maximumPressureBara upper solver pressure bound
     * @param transformedDistanceTolerance accepted distance to a bound
     * @return independently classified domain endpoint evidence
     */
    public static EndpointEvidence domainExit(String identifier, BoundaryFamily family,
        TwoToThreePhasePseudoArcLengthTracer.Result trace, double minimumTemperatureK, double maximumTemperatureK,
        double minimumPressureBara, double maximumPressureBara, double transformedDistanceTolerance) {
      if (family == null || trace == null || trace.getPoints().size() < 3 || trace.getAcceptedCorrections().isEmpty()) {
        throw new IllegalArgumentException("a terminated boundary trace with corrected points is required");
      }
      HydrocarbonWaterBoundaryEventDetector detector = new HydrocarbonWaterBoundaryEventDetector();
      HydrocarbonWaterBoundaryEventDetector.Result events = detector.attachDomainExit(
          detector.detect(trace.getPoints()), trace, minimumTemperatureK, maximumTemperatureK, minimumPressureBara,
          maximumPressureBara, transformedDistanceTolerance);
      HydrocarbonWaterBoundaryEventDetector.Event event = events.getEvent(SpecialPointType.DOMAIN_EXIT);
      if (event == null) {
        throw new IllegalArgumentException("the trace did not reach a declared PT domain bound");
      }
      EvidencePoint exact = EvidencePoint
          .from(trace.getAcceptedCorrections().get(trace.getAcceptedCorrections().size() - 1));
      return new EndpointEvidence(identifier, SpecialPointType.DOMAIN_EXIT, event.getTemperatureK(),
          event.getPressureBara(), EnumSet.of(family), 1, EvidenceSource.DOMAIN_BOUND_CLASSIFIER, exact, null);
    }

    /** Creates declared-domain evidence from a strict fixed-pressure trace that terminated at one configured bound. */
    public static EndpointEvidence domainExit(String identifier, BoundaryFamily family,
        HydrocarbonWaterRegularPressureBoundaryTracer.Result trace, double minimumTemperatureK,
        double maximumTemperatureK, double minimumPressureBara, double maximumPressureBara,
        double transformedDistanceTolerance) {
      if (family == null || trace == null
          || trace.getTerminationReason() != HydrocarbonWaterRegularPressureBoundaryTracer.TerminationReason.DOMAIN_EXIT
          || trace.getPoints().size() < 3 || !positive(transformedDistanceTolerance)) {
        throw new IllegalArgumentException("a fixed-pressure trace terminated at a declared bound is required");
      }
      HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint last = trace.getPoints().get(trace.getPoints().size() - 1);
      double temperatureDistance = Math.min(relativeTemperatureDistance(last.getTemperatureK(), minimumTemperatureK),
          relativeTemperatureDistance(last.getTemperatureK(), maximumTemperatureK));
      double pressureDistance = Math.min(logPressureDistance(last.getPressureBara(), minimumPressureBara),
          logPressureDistance(last.getPressureBara(), maximumPressureBara));
      if (Math.min(temperatureDistance, pressureDistance) > transformedDistanceTolerance) {
        throw new IllegalArgumentException("the fixed-pressure trace did not reach a declared PT domain bound");
      }
      return new EndpointEvidence(identifier, SpecialPointType.DOMAIN_EXIT, last.getTemperatureK(),
          last.getPressureBara(), EnumSet.of(family), 1, EvidenceSource.DOMAIN_BOUND_CLASSIFIER, last.toEvidencePoint(),
          null);
    }

    private static EndpointEvidence threePhasePoint(String identifier, ThreePhasePointSolver.Result endpoint,
        EvidencePoint exactBoundaryPoint, EvidenceSource source) {
      if (endpoint == null || !endpoint.isConverged()) {
        throw new IllegalArgumentException("a converged distinct simultaneous-incidence point is required");
      }
      BoundaryFamily firstFamily = family(endpoint.getMotherPhase(), endpoint.getFirstIncipientPhase());
      BoundaryFamily secondFamily = family(endpoint.getMotherPhase(), endpoint.getSecondIncipientPhase());
      return new EndpointEvidence(identifier, SpecialPointType.THREE_PHASE_POINT, endpoint.getTemperatureK(),
          endpoint.getPressureBara(), EnumSet.of(firstFamily, secondFamily), 2, source, exactBoundaryPoint, null);
    }

    private static EndpointEvidence spinodal(String identifier, BoundaryFamily family,
        HydrocarbonWaterBoundaryTerminationClassifier.Result termination) {
      if (family == null || termination == null || !termination.isPhysicalEndpoint()
          || termination.getType() != HydrocarbonWaterBoundaryTerminationClassifier.Type.RETAINED_PHASE_SPINODAL) {
        throw new IllegalArgumentException("a refined retained-phase spinodal is required");
      }
      return new EndpointEvidence(identifier, SpecialPointType.RETAINED_PHASE_SPINODAL, termination.getTemperatureK(),
          termination.getPressureBara(), EnumSet.of(family), 1, EvidenceSource.TERMINATION_CLASSIFIER,
          EvidencePoint.from(termination.getBoundaryRoot()), termination);
    }

    /** Creates a one-branch endpoint where a proven secondary stationary identity merges with the target identity. */
    public static EndpointEvidence targetBranchMerge(String identifier, BoundaryFamily family,
        HydrocarbonWaterSecondaryStationaryBranchTracker.TargetBranchMergeBracket bracket, End requestedBranchEnd) {
      if (bracket == null || bracket.getMerged() == null || bracket.getDistinct() == null
          || !bracket.getMerged().isCollapsedTargetBranch() || !bracket.getDistinct().isTracked()) {
        throw new IllegalArgumentException("a collapsed-target/distinct stationary-branch bracket is required");
      }
      return targetBranchMerge(identifier, family, EvidencePoint.from(bracket.getMerged().getBoundaryRoot()),
          EvidencePoint.from(bracket.getDistinct().getBoundaryRoot()), requestedBranchEnd);
    }

    static EndpointEvidence targetBranchMerge(String identifier, BoundaryFamily family, EvidencePoint mergedPoint,
        EvidencePoint excludedSidePoint, End requestedBranchEnd) {
      if (family == null || mergedPoint == null || excludedSidePoint == null || requestedBranchEnd == null) {
        throw new IllegalArgumentException("target-branch merge evidence and retained branch orientation are required");
      }
      return new EndpointEvidence(identifier, SpecialPointType.TARGET_BRANCH_MERGE, mergedPoint.getTemperatureK(),
          mergedPoint.getPressureBara(), EnumSet.of(family), 1, EvidenceSource.SECONDARY_STATIONARY_BRANCH_TRACKER,
          mergedPoint, null, excludedSidePoint, requestedBranchEnd);
    }

    private void recordAttachment(String branchIdentifier, BoundaryFamily family, End end) {
      attachments.add(new EndpointUse(branchIdentifier, family, end));
    }

    private EndpointEvidence copyWithoutAttachments() {
      return new EndpointEvidence(identifier, type, temperatureK, pressureBara, compatibleFamilies, requiredDegree,
          source, exactBoundaryPoint, termination, excludedSidePoint, requestedBranchEnd);
    }

    private boolean isAttachedTo(String branchIdentifier, End end) {
      for (EndpointUse attachment : attachments) {
        if (attachment.branchIdentifier.equals(branchIdentifier) && attachment.end == end) {
          return true;
        }
      }
      return false;
    }

    public String getIdentifier() {
      return identifier;
    }

    public SpecialPointType getType() {
      return type;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public Set<BoundaryFamily> getCompatibleFamilies() {
      return compatibleFamilies;
    }

    public int getRequiredDegree() {
      return requiredDegree;
    }

    public int getAttachmentCount() {
      return attachments.size();
    }

    public EvidenceSource getSource() {
      return source;
    }

    public List<EndpointUse> getAttachments() {
      return Collections.unmodifiableList(new ArrayList<EndpointUse>(attachments));
    }

    public End getRequestedBranchEnd() {
      return requestedBranchEnd;
    }
  }

  /** One physical endpoint use by a named branch end. */
  public static final class EndpointUse {
    private final String branchIdentifier;
    private final BoundaryFamily family;
    private final End end;

    private EndpointUse(String branchIdentifier, BoundaryFamily family, End end) {
      this.branchIdentifier = branchIdentifier;
      this.family = family;
      this.end = end;
    }

    public String getBranchIdentifier() {
      return branchIdentifier;
    }

    public BoundaryFamily getFamily() {
      return family;
    }

    public End getEnd() {
      return end;
    }
  }

  /** Immutable branch-to-endpoint attachment recorded in the assembled network. */
  public static final class EndpointAttachment {
    private final String endpointIdentifier;
    private final SpecialPointType type;
    private final double temperatureK;
    private final double pressureBara;
    private final EvidenceSource source;
    private final double transformedDistance;

    private EndpointAttachment(EndpointEvidence endpoint, double transformedDistance) {
      this.endpointIdentifier = endpoint.identifier;
      this.type = endpoint.type;
      this.temperatureK = endpoint.temperatureK;
      this.pressureBara = endpoint.pressureBara;
      this.source = endpoint.source;
      this.transformedDistance = transformedDistance;
    }

    public String getEndpointIdentifier() {
      return endpointIdentifier;
    }

    public SpecialPointType getType() {
      return type;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public EvidenceSource getSource() {
      return source;
    }

    public double getTransformedDistance() {
      return transformedDistance;
    }
  }

  /** One preserved branch after physical endpoint attachment. */
  public static final class NetworkBranch {
    private final AssembledBranch sourceBranch;
    private final List<EvidencePoint> points;
    private final EndpointAttachment startAttachment;
    private final EndpointAttachment endAttachment;
    private final int trimmedPointCount;

    private NetworkBranch(AssembledBranch sourceBranch, List<EvidencePoint> points, EndpointAttachment startAttachment,
        EndpointAttachment endAttachment, int trimmedPointCount) {
      this.sourceBranch = sourceBranch;
      this.points = Collections.unmodifiableList(new ArrayList<EvidencePoint>(points));
      this.startAttachment = startAttachment;
      this.endAttachment = endAttachment;
      this.trimmedPointCount = trimmedPointCount;
    }

    public String getIdentifier() {
      return sourceBranch.getIdentifier();
    }

    public BoundaryFamily getFamily() {
      return sourceBranch.getFamily();
    }

    public boolean isIsolated() {
      return sourceBranch.isIsolated();
    }

    /** @return strict globally stable points inserted between sparse discovery anchors */
    public int getInternalBridgePointCount() {
      return sourceBranch.getInternalBridgePointCount();
    }

    /** @return strict fixed-pressure seed points generated from a single discovered anchor */
    public int getGeneratedSeedPointCount() {
      return sourceBranch.getGeneratedSeedPointCount();
    }

    public List<EvidencePoint> getPoints() {
      return points;
    }

    public EndpointAttachment getStartAttachment() {
      return startAttachment;
    }

    public EndpointAttachment getEndAttachment() {
      return endAttachment;
    }

    /** @return converged but excluded points removed beyond a proven physical branch-merge endpoint */
    public int getTrimmedPointCount() {
      return trimmedPointCount;
    }

    public TwoToThreePhasePseudoArcLengthTracer.Result getBackwardTrace() {
      return sourceBranch.getBackwardTrace();
    }

    public TwoToThreePhasePseudoArcLengthTracer.Result getForwardTrace() {
      return sourceBranch.getForwardTrace();
    }

    public HydrocarbonWaterRegularPressureBoundaryTracer.Result getRegularBackwardTrace() {
      return sourceBranch.getRegularBackwardTrace();
    }

    public HydrocarbonWaterRegularPressureBoundaryTracer.Result getRegularForwardTrace() {
      return sourceBranch.getRegularForwardTrace();
    }

    public TwoToThreePhasePseudoArcLengthTracer.Result getPostRegularBackwardTrace() {
      return sourceBranch.getPostRegularBackwardTrace();
    }

    public TwoToThreePhasePseudoArcLengthTracer.Result getPostRegularForwardTrace() {
      return sourceBranch.getPostRegularForwardTrace();
    }

    public HydrocarbonWaterBoundaryTerminationClassifier.Result getBackwardTermination() {
      return sourceBranch.getBackwardTermination();
    }

    public HydrocarbonWaterBoundaryTerminationClassifier.Result getForwardTermination() {
      return sourceBranch.getForwardTermination();
    }

    /** Converts this network branch to the existing per-branch engineering gate contract. */
    public TwoToThreePhaseBoundaryQualityGate.Branch toQualityGateBranch() {
      if (isIsolated() || points.size() < 2) {
        throw new IllegalStateException("an isolated network branch cannot enter the boundary quality gate");
      }
      SpecialPointType startType = startAttachment == null ? null : startAttachment.type;
      SpecialPointType endType = endAttachment == null ? null : endAttachment.type;
      HydrocarbonWaterBoundaryTerminationClassifier.Result startTermination = startAttachment != null
          && startAttachment.type == SpecialPointType.RETAINED_PHASE_SPINODAL ? sourceBranch.getBackwardTermination()
              : null;
      HydrocarbonWaterBoundaryTerminationClassifier.Result endTermination = endAttachment != null
          && endAttachment.type == SpecialPointType.RETAINED_PHASE_SPINODAL ? sourceBranch.getForwardTermination()
              : null;
      return new TwoToThreePhaseBoundaryQualityGate.Branch(sourceBranch.getIdentifier(),
          sourceBranch.getFamily().getDefinition(), sourceBranch.getFamily().getRetainedPhaseZero(),
          sourceBranch.getFamily().getRetainedPhaseOne(), sourceBranch.getFamily().getIncipientPhase(), points,
          startType, endType, startTermination, endTermination);
    }
  }

  /** Immutable network-level result. */
  public static final class Result {
    private final List<NetworkBranch> branches;
    private final List<EndpointEvidence> endpoints;
    private final int originalEndpointCandidateCount;
    private final TwoToThreePhaseBoundaryQualityGate.EnvelopeReport qualityReport;
    private final List<String> violations;
    private final boolean engineeringEligible;

    private Result(List<NetworkBranch> branches, List<EndpointEvidence> endpoints, int originalEndpointCandidateCount,
        TwoToThreePhaseBoundaryQualityGate.EnvelopeReport qualityReport, List<String> violations,
        boolean engineeringEligible) {
      this.branches = Collections.unmodifiableList(new ArrayList<NetworkBranch>(branches));
      this.endpoints = Collections.unmodifiableList(new ArrayList<EndpointEvidence>(endpoints));
      this.originalEndpointCandidateCount = originalEndpointCandidateCount;
      this.qualityReport = qualityReport;
      this.violations = Collections.unmodifiableList(new ArrayList<String>(violations));
      this.engineeringEligible = engineeringEligible;
    }

    public List<NetworkBranch> getBranches() {
      return branches;
    }

    public List<EndpointEvidence> getEndpoints() {
      return endpoints;
    }

    public int getOriginalEndpointCandidateCount() {
      return originalEndpointCandidateCount;
    }

    public TwoToThreePhaseBoundaryQualityGate.EnvelopeReport getQualityReport() {
      return qualityReport;
    }

    public List<String> getViolations() {
      return violations;
    }

    public boolean isEngineeringEligible() {
      return engineeringEligible;
    }
  }

  private static final class MutableBranch {
    private final AssembledBranch branch;
    private final List<EvidencePoint> points;
    private EndpointAttachment startAttachment;
    private EndpointAttachment endAttachment;
    private int trimmedPointCount;

    private MutableBranch(AssembledBranch branch) {
      this.branch = branch;
      this.points = new ArrayList<EvidencePoint>(branch.getEvidencePoints());
    }

    private boolean isAttached(End end) {
      return end == End.START ? startAttachment != null : endAttachment != null;
    }

    private EvidencePoint endpointPoint(End end) {
      if (points.isEmpty()) {
        throw new IllegalStateException("a continuation branch must carry endpoint evidence");
      }
      return end == End.START ? points.get(0) : points.get(points.size() - 1);
    }

    private void attach(End end, EndpointEvidence endpoint, double distance) {
      EndpointAttachment attachment = new EndpointAttachment(endpoint, distance);
      if (end == End.START) {
        startAttachment = attachment;
        addExactEndpointPoint(end, endpoint.exactBoundaryPoint);
      } else {
        endAttachment = attachment;
        addExactEndpointPoint(end, endpoint.exactBoundaryPoint);
      }
    }

    private void applyTargetBranchMergeTrim(TargetBranchMergeTrimResult trim, EndpointEvidence endpoint,
        double distance) {
      points.clear();
      points.addAll(trim.getRetainedPoints());
      trimmedPointCount += trim.getTrimmedPointCount();
      attach(endpoint.requestedBranchEnd, endpoint, distance);
    }

    private void addExactEndpointPoint(End end, EvidencePoint exact) {
      if (exact == null || points.isEmpty()) {
        return;
      }
      EvidencePoint current = endpointPoint(end);
      if (relativeTemperatureDistance(current.getTemperatureK(), exact.getTemperatureK()) <= 1.0e-10
          && logPressureDistance(current.getPressureBara(), exact.getPressureBara()) <= 1.0e-10) {
        if (end == End.START) {
          points.set(0, exact);
        } else {
          points.set(points.size() - 1, exact);
        }
      } else if (end == End.START) {
        points.add(0, exact);
      } else {
        points.add(exact);
      }
    }

    private NetworkBranch freeze() {
      return new NetworkBranch(branch, points, startAttachment, endAttachment, trimmedPointCount);
    }
  }

  private static final class TargetBranchMergeTrial {
    private final MutableBranch branch;
    private final EndpointEvidence endpoint;
    private final TargetBranchMergeTrimResult trim;
    private final double distance;

    private TargetBranchMergeTrial(MutableBranch branch, EndpointEvidence endpoint, TargetBranchMergeTrimResult trim,
        double distance) {
      this.branch = branch;
      this.endpoint = endpoint;
      this.trim = trim;
      this.distance = distance;
    }

    private double getDistance() {
      return distance;
    }
  }

  /** Immutable result of identity-aware target-branch merge truncation. */
  public static final class TargetBranchMergeTrimResult {
    private final List<EvidencePoint> retainedPoints;
    private final int trimmedPointCount;
    private final int excludedPointIndex;
    private final double transformedDistance;
    private final boolean accepted;
    private final String failureMessage;

    private TargetBranchMergeTrimResult(List<EvidencePoint> retainedPoints, int trimmedPointCount,
        int excludedPointIndex, double transformedDistance, boolean accepted, String failureMessage) {
      this.retainedPoints = Collections.unmodifiableList(new ArrayList<EvidencePoint>(retainedPoints));
      this.trimmedPointCount = trimmedPointCount;
      this.excludedPointIndex = excludedPointIndex;
      this.transformedDistance = transformedDistance;
      this.accepted = accepted;
      this.failureMessage = failureMessage;
    }

    private static TargetBranchMergeTrimResult accepted(List<EvidencePoint> retainedPoints, int trimmedPointCount,
        int excludedPointIndex, double transformedDistance) {
      return new TargetBranchMergeTrimResult(retainedPoints, trimmedPointCount, excludedPointIndex, transformedDistance,
          true, null);
    }

    private static TargetBranchMergeTrimResult rejected(String failureMessage) {
      return new TargetBranchMergeTrimResult(Collections.<EvidencePoint>emptyList(), 0, -1, Double.POSITIVE_INFINITY,
          false, failureMessage);
    }

    public List<EvidencePoint> getRetainedPoints() {
      return retainedPoints;
    }

    public int getTrimmedPointCount() {
      return trimmedPointCount;
    }

    public int getExcludedPointIndex() {
      return excludedPointIndex;
    }

    public double getTransformedDistance() {
      return transformedDistance;
    }

    public boolean isAccepted() {
      return accepted;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }

  private static final class AttachmentTrial {
    private final MutableBranch branch;
    private final End end;
    private final EndpointEvidence endpoint;
    private final double distance;

    private AttachmentTrial(MutableBranch branch, End end, EndpointEvidence endpoint, double distance) {
      this.branch = branch;
      this.end = end;
      this.endpoint = endpoint;
      this.distance = distance;
    }

    private double getDistance() {
      return distance;
    }
  }
}
