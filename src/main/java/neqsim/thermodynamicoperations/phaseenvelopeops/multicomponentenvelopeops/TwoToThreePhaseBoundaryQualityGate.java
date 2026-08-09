package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.EndpointCandidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.BoundaryDefinition;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Phase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Global topology and numerical quality gate for two-to-three-phase PT boundary branches. */
public final class TwoToThreePhaseBoundaryQualityGate {
  private double residualTolerance = 1.0e-7;
  private double compositionSumTolerance = 1.0e-8;
  private double minimumDistinctPhaseDistance = 1.0e-5;
  private double maximumCompositionJump = 0.35;
  private double maximumLogPressureJump = 1.0;
  private double maximumRelativeTemperatureJump = 0.25;

  /** Sets thermodynamic, composition, phase-distinctness, and adjacent-point tolerances. */
  public TwoToThreePhaseBoundaryQualityGate setTolerances(double residualTolerance, double compositionSumTolerance,
      double minimumDistinctPhaseDistance, double maximumCompositionJump, double maximumLogPressureJump,
      double maximumRelativeTemperatureJump) {
    if (!Double.isFinite(residualTolerance) || residualTolerance <= 0.0 || !Double.isFinite(compositionSumTolerance)
        || compositionSumTolerance <= 0.0 || !Double.isFinite(minimumDistinctPhaseDistance)
        || minimumDistinctPhaseDistance <= 0.0 || !Double.isFinite(maximumCompositionJump)
        || maximumCompositionJump <= 0.0 || !Double.isFinite(maximumLogPressureJump) || maximumLogPressureJump <= 0.0
        || !Double.isFinite(maximumRelativeTemperatureJump) || maximumRelativeTemperatureJump <= 0.0) {
      throw new IllegalArgumentException("invalid two-to-three-phase quality-gate tolerances");
    }
    this.residualTolerance = residualTolerance;
    this.compositionSumTolerance = compositionSumTolerance;
    this.minimumDistinctPhaseDistance = minimumDistinctPhaseDistance;
    this.maximumCompositionJump = maximumCompositionJump;
    this.maximumLogPressureJump = maximumLogPressureJump;
    this.maximumRelativeTemperatureJump = maximumRelativeTemperatureJump;
    return this;
  }

  /**
   * Validates every supplied branch without dropping failed or isolated branches.
   *
   * @param branches named physical branches
   * @return one report per input branch plus an envelope-level decision
   */
  public EnvelopeReport validate(List<Branch> branches) {
    return validate(branches, Collections.<EndpointCandidate>emptyList());
  }

  /**
   * Validates classified branches and rejects every corrected special point that remains unattached to the topology.
   *
   * <p>
   * Endpoint candidates are never discarded to obtain a green envelope. Even an independently resolved special point
   * must be connected to the appropriate branch before the global topology can enter engineering use.
   * </p>
   */
  public EnvelopeReport validate(List<Branch> branches, List<EndpointCandidate> endpointCandidates) {
    if (branches == null || endpointCandidates == null || branches.isEmpty() && endpointCandidates.isEmpty()) {
      throw new IllegalArgumentException("at least one boundary branch is required");
    }
    List<BranchReport> reports = new ArrayList<BranchReport>();
    Set<String> identifiers = new HashSet<String>();
    List<String> envelopeViolations = new ArrayList<String>();
    for (Branch branch : branches) {
      if (branch == null) {
        envelopeViolations.add("NULL_BRANCH");
        continue;
      }
      if (!identifiers.add(branch.identifier)) {
        envelopeViolations.add("DUPLICATE_BRANCH_ID:" + branch.identifier);
      }
      reports.add(validateBranch(branch));
    }
    for (int candidateIndex = 0; candidateIndex < endpointCandidates.size(); candidateIndex++) {
      EndpointCandidate candidate = endpointCandidates.get(candidateIndex);
      if (candidate == null) {
        envelopeViolations.add("NULL_ENDPOINT_CANDIDATE:" + candidateIndex);
      } else {
        envelopeViolations.add("UNATTACHED_ENDPOINT_CANDIDATE:" + candidate.getFamily().name() + ":"
            + candidate.getClassification().getClassification().name());
      }
    }
    if (reports.size() != branches.size()) {
      envelopeViolations.add("BRANCH_PRESERVATION_FAILURE");
    }
    boolean accepted = envelopeViolations.isEmpty();
    for (BranchReport report : reports) {
      accepted &= report.isAccepted();
    }
    return new EnvelopeReport(reports, endpointCandidates.size(), envelopeViolations, accepted);
  }

  private BranchReport validateBranch(Branch branch) {
    List<String> violations = new ArrayList<String>();
    validateTopology(branch, violations);
    if (branch.points.size() < 2) {
      violations.add("TOO_FEW_POINTS");
    }
    double maximumResidual = 0.0;
    double maximumObservedCompositionJump = 0.0;
    double minimumObservedPhaseDistance = Double.POSITIVE_INFINITY;
    for (int pointIndex = 0; pointIndex < branch.points.size(); pointIndex++) {
      EvidencePoint point = branch.points.get(pointIndex);
      if (point == null) {
        violations.add("NULL_POINT:" + pointIndex);
        continue;
      }
      if (!Double.isFinite(point.temperatureK) || point.temperatureK < 50.0 || !Double.isFinite(point.pressureBara)
          || point.pressureBara <= 0.0 || !Double.isFinite(point.beta) || point.beta <= 0.0 || point.beta >= 1.0) {
        violations.add("NONPHYSICAL_PT_OR_BETA:" + pointIndex);
      }
      double pointMaximumResidual = maximumAbs(point.retainedFlashResidual, point.tangentPlaneDistance,
          point.stationarityResidual);
      maximumResidual = Math.max(maximumResidual, pointMaximumResidual);
      if (!Double.isFinite(pointMaximumResidual) || pointMaximumResidual > residualTolerance) {
        violations.add("RESIDUAL_GATE:" + pointIndex);
      }
      validateComposition(point.phaseZeroComposition, pointIndex, "PHASE_ZERO", violations);
      validateComposition(point.phaseOneComposition, pointIndex, "PHASE_ONE", violations);
      validateComposition(point.incipientComposition, pointIndex, "INCIPIENT", violations);
      if (sameLength(point.phaseZeroComposition, point.phaseOneComposition, point.incipientComposition)) {
        double firstSecond = compositionDistance(point.phaseZeroComposition, point.phaseOneComposition);
        double firstThird = compositionDistance(point.phaseZeroComposition, point.incipientComposition);
        double secondThird = compositionDistance(point.phaseOneComposition, point.incipientComposition);
        double minimumDistance = Math.min(firstSecond, Math.min(firstThird, secondThird));
        minimumObservedPhaseDistance = Math.min(minimumObservedPhaseDistance, minimumDistance);
        if (minimumDistance <= minimumDistinctPhaseDistance && !isCriticalCoalescenceEndpoint(branch, pointIndex)) {
          violations.add("TRIVIAL_OR_DUPLICATE_PHASE:" + pointIndex);
        }
      } else {
        violations.add("COMPONENT_COUNT_MISMATCH:" + pointIndex);
      }
      if (pointIndex > 0) {
        EvidencePoint previous = branch.points.get(pointIndex - 1);
        if (previous != null && sameLength(previous.phaseZeroComposition, point.phaseZeroComposition)
            && sameLength(previous.phaseOneComposition, point.phaseOneComposition)
            && sameLength(previous.incipientComposition, point.incipientComposition)) {
          double compositionJump = Math.max(
              compositionDistance(previous.phaseZeroComposition, point.phaseZeroComposition),
              Math.max(compositionDistance(previous.phaseOneComposition, point.phaseOneComposition),
                  compositionDistance(previous.incipientComposition, point.incipientComposition)));
          maximumObservedCompositionJump = Math.max(maximumObservedCompositionJump, compositionJump);
          double logPressureJump = Math.abs(Math.log(point.pressureBara / previous.pressureBara));
          double relativeTemperatureJump = Math.abs(point.temperatureK - previous.temperatureK)
              / Math.max(1.0, Math.min(point.temperatureK, previous.temperatureK));
          if (compositionJump > maximumCompositionJump || logPressureJump > maximumLogPressureJump
              || relativeTemperatureJump > maximumRelativeTemperatureJump) {
            violations.add("ADJACENT_CONTINUITY_GATE:" + (pointIndex - 1) + "->" + pointIndex);
          }
        }
      }
    }
    validateEndpoints(branch, minimumObservedPhaseDistance, violations);
    IntersectionSummary intersections = selfIntersections(branch.points);
    if (intersections.identityIntersectionCount > 0) {
      violations.add("SELF_INTERSECTION");
    }
    return new BranchReport(branch.identifier, branch.points.size(), maximumResidual, maximumObservedCompositionJump,
        minimumObservedPhaseDistance, intersections.projectedIntersectionCount, intersections.identityIntersectionCount,
        intersections.minimumStateDistance, violations, violations.isEmpty());
  }

  private static void validateTopology(Branch branch, List<String> violations) {
    if (branch.definition == null || branch.retainedPhaseZero == null || branch.retainedPhaseOne == null
        || branch.incipientPhase == null) {
      violations.add("MISSING_TOPOLOGY_DEFINITION");
      return;
    }
    Set<Phase> retained = new HashSet<Phase>();
    retained.add(toTopologyPhase(branch.retainedPhaseZero));
    retained.add(toTopologyPhase(branch.retainedPhaseOne));
    if (!branch.definition.getLowerPhaseRegion().getPhases().equals(retained)) {
      violations.add("LOWER_REGION_PHASE_MISMATCH");
    }
    if (branch.definition.getIncipientPhase() != toTopologyPhase(branch.incipientPhase)) {
      violations.add("INCIPIENT_PHASE_MISMATCH");
    }
    if (branch.definition.getHigherPhaseRegion().getPhases().size() != 3) {
      violations.add("HIGHER_REGION_IS_NOT_THREE_PHASE");
    }
  }

  private void validateEndpoints(Branch branch, double minimumPhaseDistance, List<String> violations) {
    boolean closed = isClosed(branch.points);
    if (!closed && (branch.startPointType == null || branch.endPointType == null)) {
      violations.add("OPEN_ENDPOINT_UNCLASSIFIED");
      return;
    }
    if (closed && (branch.startPointType != null || branch.endPointType != null)) {
      violations.add("CLOSED_BRANCH_HAS_ENDPOINT_LABEL");
    }
    validateEndpointType(branch.startPointType, "START", branch.startTermination, violations);
    validateEndpointType(branch.endPointType, "END", branch.endTermination, violations);
    if (!branch.points.isEmpty()) {
      validateTerminationLocation(branch.startPointType, branch.startTermination, branch.points.get(0), "START",
          violations);
      validateTerminationLocation(branch.endPointType, branch.endTermination,
          branch.points.get(branch.points.size() - 1), "END", violations);
    }
    if (!closed && (requiresPhaseCoalescence(branch.startPointType) || requiresPhaseCoalescence(branch.endPointType))
        && minimumPhaseDistance > 10.0 * minimumDistinctPhaseDistance) {
      violations.add("CRITICAL_ENDPOINT_WITHOUT_PHASE_COALESCENCE");
    }
  }

  private static void validateEndpointType(SpecialPointType type, String endpoint,
      HydrocarbonWaterBoundaryTerminationClassifier.Result termination, List<String> violations) {
    if (type == SpecialPointType.STABILITY_STATIONARY_FOLD) {
      violations.add("UNRESOLVED_STABILITY_FOLD_ENDPOINT:" + endpoint);
    } else if (type == SpecialPointType.RETAINED_PHASE_SPINODAL) {
      if (termination == null || !termination.isPhysicalEndpoint()
          || termination.getType() != HydrocarbonWaterBoundaryTerminationClassifier.Type.RETAINED_PHASE_SPINODAL) {
        violations.add("UNRESOLVED_RETAINED_PHASE_SPINODAL_ENDPOINT:" + endpoint);
      }
    } else if (type == SpecialPointType.ALTERNATE_STATIONARY_BRANCH_ONSET) {
      violations.add("UNRESOLVED_ALTERNATE_STATIONARY_BRANCH_ONSET_ENDPOINT:" + endpoint);
    } else if (type == SpecialPointType.PHASE_COALESCENCE_CANDIDATE) {
      violations.add("UNRESOLVED_PHASE_COALESCENCE_ENDPOINT:" + endpoint);
    } else if (type == SpecialPointType.CRICONDENBAR || type == SpecialPointType.CRICONDENTHERM) {
      violations.add("NONTERMINATING_SPECIAL_POINT_AS_ENDPOINT:" + endpoint);
    }
    if (termination != null && termination.isPhysicalEndpoint() && type != SpecialPointType.RETAINED_PHASE_SPINODAL) {
      violations.add("PHYSICAL_TERMINATION_TYPE_MISMATCH:" + endpoint);
    }
  }

  private static void validateTerminationLocation(SpecialPointType type,
      HydrocarbonWaterBoundaryTerminationClassifier.Result termination, EvidencePoint endpointPoint, String endpoint,
      List<String> violations) {
    if (type != SpecialPointType.RETAINED_PHASE_SPINODAL || termination == null || !termination.isPhysicalEndpoint()) {
      return;
    }
    double relativeTemperature = Math.abs(endpointPoint.temperatureK - termination.getTemperatureK())
        / Math.max(1.0, termination.getTemperatureK());
    double logPressure = Math.abs(Math.log(endpointPoint.pressureBara / termination.getPressureBara()));
    if (relativeTemperature > 1.0e-7 || logPressure > 1.0e-7) {
      violations.add("PHYSICAL_TERMINATION_POINT_MISMATCH:" + endpoint);
    }
  }

  private static boolean requiresPhaseCoalescence(SpecialPointType type) {
    return type == SpecialPointType.CRITICAL_POINT || type == SpecialPointType.CRITICAL_END_POINT;
  }

  private static boolean isCriticalCoalescenceEndpoint(Branch branch, int pointIndex) {
    return pointIndex == 0 && requiresPhaseCoalescence(branch.startPointType)
        || pointIndex == branch.points.size() - 1 && requiresPhaseCoalescence(branch.endPointType);
  }

  private static boolean isClosed(List<EvidencePoint> points) {
    if (points.size() < 3 || points.get(0) == null || points.get(points.size() - 1) == null) {
      return false;
    }
    EvidencePoint first = points.get(0);
    EvidencePoint last = points.get(points.size() - 1);
    return Math.abs(first.temperatureK - last.temperatureK) <= 1.0e-7
        && Math.abs(first.pressureBara - last.pressureBara) <= 1.0e-7;
  }

  private static IntersectionSummary selfIntersections(List<EvidencePoint> points) {
    int projectedIntersectionCount = 0;
    int identityIntersectionCount = 0;
    double minimumStateDistance = Double.POSITIVE_INFINITY;
    for (int firstSegment = 0; firstSegment + 1 < points.size(); firstSegment++) {
      EvidencePoint a = points.get(firstSegment);
      EvidencePoint b = points.get(firstSegment + 1);
      if (a == null || b == null) {
        continue;
      }
      for (int secondSegment = firstSegment + 2; secondSegment + 1 < points.size(); secondSegment++) {
        if (firstSegment == 0 && secondSegment + 1 == points.size() - 1 && isClosed(points)) {
          continue;
        }
        EvidencePoint c = points.get(secondSegment);
        EvidencePoint d = points.get(secondSegment + 1);
        if (c != null && d != null && segmentsProperlyIntersect(a, b, c, d)) {
          projectedIntersectionCount++;
          double stateDistance = interpolatedIntersectionStateDistance(a, b, c, d);
          minimumStateDistance = Math.min(minimumStateDistance, stateDistance);
          if (stateDistance <= 1.0e-7) {
            identityIntersectionCount++;
          }
        }
      }
    }
    return new IntersectionSummary(projectedIntersectionCount, identityIntersectionCount, minimumStateDistance);
  }

  private static double interpolatedIntersectionStateDistance(EvidencePoint a, EvidencePoint b, EvidencePoint c,
      EvidencePoint d) {
    double firstTemperatureStep = b.temperatureK - a.temperatureK;
    double firstPressureStep = b.pressureBara - a.pressureBara;
    double secondTemperatureStep = d.temperatureK - c.temperatureK;
    double secondPressureStep = d.pressureBara - c.pressureBara;
    double denominator = cross(firstTemperatureStep, firstPressureStep, secondTemperatureStep, secondPressureStep);
    if (Math.abs(denominator) <= 1.0e-20) {
      return Double.POSITIVE_INFINITY;
    }
    double offsetTemperature = c.temperatureK - a.temperatureK;
    double offsetPressure = c.pressureBara - a.pressureBara;
    double firstFraction = cross(offsetTemperature, offsetPressure, secondTemperatureStep, secondPressureStep)
        / denominator;
    double secondFraction = cross(offsetTemperature, offsetPressure, firstTemperatureStep, firstPressureStep)
        / denominator;
    double compositionDistance = Math.max(
        interpolatedCompositionDistance(a.phaseZeroComposition, b.phaseZeroComposition, firstFraction,
            c.phaseZeroComposition, d.phaseZeroComposition, secondFraction),
        Math.max(
            interpolatedCompositionDistance(a.phaseOneComposition, b.phaseOneComposition, firstFraction,
                c.phaseOneComposition, d.phaseOneComposition, secondFraction),
            interpolatedCompositionDistance(a.incipientComposition, b.incipientComposition, firstFraction,
                c.incipientComposition, d.incipientComposition, secondFraction)));
    double betaDistance = Math
        .abs(interpolate(a.beta, b.beta, firstFraction) - interpolate(c.beta, d.beta, secondFraction));
    return Math.max(compositionDistance, betaDistance);
  }

  private static double interpolatedCompositionDistance(double[] firstStart, double[] firstEnd, double firstFraction,
      double[] secondStart, double[] secondEnd, double secondFraction) {
    if (!sameLength(firstStart, firstEnd, secondStart, secondEnd)) {
      return Double.POSITIVE_INFINITY;
    }
    double distance = 0.0;
    for (int index = 0; index < firstStart.length; index++) {
      distance += Math.abs(interpolate(firstStart[index], firstEnd[index], firstFraction)
          - interpolate(secondStart[index], secondEnd[index], secondFraction));
    }
    return distance;
  }

  private static double interpolate(double start, double end, double fraction) {
    return start + fraction * (end - start);
  }

  private static double cross(double firstX, double firstY, double secondX, double secondY) {
    return firstX * secondY - firstY * secondX;
  }

  private static boolean segmentsProperlyIntersect(EvidencePoint a, EvidencePoint b, EvidencePoint c, EvidencePoint d) {
    double first = orientation(a, b, c);
    double second = orientation(a, b, d);
    double third = orientation(c, d, a);
    double fourth = orientation(c, d, b);
    double tolerance = 1.0e-12;
    return first * second < -tolerance && third * fourth < -tolerance;
  }

  private static double orientation(EvidencePoint a, EvidencePoint b, EvidencePoint c) {
    return (b.temperatureK - a.temperatureK) * (c.pressureBara - a.pressureBara)
        - (b.pressureBara - a.pressureBara) * (c.temperatureK - a.temperatureK);
  }

  private static final class IntersectionSummary {
    private final int projectedIntersectionCount;
    private final int identityIntersectionCount;
    private final double minimumStateDistance;

    private IntersectionSummary(int projectedIntersectionCount, int identityIntersectionCount,
        double minimumStateDistance) {
      this.projectedIntersectionCount = projectedIntersectionCount;
      this.identityIntersectionCount = identityIntersectionCount;
      this.minimumStateDistance = minimumStateDistance;
    }
  }

  private void validateComposition(double[] composition, int pointIndex, String label, List<String> violations) {
    if (composition == null || composition.length == 0) {
      violations.add("EMPTY_" + label + ":" + pointIndex);
      return;
    }
    double sum = 0.0;
    for (double value : composition) {
      if (!Double.isFinite(value) || value < 0.0 || value > 1.0 + compositionSumTolerance) {
        violations.add("INVALID_" + label + ":" + pointIndex);
        return;
      }
      sum += value;
    }
    if (Math.abs(sum - 1.0) > compositionSumTolerance) {
      violations.add("UNNORMALIZED_" + label + ":" + pointIndex);
    }
  }

  private static double maximumAbs(double... values) {
    double maximum = 0.0;
    for (double value : values) {
      if (!Double.isFinite(value)) {
        return Double.POSITIVE_INFINITY;
      }
      maximum = Math.max(maximum, Math.abs(value));
    }
    return maximum;
  }

  private static boolean sameLength(double[]... compositions) {
    if (compositions.length == 0 || compositions[0] == null) {
      return false;
    }
    int length = compositions[0].length;
    for (double[] composition : compositions) {
      if (composition == null || composition.length != length) {
        return false;
      }
    }
    return true;
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

  private static Phase toTopologyPhase(CandidatePhase phase) {
    switch (phase) {
    case GAS:
      return Phase.GAS;
    case OIL:
      return Phase.OIL;
    case AQUEOUS:
      return Phase.AQUEOUS;
    default:
      throw new IllegalArgumentException("unsupported phase family " + phase);
    }
  }

  /** Immutable evidence for one corrected PT boundary point. */
  public static final class EvidencePoint {
    private final double temperatureK;
    private final double pressureBara;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double[] incipientComposition;
    private final double retainedFlashResidual;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;

    /** Creates a point from explicit numerical evidence. */
    public EvidencePoint(double temperatureK, double pressureBara, double beta, double[] phaseZeroComposition,
        double[] phaseOneComposition, double[] incipientComposition, double retainedFlashResidual,
        double tangentPlaneDistance, double stationarityResidual) {
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.beta = beta;
      this.phaseZeroComposition = cloneOrEmpty(phaseZeroComposition);
      this.phaseOneComposition = cloneOrEmpty(phaseOneComposition);
      this.incipientComposition = cloneOrEmpty(incipientComposition);
      this.retainedFlashResidual = retainedFlashResidual;
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.stationarityResidual = stationarityResidual;
    }

    /** Converts one accepted regular-pressure continuation point. */
    public static EvidencePoint from(TwoToThreePhaseBoundaryTracer.Point point) {
      if (point == null) {
        throw new IllegalArgumentException("continuation point must not be null");
      }
      return new EvidencePoint(point.getTemperatureK(), point.getPressureBara(), point.getBeta(),
          point.getPhaseZeroComposition(), point.getPhaseOneComposition(), point.getIncipientComposition(),
          point.getRetainedFlashResidual(), point.getTangentPlaneDistance(), point.getStationarityResidual());
    }

    /**
     * Converts one accepted pseudo-arclength Newton correction using its conservative maximum thermodynamic residual.
     */
    public static EvidencePoint from(TwoToThreePhaseArcLengthCorrector.Result correction) {
      if (correction == null || !correction.isConverged()) {
        throw new IllegalArgumentException("accepted pseudo-arclength correction must be converged");
      }
      TwoToThreePhaseArcLengthCorrector.State state = correction.getState();
      double residual = correction.getThermodynamicMaximumResidual();
      return new EvidencePoint(state.getTemperatureK(), state.getPressureBara(), state.getBeta(),
          state.getPhaseZeroComposition(), state.getPhaseOneComposition(), state.getIncipientComposition(), residual,
          residual, residual);
    }

    /** Converts one strictly corrected two-to-three-phase boundary root. */
    public static EvidencePoint from(TwoToThreePhaseBoundaryPointSolver.Result root) {
      if (root == null || !root.isConverged()) {
        throw new IllegalArgumentException("corrected boundary root must be converged");
      }
      return new EvidencePoint(root.getTemperatureK(), root.getPressureBara(), root.getBeta(),
          root.getPhaseZeroComposition(), root.getPhaseOneComposition(), root.getIncipientComposition(),
          root.getRetainedFlashResidual(), root.getTangentPlaneDistance(), root.getStationarityResidual());
    }

    private static double[] cloneOrEmpty(double[] composition) {
      return composition == null ? new double[0] : composition.clone();
    }

    /** @return corrected temperature in kelvin */
    public double getTemperatureK() {
      return temperatureK;
    }

    /** @return corrected absolute pressure in bara */
    public double getPressureBara() {
      return pressureBara;
    }

    /** @return retained phase-zero mole fraction */
    public double getBeta() {
      return beta;
    }

    /** @return defensive copy of the first retained-phase composition */
    public double[] getPhaseZeroComposition() {
      return phaseZeroComposition.clone();
    }

    /** @return defensive copy of the second retained-phase composition */
    public double[] getPhaseOneComposition() {
      return phaseOneComposition.clone();
    }

    /** @return defensive copy of the incipient-phase composition */
    public double[] getIncipientComposition() {
      return incipientComposition.clone();
    }

    /** @return specified retained-flash residual */
    public double getRetainedFlashResidual() {
      return retainedFlashResidual;
    }

    /** @return incipient tangent-plane distance */
    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    /** @return incipient stationarity residual */
    public double getStationarityResidual() {
      return stationarityResidual;
    }
  }

  /** One named physical boundary branch, including explicit endpoint classifications. */
  public static final class Branch {
    private final String identifier;
    private final BoundaryDefinition definition;
    private final CandidatePhase retainedPhaseZero;
    private final CandidatePhase retainedPhaseOne;
    private final CandidatePhase incipientPhase;
    private final List<EvidencePoint> points;
    private final SpecialPointType startPointType;
    private final SpecialPointType endPointType;
    private final HydrocarbonWaterBoundaryTerminationClassifier.Result startTermination;
    private final HydrocarbonWaterBoundaryTerminationClassifier.Result endTermination;

    public Branch(String identifier, BoundaryDefinition definition, CandidatePhase retainedPhaseZero,
        CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase, List<EvidencePoint> points,
        SpecialPointType startPointType, SpecialPointType endPointType) {
      this(identifier, definition, retainedPhaseZero, retainedPhaseOne, incipientPhase, points, startPointType,
          endPointType, null, null);
    }

    /**
     * Creates a branch with explicit independently refined endpoint evidence.
     *
     * @param identifier unique branch identifier
     * @param definition expected lower/higher phase topology
     * @param retainedPhaseZero first retained phase
     * @param retainedPhaseOne second retained phase
     * @param incipientPhase zero-fraction third phase
     * @param points ordered corrected boundary evidence
     * @param startPointType classified start endpoint
     * @param endPointType classified end endpoint
     * @param startTermination refined start termination evidence, when applicable
     * @param endTermination refined end termination evidence, when applicable
     */
    public Branch(String identifier, BoundaryDefinition definition, CandidatePhase retainedPhaseZero,
        CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase, List<EvidencePoint> points,
        SpecialPointType startPointType, SpecialPointType endPointType,
        HydrocarbonWaterBoundaryTerminationClassifier.Result startTermination,
        HydrocarbonWaterBoundaryTerminationClassifier.Result endTermination) {
      if (identifier == null || identifier.trim().isEmpty() || points == null) {
        throw new IllegalArgumentException("branch identifier and point list are required");
      }
      this.identifier = identifier;
      this.definition = definition;
      this.retainedPhaseZero = retainedPhaseZero;
      this.retainedPhaseOne = retainedPhaseOne;
      this.incipientPhase = incipientPhase;
      this.points = Collections.unmodifiableList(new ArrayList<EvidencePoint>(points));
      this.startPointType = startPointType;
      this.endPointType = endPointType;
      this.startTermination = startTermination;
      this.endTermination = endTermination;
    }
  }

  /** Immutable per-branch gate result. */
  public static final class BranchReport {
    private final String identifier;
    private final int pointCount;
    private final double maximumResidual;
    private final double maximumCompositionJump;
    private final double minimumPhaseDistance;
    private final int projectedIntersectionCount;
    private final int identityIntersectionCount;
    private final double minimumProjectedIntersectionStateDistance;
    private final List<String> violations;
    private final boolean accepted;

    private BranchReport(String identifier, int pointCount, double maximumResidual, double maximumCompositionJump,
        double minimumPhaseDistance, int projectedIntersectionCount, int identityIntersectionCount,
        double minimumProjectedIntersectionStateDistance, List<String> violations, boolean accepted) {
      this.identifier = identifier;
      this.pointCount = pointCount;
      this.maximumResidual = maximumResidual;
      this.maximumCompositionJump = maximumCompositionJump;
      this.minimumPhaseDistance = minimumPhaseDistance;
      this.projectedIntersectionCount = projectedIntersectionCount;
      this.identityIntersectionCount = identityIntersectionCount;
      this.minimumProjectedIntersectionStateDistance = minimumProjectedIntersectionStateDistance;
      this.violations = Collections.unmodifiableList(new ArrayList<String>(violations));
      this.accepted = accepted;
    }

    public String getIdentifier() {
      return identifier;
    }

    public int getPointCount() {
      return pointCount;
    }

    public double getMaximumResidual() {
      return maximumResidual;
    }

    public double getMaximumCompositionJump() {
      return maximumCompositionJump;
    }

    public double getMinimumPhaseDistance() {
      return minimumPhaseDistance;
    }

    public int getProjectedIntersectionCount() {
      return projectedIntersectionCount;
    }

    public int getIdentityIntersectionCount() {
      return identityIntersectionCount;
    }

    public double getMinimumProjectedIntersectionStateDistance() {
      return minimumProjectedIntersectionStateDistance;
    }

    public List<String> getViolations() {
      return violations;
    }

    public boolean isAccepted() {
      return accepted;
    }
  }

  /** Immutable envelope-level decision retaining every per-branch report. */
  public static final class EnvelopeReport {
    private final List<BranchReport> branchReports;
    private final int endpointCandidateCount;
    private final List<String> violations;
    private final boolean accepted;

    private EnvelopeReport(List<BranchReport> branchReports, int endpointCandidateCount, List<String> violations,
        boolean accepted) {
      this.branchReports = Collections.unmodifiableList(new ArrayList<BranchReport>(branchReports));
      this.endpointCandidateCount = endpointCandidateCount;
      this.violations = Collections.unmodifiableList(new ArrayList<String>(violations));
      this.accepted = accepted;
    }

    public List<BranchReport> getBranchReports() {
      return branchReports;
    }

    public List<String> getViolations() {
      return violations;
    }

    public int getEndpointCandidateCount() {
      return endpointCandidateCount;
    }

    public boolean isAccepted() {
      return accepted;
    }
  }
}
