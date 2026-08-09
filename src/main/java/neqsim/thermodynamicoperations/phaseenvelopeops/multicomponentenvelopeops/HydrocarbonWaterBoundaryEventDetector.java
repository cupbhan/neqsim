package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryGlobalStabilityGate.RetainedPhaseLocalStability;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Detects global PT extrema and endpoint phase-coalescence events on a corrected boundary branch. */
public final class HydrocarbonWaterBoundaryEventDetector {
  private double phaseCoalescenceDistance = 1.0e-4;
  private double endpointWindowFraction = 0.1;
  private double stationaryFoldEigenvalueTolerance = 1.0e-7;
  private double stationaryFoldTpdTolerance = 1.0e-7;
  private double stationaryFoldGradientTolerance = 1.0e-6;
  private double stationaryFoldConditionNumber = 1.0e7;

  /** Sets the L1 phase-distance and endpoint-window controls for critical-endpoint candidates. */
  public HydrocarbonWaterBoundaryEventDetector setEndpointControls(double phaseCoalescenceDistance,
      double endpointWindowFraction) {
    if (!Double.isFinite(phaseCoalescenceDistance) || phaseCoalescenceDistance <= 0.0
        || !Double.isFinite(endpointWindowFraction) || endpointWindowFraction <= 0.0 || endpointWindowFraction > 0.5) {
      throw new IllegalArgumentException("invalid boundary event detector controls");
    }
    this.phaseCoalescenceDistance = phaseCoalescenceDistance;
    this.endpointWindowFraction = endpointWindowFraction;
    return this;
  }

  /** Sets the numerical gates for a stability-stationary-fold continuation endpoint. */
  public HydrocarbonWaterBoundaryEventDetector setStationaryFoldControls(double eigenvalueTolerance,
      double tangentPlaneTolerance, double gradientTolerance, double minimumJacobianConditionNumber) {
    if (!Double.isFinite(eigenvalueTolerance) || eigenvalueTolerance <= 0.0 || !Double.isFinite(tangentPlaneTolerance)
        || tangentPlaneTolerance <= 0.0 || !Double.isFinite(gradientTolerance) || gradientTolerance <= 0.0
        || !Double.isFinite(minimumJacobianConditionNumber) || minimumJacobianConditionNumber <= 1.0) {
      throw new IllegalArgumentException("invalid stationary-fold endpoint controls");
    }
    this.stationaryFoldEigenvalueTolerance = eigenvalueTolerance;
    this.stationaryFoldTpdTolerance = tangentPlaneTolerance;
    this.stationaryFoldGradientTolerance = gradientTolerance;
    this.stationaryFoldConditionNumber = minimumJacobianConditionNumber;
    return this;
  }

  /**
   * Detects special points from a pseudo-arclength ordered branch.
   *
   * @param states corrected states in continuation order
   * @return immutable events, including global interior T/P maxima and endpoint coalescence candidates
   */
  public Result detect(List<TwoToThreePhaseArcLengthCorrector.State> states) {
    if (states == null || states.size() < 3) {
      throw new IllegalArgumentException("at least three corrected states are required");
    }
    for (TwoToThreePhaseArcLengthCorrector.State state : states) {
      if (state == null || !Double.isFinite(state.getTemperatureK()) || !Double.isFinite(state.getPressureBara())) {
        throw new IllegalArgumentException("event detection requires finite corrected states");
      }
    }
    List<Event> events = new ArrayList<Event>();
    int maximumTemperatureIndex = globalMaximumTemperatureIndex(states);
    if (maximumTemperatureIndex > 0 && maximumTemperatureIndex < states.size() - 1) {
      events.add(interpolateExtremum(states, maximumTemperatureIndex, true, SpecialPointType.CRICONDENTHERM));
    }
    int maximumPressureIndex = globalMaximumPressureIndex(states);
    if (maximumPressureIndex > 0 && maximumPressureIndex < states.size() - 1) {
      events.add(interpolateExtremum(states, maximumPressureIndex, false, SpecialPointType.CRICONDENBAR));
    }
    int endpointWindow = Math.max(1, (int) Math.ceil(endpointWindowFraction * (states.size() - 1)));
    Event startCoalescence = endpointCoalescence(states, 0, endpointWindow, 1);
    if (startCoalescence != null) {
      events.add(startCoalescence);
    }
    Event endCoalescence = endpointCoalescence(states, states.size() - 1, endpointWindow, -1);
    if (endCoalescence != null) {
      events.add(endCoalescence);
    }
    return new Result(events, maximumTemperatureIndex, maximumPressureIndex);
  }

  /**
   * Detects branch events and classifies a rigorously gated stationary-fold termination.
   *
   * @param trace terminated or completed pseudo-arclength trace
   * @param curvatureAnalyzer incipient-phase TPD Hessian analyzer
   * @return detected events including a stability fold only when every rank-loss gate passes
   */
  public Result detect(TwoToThreePhasePseudoArcLengthTracer.Result trace,
      IncipientPhaseCurvatureAnalyzer curvatureAnalyzer) {
    if (trace == null || curvatureAnalyzer == null) {
      throw new IllegalArgumentException("trace and curvature analyzer are required");
    }
    Result base = detect(trace.getPoints());
    List<Event> events = new ArrayList<Event>(base.events);
    if (trace.getTerminationReason() == TwoToThreePhasePseudoArcLengthTracer.TerminationReason.MINIMUM_ARC_STEP
        && trace.getTerminalCorrection() != null && !trace.getPoints().isEmpty()) {
      TwoToThreePhaseArcLengthCorrector.Result terminalCorrection = trace.getTerminalCorrection();
      TwoToThreePhaseArcLengthCorrector.State endpoint = trace.getPoints().get(trace.getPoints().size() - 1);
      IncipientPhaseCurvatureAnalyzer.Result curvature = curvatureAnalyzer.analyze(endpoint);
      boolean conditionGate = Double.isFinite(terminalCorrection.getJacobianConditionNumber())
          && terminalCorrection.getJacobianConditionNumber() >= stationaryFoldConditionNumber;
      boolean tpdGate = Math.abs(curvature.getTangentPlaneDistance()) <= stationaryFoldTpdTolerance;
      boolean gradientGate = curvature.getMaximumGradient() <= stationaryFoldGradientTolerance;
      boolean eigenvalueGate = Math.abs(curvature.getMinimumEigenvalue()) <= stationaryFoldEigenvalueTolerance
          && Math.abs(curvature.getSecondEigenvalue()) > 10.0 * Math.abs(curvature.getMinimumEigenvalue());
      if (conditionGate && tpdGate && gradientGate && eigenvalueGate) {
        int endpointIndex = trace.getPoints().size() - 1;
        events.add(new Event(SpecialPointType.STABILITY_STATIONARY_FOLD, endpoint.getTemperatureK(),
            endpoint.getPressureBara(), Math.max(0, endpointIndex - 1), endpointIndex,
            Math.abs(curvature.getMinimumEigenvalue()), false));
      }
    }
    return new Result(events, base.maximumTemperaturePointIndex, base.maximumPressurePointIndex);
  }

  /** Adds an independently solved critical endpoint; rejected mathematical candidates add no event. */
  public Result attachCriticalEndpoint(Result detected, HydrocarbonWaterCriticalEndpointSolver.Result endpoint) {
    if (detected == null || endpoint == null) {
      throw new IllegalArgumentException("detected events and critical endpoint result are required");
    }
    if (!endpoint.isPhysicalEndpoint()) {
      return detected;
    }
    List<Event> events = new ArrayList<Event>();
    for (Event event : detected.events) {
      if (event.type != SpecialPointType.CRITICAL_END_POINT) {
        events.add(event);
      }
    }
    double quality = maximumAbs(endpoint.getMinimumEigenvalue(), endpoint.getThirdDirectionalDerivative(),
        endpoint.getHomogeneousTpdGradient(), endpoint.getPhaseDistance(), endpoint.getIncipientTangentPlaneDistance(),
        endpoint.getIncipientStationarityResidual());
    events.add(new Event(SpecialPointType.CRITICAL_END_POINT, endpoint.getTemperatureK(), endpoint.getPressureBara(),
        -1, -1, quality, false));
    return new Result(events, detected.maximumTemperaturePointIndex, detected.maximumPressurePointIndex);
  }

  /**
   * Adds a retained-phase spinodal or alternate-stationary-branch onset when the next root fails global stability.
   *
   * <p>
   * A negative retained-phase stationary-Jacobian eigenvalue takes precedence because it locates loss of local
   * composition stability even before a generic TPD seed can find the emerging non-trivial branch. Otherwise the method
   * records a lower, negative-TPD stationary composition. Neither event is mislabeled as a critical endpoint or fourth
   * phase.
   * </p>
   *
   * @param detected already detected branch events
   * @param trace terminated pseudo-arclength trace retaining all global-stability trials
   * @return events with an auditable alternate-branch onset when the evidence gates pass
   */
  public Result attachGlobalStabilityLimit(Result detected, TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    if (detected == null || trace == null) {
      throw new IllegalArgumentException("detected events and pseudo-arclength trace are required");
    }
    if (trace.hasCompletedRequestedPoints() || trace.getPoints().isEmpty()
        || trace.getGlobalStabilityEvidence().isEmpty()) {
      return detected;
    }
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejected = trace.getGlobalStabilityEvidence()
        .get(trace.getGlobalStabilityEvidence().size() - 1);
    TwoToThreePhaseBoundaryPointSolver.Result root = rejected.getBoundaryRoot();
    RetainedPhaseLocalStability unstableRetainedPhase = null;
    for (RetainedPhaseLocalStability retainedPhase : rejected.getRetainedPhaseStability()) {
      if (!retainedPhase.isAccepted() && Double.isFinite(retainedPhase.getMinimumEigenvalue())
          && retainedPhase.getMinimumEigenvalue() < 0.0 && (unstableRetainedPhase == null
              || retainedPhase.getMinimumEigenvalue() < unstableRetainedPhase.getMinimumEigenvalue())) {
        unstableRetainedPhase = retainedPhase;
      }
    }
    if (!rejected.isAccepted() && rejected.isTargetMatched() && root != null && unstableRetainedPhase != null) {
      List<Event> events = new ArrayList<Event>();
      for (Event event : detected.events) {
        if (event.type != SpecialPointType.RETAINED_PHASE_SPINODAL
            && event.type != SpecialPointType.ALTERNATE_STATIONARY_BRANCH_ONSET) {
          events.add(event);
        }
      }
      int endpointIndex = trace.getPoints().size() - 1;
      double eigenvalue = unstableRetainedPhase.getMinimumEigenvalue();
      String diagnostic = "retained " + unstableRetainedPhase.getPhase()
          + " stationary Jacobian crossed into local instability: eigenvalue=" + eigenvalue
          + (rejected.getFailureMessage() == null ? "" : "; " + rejected.getFailureMessage());
      events.add(new Event(SpecialPointType.RETAINED_PHASE_SPINODAL, root.getTemperatureK(), root.getPressureBara(),
          endpointIndex, -1, Math.abs(eigenvalue), false, unstableRetainedPhase.getPhase(), diagnostic));
      return new Result(events, detected.maximumTemperaturePointIndex, detected.maximumPressurePointIndex);
    }

    Candidate destabilizing = rejected.getMostUnstableCandidate();
    double minimumTpd = rejected.getMinimumNonTrivialTangentPlaneDistance();
    if (rejected.isAccepted() || !rejected.isTargetMatched() || root == null || destabilizing == null
        || !destabilizing.isConverged() || destabilizing.isTrivial() || !Double.isFinite(minimumTpd)
        || minimumTpd >= 0.0) {
      return detected;
    }
    List<Event> events = new ArrayList<Event>();
    for (Event event : detected.events) {
      if (event.type != SpecialPointType.ALTERNATE_STATIONARY_BRANCH_ONSET) {
        events.add(event);
      }
    }
    int endpointIndex = trace.getPoints().size() - 1;
    String diagnostic = "lower-TPD " + destabilizing.getPhase() + " stationary branch: minimumTPD=" + minimumTpd
        + (rejected.getFailureMessage() == null ? "" : "; " + rejected.getFailureMessage());
    events.add(new Event(SpecialPointType.ALTERNATE_STATIONARY_BRANCH_ONSET, root.getTemperatureK(),
        root.getPressureBara(), endpointIndex, -1, Math.abs(minimumTpd), false, destabilizing.getPhase(), diagnostic));
    return new Result(events, detected.maximumTemperaturePointIndex, detected.maximumPressurePointIndex);
  }

  /**
   * Adds a domain-exit endpoint only when continuation actually reaches or crosses a declared PT bound.
   *
   * <p>
   * A minimum-step failure inside the domain is not a domain exit. This distinction prevents a numerical fold from
   * being accepted merely by labeling it as a scan limit.
   * </p>
   *
   * @param detected already detected branch events
   * @param trace terminated pseudo-arclength trace
   * @param minimumTemperatureK lower solver temperature bound
   * @param maximumTemperatureK upper solver temperature bound
   * @param minimumPressureBara lower solver pressure bound
   * @param maximumPressureBara upper solver pressure bound
   * @param transformedDistanceTolerance maximum relative-T or log-P distance to a bound
   * @return events with a rigorously gated domain exit when applicable
   */
  public Result attachDomainExit(Result detected, TwoToThreePhasePseudoArcLengthTracer.Result trace,
      double minimumTemperatureK, double maximumTemperatureK, double minimumPressureBara, double maximumPressureBara,
      double transformedDistanceTolerance) {
    if (detected == null || trace == null || trace.getPoints().isEmpty() || !Double.isFinite(minimumTemperatureK)
        || !Double.isFinite(maximumTemperatureK) || maximumTemperatureK <= minimumTemperatureK
        || !Double.isFinite(minimumPressureBara) || minimumPressureBara <= 0.0 || !Double.isFinite(maximumPressureBara)
        || maximumPressureBara <= minimumPressureBara || !Double.isFinite(transformedDistanceTolerance)
        || transformedDistanceTolerance <= 0.0) {
      throw new IllegalArgumentException("invalid trace or domain-exit controls");
    }
    TwoToThreePhaseArcLengthCorrector.State endpoint = trace.getPoints().get(trace.getPoints().size() - 1);
    double temperatureDistance = Math.min(
        Math.abs(endpoint.getTemperatureK() - minimumTemperatureK) / Math.max(1.0, minimumTemperatureK),
        Math.abs(endpoint.getTemperatureK() - maximumTemperatureK) / Math.max(1.0, maximumTemperatureK));
    double pressureDistance = Math.min(Math.abs(Math.log(endpoint.getPressureBara() / minimumPressureBara)),
        Math.abs(Math.log(endpoint.getPressureBara() / maximumPressureBara)));
    double distance = Math.min(temperatureDistance, pressureDistance);
    boolean endpointOutside = !insideDomain(endpoint, minimumTemperatureK, maximumTemperatureK, minimumPressureBara,
        maximumPressureBara);
    boolean earlierPointInside = false;
    for (int index = 0; index < trace.getPoints().size() - 1; index++) {
      if (insideDomain(trace.getPoints().get(index), minimumTemperatureK, maximumTemperatureK, minimumPressureBara,
          maximumPressureBara)) {
        earlierPointInside = true;
        break;
      }
    }
    if (endpointOutside && earlierPointInside) {
      List<Event> events = withoutDomainExit(detected);
      int endpointIndex = trace.getPoints().size() - 1;
      events.add(new Event(SpecialPointType.DOMAIN_EXIT, endpoint.getTemperatureK(), endpoint.getPressureBara(),
          endpointIndex, endpointIndex, distance, false, null, "continuation crossed a declared PT domain bound"));
      return new Result(events, detected.maximumTemperaturePointIndex, detected.maximumPressurePointIndex);
    }
    if (trace.hasCompletedRequestedPoints()) {
      return detected;
    }
    if (!Double.isFinite(distance) || distance > transformedDistanceTolerance) {
      return detected;
    }
    List<Event> events = withoutDomainExit(detected);
    int endpointIndex = trace.getPoints().size() - 1;
    events.add(new Event(SpecialPointType.DOMAIN_EXIT, endpoint.getTemperatureK(), endpoint.getPressureBara(),
        endpointIndex, endpointIndex, distance, false));
    return new Result(events, detected.maximumTemperaturePointIndex, detected.maximumPressurePointIndex);
  }

  private static boolean insideDomain(TwoToThreePhaseArcLengthCorrector.State point, double minimumTemperatureK,
      double maximumTemperatureK, double minimumPressureBara, double maximumPressureBara) {
    return point.getTemperatureK() >= minimumTemperatureK && point.getTemperatureK() <= maximumTemperatureK
        && point.getPressureBara() >= minimumPressureBara && point.getPressureBara() <= maximumPressureBara;
  }

  private static List<Event> withoutDomainExit(Result detected) {
    List<Event> events = new ArrayList<Event>();
    for (Event event : detected.events) {
      if (event.type != SpecialPointType.DOMAIN_EXIT) {
        events.add(event);
      }
    }
    return events;
  }

  private Event endpointCoalescence(List<TwoToThreePhaseArcLengthCorrector.State> states, int endpointIndex, int window,
      int direction) {
    int bestIndex = endpointIndex;
    double bestDistance = minimumPhaseDistance(states.get(endpointIndex));
    for (int offset = 1; offset <= window; offset++) {
      int index = endpointIndex + direction * offset;
      if (index < 0 || index >= states.size()) {
        break;
      }
      double distance = minimumPhaseDistance(states.get(index));
      if (distance < bestDistance) {
        bestDistance = distance;
        bestIndex = index;
      }
    }
    if (bestIndex != endpointIndex || bestDistance > phaseCoalescenceDistance) {
      return null;
    }
    TwoToThreePhaseArcLengthCorrector.State endpoint = states.get(endpointIndex);
    return new Event(SpecialPointType.PHASE_COALESCENCE_CANDIDATE, endpoint.getTemperatureK(),
        endpoint.getPressureBara(), endpointIndex, endpointIndex, bestDistance, false);
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

  private static Event interpolateExtremum(List<TwoToThreePhaseArcLengthCorrector.State> states, int centerIndex,
      boolean temperatureExtremum, SpecialPointType type) {
    TwoToThreePhaseArcLengthCorrector.State previous = states.get(centerIndex - 1);
    TwoToThreePhaseArcLengthCorrector.State center = states.get(centerIndex);
    TwoToThreePhaseArcLengthCorrector.State next = states.get(centerIndex + 1);
    double previousValue = temperatureExtremum ? previous.getTemperatureK() : previous.getPressureBara();
    double centerValue = temperatureExtremum ? center.getTemperatureK() : center.getPressureBara();
    double nextValue = temperatureExtremum ? next.getTemperatureK() : next.getPressureBara();
    double previousDistance = transformedDistance(previous, center);
    double nextDistance = transformedDistance(center, next);
    double leftSlope = (centerValue - previousValue) / previousDistance;
    double rightSlope = (nextValue - centerValue) / nextDistance;
    double quadratic = (rightSlope - leftSlope) / (previousDistance + nextDistance);
    double linear = leftSlope + quadratic * previousDistance;
    double offset = Math.abs(quadratic) > 1.0e-20 ? -linear / (2.0 * quadratic) : 0.0;
    offset = Math.max(-previousDistance, Math.min(nextDistance, offset));
    double extremumValue = centerValue + linear * offset + quadratic * offset * offset;
    double otherPrevious = temperatureExtremum ? previous.getPressureBara() : previous.getTemperatureK();
    double otherCenter = temperatureExtremum ? center.getPressureBara() : center.getTemperatureK();
    double otherNext = temperatureExtremum ? next.getPressureBara() : next.getTemperatureK();
    double otherLeftSlope = (otherCenter - otherPrevious) / previousDistance;
    double otherRightSlope = (otherNext - otherCenter) / nextDistance;
    double otherQuadratic = (otherRightSlope - otherLeftSlope) / (previousDistance + nextDistance);
    double otherLinear = otherLeftSlope + otherQuadratic * previousDistance;
    double otherValue = otherCenter + otherLinear * offset + otherQuadratic * offset * offset;
    double temperatureK = temperatureExtremum ? extremumValue : otherValue;
    double pressureBara = temperatureExtremum ? otherValue : extremumValue;
    double bracketWidth = Math.hypot(next.getTemperatureK() - previous.getTemperatureK(),
        next.getPressureBara() - previous.getPressureBara());
    return new Event(type, temperatureK, pressureBara, centerIndex - 1, centerIndex + 1, bracketWidth, true);
  }

  private static double transformedDistance(TwoToThreePhaseArcLengthCorrector.State first,
      TwoToThreePhaseArcLengthCorrector.State second) {
    double squaredDistance = square(Math.log(second.getTemperatureK() / first.getTemperatureK()))
        + square(Math.log(second.getPressureBara() / first.getPressureBara()))
        + square(second.getBeta() - first.getBeta());
    double[] firstPhaseZero = first.getPhaseZeroComposition();
    double[] firstPhaseOne = first.getPhaseOneComposition();
    double[] firstIncipient = first.getIncipientComposition();
    double[] secondPhaseZero = second.getPhaseZeroComposition();
    double[] secondPhaseOne = second.getPhaseOneComposition();
    double[] secondIncipient = second.getIncipientComposition();
    for (int componentIndex = 0; componentIndex < firstPhaseZero.length; componentIndex++) {
      double firstLogK = Math
          .log(Math.max(firstPhaseZero[componentIndex], 1.0e-100) / Math.max(firstPhaseOne[componentIndex], 1.0e-100));
      double secondLogK = Math.log(
          Math.max(secondPhaseZero[componentIndex], 1.0e-100) / Math.max(secondPhaseOne[componentIndex], 1.0e-100));
      squaredDistance += square(secondLogK - firstLogK);
      squaredDistance += square(Math.log(Math.max(secondIncipient[componentIndex], 1.0e-100))
          - Math.log(Math.max(firstIncipient[componentIndex], 1.0e-100)));
    }
    return Math.max(1.0e-12, Math.sqrt(squaredDistance));
  }

  private static double square(double value) {
    return value * value;
  }

  private static int globalMaximumTemperatureIndex(List<TwoToThreePhaseArcLengthCorrector.State> states) {
    int index = 0;
    for (int candidate = 1; candidate < states.size(); candidate++) {
      if (states.get(candidate).getTemperatureK() > states.get(index).getTemperatureK()) {
        index = candidate;
      }
    }
    return index;
  }

  private static int globalMaximumPressureIndex(List<TwoToThreePhaseArcLengthCorrector.State> states) {
    int index = 0;
    for (int candidate = 1; candidate < states.size(); candidate++) {
      if (states.get(candidate).getPressureBara() > states.get(index).getPressureBara()) {
        index = candidate;
      }
    }
    return index;
  }

  private static double minimumPhaseDistance(TwoToThreePhaseArcLengthCorrector.State state) {
    double[] first = state.getPhaseZeroComposition();
    double[] second = state.getPhaseOneComposition();
    double[] third = state.getIncipientComposition();
    return Math.min(compositionDistance(first, second),
        Math.min(compositionDistance(first, third), compositionDistance(second, third)));
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

  /** One detected and optionally locally interpolated special point. */
  public static final class Event {
    private final SpecialPointType type;
    private final double temperatureK;
    private final double pressureBara;
    private final int lowerPointIndex;
    private final int upperPointIndex;
    private final double qualityMeasure;
    private final boolean interpolated;
    private final CandidatePhase destabilizingPhase;
    private final String diagnostic;

    private Event(SpecialPointType type, double temperatureK, double pressureBara, int lowerPointIndex,
        int upperPointIndex, double qualityMeasure, boolean interpolated) {
      this(type, temperatureK, pressureBara, lowerPointIndex, upperPointIndex, qualityMeasure, interpolated, null,
          null);
    }

    private Event(SpecialPointType type, double temperatureK, double pressureBara, int lowerPointIndex,
        int upperPointIndex, double qualityMeasure, boolean interpolated, CandidatePhase destabilizingPhase,
        String diagnostic) {
      this.type = type;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.lowerPointIndex = lowerPointIndex;
      this.upperPointIndex = upperPointIndex;
      this.qualityMeasure = qualityMeasure;
      this.interpolated = interpolated;
      this.destabilizingPhase = destabilizingPhase;
      this.diagnostic = diagnostic;
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

    public int getLowerPointIndex() {
      return lowerPointIndex;
    }

    public int getUpperPointIndex() {
      return upperPointIndex;
    }

    public double getQualityMeasure() {
      return qualityMeasure;
    }

    public boolean isInterpolated() {
      return interpolated;
    }

    /** @return phase family of the lower-TPD stationary point, or null for ordinary events */
    public CandidatePhase getDestabilizingPhase() {
      return destabilizingPhase;
    }

    /** @return immutable numerical diagnostic, or null for ordinary events */
    public String getDiagnostic() {
      return diagnostic;
    }
  }

  /** Immutable detected-event set and raw global-extremum indices. */
  public static final class Result {
    private final List<Event> events;
    private final int maximumTemperaturePointIndex;
    private final int maximumPressurePointIndex;

    private Result(List<Event> events, int maximumTemperaturePointIndex, int maximumPressurePointIndex) {
      this.events = Collections.unmodifiableList(new ArrayList<Event>(events));
      this.maximumTemperaturePointIndex = maximumTemperaturePointIndex;
      this.maximumPressurePointIndex = maximumPressurePointIndex;
    }

    public List<Event> getEvents() {
      return events;
    }

    public Event getEvent(SpecialPointType type) {
      for (Event event : events) {
        if (event.type == type) {
          return event;
        }
      }
      return null;
    }

    public int getMaximumTemperaturePointIndex() {
      return maximumTemperaturePointIndex;
    }

    public int getMaximumPressurePointIndex() {
      return maximumPressurePointIndex;
    }
  }
}
