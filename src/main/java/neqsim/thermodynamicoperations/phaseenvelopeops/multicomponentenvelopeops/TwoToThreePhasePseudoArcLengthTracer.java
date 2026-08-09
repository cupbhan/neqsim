package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryGlobalStabilityGate.StationarySeed;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Adaptive pseudo-arclength continuation driver for a specified two-to-three-phase PT boundary. */
public final class TwoToThreePhasePseudoArcLengthTracer {
  private final TwoToThreePhaseArcLengthCorrector corrector;
  private final HydrocarbonWaterBoundaryGlobalStabilityGate globalStabilityGate;
  private double initialArcStep = 0.5;
  private double minimumArcStep = 0.02;
  private double maximumArcStep = 1.0;
  private int maximumRetriesPerPoint = 8;
  private double maximumCompositionJump = 0.35;
  private double minimumForwardCosine = 0.0;
  private boolean domainBoundsEnabled;
  private double minimumTemperatureK;
  private double maximumTemperatureK;
  private double minimumPressureBara;
  private double maximumPressureBara;

  /** Creates a pseudo-arclength driver for three distinct physical phase families. */
  public TwoToThreePhasePseudoArcLengthTracer(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase) {
    corrector = new TwoToThreePhaseArcLengthCorrector(template, retainedPhaseZero, retainedPhaseOne, incipientPhase);
    globalStabilityGate = new HydrocarbonWaterBoundaryGlobalStabilityGate(template);
  }

  /** Sets the underlying Newton and numerical-Jacobian controls. */
  public TwoToThreePhasePseudoArcLengthTracer setCorrectorControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    corrector.setNumericalControls(maximumIterations, residualTolerance, finiteDifferenceStep);
    return this;
  }

  /** Sets adaptive absolute pseudo-arclength step controls in transformed-variable space. */
  public TwoToThreePhasePseudoArcLengthTracer setStepControls(double initialArcStep, double minimumArcStep,
      double maximumArcStep, int maximumRetriesPerPoint) {
    if (!Double.isFinite(initialArcStep) || !Double.isFinite(minimumArcStep) || !Double.isFinite(maximumArcStep)
        || minimumArcStep <= 0.0 || maximumArcStep < minimumArcStep || initialArcStep < minimumArcStep
        || initialArcStep > maximumArcStep || maximumRetriesPerPoint < 1) {
      throw new IllegalArgumentException("invalid pseudo-arclength step controls");
    }
    this.initialArcStep = initialArcStep;
    this.minimumArcStep = minimumArcStep;
    this.maximumArcStep = maximumArcStep;
    this.maximumRetriesPerPoint = maximumRetriesPerPoint;
    return this;
  }

  /** Sets the maximum L1 composition change between adjacent corrected states. */
  public TwoToThreePhasePseudoArcLengthTracer setMaximumCompositionJump(double maximumCompositionJump) {
    if (!Double.isFinite(maximumCompositionJump) || maximumCompositionJump <= 0.0) {
      throw new IllegalArgumentException("maximumCompositionJump must be positive");
    }
    this.maximumCompositionJump = maximumCompositionJump;
    return this;
  }

  /**
   * Sets the minimum scaled-state cosine between the incoming and newly corrected continuation chords.
   *
   * <p>
   * A smooth one-dimensional branch has a forward chord for a sufficiently small pseudo-arclength step. A non-positive
   * cosine means that Newton selected the already traversed side of a fold or branch point. Such a correction must be
   * retried with a shorter step rather than accepted as a new boundary state.
   * </p>
   */
  public TwoToThreePhasePseudoArcLengthTracer setMinimumForwardCosine(double minimumForwardCosine) {
    if (!Double.isFinite(minimumForwardCosine) || minimumForwardCosine < -1.0 || minimumForwardCosine >= 1.0) {
      throw new IllegalArgumentException("minimumForwardCosine must be finite and in [-1, 1)");
    }
    this.minimumForwardCosine = minimumForwardCosine;
    return this;
  }

  /** Stops after accepting the first strict point that crosses one declared PT domain bound. */
  public TwoToThreePhasePseudoArcLengthTracer setDomainBounds(double minimumTemperatureK, double maximumTemperatureK,
      double minimumPressureBara, double maximumPressureBara) {
    if (!Double.isFinite(minimumTemperatureK) || !Double.isFinite(maximumTemperatureK)
        || maximumTemperatureK <= minimumTemperatureK || !Double.isFinite(minimumPressureBara)
        || minimumPressureBara <= 0.0 || !Double.isFinite(maximumPressureBara)
        || maximumPressureBara <= minimumPressureBara) {
      throw new IllegalArgumentException("invalid pseudo-arclength PT domain bounds");
    }
    this.minimumTemperatureK = minimumTemperatureK;
    this.maximumTemperatureK = maximumTemperatureK;
    this.minimumPressureBara = minimumPressureBara;
    this.maximumPressureBara = maximumPressureBara;
    this.domainBoundsEnabled = true;
    return this;
  }

  /**
   * Continues from two accepted seed states.
   *
   * @param previous first seed
   * @param current second seed, defining the initial orientation
   * @param additionalPointCount requested number of newly corrected states
   * @return immutable accepted branch and termination diagnostics
   */
  public Result trace(TwoToThreePhaseArcLengthCorrector.State previous, TwoToThreePhaseArcLengthCorrector.State current,
      int additionalPointCount) {
    if (previous == null || current == null || additionalPointCount < 1) {
      throw new IllegalArgumentException("two seeds and a positive additional point count are required");
    }
    List<TwoToThreePhaseArcLengthCorrector.State> points = new ArrayList<TwoToThreePhaseArcLengthCorrector.State>();
    List<TwoToThreePhaseArcLengthCorrector.Result> acceptedCorrections = new ArrayList<TwoToThreePhaseArcLengthCorrector.Result>();
    List<HydrocarbonWaterBoundaryGlobalStabilityGate.Result> globalStabilityEvidence = new ArrayList<HydrocarbonWaterBoundaryGlobalStabilityGate.Result>();
    List<StationarySeed> stationarySeeds = new ArrayList<StationarySeed>();
    points.add(previous);
    points.add(current);
    double arcStep = initialArcStep;
    int attempts = 0;
    int rejections = 0;
    int localTangentRescueAttempts = 0;
    int acceptedLocalTangentRescues = 0;
    String failureMessage = null;
    TwoToThreePhaseArcLengthCorrector.Result terminalCorrection = null;
    TerminationReason terminationReason = null;

    while (points.size() < additionalPointCount + 2) {
      TwoToThreePhaseArcLengthCorrector.State first = points.get(points.size() - 2);
      TwoToThreePhaseArcLengthCorrector.State second = points.get(points.size() - 1);
      boolean accepted = false;
      String rejectionReason = null;
      for (int retry = 0; retry < maximumRetriesPerPoint; retry++) {
        attempts++;
        TwoToThreePhaseArcLengthCorrector.Result correction = corrector.correctWithArcStep(first, second, arcStep);
        if (correction.isConverged()) {
          double compositionJump = compositionJump(second, correction.getState());
          double physicalDistance = physicalDistance(second, correction.getState());
          double forwardCosine = forwardCosine(first, second, correction.getState());
          boolean rescued = false;
          if (compositionJump <= maximumCompositionJump && physicalDistance > 1.0e-8
              && (!Double.isFinite(forwardCosine) || forwardCosine <= minimumForwardCosine)) {
            LocalTangentRescue rescue = localTangentRescue(first, second, arcStep);
            attempts += rescue.attemptCount;
            localTangentRescueAttempts += rescue.attemptCount;
            if (rescue.correction != null) {
              correction = rescue.correction;
              rescued = true;
              rejections++;
              compositionJump = compositionJump(second, correction.getState());
              physicalDistance = physicalDistance(second, correction.getState());
              forwardCosine = forwardCosine(first, second, correction.getState());
            }
          }
          if (geometricallyAccepted(compositionJump, physicalDistance, forwardCosine)) {
            HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStability = globalStabilityGate
                .evaluate(TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(correction.getState(),
                    correction.getThermodynamicMaximumResidual()), stationarySeeds);
            globalStabilityEvidence.add(globalStability);
            updateStationarySeeds(stationarySeeds, globalStability, correction.getState().getIncipientComposition());
            if (globalStability.isAccepted()) {
              points.add(correction.getState());
              acceptedCorrections.add(correction);
              accepted = true;
              if (rescued) {
                acceptedLocalTangentRescues++;
                arcStep = correction.getArcStep();
              }
              if (correction.getIterations() <= 5 && compositionJump < 0.25 * maximumCompositionJump) {
                arcStep = Math.min(maximumArcStep, 1.2 * arcStep);
              } else if (correction.getIterations() > 12) {
                arcStep = Math.max(minimumArcStep, 0.7 * arcStep);
              }
              break;
            }
            rejectionReason = "global stability gate rejected pseudo-arclength correction: "
                + globalStability.getFailureMessage();
          } else {
            rejectionReason = geometryRejectionReason(compositionJump, physicalDistance, forwardCosine);
          }
        } else {
          rejectionReason = correction.getFailureMessage();
        }
        terminalCorrection = correction;
        rejections++;
        arcStep *= 0.5;
        if (arcStep < minimumArcStep) {
          break;
        }
      }
      if (!accepted) {
        failureMessage = rejectionReason == null ? "pseudo-arclength correction failed" : rejectionReason;
        terminationReason = arcStep < minimumArcStep ? TerminationReason.MINIMUM_ARC_STEP
            : TerminationReason.CORRECTION_RETRY_LIMIT;
        break;
      }
      if (domainBoundsEnabled && outsideDomain(points.get(points.size() - 1))) {
        terminationReason = TerminationReason.DOMAIN_EXIT;
        break;
      }
    }
    boolean completedRequestedPoints = points.size() == additionalPointCount + 2;
    if (completedRequestedPoints) {
      terminationReason = TerminationReason.REQUESTED_POINT_COUNT;
      terminalCorrection = null;
    }
    return new Result(points, acceptedCorrections, globalStabilityEvidence, stationarySeeds, attempts, rejections,
        localTangentRescueAttempts, acceptedLocalTangentRescues, completedRequestedPoints, arcStep, terminationReason,
        terminalCorrection, failureMessage);
  }

  private LocalTangentRescue localTangentRescue(TwoToThreePhaseArcLengthCorrector.State previous,
      TwoToThreePhaseArcLengthCorrector.State current, double arcStep) {
    TwoToThreePhaseArcLengthCorrector.Result selected = null;
    double selectedDistance = Double.NEGATIVE_INFINITY;
    double selectedCosine = Double.NEGATIVE_INFINITY;
    int attemptCount = 0;
    double probeStep = Math.min(maximumArcStep, Math.max(arcStep, initialArcStep));
    while (true) {
      for (int orientation : new int[] { 1, -1 }) {
        attemptCount++;
        TwoToThreePhaseArcLengthCorrector.Result candidate = corrector.correctFromLocalTangent(current, probeStep,
            orientation);
        if (!candidate.isConverged()) {
          continue;
        }
        double compositionJump = compositionJump(current, candidate.getState());
        double physicalDistance = physicalDistance(current, candidate.getState());
        double cosine = forwardCosine(previous, current, candidate.getState());
        if (geometricallyAccepted(compositionJump, physicalDistance, cosine) && (physicalDistance > selectedDistance
            || physicalDistance == selectedDistance && cosine > selectedCosine)) {
          selected = candidate;
          selectedDistance = physicalDistance;
          selectedCosine = cosine;
        }
      }
      if (probeStep >= maximumArcStep) {
        break;
      }
      probeStep = Math.min(maximumArcStep, 2.0 * probeStep);
    }
    return new LocalTangentRescue(selected, attemptCount);
  }

  private boolean geometricallyAccepted(double compositionJump, double physicalDistance, double forwardCosine) {
    return compositionJump <= maximumCompositionJump && physicalDistance > 1.0e-8 && Double.isFinite(forwardCosine)
        && forwardCosine > minimumForwardCosine;
  }

  private String geometryRejectionReason(double compositionJump, double physicalDistance, double forwardCosine) {
    if (compositionJump > maximumCompositionJump) {
      return "composition continuity gate rejected pseudo-arclength correction: jump=" + compositionJump;
    }
    if (!(physicalDistance > 1.0e-8)) {
      return "pseudo-arclength correction returned a duplicate physical state";
    }
    return "forward-orientation gate rejected pseudo-arclength correction: cosine=" + forwardCosine + "; minimum="
        + minimumForwardCosine;
  }

  private static void updateStationarySeeds(List<StationarySeed> seeds,
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result evidence, double[] targetComposition) {
    for (Candidate candidate : evidence.getTrials()) {
      if (!candidate.isConverged() || candidate.isTrivial() || !Double.isFinite(candidate.getTangentPlaneDistance())
          || compositionDistance(candidate.getComposition(), targetComposition) <= 1.0e-4) {
        continue;
      }
      StationarySeed replacement = StationarySeed.from(candidate);
      int nearestIndex = -1;
      double nearestDistance = Double.POSITIVE_INFINITY;
      for (int seedIndex = 0; seedIndex < seeds.size(); seedIndex++) {
        StationarySeed seed = seeds.get(seedIndex);
        if (seed.getPhase() != replacement.getPhase()) {
          continue;
        }
        double distance = compositionDistance(seed.getComposition(), replacement.getComposition());
        if (distance < nearestDistance) {
          nearestDistance = distance;
          nearestIndex = seedIndex;
        }
      }
      if (nearestIndex >= 0 && nearestDistance <= 0.1) {
        seeds.set(nearestIndex, replacement);
      } else if (seeds.size() < 12) {
        seeds.add(replacement);
      }
    }
  }

  private boolean outsideDomain(TwoToThreePhaseArcLengthCorrector.State state) {
    return state.getTemperatureK() < minimumTemperatureK || state.getTemperatureK() > maximumTemperatureK
        || state.getPressureBara() < minimumPressureBara || state.getPressureBara() > maximumPressureBara;
  }

  private static double compositionJump(TwoToThreePhaseArcLengthCorrector.State first,
      TwoToThreePhaseArcLengthCorrector.State second) {
    return Math.max(compositionDistance(first.getPhaseZeroComposition(), second.getPhaseZeroComposition()),
        Math.max(compositionDistance(first.getPhaseOneComposition(), second.getPhaseOneComposition()),
            compositionDistance(first.getIncipientComposition(), second.getIncipientComposition())));
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

  private static double physicalDistance(TwoToThreePhaseArcLengthCorrector.State first,
      TwoToThreePhaseArcLengthCorrector.State second) {
    double temperatureScale = Math.max(1.0, Math.abs(first.getTemperatureK()));
    double pressureScale = Math.max(1.0, Math.abs(first.getPressureBara()));
    double temperatureDistance = (second.getTemperatureK() - first.getTemperatureK()) / temperatureScale;
    double pressureDistance = (second.getPressureBara() - first.getPressureBara()) / pressureScale;
    double compositionDistance = compositionJump(first, second);
    return Math.sqrt(temperatureDistance * temperatureDistance + pressureDistance * pressureDistance
        + compositionDistance * compositionDistance);
  }

  /**
   * Returns the cosine between two consecutive chords in the corrector's scaled transformed-state metric.
   *
   * @param previous state before the current continuation point
   * @param current current continuation point
   * @param candidate newly corrected candidate point
   * @return finite cosine in {@code [-1, 1]}, or {@link Double#NaN} for a degenerate chord
   */
  public static double forwardCosine(TwoToThreePhaseArcLengthCorrector.State previous,
      TwoToThreePhaseArcLengthCorrector.State current, TwoToThreePhaseArcLengthCorrector.State candidate) {
    if (previous == null || current == null || candidate == null) {
      return Double.NaN;
    }
    double[] previousVariables = TwoToThreePhaseArcLengthCorrector.transformedVariables(previous);
    double[] currentVariables = TwoToThreePhaseArcLengthCorrector.transformedVariables(current);
    double[] candidateVariables = TwoToThreePhaseArcLengthCorrector.transformedVariables(candidate);
    if (previousVariables.length != currentVariables.length || currentVariables.length != candidateVariables.length) {
      return Double.NaN;
    }
    double[] weights = TwoToThreePhaseArcLengthCorrector.continuationMetricWeights(current);
    double dot = 0.0;
    double incomingNorm = 0.0;
    double outgoingNorm = 0.0;
    for (int index = 0; index < weights.length; index++) {
      double incoming = weights[index] * (currentVariables[index] - previousVariables[index]);
      double outgoing = weights[index] * (candidateVariables[index] - currentVariables[index]);
      dot += incoming * outgoing;
      incomingNorm += incoming * incoming;
      outgoingNorm += outgoing * outgoing;
    }
    double denominator = Math.sqrt(incomingNorm * outgoingNorm);
    if (!(denominator > 1.0e-20) || !Double.isFinite(denominator)) {
      return Double.NaN;
    }
    return Math.max(-1.0, Math.min(1.0, dot / denominator));
  }

  private static final class LocalTangentRescue {
    private final TwoToThreePhaseArcLengthCorrector.Result correction;
    private final int attemptCount;

    private LocalTangentRescue(TwoToThreePhaseArcLengthCorrector.Result correction, int attemptCount) {
      this.correction = correction;
      this.attemptCount = attemptCount;
    }
  }

  /** Immutable pseudo-arclength branch result. */
  public static final class Result {
    private final List<TwoToThreePhaseArcLengthCorrector.State> points;
    private final List<TwoToThreePhaseArcLengthCorrector.Result> acceptedCorrections;
    private final List<HydrocarbonWaterBoundaryGlobalStabilityGate.Result> globalStabilityEvidence;
    private final List<StationarySeed> stationarySeeds;
    private final int attemptedCorrections;
    private final int rejectedCorrections;
    private final int localTangentRescueAttempts;
    private final int acceptedLocalTangentRescues;
    private final boolean completedRequestedPoints;
    private final double finalArcStep;
    private final TerminationReason terminationReason;
    private final TwoToThreePhaseArcLengthCorrector.Result terminalCorrection;
    private final String failureMessage;

    private Result(List<TwoToThreePhaseArcLengthCorrector.State> points,
        List<TwoToThreePhaseArcLengthCorrector.Result> acceptedCorrections,
        List<HydrocarbonWaterBoundaryGlobalStabilityGate.Result> globalStabilityEvidence,
        List<StationarySeed> stationarySeeds, int attemptedCorrections, int rejectedCorrections,
        int localTangentRescueAttempts, int acceptedLocalTangentRescues, boolean completedRequestedPoints,
        double finalArcStep, TerminationReason terminationReason,
        TwoToThreePhaseArcLengthCorrector.Result terminalCorrection, String failureMessage) {
      this.points = Collections.unmodifiableList(new ArrayList<TwoToThreePhaseArcLengthCorrector.State>(points));
      this.acceptedCorrections = Collections
          .unmodifiableList(new ArrayList<TwoToThreePhaseArcLengthCorrector.Result>(acceptedCorrections));
      this.globalStabilityEvidence = Collections
          .unmodifiableList(new ArrayList<HydrocarbonWaterBoundaryGlobalStabilityGate.Result>(globalStabilityEvidence));
      this.stationarySeeds = Collections.unmodifiableList(new ArrayList<StationarySeed>(stationarySeeds));
      this.attemptedCorrections = attemptedCorrections;
      this.rejectedCorrections = rejectedCorrections;
      this.localTangentRescueAttempts = localTangentRescueAttempts;
      this.acceptedLocalTangentRescues = acceptedLocalTangentRescues;
      this.completedRequestedPoints = completedRequestedPoints;
      this.finalArcStep = finalArcStep;
      this.terminationReason = terminationReason;
      this.terminalCorrection = terminalCorrection;
      this.failureMessage = failureMessage;
    }

    public List<TwoToThreePhaseArcLengthCorrector.State> getPoints() {
      return points;
    }

    /** @return strict Newton evidence for every newly accepted point after the two input seeds */
    public List<TwoToThreePhaseArcLengthCorrector.Result> getAcceptedCorrections() {
      return acceptedCorrections;
    }

    /** @return global TPD evidence for every residual-converged correction considered by the tracer */
    public List<HydrocarbonWaterBoundaryGlobalStabilityGate.Result> getGlobalStabilityEvidence() {
      return globalStabilityEvidence;
    }

    /** @return distinct non-target stationary branches retained for identity-preserving stability checks */
    public List<StationarySeed> getStationarySeeds() {
      return stationarySeeds;
    }

    public int getAttemptedCorrections() {
      return attemptedCorrections;
    }

    public int getRejectedCorrections() {
      return rejectedCorrections;
    }

    /** @return number of two-sided local-null-tangent correction attempts */
    public int getLocalTangentRescueAttempts() {
      return localTangentRescueAttempts;
    }

    /** @return number of forward, globally stable local-tangent candidates accepted into the trace */
    public int getAcceptedLocalTangentRescues() {
      return acceptedLocalTangentRescues;
    }

    public boolean hasCompletedRequestedPoints() {
      return completedRequestedPoints;
    }

    public double getFinalArcStep() {
      return finalArcStep;
    }

    public TerminationReason getTerminationReason() {
      return terminationReason;
    }

    public TwoToThreePhaseArcLengthCorrector.Result getTerminalCorrection() {
      return terminalCorrection;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }

  /** Explicit reason why pseudo-arclength continuation stopped. */
  public enum TerminationReason {
    REQUESTED_POINT_COUNT, DOMAIN_EXIT, MINIMUM_ARC_STEP, CORRECTION_RETRY_LIMIT
  }
}
