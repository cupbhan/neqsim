package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;

/** Closes two identity-compatible paths that share an ordinary state and one target-branch merge. */
public final class HydrocarbonWaterClosedBoundaryLoopAssembler {
  private double maximumEndpointPtDistance = 1.0e-4;
  private double maximumEndpointCompositionDistance = 0.01;

  /** Sets the strict transformed-PT and phase-composition gates for both shared path endpoints. */
  public HydrocarbonWaterClosedBoundaryLoopAssembler setEndpointTolerances(double maximumEndpointPtDistance,
      double maximumEndpointCompositionDistance) {
    if (!positive(maximumEndpointPtDistance) || !positive(maximumEndpointCompositionDistance)) {
      throw new IllegalArgumentException("closed-loop endpoint tolerances must be positive");
    }
    this.maximumEndpointPtDistance = maximumEndpointPtDistance;
    this.maximumEndpointCompositionDistance = maximumEndpointCompositionDistance;
    return this;
  }

  /**
   * Joins two independently traced paths into one ordered closed loop without duplicating their shared endpoints.
   *
   * @param firstPath first path from the shared ordinary state to the merge state
   * @param secondPath second path from the same ordinary state to the same merge state
   * @return immutable closed state sequence and closure diagnostics
   */
  public Result assemble(List<State> firstPath, List<State> secondPath) {
    return assemble(firstPath, secondPath, SpecialPointType.TARGET_BRANCH_MERGE);
  }

  /**
   * Joins two paths at two independently matched ordinary states without inventing a physical merge event.
   *
   * <p>
   * {@link SpecialPointType#CLOSED_LOOP_SEAM} is bookkeeping only: it identifies where the serialized state sequence
   * was joined. It is not an endpoint and must not be promoted as a critical point or target-branch merge.
   * </p>
   */
  public Result assembleOrdinaryClosure(List<State> firstPath, List<State> secondPath) {
    return assemble(firstPath, secondPath, SpecialPointType.CLOSED_LOOP_SEAM);
  }

  private Result assemble(List<State> firstPath, List<State> secondPath, SpecialPointType closureType) {
    if (firstPath == null || secondPath == null || firstPath.size() < 3 || secondPath.size() < 3) {
      return Result.rejected("two paths with at least three states are required");
    }
    List<State> first = orientedToCommonStart(firstPath, secondPath);
    List<State> second = orientedToCommonStart(secondPath, first);
    if (first == null || second == null || !sameTopology(first.get(0), second.get(0))) {
      return Result.rejected("paths do not share one topology-compatible ordinary endpoint");
    }
    State firstStart = first.get(0);
    State secondStart = second.get(0);
    State firstMerge = first.get(first.size() - 1);
    State secondMerge = second.get(second.size() - 1);
    double startPtDistance = ptDistance(firstStart, secondStart);
    double startCompositionDistance = compositionDistance(firstStart, secondStart);
    double mergePtDistance = ptDistance(firstMerge, secondMerge);
    double mergeCompositionDistance = compositionDistance(firstMerge, secondMerge);
    if (startPtDistance > maximumEndpointPtDistance || startCompositionDistance > maximumEndpointCompositionDistance) {
      return Result.rejected("the two paths do not share the same ordinary endpoint identity");
    }
    if (mergePtDistance > maximumEndpointPtDistance || mergeCompositionDistance > maximumEndpointCompositionDistance) {
      return Result.rejected("the two paths do not share the same target-branch merge identity");
    }

    List<State> loop = new ArrayList<State>(first);
    for (int index = second.size() - 2; index >= 1; index--) {
      loop.add(second.get(index));
    }
    loop.add(firstStart);
    int mergePointIndex = first.size() - 1;
    return Result.accepted(loop, mergePointIndex, startPtDistance, startCompositionDistance, mergePtDistance,
        mergeCompositionDistance, closureType);
  }

  private List<State> orientedToCommonStart(List<State> path, List<State> other) {
    List<State> forward = new ArrayList<State>(path);
    State otherStart = other.get(0);
    State otherEnd = other.get(other.size() - 1);
    double forwardScore = endpointIdentityDistance(forward.get(0), otherStart)
        + endpointIdentityDistance(forward.get(forward.size() - 1), otherEnd);
    Collections.reverse(forward);
    double reverseScore = endpointIdentityDistance(forward.get(0), otherStart)
        + endpointIdentityDistance(forward.get(forward.size() - 1), otherEnd);
    if (!Double.isFinite(forwardScore) && !Double.isFinite(reverseScore)) {
      return null;
    }
    if (forwardScore <= reverseScore) {
      Collections.reverse(forward);
    }
    return forward;
  }

  private static double endpointIdentityDistance(State first, State second) {
    return sameTopology(first, second) ? ptDistance(first, second) + compositionDistance(first, second)
        : Double.POSITIVE_INFINITY;
  }

  private static boolean sameTopology(State first, State second) {
    return first.getRetainedPhaseZero() == second.getRetainedPhaseZero()
        && first.getRetainedPhaseOne() == second.getRetainedPhaseOne()
        && first.getIncipientPhase() == second.getIncipientPhase();
  }

  private static double ptDistance(State first, State second) {
    return Math.abs(Math.log(first.getTemperatureK() / second.getTemperatureK()))
        + Math.abs(Math.log(first.getPressureBara() / second.getPressureBara()));
  }

  private static double compositionDistance(State first, State second) {
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

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Immutable state-level closed loop with the internal merge-event position. */
  public static final class Result {
    private final List<State> states;
    private final int mergePointIndex;
    private final SpecialPointType mergePointType;
    private final double ordinaryEndpointPtDistance;
    private final double ordinaryEndpointCompositionDistance;
    private final double mergeEndpointPtDistance;
    private final double mergeEndpointCompositionDistance;
    private final boolean accepted;
    private final String failureMessage;

    private Result(List<State> states, int mergePointIndex, SpecialPointType mergePointType,
        double ordinaryEndpointPtDistance, double ordinaryEndpointCompositionDistance, double mergeEndpointPtDistance,
        double mergeEndpointCompositionDistance, boolean accepted, String failureMessage) {
      this.states = Collections.unmodifiableList(new ArrayList<State>(states));
      this.mergePointIndex = mergePointIndex;
      this.mergePointType = mergePointType;
      this.ordinaryEndpointPtDistance = ordinaryEndpointPtDistance;
      this.ordinaryEndpointCompositionDistance = ordinaryEndpointCompositionDistance;
      this.mergeEndpointPtDistance = mergeEndpointPtDistance;
      this.mergeEndpointCompositionDistance = mergeEndpointCompositionDistance;
      this.accepted = accepted;
      this.failureMessage = failureMessage;
    }

    private static Result accepted(List<State> states, int mergePointIndex, double ordinaryEndpointPtDistance,
        double ordinaryEndpointCompositionDistance, double mergeEndpointPtDistance,
        double mergeEndpointCompositionDistance, SpecialPointType closureType) {
      return new Result(states, mergePointIndex, closureType, ordinaryEndpointPtDistance,
          ordinaryEndpointCompositionDistance, mergeEndpointPtDistance, mergeEndpointCompositionDistance, true, null);
    }

    private static Result rejected(String failureMessage) {
      return new Result(Collections.<State>emptyList(), -1, null, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
          Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, false, failureMessage);
    }

    public List<State> getStates() {
      return states;
    }

    public int getMergePointIndex() {
      return mergePointIndex;
    }

    public SpecialPointType getMergePointType() {
      return mergePointType;
    }

    public double getOrdinaryEndpointPtDistance() {
      return ordinaryEndpointPtDistance;
    }

    public double getOrdinaryEndpointCompositionDistance() {
      return ordinaryEndpointCompositionDistance;
    }

    public double getMergeEndpointPtDistance() {
      return mergeEndpointPtDistance;
    }

    public double getMergeEndpointCompositionDistance() {
      return mergeEndpointCompositionDistance;
    }

    public boolean isAccepted() {
      return accepted;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
