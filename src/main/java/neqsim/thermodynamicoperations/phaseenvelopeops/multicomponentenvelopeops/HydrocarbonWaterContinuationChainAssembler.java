package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;

/** Joins serialized continuation checkpoints only through identity-matching overlapping states. */
public final class HydrocarbonWaterContinuationChainAssembler {
  private double maximumOverlapPtDistance = 1.0e-8;
  private double maximumOverlapCompositionDistance = 1.0e-8;

  /** Sets strict transformed-PT and phase-composition tolerances for checkpoint identity matching. */
  public HydrocarbonWaterContinuationChainAssembler setOverlapTolerances(double maximumOverlapPtDistance,
      double maximumOverlapCompositionDistance) {
    if (!positive(maximumOverlapPtDistance) || !positive(maximumOverlapCompositionDistance)) {
      throw new IllegalArgumentException("checkpoint overlap tolerances must be positive");
    }
    this.maximumOverlapPtDistance = maximumOverlapPtDistance;
    this.maximumOverlapCompositionDistance = maximumOverlapCompositionDistance;
    return this;
  }

  /**
   * Builds every identity-connected chain without inventing points across an unresolved gap.
   *
   * @param segments named serialized checkpoint segments
   * @return all joined and still-disconnected chains with an explicit connection inventory
   */
  public Result assemble(List<Segment> segments) {
    if (segments == null || segments.isEmpty()) {
      throw new IllegalArgumentException("at least one continuation checkpoint segment is required");
    }
    List<MutableChain> chains = new ArrayList<MutableChain>();
    Set<String> identifiers = new HashSet<String>();
    for (Segment segment : segments) {
      if (segment == null || !identifiers.add(segment.getIdentifier())) {
        throw new IllegalArgumentException("checkpoint segment identifiers must be non-null and unique");
      }
      chains.add(new MutableChain(segment));
    }

    List<Join> joins = new ArrayList<Join>();
    while (true) {
      Match best = bestMatch(chains);
      if (best == null) {
        break;
      }
      MutableChain merged = merge(best, joins.size());
      int highIndex = Math.max(best.firstIndex, best.secondIndex);
      int lowIndex = Math.min(best.firstIndex, best.secondIndex);
      chains.remove(highIndex);
      chains.remove(lowIndex);
      chains.add(merged);
      joins.add(merged.lastJoin);
    }

    List<Chain> outputChains = new ArrayList<Chain>();
    for (MutableChain chain : chains) {
      outputChains.add(chain.freeze());
    }
    List<String> violations = new ArrayList<String>();
    if (outputChains.size() > 1) {
      violations.add("DISCONNECTED_CHECKPOINT_CHAINS:" + outputChains.size());
    }
    return new Result(outputChains, joins, violations, outputChains.size() == 1);
  }

  private Match bestMatch(List<MutableChain> chains) {
    Match best = null;
    for (int firstIndex = 0; firstIndex < chains.size(); firstIndex++) {
      for (int secondIndex = firstIndex + 1; secondIndex < chains.size(); secondIndex++) {
        MutableChain first = chains.get(firstIndex);
        MutableChain second = chains.get(secondIndex);
        for (boolean reverseFirst : new boolean[] {false, true}) {
          for (boolean reverseSecond : new boolean[] {false, true}) {
            List<State> orientedFirst = oriented(first.states, reverseFirst);
            List<State> orientedSecond = oriented(second.states, reverseSecond);
            Overlap overlap = overlap(orientedFirst, orientedSecond);
            if (overlap.count == 0) {
              continue;
            }
            Match candidate = new Match(firstIndex, secondIndex, first, second, reverseFirst, reverseSecond, overlap);
            if (best == null || candidate.betterThan(best)) {
              best = candidate;
            }
          }
        }
      }
    }
    return best;
  }

  private Overlap overlap(List<State> first, List<State> second) {
    int maximum = Math.min(first.size(), second.size());
    for (int count = maximum; count >= 1; count--) {
      double maximumPt = 0.0;
      double maximumComposition = 0.0;
      boolean matches = true;
      for (int offset = 0; offset < count; offset++) {
        State left = first.get(first.size() - count + offset);
        State right = second.get(offset);
        if (!sameTopology(left, right)) {
          matches = false;
          break;
        }
        double pt = ptDistance(left, right);
        double composition = compositionDistance(left, right);
        if (pt > maximumOverlapPtDistance || composition > maximumOverlapCompositionDistance) {
          matches = false;
          break;
        }
        maximumPt = Math.max(maximumPt, pt);
        maximumComposition = Math.max(maximumComposition, composition);
      }
      if (matches) {
        return new Overlap(count, maximumPt, maximumComposition);
      }
    }
    return new Overlap(0, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
  }

  private MutableChain merge(Match match, int joinIndex) {
    List<State> first = oriented(match.first.states, match.reverseFirst);
    List<State> second = oriented(match.second.states, match.reverseSecond);
    List<State> states = new ArrayList<State>(first);
    states.addAll(second.subList(match.overlap.count, second.size()));
    List<String> firstSources = orientedStrings(match.first.sourceIdentifiers, match.reverseFirst);
    List<String> secondSources = orientedStrings(match.second.sourceIdentifiers, match.reverseSecond);
    List<String> sources = new ArrayList<String>(firstSources);
    sources.addAll(secondSources);
    Join join = new Join("checkpoint-join-" + joinIndex, firstSources.get(firstSources.size() - 1),
        secondSources.get(0), match.overlap.count, match.overlap.maximumPtDistance,
        match.overlap.maximumCompositionDistance, match.reverseFirst, match.reverseSecond);
    return new MutableChain(states, sources, join);
  }

  private List<State> normalizedStates(Segment segment) {
    List<State> normalized = new ArrayList<State>();
    State topologyReference = null;
    for (State state : segment.getStates()) {
      if (state == null || topologyReference != null && !sameTopology(topologyReference, state)) {
        throw new IllegalArgumentException("every checkpoint segment must contain one non-null topology");
      }
      if (normalized.isEmpty() || ptDistance(normalized.get(normalized.size() - 1), state) > maximumOverlapPtDistance
          || compositionDistance(normalized.get(normalized.size() - 1), state) > maximumOverlapCompositionDistance) {
        normalized.add(state);
      }
      topologyReference = state;
    }
    if (normalized.isEmpty()) {
      throw new IllegalArgumentException("checkpoint segment must contain at least one state");
    }
    return normalized;
  }

  private static List<State> oriented(List<State> states, boolean reverse) {
    List<State> result = new ArrayList<State>(states);
    if (reverse) {
      Collections.reverse(result);
    }
    return result;
  }

  private static List<String> orientedStrings(List<String> values, boolean reverse) {
    List<String> result = new ArrayList<String>(values);
    if (reverse) {
      Collections.reverse(result);
    }
    return result;
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

  /** One named, ordered checkpoint segment. */
  public static final class Segment {
    private final String identifier;
    private final List<State> states;

    public Segment(String identifier, List<State> states) {
      if (identifier == null || identifier.trim().isEmpty() || states == null || states.isEmpty()) {
        throw new IllegalArgumentException("checkpoint segment identifier and states are required");
      }
      this.identifier = identifier;
      this.states = Collections.unmodifiableList(new ArrayList<State>(states));
    }

    public String getIdentifier() {
      return identifier;
    }

    public List<State> getStates() {
      return states;
    }
  }

  /** One identity-proven overlap removed while joining two checkpoint chains. */
  public static final class Join {
    private final String identifier;
    private final String firstSegmentIdentifier;
    private final String secondSegmentIdentifier;
    private final int overlapStateCount;
    private final double maximumPtDistance;
    private final double maximumCompositionDistance;
    private final boolean reversedFirst;
    private final boolean reversedSecond;

    private Join(String identifier, String firstSegmentIdentifier, String secondSegmentIdentifier,
        int overlapStateCount, double maximumPtDistance, double maximumCompositionDistance, boolean reversedFirst,
        boolean reversedSecond) {
      this.identifier = identifier;
      this.firstSegmentIdentifier = firstSegmentIdentifier;
      this.secondSegmentIdentifier = secondSegmentIdentifier;
      this.overlapStateCount = overlapStateCount;
      this.maximumPtDistance = maximumPtDistance;
      this.maximumCompositionDistance = maximumCompositionDistance;
      this.reversedFirst = reversedFirst;
      this.reversedSecond = reversedSecond;
    }

    public String getIdentifier() {
      return identifier;
    }

    public String getFirstSegmentIdentifier() {
      return firstSegmentIdentifier;
    }

    public String getSecondSegmentIdentifier() {
      return secondSegmentIdentifier;
    }

    public int getOverlapStateCount() {
      return overlapStateCount;
    }

    public double getMaximumPtDistance() {
      return maximumPtDistance;
    }

    public double getMaximumCompositionDistance() {
      return maximumCompositionDistance;
    }

    public boolean isReversedFirst() {
      return reversedFirst;
    }

    public boolean isReversedSecond() {
      return reversedSecond;
    }
  }

  /** One assembled or still isolated chain. */
  public static final class Chain {
    private final List<State> states;
    private final List<String> sourceIdentifiers;

    private Chain(List<State> states, List<String> sourceIdentifiers) {
      this.states = Collections.unmodifiableList(new ArrayList<State>(states));
      this.sourceIdentifiers = Collections.unmodifiableList(new ArrayList<String>(sourceIdentifiers));
    }

    public List<State> getStates() {
      return states;
    }

    public List<String> getSourceIdentifiers() {
      return sourceIdentifiers;
    }
  }

  /** Immutable full checkpoint-chain inventory. */
  public static final class Result {
    private final List<Chain> chains;
    private final List<Join> joins;
    private final List<String> violations;
    private final boolean singleConnectedChain;

    private Result(List<Chain> chains, List<Join> joins, List<String> violations, boolean singleConnectedChain) {
      this.chains = Collections.unmodifiableList(new ArrayList<Chain>(chains));
      this.joins = Collections.unmodifiableList(new ArrayList<Join>(joins));
      this.violations = Collections.unmodifiableList(new ArrayList<String>(violations));
      this.singleConnectedChain = singleConnectedChain;
    }

    public List<Chain> getChains() {
      return chains;
    }

    public List<Join> getJoins() {
      return joins;
    }

    public List<String> getViolations() {
      return violations;
    }

    public boolean isSingleConnectedChain() {
      return singleConnectedChain;
    }
  }

  private final class MutableChain {
    private final List<State> states;
    private final List<String> sourceIdentifiers;
    private final Join lastJoin;

    private MutableChain(Segment segment) {
      this.states = normalizedStates(segment);
      this.sourceIdentifiers = new ArrayList<String>();
      this.sourceIdentifiers.add(segment.getIdentifier());
      this.lastJoin = null;
    }

    private MutableChain(List<State> states, List<String> sourceIdentifiers, Join lastJoin) {
      this.states = states;
      this.sourceIdentifiers = sourceIdentifiers;
      this.lastJoin = lastJoin;
    }

    private Chain freeze() {
      return new Chain(states, sourceIdentifiers);
    }
  }

  private static final class Overlap {
    private final int count;
    private final double maximumPtDistance;
    private final double maximumCompositionDistance;

    private Overlap(int count, double maximumPtDistance, double maximumCompositionDistance) {
      this.count = count;
      this.maximumPtDistance = maximumPtDistance;
      this.maximumCompositionDistance = maximumCompositionDistance;
    }
  }

  private static final class Match {
    private final int firstIndex;
    private final int secondIndex;
    private final MutableChain first;
    private final MutableChain second;
    private final boolean reverseFirst;
    private final boolean reverseSecond;
    private final Overlap overlap;

    private Match(int firstIndex, int secondIndex, MutableChain first, MutableChain second, boolean reverseFirst,
        boolean reverseSecond, Overlap overlap) {
      this.firstIndex = firstIndex;
      this.secondIndex = secondIndex;
      this.first = first;
      this.second = second;
      this.reverseFirst = reverseFirst;
      this.reverseSecond = reverseSecond;
      this.overlap = overlap;
    }

    private boolean betterThan(Match other) {
      if (overlap.count != other.overlap.count) {
        return overlap.count > other.overlap.count;
      }
      return overlap.maximumPtDistance + overlap.maximumCompositionDistance < other.overlap.maximumPtDistance
          + other.overlap.maximumCompositionDistance;
    }
  }
}
