package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Phase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Finds stable 2P/GOW topology brackets before any zero-TPD boundary correction is attempted. */
public final class HydrocarbonWaterStableRegionTransitionScanner {
  private final SystemInterface template;
  private int maximumBisections = 80;
  private double temperatureToleranceK = 1.0e-4;
  private double minimumPhaseFraction = 1.0e-9;

  /** Creates a non-destructive stable-region transition scanner. */
  public HydrocarbonWaterStableRegionTransitionScanner(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets topology-bracket bisection, temperature, and active-phase-fraction controls. */
  public HydrocarbonWaterStableRegionTransitionScanner setNumericalControls(int maximumBisections,
      double temperatureToleranceK, double minimumPhaseFraction) {
    if (maximumBisections < 1 || !positive(temperatureToleranceK) || !positive(minimumPhaseFraction)
        || minimumPhaseFraction >= 0.1) {
      throw new IllegalArgumentException("invalid stable-region transition controls");
    }
    this.maximumBisections = maximumBisections;
    this.temperatureToleranceK = temperatureToleranceK;
    this.minimumPhaseFraction = minimumPhaseFraction;
    return this;
  }

  /**
   * Scans ordered temperature probes at every pressure and refines every adjacent stable 2P/GOW transition.
   *
   * <p>
   * A transition is accepted as a bracket only when the two stable regions differ by the requested third phase. A
   * direct jump over an intermediate region is preserved as an unresolved diagnostic and is never rendered as a
   * boundary.
   * </p>
   */
  public Result scan(double[] temperaturesK, double[] pressuresBara) {
    validateGrid(temperaturesK, pressuresBara);
    double[] sortedTemperatures = sortedUnique(temperaturesK);
    double[] sortedPressures = sortedUnique(pressuresBara);
    List<StableState> states = new ArrayList<StableState>();
    List<TransitionBracket> brackets = new ArrayList<TransitionBracket>();
    List<Failure> failures = new ArrayList<Failure>();
    int flashEvaluations = 0;
    for (double pressureBara : sortedPressures) {
      StableState previous = null;
      for (double temperatureK : sortedTemperatures) {
        Evaluation evaluation = evaluate(temperatureK, pressureBara);
        flashEvaluations++;
        if (!evaluation.isSuccessful()) {
          failures.add(new Failure(temperatureK, pressureBara, evaluation.failureMessage));
          previous = null;
          continue;
        }
        StableState current = evaluation.state;
        states.add(current);
        if (previous != null && previous.region != current.region) {
          BoundaryFamily family = familyFor(previous.region, current.region);
          if (family == null) {
            failures.add(new Failure(0.5 * (previous.temperatureK + current.temperatureK), pressureBara,
                "OUT_OF_SCOPE_STABLE_REGION_TRANSITION:" + previous.region + "->" + current.region));
          } else {
            Refined refined = refine(previous, current, family);
            flashEvaluations += refined.flashEvaluations;
            if (refined.bracket != null) {
              addDistinct(brackets, refined.bracket);
            } else {
              failures.add(new Failure(0.5 * (previous.temperatureK + current.temperatureK), pressureBara,
                  refined.failureMessage));
            }
          }
        }
        previous = current;
      }
    }
    brackets.sort(Comparator.comparingDouble(TransitionBracket::getPressureBara)
        .thenComparingDouble(TransitionBracket::getMidpointTemperatureK));
    return new Result(states, brackets, failures, flashEvaluations);
  }

  private Refined refine(StableState first, StableState second, BoundaryFamily family) {
    StableState lowerTemperature = first.temperatureK < second.temperatureK ? first : second;
    StableState upperTemperature = first.temperatureK < second.temperatureK ? second : first;
    Region firstRegion = lowerTemperature.region;
    Region secondRegion = upperTemperature.region;
    int evaluations = 0;
    int bisections = 0;
    while (bisections < maximumBisections
        && upperTemperature.temperatureK - lowerTemperature.temperatureK > temperatureToleranceK) {
      bisections++;
      double midpoint = 0.5 * (lowerTemperature.temperatureK + upperTemperature.temperatureK);
      Evaluation evaluation = evaluate(midpoint, first.pressureBara);
      evaluations++;
      if (!evaluation.isSuccessful()) {
        return Refined.failure(evaluations,
            "TP_FLASH_FAILED_DURING_STABLE_REGION_REFINEMENT:" + evaluation.failureMessage);
      }
      StableState middle = evaluation.state;
      if (middle.region == firstRegion) {
        lowerTemperature = middle;
      } else if (middle.region == secondRegion) {
        upperTemperature = middle;
      } else {
        return Refined.failure(evaluations,
            "INTERMEDIATE_STABLE_REGION_DISCOVERED:" + firstRegion + "->" + middle.region + "->" + secondRegion);
      }
    }
    if (upperTemperature.temperatureK - lowerTemperature.temperatureK > temperatureToleranceK) {
      return Refined.failure(evaluations, "MAXIMUM_STABLE_REGION_BISECTION_COUNT_REACHED");
    }
    StableState twoPhase = lowerTemperature.region == family.getDefinition().getLowerPhaseRegion() ? lowerTemperature
        : upperTemperature;
    StableState threePhase = lowerTemperature.region == Region.GAS_OIL_AQUEOUS ? lowerTemperature : upperTemperature;
    if (twoPhase.region != family.getDefinition().getLowerPhaseRegion()
        || threePhase.region != Region.GAS_OIL_AQUEOUS) {
      return Refined.failure(evaluations, "REFINED_BRACKET_LOST_REQUESTED_TOPOLOGY");
    }
    return Refined.success(
        new TransitionBracket(family, lowerTemperature, upperTemperature, twoPhase, threePhase, bisections),
        evaluations);
  }

  private Evaluation evaluate(double temperatureK, double pressureBara) {
    SystemInterface working = template.clone();
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    working.setMultiPhaseCheck(true);
    working.setMaxNumberOfPhases(3);
    try {
      new ThermodynamicOperations(working).TPflash();
      working.init(1);
      double[][] compositions = new double[CandidatePhase.values().length][];
      double[] fractions = new double[CandidatePhase.values().length];
      Set<Phase> activePhases = EnumSet.noneOf(Phase.class);
      for (int phaseIndex = 0; phaseIndex < working.getNumberOfPhases(); phaseIndex++) {
        double fraction = working.getBeta(phaseIndex);
        if (!Double.isFinite(fraction) || fraction < 0.0) {
          return Evaluation.failure("TP_FLASH_RETURNED_INVALID_PHASE_FRACTION");
        }
        CandidatePhase phase = candidatePhase(working.getPhase(phaseIndex).getType());
        if (fraction <= minimumPhaseFraction) {
          continue;
        }
        if (compositions[phase.ordinal()] != null) {
          return Evaluation.failure("TP_FLASH_RETURNED_DUPLICATE_PHASE_FAMILY:" + phase);
        }
        compositions[phase.ordinal()] = composition(working, phaseIndex);
        fractions[phase.ordinal()] = fraction;
        activePhases.add(topologyPhase(phase));
      }
      Region region = region(activePhases);
      return Evaluation.success(
          new StableState(temperatureK, pressureBara, region, compositions, fractions, working.getNumberOfPhases()));
    } catch (RuntimeException error) {
      return Evaluation.failure(error.getMessage() == null ? error.getClass().getSimpleName() : error.getMessage());
    }
  }

  private static BoundaryFamily familyFor(Region first, Region second) {
    for (BoundaryFamily family : BoundaryFamily.values()) {
      Region lower = family.getDefinition().getLowerPhaseRegion();
      Region higher = family.getDefinition().getHigherPhaseRegion();
      if (first == lower && second == higher || first == higher && second == lower) {
        return family;
      }
    }
    return null;
  }

  private static Region region(Set<Phase> phases) {
    for (Region candidate : Region.values()) {
      if (candidate.getPhases().equals(phases)) {
        return candidate;
      }
    }
    throw new IllegalStateException("unsupported stable phase set " + phases);
  }

  private static CandidatePhase candidatePhase(PhaseType phaseType) {
    if (phaseType == PhaseType.GAS) {
      return CandidatePhase.GAS;
    }
    if (phaseType == PhaseType.AQUEOUS) {
      return CandidatePhase.AQUEOUS;
    }
    return CandidatePhase.OIL;
  }

  private static Phase topologyPhase(CandidatePhase phase) {
    if (phase == CandidatePhase.GAS) {
      return Phase.GAS;
    }
    if (phase == CandidatePhase.AQUEOUS) {
      return Phase.AQUEOUS;
    }
    return Phase.OIL;
  }

  private static double[] composition(SystemInterface system, int phaseIndex) {
    int componentCount = system.getPhase(phaseIndex).getNumberOfComponents();
    double[] composition = new double[componentCount];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = Math.max(0.0, system.getPhase(phaseIndex).getComponent(componentIndex).getx());
      total += composition[componentIndex];
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalStateException("TP flash returned a non-normalizable phase composition");
    }
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }

  private static void addDistinct(List<TransitionBracket> brackets, TransitionBracket candidate) {
    for (TransitionBracket bracket : brackets) {
      if (bracket.family == candidate.family
          && Math.abs(bracket.pressureBara - candidate.pressureBara) <= 1.0e-10 * Math.max(1.0, candidate.pressureBara)
          && Math.abs(bracket.getMidpointTemperatureK() - candidate.getMidpointTemperatureK()) <= 1.0e-4) {
        return;
      }
    }
    brackets.add(candidate);
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

  private static double[] sortedUnique(double[] values) {
    double[] sorted = values.clone();
    java.util.Arrays.sort(sorted);
    List<Double> unique = new ArrayList<Double>();
    for (double value : sorted) {
      if (unique.isEmpty() || Math.abs(value - unique.get(unique.size() - 1)) > 1.0e-12 * Math.max(1.0, value)) {
        unique.add(value);
      }
    }
    double[] result = new double[unique.size()];
    for (int index = 0; index < result.length; index++) {
      result[index] = unique.get(index);
    }
    return result;
  }

  private static void validateGrid(double[] temperaturesK, double[] pressuresBara) {
    if (temperaturesK == null || temperaturesK.length < 2 || pressuresBara == null || pressuresBara.length == 0) {
      throw new IllegalArgumentException("at least two temperatures and one pressure are required");
    }
    for (double temperatureK : temperaturesK) {
      if (!Double.isFinite(temperatureK) || temperatureK < 50.0 || temperatureK > 2500.0) {
        throw new IllegalArgumentException("temperature grid contains a nonphysical value");
      }
    }
    for (double pressureBara : pressuresBara) {
      if (!Double.isFinite(pressureBara) || pressureBara <= 0.0 || pressureBara > 1.0e6) {
        throw new IllegalArgumentException("pressure grid contains a nonphysical value");
      }
    }
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** One stable TP-flash state with phase identities, compositions, and fractions preserved. */
  public static final class StableState {
    private final double temperatureK;
    private final double pressureBara;
    private final Region region;
    private final double[][] phaseCompositions;
    private final double[] phaseFractions;
    private final int rawPhaseCount;

    private StableState(double temperatureK, double pressureBara, Region region, double[][] phaseCompositions,
        double[] phaseFractions, int rawPhaseCount) {
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.region = region;
      this.phaseCompositions = new double[phaseCompositions.length][];
      for (int phaseIndex = 0; phaseIndex < phaseCompositions.length; phaseIndex++) {
        this.phaseCompositions[phaseIndex] = phaseCompositions[phaseIndex] == null ? new double[0]
            : phaseCompositions[phaseIndex].clone();
      }
      this.phaseFractions = phaseFractions.clone();
      this.rawPhaseCount = rawPhaseCount;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public Region getRegion() {
      return region;
    }

    public boolean hasPhase(CandidatePhase phase) {
      return phase != null && phaseCompositions[phase.ordinal()].length > 0;
    }

    public double[] getPhaseComposition(CandidatePhase phase) {
      if (phase == null) {
        throw new IllegalArgumentException("phase is required");
      }
      return phaseCompositions[phase.ordinal()].clone();
    }

    public double getPhaseFraction(CandidatePhase phase) {
      if (phase == null) {
        throw new IllegalArgumentException("phase is required");
      }
      return phaseFractions[phase.ordinal()];
    }

    public int getRawPhaseCount() {
      return rawPhaseCount;
    }
  }

  /** One stable adjacent-region bracket suitable for a subsequent strict zero-TPD correction. */
  public static final class TransitionBracket {
    private final BoundaryFamily family;
    private final StableState lowerTemperatureState;
    private final StableState upperTemperatureState;
    private final StableState twoPhaseState;
    private final StableState threePhaseState;
    private final double pressureBara;
    private final int bisections;

    private TransitionBracket(BoundaryFamily family, StableState lowerTemperatureState,
        StableState upperTemperatureState, StableState twoPhaseState, StableState threePhaseState, int bisections) {
      this.family = family;
      this.lowerTemperatureState = lowerTemperatureState;
      this.upperTemperatureState = upperTemperatureState;
      this.twoPhaseState = twoPhaseState;
      this.threePhaseState = threePhaseState;
      this.pressureBara = lowerTemperatureState.pressureBara;
      this.bisections = bisections;
    }

    public BoundaryFamily getFamily() {
      return family;
    }

    public StableState getLowerTemperatureState() {
      return lowerTemperatureState;
    }

    public StableState getUpperTemperatureState() {
      return upperTemperatureState;
    }

    public StableState getTwoPhaseState() {
      return twoPhaseState;
    }

    public StableState getThreePhaseState() {
      return threePhaseState;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double getMidpointTemperatureK() {
      return 0.5 * (lowerTemperatureState.temperatureK + upperTemperatureState.temperatureK);
    }

    public double getTemperatureWidthK() {
      return upperTemperatureState.temperatureK - lowerTemperatureState.temperatureK;
    }

    public int getBisections() {
      return bisections;
    }

    /** Fraction of the appearing phase on the GOW side of the refined phase-count switch. */
    public double getIncipientPhaseFractionOnThreePhaseSide() {
      return threePhaseState.getPhaseFraction(family.getIncipientPhase());
    }

    /** Normalized retained-phase-zero split on the two-phase side. */
    public double getRetainedPhaseZeroSplitOnTwoPhaseSide() {
      return normalizedRetainedPhaseZeroSplit(twoPhaseState);
    }

    /** Normalized retained-phase-zero split after removing the incipient phase on the GOW side. */
    public double getRetainedPhaseZeroSplitOnThreePhaseSide() {
      return normalizedRetainedPhaseZeroSplit(threePhaseState);
    }

    public double getRetainedPhaseSplitJump() {
      return Math.abs(getRetainedPhaseZeroSplitOnTwoPhaseSide() - getRetainedPhaseZeroSplitOnThreePhaseSide());
    }

    public double getRetainedPhaseZeroCompositionJump() {
      return compositionDistance(twoPhaseState.getPhaseComposition(family.getRetainedPhaseZero()),
          threePhaseState.getPhaseComposition(family.getRetainedPhaseZero()));
    }

    public double getRetainedPhaseOneCompositionJump() {
      return compositionDistance(twoPhaseState.getPhaseComposition(family.getRetainedPhaseOne()),
          threePhaseState.getPhaseComposition(family.getRetainedPhaseOne()));
    }

    public double getMaximumRetainedPhaseCompositionJump() {
      return Math.max(getRetainedPhaseZeroCompositionJump(), getRetainedPhaseOneCompositionJump());
    }

    private double normalizedRetainedPhaseZeroSplit(StableState state) {
      double phaseZero = state.getPhaseFraction(family.getRetainedPhaseZero());
      double retainedTotal = phaseZero + state.getPhaseFraction(family.getRetainedPhaseOne());
      return phaseZero / retainedTotal;
    }
  }

  /** Explicit TP-flash or non-adjacent-region diagnostic. */
  public static final class Failure {
    private final double temperatureK;
    private final double pressureBara;
    private final String reason;

    private Failure(double temperatureK, double pressureBara, String reason) {
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.reason = reason;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public String getReason() {
      return reason;
    }
  }

  /** Immutable stable-region inventory retaining every grid state, bracket, and failure. */
  public static final class Result {
    private final List<StableState> states;
    private final List<TransitionBracket> brackets;
    private final List<Failure> failures;
    private final int flashEvaluations;

    private Result(List<StableState> states, List<TransitionBracket> brackets, List<Failure> failures,
        int flashEvaluations) {
      this.states = Collections.unmodifiableList(new ArrayList<StableState>(states));
      this.brackets = Collections.unmodifiableList(new ArrayList<TransitionBracket>(brackets));
      this.failures = Collections.unmodifiableList(new ArrayList<Failure>(failures));
      this.flashEvaluations = flashEvaluations;
    }

    public List<StableState> getStates() {
      return states;
    }

    public List<TransitionBracket> getBrackets() {
      return brackets;
    }

    public List<Failure> getFailures() {
      return failures;
    }

    public int getFlashEvaluations() {
      return flashEvaluations;
    }
  }

  private static final class Evaluation {
    private final StableState state;
    private final String failureMessage;

    private Evaluation(StableState state, String failureMessage) {
      this.state = state;
      this.failureMessage = failureMessage;
    }

    private static Evaluation success(StableState state) {
      return new Evaluation(state, null);
    }

    private static Evaluation failure(String message) {
      return new Evaluation(null, message);
    }

    private boolean isSuccessful() {
      return state != null && failureMessage == null;
    }
  }

  private static final class Refined {
    private final TransitionBracket bracket;
    private final int flashEvaluations;
    private final String failureMessage;

    private Refined(TransitionBracket bracket, int flashEvaluations, String failureMessage) {
      this.bracket = bracket;
      this.flashEvaluations = flashEvaluations;
      this.failureMessage = failureMessage;
    }

    private static Refined success(TransitionBracket bracket, int evaluations) {
      return new Refined(bracket, evaluations, null);
    }

    private static Refined failure(int evaluations, String message) {
      return new Refined(null, evaluations, message);
    }
  }
}
