package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Traces a two-hydrocarbon-phase PT boundary with the Michelsen predictor-corrector continuation solver.
 *
 * <p>
 * The supplied thermodynamic system is never mutated. The solver works on a private clone, traces the incipient
 * two-phase boundary, preserves disjoint branch segments, classifies the physical dew side from the cricondentherm, and
 * optionally checks the traced boundary for multiphase stability. If native continuation cannot initialize, an optional
 * stable-phase topology fallback uses marching squares and edge refinement. The fallback preserves every disconnected
 * contour and never joins boundary points merely because they share a temperature row.
 * </p>
 */
public final class TwoHydrocarbonPhaseEnvelopeSolver {
  private static final double MINIMUM_PHASE_FRACTION = 1.0e-10;
  /** Factor the pressure ceiling grows by on each extension attempt. */
  private static final double PRESSURE_EXTENSION_FACTOR = 4.0;
  /** Hard bound on how many times the pressure ceiling may be raised for one start side. */
  private static final int MAXIMUM_PRESSURE_EXTENSIONS = 4;

  private final SystemInterface sourceSystem;
  private double pressureCeilingExtensionLimitBara = 0.0;
  private double minimumPressureBara = 0.05;
  private double maximumPressureBara = 10000.0;
  private double maximumTemperatureStepK = 5.0;
  private double maximumPressureStepBara = 5.0;
  private double minimumTemperatureK = 273.15 - 60.0;
  private double maximumTemperatureK = 273.15 + 500.0;
  private int topologyTemperaturePoints = 41;
  private int topologyPressurePoints = 25;
  private int topologyBisectionIterations = 12;
  private int maximumContinuationIterations = 1200;
  private int threePhaseRefinementIterations = 0;
  private double threePhaseTemperatureToleranceK = 0.02;
  private double threePhasePressureToleranceBara = 0.02;
  private boolean stabilityAnalysisEnabled = true;
  private boolean topologyFallbackEnabled = true;

  /**
   * Creates a continuation solver for a thermodynamic system.
   *
   * @param sourceSystem fully configured thermodynamic system
   */
  public TwoHydrocarbonPhaseEnvelopeSolver(SystemInterface sourceSystem) {
    if (sourceSystem == null) {
      throw new IllegalArgumentException("sourceSystem must not be null");
    }
    this.sourceSystem = sourceSystem;
  }

  /**
   * Sets the pressure interval used by the continuation solver.
   *
   * @param minimumPressureBara lower pressure and initial tracing pressure in bara
   * @param maximumPressureBara upper pressure termination limit in bara
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setPressureRange(double minimumPressureBara, double maximumPressureBara) {
    if (!Double.isFinite(minimumPressureBara) || !Double.isFinite(maximumPressureBara) || minimumPressureBara <= 0.0
        || maximumPressureBara <= minimumPressureBara) {
      throw new IllegalArgumentException("pressure range must be finite, positive, and increasing");
    }
    this.minimumPressureBara = minimumPressureBara;
    this.maximumPressureBara = maximumPressureBara;
    return this;
  }

  /**
   * Sets the maximum continuation steps.
   *
   * @param maximumTemperatureStepK maximum temperature step in kelvin
   * @param maximumPressureStepBara maximum pressure step in bara
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setMaximumSteps(double maximumTemperatureStepK,
      double maximumPressureStepBara) {
    if (!Double.isFinite(maximumTemperatureStepK) || !Double.isFinite(maximumPressureStepBara)
        || maximumTemperatureStepK <= 0.0 || maximumPressureStepBara <= 0.0) {
      throw new IllegalArgumentException("continuation steps must be finite and positive");
    }
    this.maximumTemperatureStepK = maximumTemperatureStepK;
    this.maximumPressureStepBara = maximumPressureStepBara;
    return this;
  }

  /**
   * Sets the temperature interval used by the stability-contour fallback.
   *
   * @param minimumTemperatureK minimum temperature in kelvin
   * @param maximumTemperatureK maximum temperature in kelvin
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setTemperatureRange(double minimumTemperatureK, double maximumTemperatureK) {
    if (!Double.isFinite(minimumTemperatureK) || !Double.isFinite(maximumTemperatureK) || minimumTemperatureK < 50.0
        || maximumTemperatureK <= minimumTemperatureK) {
      throw new IllegalArgumentException("temperature range must be finite, physical, and increasing");
    }
    this.minimumTemperatureK = minimumTemperatureK;
    this.maximumTemperatureK = maximumTemperatureK;
    return this;
  }

  /**
   * Sets the seed grid used only when direct Michelsen continuation cannot initialize.
   *
   * @param temperaturePoints number of temperature nodes
   * @param pressurePoints number of logarithmic-pressure nodes
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setTopologySeedGrid(int temperaturePoints, int pressurePoints) {
    if (temperaturePoints < 11 || pressurePoints < 9) {
      throw new IllegalArgumentException("topology seed grid must contain at least 11 by 9 nodes");
    }
    this.topologyTemperaturePoints = temperaturePoints;
    this.topologyPressurePoints = pressurePoints;
    return this;
  }

  /**
   * Enables or disables post-trace multiphase stability checks.
   *
   * @param enabled true to check sampled boundary points with multiphase TP flashes
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setStabilityAnalysisEnabled(boolean enabled) {
    this.stabilityAnalysisEnabled = enabled;
    return this;
  }

  /**
   * Enables fixed-boundary refinement of every detected two-phase/three-phase transition bracket.
   *
   * @param maximumBisections maximum bisections per coarse bracket; zero disables refinement
   * @param temperatureToleranceK terminal bracket width in kelvin
   * @param pressureToleranceBara terminal bracket width in bara
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setThreePhaseRefinement(int maximumBisections, double temperatureToleranceK,
      double pressureToleranceBara) {
    if (maximumBisections < 0 || !Double.isFinite(temperatureToleranceK) || temperatureToleranceK <= 0.0
        || !Double.isFinite(pressureToleranceBara) || pressureToleranceBara <= 0.0) {
      throw new IllegalArgumentException("invalid three-phase refinement configuration");
    }
    this.threePhaseRefinementIterations = maximumBisections;
    this.threePhaseTemperatureToleranceK = temperatureToleranceK;
    this.threePhasePressureToleranceBara = pressureToleranceBara;
    return this;
  }

  /**
   * Sets the hard iteration limit for each Michelsen continuation pass.
   *
   * @param maximumIterations maximum iterations per pass
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setMaximumContinuationIterations(int maximumIterations) {
    if (maximumIterations < 10) {
      throw new IllegalArgumentException("maximumIterations must be at least 10");
    }
    this.maximumContinuationIterations = maximumIterations;
    return this;
  }

  /**
   * Enables or disables the stable-phase topology fallback.
   *
   * @param enabled true to run the fallback when continuation is incomplete
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setTopologyFallbackEnabled(boolean enabled) {
    this.topologyFallbackEnabled = enabled;
    return this;
  }

  /**
   * Traces the phase boundary. A dew-side start is attempted first; a bubble-side start is used only when it produces a
   * more complete physical envelope.
   *
   * @return immutable result with contiguous segments and solver diagnostics
   */
  public Result solve() {
    Attempt dewFirst = trace(false);
    if (dewFirst.isEngineeringComplete()) {
      return dewFirst.toResult();
    }
    Attempt bubbleFirst = trace(true);
    if (bubbleFirst.isEngineeringComplete()) {
      return bubbleFirst.toResult();
    }
    if (!topologyFallbackEnabled) {
      return bubbleFirst.score() > dewFirst.score() ? bubbleFirst.toResult() : dewFirst.toResult();
    }
    TopologyAttempt topology = traceStablePhaseTopology();
    if (topology.isEngineeringComplete()) {
      return topology.toResult();
    }
    Attempt bestContinuation = bubbleFirst.score() > dewFirst.score() ? bubbleFirst : dewFirst;
    return topology.score() > bestContinuation.score() ? topology.toResult() : bestContinuation.toResult();
  }

  /**
   * Traces one start side, raising the pressure ceiling while the boundary is still climbing when it is reached.
   *
   * <p>
   * A caller's pressure range is a plotting window, not a statement about where the boundary ends. Stopping at that
   * window produces a branch that looks like a traced envelope while covering only the part of it below the window -
   * and, because points then exist on both sides, it can pass a naive closure check. Extending the ceiling until the
   * continuation stops on its own keeps the traced object the whole envelope.
   * </p>
   *
   * @param bubblePointFirst true to start from the bubble side
   * @return the best attempt reachable within the extension limit
   */
  /**
   * Allows the traced pressure ceiling to be raised above the configured range when the boundary is still climbing.
   *
   * @param limitBara hard upper bound the ceiling may be raised to; zero or negative disables extension
   * @return this solver
   */
  public TwoHydrocarbonPhaseEnvelopeSolver setPressureCeilingExtensionLimit(double limitBara) {
    this.pressureCeilingExtensionLimitBara = limitBara;
    return this;
  }

  private Attempt trace(boolean bubblePointFirst) {
    double ceiling = maximumPressureBara;
    Attempt attempt = traceAtCeiling(bubblePointFirst, ceiling);
    for (int extension = 0; extension < MAXIMUM_PRESSURE_EXTENSIONS && attempt.stoppedAtPressureCeiling()
        && ceiling < pressureCeilingExtensionLimitBara; extension++) {
      ceiling = Math.min(ceiling * PRESSURE_EXTENSION_FACTOR, pressureCeilingExtensionLimitBara);
      Attempt extended = traceAtCeiling(bubblePointFirst, ceiling);
      if (extended.pointCount() < attempt.pointCount()) {
        break;
      }
      attempt = extended;
    }
    return attempt;
  }

  private Attempt traceAtCeiling(boolean bubblePointFirst, double ceilingBara) {
    SystemInterface workingSystem = sourceSystem.clone();
    workingSystem.setMultiPhaseCheck(false);
    workingSystem.setMaxNumberOfPhases(2);
    double phaseFraction = bubblePointFirst ? MINIMUM_PHASE_FRACTION : 1.0 - MINIMUM_PHASE_FRACTION;
    PTPhaseEnvelopeMichelsen operation = new PTPhaseEnvelopeMichelsen(workingSystem, "", phaseFraction,
        minimumPressureBara, bubblePointFirst);
    operation.setMinPressure(minimumPressureBara);
    operation.setMaxPressure(ceilingBara);
    operation.setDTmax(maximumTemperatureStepK);
    operation.setDPmax(maximumPressureStepBara);
    operation.setMaximumEnvelopeIterations(maximumContinuationIterations);
    try {
      operation.run();
      if (stabilityAnalysisEnabled && !operation.getSegments().isEmpty()) {
        operation.checkStabilityAlongEnvelope();
        if (threePhaseRefinementIterations > 0) {
          operation.refineSecondaryStabilityBrackets(threePhaseRefinementIterations, threePhaseTemperatureToleranceK,
              threePhasePressureToleranceBara);
        }
      }
      return Attempt.success(operation, bubblePointFirst);
    } catch (Exception error) {
      return Attempt.failure(operation, bubblePointFirst, error.getMessage());
    }
  }

  private TopologyAttempt traceStablePhaseTopology() {
    return new PhaseTopologyTracer(sourceSystem, minimumTemperatureK, maximumTemperatureK, minimumPressureBara,
        maximumPressureBara, topologyTemperaturePoints, topologyPressurePoints, topologyBisectionIterations).trace();
  }

  /** Immutable continuation result. */
  public static final class Result {
    private final List<Segment> segments;
    private final boolean envelopeClosed;
    private final boolean finitePhysicalValues;
    private final int criticalPointCount;
    private final double[] criticalPoint;
    private final double[] cricondenbar;
    private final double[] cricondentherm;
    private final int threePhaseStabilityPointCount;
    private final boolean bubblePointFirst;
    private final String failureMessage;
    private final String method;
    private final int flashCount;
    private final int failedFlashCount;
    private final boolean iterationLimitReached;
    private final List<PTPhaseEnvelopeMichelsen.SecondaryStabilitySample> secondaryStabilitySamples;
    private final List<PTPhaseEnvelopeMichelsen.SecondaryStabilityBracket> secondaryStabilityBrackets;
    private final List<PTPhaseEnvelopeMichelsen.ThreePhasePointCandidate> threePhasePointCandidates;
    private final String terminationSummary;

    private Result(List<Segment> segments, boolean envelopeClosed, boolean finitePhysicalValues, int criticalPointCount,
        double[] criticalPoint, double[] cricondenbar, double[] cricondentherm, int threePhaseStabilityPointCount,
        boolean bubblePointFirst, String failureMessage, String method, int flashCount, int failedFlashCount) {
      this(segments, envelopeClosed, finitePhysicalValues, criticalPointCount, criticalPoint, cricondenbar,
          cricondentherm, threePhaseStabilityPointCount, bubblePointFirst, failureMessage, method, flashCount,
          failedFlashCount, false, Collections.<PTPhaseEnvelopeMichelsen.SecondaryStabilitySample>emptyList(),
          Collections.<PTPhaseEnvelopeMichelsen.SecondaryStabilityBracket>emptyList(),
          Collections.<PTPhaseEnvelopeMichelsen.ThreePhasePointCandidate>emptyList(), "NOT_APPLICABLE");
    }

    private Result(List<Segment> segments, boolean envelopeClosed, boolean finitePhysicalValues, int criticalPointCount,
        double[] criticalPoint, double[] cricondenbar, double[] cricondentherm, int threePhaseStabilityPointCount,
        boolean bubblePointFirst, String failureMessage, String method, int flashCount, int failedFlashCount,
        boolean iterationLimitReached) {
      this(segments, envelopeClosed, finitePhysicalValues, criticalPointCount, criticalPoint, cricondenbar,
          cricondentherm, threePhaseStabilityPointCount, bubblePointFirst, failureMessage, method, flashCount,
          failedFlashCount, iterationLimitReached,
          Collections.<PTPhaseEnvelopeMichelsen.SecondaryStabilitySample>emptyList(),
          Collections.<PTPhaseEnvelopeMichelsen.SecondaryStabilityBracket>emptyList(),
          Collections.<PTPhaseEnvelopeMichelsen.ThreePhasePointCandidate>emptyList(), "NOT_APPLICABLE");
    }

    private Result(List<Segment> segments, boolean envelopeClosed, boolean finitePhysicalValues, int criticalPointCount,
        double[] criticalPoint, double[] cricondenbar, double[] cricondentherm, int threePhaseStabilityPointCount,
        boolean bubblePointFirst, String failureMessage, String method, int flashCount, int failedFlashCount,
        boolean iterationLimitReached,
        List<PTPhaseEnvelopeMichelsen.SecondaryStabilitySample> secondaryStabilitySamples,
        List<PTPhaseEnvelopeMichelsen.SecondaryStabilityBracket> secondaryStabilityBrackets,
        List<PTPhaseEnvelopeMichelsen.ThreePhasePointCandidate> threePhasePointCandidates, String terminationSummary) {
      this.terminationSummary = terminationSummary;
      this.segments = Collections.unmodifiableList(new ArrayList<Segment>(segments));
      this.envelopeClosed = envelopeClosed;
      this.finitePhysicalValues = finitePhysicalValues;
      this.criticalPointCount = criticalPointCount;
      this.criticalPoint = criticalPoint.clone();
      this.cricondenbar = cricondenbar.clone();
      this.cricondentherm = cricondentherm.clone();
      this.threePhaseStabilityPointCount = threePhaseStabilityPointCount;
      this.bubblePointFirst = bubblePointFirst;
      this.failureMessage = failureMessage;
      this.method = method;
      this.flashCount = flashCount;
      this.failedFlashCount = failedFlashCount;
      this.iterationLimitReached = iterationLimitReached;
      this.secondaryStabilitySamples = Collections.unmodifiableList(
          new ArrayList<PTPhaseEnvelopeMichelsen.SecondaryStabilitySample>(secondaryStabilitySamples));
      this.secondaryStabilityBrackets = Collections.unmodifiableList(
          new ArrayList<PTPhaseEnvelopeMichelsen.SecondaryStabilityBracket>(secondaryStabilityBrackets));
      this.threePhasePointCandidates = Collections.unmodifiableList(
          new ArrayList<PTPhaseEnvelopeMichelsen.ThreePhasePointCandidate>(threePhasePointCandidates));
    }

    /** @return immutable contiguous envelope segments */
    public List<Segment> getSegments() {
      return segments;
    }

    /**
     * Getter for why each continuation pass stopped, joined with {@code +} in pass order.
     *
     * <p>
     * Anything other than {@code PRESSURE_FLOOR} or {@code CONTINUATION_END} means the traced curve is a truncated
     * piece of the real boundary.
     * </p>
     *
     * @return termination summary, never null
     */
    public String getTerminationSummary() {
      return terminationSummary;
    }

    /** @return total number of finite boundary points */
    public int getPointCount() {
      int count = 0;
      for (Segment segment : segments) {
        count += segment.size();
      }
      return count;
    }

    /** @return true when both physical sides contain sufficient points */
    public boolean isEnvelopeClosed() {
      return envelopeClosed;
    }

    /** @return true when all exported values are finite and physically bounded */
    public boolean hasFinitePhysicalValues() {
      return finitePhysicalValues;
    }

    /** @return number of detected critical points */
    public int getCriticalPointCount() {
      return criticalPointCount;
    }

    /** @return first critical point as {@code [temperature_K, pressure_bara]} */
    public double[] getCriticalPoint() {
      return criticalPoint.clone();
    }

    /** @return cricondenbar as {@code [temperature_K, pressure_bara, 0]} */
    public double[] getCricondenbar() {
      return cricondenbar.clone();
    }

    /** @return cricondentherm as {@code [temperature_K, pressure_bara, 0]} */
    public double[] getCricondentherm() {
      return cricondentherm.clone();
    }

    /** @return sampled boundary points where a multiphase TP flash found three or more phases */
    public int getThreePhaseStabilityPointCount() {
      return threePhaseStabilityPointCount;
    }

    /** @return tangent-plane stability diagnostics sampled along the native two-phase boundary */
    public List<PTPhaseEnvelopeMichelsen.SecondaryStabilitySample> getSecondaryStabilitySamples() {
      return secondaryStabilitySamples;
    }

    /** @return coarse brackets around two-phase/three-phase topology transitions */
    public List<PTPhaseEnvelopeMichelsen.SecondaryStabilityBracket> getSecondaryStabilityBrackets() {
      return secondaryStabilityBrackets;
    }

    /** @return refined two-phase/three-phase transition candidates */
    public List<PTPhaseEnvelopeMichelsen.ThreePhasePointCandidate> getThreePhasePointCandidates() {
      return threePhasePointCandidates;
    }

    /** @return true when the selected attempt started from the bubble side */
    public boolean isBubblePointFirst() {
      return bubblePointFirst;
    }

    /** @return algorithm that produced the exported boundary */
    public String getMethod() {
      return method;
    }

    /** @return number of stability TP flashes used by the topology fallback */
    public int getFlashCount() {
      return flashCount;
    }

    /** @return failed stability TP flashes */
    public int getFailedFlashCount() {
      return failedFlashCount;
    }

    /** @return true when native continuation stopped at its configured hard iteration limit */
    public boolean isIterationLimitReached() {
      return iterationLimitReached;
    }

    /** @return failure message, or {@code null} when tracing did not throw */
    public String getFailureMessage() {
      return failureMessage;
    }
  }

  /** One contiguous, physically classified envelope branch. */
  public static final class Segment {
    /** Physical branch type. */
    public enum PhysicalType {
      /** Incipient liquid from the vapor side. */
      DEW,
      /** Incipient vapor from the liquid side. */
      BUBBLE,
      /** Stable gas-and-hydrocarbon-liquid coexistence boundary. */
      TWO_HC
    }

    private final PhysicalType physicalType;
    private final double[] temperaturesK;
    private final double[] pressuresBara;

    private Segment(PhysicalType physicalType, double[] temperaturesK, double[] pressuresBara) {
      this.physicalType = physicalType;
      this.temperaturesK = temperaturesK.clone();
      this.pressuresBara = pressuresBara.clone();
    }

    /** @return physical dew or bubble classification */
    public PhysicalType getPhysicalType() {
      return physicalType;
    }

    /** @return temperatures in kelvin */
    public double[] getTemperaturesK() {
      return temperaturesK.clone();
    }

    /** @return pressures in bara */
    public double[] getPressuresBara() {
      return pressuresBara.clone();
    }

    /** @return number of points */
    public int size() {
      return temperaturesK.length;
    }
  }

  private static final class Attempt {
    private final PTPhaseEnvelopeMichelsen operation;
    private final boolean bubblePointFirst;
    private final String failureMessage;
    private final List<Segment> segments;
    private final boolean finitePhysicalValues;

    private Attempt(PTPhaseEnvelopeMichelsen operation, boolean bubblePointFirst, String failureMessage,
        List<Segment> segments, boolean finitePhysicalValues) {
      this.operation = operation;
      this.bubblePointFirst = bubblePointFirst;
      this.failureMessage = failureMessage;
      this.segments = segments;
      this.finitePhysicalValues = finitePhysicalValues;
    }

    private static Attempt success(PTPhaseEnvelopeMichelsen operation, boolean bubblePointFirst) {
      List<EnvelopeSegment> rawSegments = operation.getSegments();
      EnvelopeSegment.PhaseType physicalDewStoredType = storedTypeAtMaximumTemperature(rawSegments);
      List<Segment> output = new ArrayList<Segment>();
      boolean physical = !rawSegments.isEmpty();
      for (EnvelopeSegment raw : rawSegments) {
        double[] temperatures = raw.getTemperatures();
        double[] pressures = raw.getPressures();
        physical = physical && finitePhysicalValues(temperatures, pressures);
        Segment.PhysicalType type = raw.getPhaseType() == physicalDewStoredType ? Segment.PhysicalType.DEW
            : Segment.PhysicalType.BUBBLE;
        output.add(new Segment(type, temperatures, pressures));
      }
      return new Attempt(operation, bubblePointFirst, null, output, physical);
    }

    private static Attempt failure(PTPhaseEnvelopeMichelsen operation, boolean bubblePointFirst,
        String failureMessage) {
      return new Attempt(operation, bubblePointFirst, failureMessage, Collections.<Segment>emptyList(), false);
    }

    private int pointCount() {
      int count = 0;
      for (Segment segment : segments) {
        count += segment.size();
      }
      return count;
    }

    private boolean stoppedAtPressureCeiling() {
      return failureMessage == null
          && operation.getTerminationReasons().contains(PTPhaseEnvelopeMichelsen.TerminationReason.PRESSURE_CEILING);
    }

    private boolean hasBothPhysicalSides() {
      int dewPoints = 0;
      int bubblePoints = 0;
      for (Segment segment : segments) {
        if (segment.getPhysicalType() == Segment.PhysicalType.DEW) {
          dewPoints += segment.size();
        } else {
          bubblePoints += segment.size();
        }
      }
      return dewPoints >= 3 && bubblePoints >= 3;
    }

    private boolean isEngineeringComplete() {
      return failureMessage == null && finitePhysicalValues && hasBothPhysicalSides() && pointCount() >= 16
          && operation.getNumberOfCriticalPoints() >= 1 && !operation.isTruncated();
    }

    private int score() {
      int value = pointCount();
      if (failureMessage == null) {
        value += 100;
      }
      if (finitePhysicalValues) {
        value += 100;
      }
      if (hasBothPhysicalSides()) {
        value += 200;
      }
      value += operation.getNumberOfCriticalPoints() * 100;
      if (operation.isTruncated()) {
        // A branch cut short at the pressure ceiling can accumulate many points while covering only a fraction of the
        // real boundary. Without this penalty it outscores a shorter attempt that actually closed.
        value -= 300;
      }
      return value;
    }

    private Result toResult() {
      double[] critical = operation.get("criticalPoint1");
      if (critical == null || critical.length < 2) {
        critical = new double[] {0.0, 0.0};
      }
      double[] threePhaseTemperatures = operation.get("threePhaseT");
      return new Result(segments, operation.isEnvelopeClosed(), finitePhysicalValues,
          operation.getNumberOfCriticalPoints(), critical, safeArray(operation.getCricondenBar()),
          safeArray(operation.getCricondenTherm()), threePhaseTemperatures == null ? 0 : threePhaseTemperatures.length,
          bubblePointFirst, failureMessage, "michelsen-predictor-corrector", 0, 0, operation.isIterationLimitReached(),
          operation.getSecondaryStabilitySamples(), operation.getSecondaryStabilityBrackets(),
          operation.getThreePhasePointCandidates(), terminationSummary(operation));
    }

    private static String terminationSummary(PTPhaseEnvelopeMichelsen operation) {
      List<PTPhaseEnvelopeMichelsen.TerminationReason> reasons = operation.getTerminationReasons();
      if (reasons.isEmpty()) {
        return "NOT_STARTED";
      }
      StringBuilder text = new StringBuilder();
      for (PTPhaseEnvelopeMichelsen.TerminationReason reason : reasons) {
        if (text.length() > 0) {
          text.append('+');
        }
        text.append(reason);
      }
      return text.toString();
    }
  }

  private static final class TopologyAttempt {
    private final List<Segment> segments;
    private final boolean closed;
    private final boolean finitePhysicalValues;
    private final double[] criticalPoint;
    private final double[] cricondenbar;
    private final double[] cricondentherm;
    private final int threePhasePointCount;
    private final int flashCount;
    private final int failedFlashCount;
    private final String failureMessage;

    private TopologyAttempt(List<Segment> segments, boolean closed, boolean finitePhysicalValues,
        double[] criticalPoint, double[] cricondenbar, double[] cricondentherm, int threePhasePointCount,
        int flashCount, int failedFlashCount, String failureMessage) {
      this.segments = segments;
      this.closed = closed;
      this.finitePhysicalValues = finitePhysicalValues;
      this.criticalPoint = criticalPoint;
      this.cricondenbar = cricondenbar;
      this.cricondentherm = cricondentherm;
      this.threePhasePointCount = threePhasePointCount;
      this.flashCount = flashCount;
      this.failedFlashCount = failedFlashCount;
      this.failureMessage = failureMessage;
    }

    private int pointCount() {
      int count = 0;
      for (Segment segment : segments) {
        count += segment.size();
      }
      return count;
    }

    private boolean isEngineeringComplete() {
      return failureMessage == null && finitePhysicalValues && pointCount() >= 16
          && failedFlashCount <= Math.max(3, flashCount / 100);
    }

    private int score() {
      int value = pointCount();
      if (failureMessage == null) {
        value += 100;
      }
      if (finitePhysicalValues) {
        value += 100;
      }
      if (closed) {
        value += 200;
      }
      if (criticalPoint[0] > 0.0) {
        value += 100;
      }
      return value;
    }

    private Result toResult() {
      return new Result(segments, closed, finitePhysicalValues, criticalPoint[0] > 0.0 ? 1 : 0, criticalPoint,
          cricondenbar, cricondentherm, threePhasePointCount, false, failureMessage,
          "multiphase-stability-marching-squares", flashCount, failedFlashCount);
    }
  }

  private static final class PhaseTopologyTracer {
    private final SystemInterface template;
    private final double minimumTemperatureK;
    private final double maximumTemperatureK;
    private final double minimumPressureBara;
    private final double maximumPressureBara;
    private final int temperaturePointCount;
    private final int pressurePointCount;
    private final int bisectionIterations;
    private int flashCount;
    private int failedFlashCount;

    private PhaseTopologyTracer(SystemInterface source, double minimumTemperatureK, double maximumTemperatureK,
        double minimumPressureBara, double maximumPressureBara, int temperaturePointCount, int pressurePointCount,
        int bisectionIterations) {
      this.template = source.clone();
      this.template.setMultiPhaseCheck(true);
      this.template.setMaxNumberOfPhases(3);
      this.minimumTemperatureK = minimumTemperatureK;
      this.maximumTemperatureK = maximumTemperatureK;
      this.minimumPressureBara = minimumPressureBara;
      this.maximumPressureBara = maximumPressureBara;
      this.temperaturePointCount = temperaturePointCount;
      this.pressurePointCount = pressurePointCount;
      this.bisectionIterations = bisectionIterations;
    }

    private TopologyAttempt trace() {
      double[] temperatures = linearValues(minimumTemperatureK, maximumTemperatureK, temperaturePointCount);
      double[] pressures = logarithmicValues(minimumPressureBara, maximumPressureBara, pressurePointCount);
      PhaseState[][] states = new PhaseState[temperaturePointCount][pressurePointCount];
      for (int temperatureIndex = 0; temperatureIndex < temperaturePointCount; temperatureIndex++) {
        for (int pressureIndex = 0; pressureIndex < pressurePointCount; pressureIndex++) {
          states[temperatureIndex][pressureIndex] = evaluate(temperatures[temperatureIndex], pressures[pressureIndex]);
        }
      }

      Map<String, ContourPoint> points = new LinkedHashMap<String, ContourPoint>();
      Map<String, Set<String>> graph = new LinkedHashMap<String, Set<String>>();
      for (int temperatureIndex = 0; temperatureIndex < temperaturePointCount - 1; temperatureIndex++) {
        for (int pressureIndex = 0; pressureIndex < pressurePointCount - 1; pressureIndex++) {
          addCellContours(states, temperatureIndex, pressureIndex, points, graph);
        }
      }

      List<Polyline> polylines = tracePolylines(points, graph);
      List<Segment> segments = new ArrayList<Segment>();
      boolean closed = false;
      boolean finite = !polylines.isEmpty();
      double maximumPressure = Double.NEGATIVE_INFINITY;
      double temperatureAtMaximumPressure = 0.0;
      double maximumTemperature = Double.NEGATIVE_INFINITY;
      double pressureAtMaximumTemperature = 0.0;
      double minimumCompositionDistance = Double.POSITIVE_INFINITY;
      double[] criticalPoint = new double[] {0.0, 0.0};
      int threePhasePointCount = 0;
      for (Polyline polyline : polylines) {
        if (polyline.points.size() < 3) {
          continue;
        }
        double[] segmentTemperatures = new double[polyline.points.size()];
        double[] segmentPressures = new double[polyline.points.size()];
        for (int index = 0; index < polyline.points.size(); index++) {
          ContourPoint point = polyline.points.get(index);
          segmentTemperatures[index] = point.temperatureK;
          segmentPressures[index] = point.pressureBara;
          if (point.pressureBara > maximumPressure) {
            maximumPressure = point.pressureBara;
            temperatureAtMaximumPressure = point.temperatureK;
          }
          if (point.temperatureK > maximumTemperature) {
            maximumTemperature = point.temperatureK;
            pressureAtMaximumTemperature = point.pressureBara;
          }
          if (point.compositionDistance < minimumCompositionDistance) {
            minimumCompositionDistance = point.compositionDistance;
            criticalPoint = new double[] {point.temperatureK, point.pressureBara};
          }
          if (point.threePhase) {
            threePhasePointCount++;
          }
        }
        finite = finite && finitePhysicalValues(segmentTemperatures, segmentPressures);
        closed = closed || polyline.closed;
        segments.add(new Segment(Segment.PhysicalType.TWO_HC, segmentTemperatures, segmentPressures));
      }

      if (segments.isEmpty()) {
        return new TopologyAttempt(segments, false, false, new double[] {0.0, 0.0}, new double[] {0.0, 0.0, 0.0},
            new double[] {0.0, 0.0, 0.0}, 0, flashCount, failedFlashCount,
            "No stable gas-and-hydrocarbon-liquid coexistence boundary was found in the automatic domain");
      }
      if (!Double.isFinite(minimumCompositionDistance)) {
        criticalPoint = new double[] {0.0, 0.0};
      }
      return new TopologyAttempt(segments, closed, finite, criticalPoint,
          new double[] {temperatureAtMaximumPressure, maximumPressure, 0.0},
          new double[] {maximumTemperature, pressureAtMaximumTemperature, 0.0}, threePhasePointCount, flashCount,
          failedFlashCount, null);
    }

    private void addCellContours(PhaseState[][] states, int temperatureIndex, int pressureIndex,
        Map<String, ContourPoint> points, Map<String, Set<String>> graph) {
      PhaseState[] corners = new PhaseState[] {states[temperatureIndex][pressureIndex],
          states[temperatureIndex + 1][pressureIndex], states[temperatureIndex + 1][pressureIndex + 1],
          states[temperatureIndex][pressureIndex + 1]};
      if (!corners[0].valid || !corners[1].valid || !corners[2].valid || !corners[3].valid) {
        return;
      }
      int mask = (corners[0].hydrocarbonTwoPhase ? 1 : 0) | (corners[1].hydrocarbonTwoPhase ? 2 : 0)
          | (corners[2].hydrocarbonTwoPhase ? 4 : 0) | (corners[3].hydrocarbonTwoPhase ? 8 : 0);
      if (mask == 0 || mask == 15) {
        return;
      }
      List<int[]> edgePairs = contourEdgePairs(mask, corners);
      for (int[] pair : edgePairs) {
        String firstKey = edgeKey(temperatureIndex, pressureIndex, pair[0]);
        String secondKey = edgeKey(temperatureIndex, pressureIndex, pair[1]);
        if (!points.containsKey(firstKey)) {
          points.put(firstKey, refineEdge(corners, pair[0]));
        }
        if (!points.containsKey(secondKey)) {
          points.put(secondKey, refineEdge(corners, pair[1]));
        }
        connect(graph, firstKey, secondKey);
      }
    }

    private List<int[]> contourEdgePairs(int mask, PhaseState[] corners) {
      List<int[]> pairs = new ArrayList<int[]>();
      switch (mask) {
      case 1:
        pairs.add(new int[] {3, 0});
        break;
      case 2:
        pairs.add(new int[] {0, 1});
        break;
      case 3:
        pairs.add(new int[] {3, 1});
        break;
      case 4:
        pairs.add(new int[] {1, 2});
        break;
      case 5:
        addAmbiguousPairs(pairs, corners, true);
        break;
      case 6:
        pairs.add(new int[] {0, 2});
        break;
      case 7:
        pairs.add(new int[] {3, 2});
        break;
      case 8:
        pairs.add(new int[] {2, 3});
        break;
      case 9:
        pairs.add(new int[] {0, 2});
        break;
      case 10:
        addAmbiguousPairs(pairs, corners, false);
        break;
      case 11:
        pairs.add(new int[] {1, 2});
        break;
      case 12:
        pairs.add(new int[] {1, 3});
        break;
      case 13:
        pairs.add(new int[] {0, 1});
        break;
      case 14:
        pairs.add(new int[] {3, 0});
        break;
      default:
        break;
      }
      return pairs;
    }

    private void addAmbiguousPairs(List<int[]> pairs, PhaseState[] corners, boolean diagonalZeroTwoInside) {
      double centerTemperature = 0.5 * (corners[0].temperatureK + corners[2].temperatureK);
      double centerPressure = Math.sqrt(corners[0].pressureBara * corners[2].pressureBara);
      PhaseState center = evaluate(centerTemperature, centerPressure);
      boolean centerInside = center.valid && center.hydrocarbonTwoPhase;
      if (centerInside == diagonalZeroTwoInside) {
        pairs.add(new int[] {0, 1});
        pairs.add(new int[] {2, 3});
      } else {
        pairs.add(new int[] {3, 0});
        pairs.add(new int[] {1, 2});
      }
    }

    private ContourPoint refineEdge(PhaseState[] corners, int edge) {
      PhaseState first;
      PhaseState second;
      if (edge == 0) {
        first = corners[0];
        second = corners[1];
      } else if (edge == 1) {
        first = corners[1];
        second = corners[2];
      } else if (edge == 2) {
        first = corners[3];
        second = corners[2];
      } else {
        first = corners[0];
        second = corners[3];
      }
      PhaseState inside = first.hydrocarbonTwoPhase ? first : second;
      PhaseState outside = first.hydrocarbonTwoPhase ? second : first;
      for (int iteration = 0; iteration < bisectionIterations; iteration++) {
        double midpointTemperature = 0.5 * (inside.temperatureK + outside.temperatureK);
        double midpointPressure = Math.sqrt(inside.pressureBara * outside.pressureBara);
        PhaseState midpoint = evaluate(midpointTemperature, midpointPressure);
        if (!midpoint.valid) {
          break;
        }
        if (midpoint.hydrocarbonTwoPhase) {
          inside = midpoint;
        } else {
          outside = midpoint;
        }
      }
      return new ContourPoint(0.5 * (inside.temperatureK + outside.temperatureK),
          Math.sqrt(inside.pressureBara * outside.pressureBara), inside.compositionDistance, inside.phaseCount >= 3);
    }

    private PhaseState evaluate(double temperatureK, double pressureBara) {
      flashCount++;
      try {
        SystemInterface point = template.clone();
        point.setTemperature(temperatureK);
        point.setPressure(pressureBara);
        point.setMultiPhaseCheck(true);
        point.setMaxNumberOfPhases(3);
        new ThermodynamicOperations(point).TPflash();
        PhaseInterface gas = null;
        PhaseInterface oil = null;
        int materialPhaseCount = 0;
        for (int phaseIndex = 0; phaseIndex < point.getNumberOfPhases(); phaseIndex++) {
          if (point.getBeta(phaseIndex) <= MINIMUM_PHASE_FRACTION) {
            continue;
          }
          materialPhaseCount++;
          PhaseInterface phase = point.getPhase(phaseIndex);
          String phaseName = phase.getPhaseTypeName().toLowerCase(Locale.ROOT);
          if (phaseName.contains("gas") || phaseName.contains("vapour") || phaseName.contains("vapor")) {
            gas = phase;
          } else if (phaseName.contains("oil")
              || (phaseName.contains("liquid") && !phaseName.contains("aqueous") && !phaseName.contains("water"))) {
            oil = phase;
          }
        }
        double compositionDistance = gas == null || oil == null ? Double.POSITIVE_INFINITY
            : phaseCompositionDistance(gas, oil);
        return new PhaseState(temperatureK, pressureBara, true, gas != null && oil != null, compositionDistance,
            materialPhaseCount);
      } catch (Exception error) {
        failedFlashCount++;
        return PhaseState.failed(temperatureK, pressureBara);
      }
    }
  }

  private static final class PhaseState {
    private final double temperatureK;
    private final double pressureBara;
    private final boolean valid;
    private final boolean hydrocarbonTwoPhase;
    private final double compositionDistance;
    private final int phaseCount;

    private PhaseState(double temperatureK, double pressureBara, boolean valid, boolean hydrocarbonTwoPhase,
        double compositionDistance, int phaseCount) {
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.valid = valid;
      this.hydrocarbonTwoPhase = hydrocarbonTwoPhase;
      this.compositionDistance = compositionDistance;
      this.phaseCount = phaseCount;
    }

    private static PhaseState failed(double temperatureK, double pressureBara) {
      return new PhaseState(temperatureK, pressureBara, false, false, Double.POSITIVE_INFINITY, 0);
    }
  }

  private static final class ContourPoint {
    private final double temperatureK;
    private final double pressureBara;
    private final double compositionDistance;
    private final boolean threePhase;

    private ContourPoint(double temperatureK, double pressureBara, double compositionDistance, boolean threePhase) {
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.compositionDistance = compositionDistance;
      this.threePhase = threePhase;
    }
  }

  private static final class Polyline {
    private final List<ContourPoint> points;
    private final boolean closed;

    private Polyline(List<ContourPoint> points, boolean closed) {
      this.points = points;
      this.closed = closed;
    }
  }

  private static List<Polyline> tracePolylines(Map<String, ContourPoint> points, Map<String, Set<String>> graph) {
    List<Polyline> output = new ArrayList<Polyline>();
    Set<String> visitedEdges = new HashSet<String>();
    for (Map.Entry<String, Set<String>> entry : graph.entrySet()) {
      if (entry.getValue().size() == 2) {
        continue;
      }
      for (String neighbor : entry.getValue()) {
        if (!visitedEdges.contains(undirectedEdgeKey(entry.getKey(), neighbor))) {
          output.add(walkPolyline(entry.getKey(), neighbor, points, graph, visitedEdges));
        }
      }
    }
    for (Map.Entry<String, Set<String>> entry : graph.entrySet()) {
      for (String neighbor : entry.getValue()) {
        if (!visitedEdges.contains(undirectedEdgeKey(entry.getKey(), neighbor))) {
          output.add(walkPolyline(entry.getKey(), neighbor, points, graph, visitedEdges));
        }
      }
    }
    return output;
  }

  private static Polyline walkPolyline(String start, String next, Map<String, ContourPoint> points,
      Map<String, Set<String>> graph, Set<String> visitedEdges) {
    List<ContourPoint> output = new ArrayList<ContourPoint>();
    output.add(points.get(start));
    String previous = start;
    String current = next;
    visitedEdges.add(undirectedEdgeKey(previous, current));
    output.add(points.get(current));
    while (!current.equals(start)) {
      String candidate = null;
      for (String neighbor : graph.get(current)) {
        String key = undirectedEdgeKey(current, neighbor);
        if (!visitedEdges.contains(key) && !neighbor.equals(previous)) {
          candidate = neighbor;
          break;
        }
      }
      if (candidate == null) {
        break;
      }
      previous = current;
      current = candidate;
      visitedEdges.add(undirectedEdgeKey(previous, current));
      output.add(points.get(current));
    }
    return new Polyline(output, current.equals(start));
  }

  private static void connect(Map<String, Set<String>> graph, String first, String second) {
    if (!graph.containsKey(first)) {
      graph.put(first, new HashSet<String>());
    }
    if (!graph.containsKey(second)) {
      graph.put(second, new HashSet<String>());
    }
    graph.get(first).add(second);
    graph.get(second).add(first);
  }

  private static String edgeKey(int temperatureIndex, int pressureIndex, int edge) {
    if (edge == 0) {
      return "H:" + temperatureIndex + ":" + pressureIndex;
    }
    if (edge == 1) {
      return "V:" + (temperatureIndex + 1) + ":" + pressureIndex;
    }
    if (edge == 2) {
      return "H:" + temperatureIndex + ":" + (pressureIndex + 1);
    }
    return "V:" + temperatureIndex + ":" + pressureIndex;
  }

  private static String undirectedEdgeKey(String first, String second) {
    return first.compareTo(second) < 0 ? first + "|" + second : second + "|" + first;
  }

  private static double[] linearValues(double minimum, double maximum, int count) {
    double[] values = new double[count];
    for (int index = 0; index < count; index++) {
      values[index] = minimum + (maximum - minimum) * index / (count - 1.0);
    }
    return values;
  }

  private static double[] logarithmicValues(double minimum, double maximum, int count) {
    double[] values = new double[count];
    double logarithmicSpan = Math.log(maximum / minimum);
    for (int index = 0; index < count; index++) {
      values[index] = minimum * Math.exp(logarithmicSpan * index / (count - 1.0));
    }
    return values;
  }

  private static double phaseCompositionDistance(PhaseInterface first, PhaseInterface second) {
    int componentCount = Math.min(first.getNumberOfComponents(), second.getNumberOfComponents());
    double squaredDistance = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double difference = first.getComponent(componentIndex).getx() - second.getComponent(componentIndex).getx();
      squaredDistance += difference * difference;
    }
    return Math.sqrt(squaredDistance);
  }

  private static EnvelopeSegment.PhaseType storedTypeAtMaximumTemperature(List<EnvelopeSegment> segments) {
    EnvelopeSegment.PhaseType selected = EnvelopeSegment.PhaseType.DEW;
    double maximumTemperature = Double.NEGATIVE_INFINITY;
    for (EnvelopeSegment segment : segments) {
      for (double temperature : segment.getTemperatures()) {
        if (Double.isFinite(temperature) && temperature > maximumTemperature) {
          maximumTemperature = temperature;
          selected = segment.getPhaseType();
        }
      }
    }
    return selected;
  }

  private static boolean finitePhysicalValues(double[] temperatures, double[] pressures) {
    if (temperatures.length == 0 || temperatures.length != pressures.length) {
      return false;
    }
    for (int index = 0; index < temperatures.length; index++) {
      if (!Double.isFinite(temperatures[index]) || !Double.isFinite(pressures[index]) || temperatures[index] < 50.0
          || temperatures[index] > 2000.0 || pressures[index] <= 0.0 || pressures[index] > 100000.0) {
        return false;
      }
    }
    return true;
  }

  private static double[] safeArray(double[] values) {
    return values == null ? new double[] {0.0, 0.0, 0.0} : values;
  }
}
