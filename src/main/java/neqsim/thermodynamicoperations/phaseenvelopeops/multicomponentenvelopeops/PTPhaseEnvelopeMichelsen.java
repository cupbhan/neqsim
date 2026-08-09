/*
 * PTPhaseEnvelopeMichelsen.java
 *
 * Unified phase envelope calculation using the Michelsen continuation method. Replaces the
 * fragmented implementations (PTphaseEnvelope1, PTphaseEnvelopeMay, PTphaseEnvelopeNew,
 * PTphaseEnvelopeNew2) with a single robust algorithm. Note: PTphaseEnvelopeNew3 (grid-based) and
 * PTphaseEnvelope (base class for CricondenBarFlash/CricondenThermFlash) are retained separately.
 *
 * Algorithm: Michelsen (1980) natural parameter continuation with: - 3rd-order polynomial predictor
 * through last 4 converged points - Adaptive step size control based on Newton iteration count -
 * Automatic specification variable selection (most sensitive variable) - Iterative (non-recursive)
 * restart with automatic branch switching - Critical point detection using K-value convergence
 * monitoring - Dynamic ArrayList storage for arbitrary-length envelopes
 *
 * References: - Michelsen, M.L. (1980). "Calculation of phase envelopes and critical points for
 * multicomponent mixtures." Fluid Phase Equilibria, 4, 1-10. - Michelsen, M.L. and Mollerup, J.M.
 * (2004). "Thermodynamic Models: Fundamentals &amp; Computational Aspects." Tie-Line Publications.
 */

package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.BaseOperation;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Unified PT phase envelope calculation using the Michelsen continuation method.
 *
 * <p>
 * This implementation traces the two-phase boundary (bubble and dew point curves) in pressure-temperature space for
 * multicomponent mixtures. The algorithm uses a predictor-corrector approach with natural parameter continuation to
 * follow the envelope through the critical point region.
 * </p>
 *
 * <p>
 * The algorithm traces from one end of the envelope (dew or bubble side at low pressure), through the cricondenbar and
 * cricondentherm, across the critical point, and back down the other side. If the trace crashes before completing the
 * full envelope, it automatically restarts from the opposite side to fill in the gap.
 * </p>
 *
 * <p>
 * Key improvements over legacy implementations:
 * </p>
 * <ul>
 * <li>Non-recursive restart (no stack overflow risk)</li>
 * <li>Dynamic ArrayList storage (no fixed 10,000-point limit)</li>
 * <li>Configurable step limits and pressure bounds</li>
 * <li>K-value reset and Tmin stopping criterion for restart branch</li>
 * <li>Clean separated data output for bubble and dew point branches</li>
 * </ul>
 *
 * @author asmund
 * @version 2.0
 */
public class PTPhaseEnvelopeMichelsen extends BaseOperation {
  /** Serialization version UID. */
  private static final long serialVersionUID = 1000;
  /** Logger object for class. */
  private static final Logger logger = LogManager.getLogger(PTPhaseEnvelopeMichelsen.class);

  /** Constant in the Wilson K-value correlation (5.373). */
  private static final double WILSON_CONST = 5.373;
  /** Overall mole fractions below this value are excluded from phase-envelope tracing. */
  private static final double ZERO_FRACTION_THRESHOLD = 1.0e-12;
  /** Maximum Newton iterations for Wilson temperature estimate. */
  private static final int MAX_WILSON_ITERATIONS = 1000;
  /** Maximum number of attempts when estimating the first point. */
  private static final int FIRST_POINT_ATTEMPTS = 5;
  /** Temperature step (K) used when searching for the first point. */
  private static final double FIRST_POINT_STEP = 2.0;
  private static final double FIRST_POINT_SCAN_MINIMUM_K = 200.0;
  private static final double FIRST_POINT_SCAN_MAXIMUM_K = 1000.0;
  private static final double FIRST_POINT_SCAN_STEP_K = 20.0;
  private static final double MINIMUM_NON_TRIVIAL_LOG_K = 0.5;
  /** Maximum envelope points per branch to prevent infinite loops. */
  private static final int MAX_ENVELOPE_ITERATIONS = 9980;
  /** Maximum points per quality line. */
  private static final int MAX_QUALITY_LINE_POINTS = 5000;
  /**
   * Multiple of the maximum continuation step ({@code dTmax}/{@code dPmax}) above which two consecutive points are
   * treated as belonging to disjoint branch segments. A NaN break is inserted between them so plotters do not draw a
   * spurious straight chord across a discontinuity (e.g. between a near-critical sub-branch and the low-pressure tail).
   */
  private static final double JUMP_BREAK_FACTOR = 3.0;

  // --- Configuration ---
  private double maxPressure = 1000.0;
  private double minPressure = 1.0;
  private double lowPres = 1.0;
  private double dTmax = 10.0;
  private double dPmax = 10.0;
  private double phaseFraction = 1e-10;
  private boolean bubblePointFirst = true;
  private int maximumEnvelopeIterations = MAX_ENVELOPE_ITERATIONS;
  private boolean iterationLimitReached = false;

  /**
   * Why a continuation pass stopped adding points.
   *
   * <p>
   * Only {@link #PRESSURE_FLOOR} and {@link #CONTINUATION_END} are physical ends of the traced branch. Every other
   * value means the branch was cut short by a solver or domain limit, so the traced curve is a truncated piece of the
   * real boundary and must not be reported as a closed envelope.
   * </p>
   */
  public enum TerminationReason {
    /** The branch came back down below the minimum pressure after passing the cricondentherm. */
    PRESSURE_FLOOR,
    /** The branch ran into the configured maximum pressure while still climbing. */
    PRESSURE_CEILING,
    /** The restart pass reached the temperature where the first pass had already traced the boundary. */
    RESTART_OVERLAP,
    /** The continuation used up its iteration budget. */
    ITERATION_LIMIT,
    /** The continuation stopped on its own without hitting a configured limit. */
    CONTINUATION_END
  }

  private final List<TerminationReason> terminationReasons = new ArrayList<TerminationReason>();

  // --- System reference ---
  private SystemInterface system;

  // --- State tracking ---
  /** True while adding points to dew lists; false for bubble lists. */
  private boolean isDewPhase = true;
  /** Starting component index for Wilson estimate. */
  private int speceq = 0;
  /** True once K-values have diverged sufficiently after a CP (for multiple CP detection). */
  private boolean kValuesDivergedAfterCP = false;

  // --- Results: dynamic storage ---
  private ArrayList<Double> dewPointTemperatures = new ArrayList<Double>();
  private ArrayList<Double> dewPointPressures = new ArrayList<Double>();
  private ArrayList<Double> bubblePointTemperatures = new ArrayList<Double>();
  private ArrayList<Double> bubblePointPressures = new ArrayList<Double>();
  private ArrayList<Double> dewPointEnthalpies = new ArrayList<Double>();
  private ArrayList<Double> bubblePointEnthalpies = new ArrayList<Double>();
  private ArrayList<Double> dewPointDensities = new ArrayList<Double>();
  private ArrayList<Double> bubblePointDensities = new ArrayList<Double>();
  private ArrayList<Double> dewPointEntropies = new ArrayList<Double>();
  private ArrayList<Double> bubblePointEntropies = new ArrayList<Double>();
  /** Exact continuation states retained for secondary-phase stability analysis. */
  private ArrayList<BoundaryStateSeed> boundaryStateSeeds = new ArrayList<BoundaryStateSeed>();

  // --- Output arrays (built after run) ---
  private double[] dewTempArray = new double[0];
  private double[] dewPresArray = new double[0];
  private double[] dewEnthalpyArray = new double[0];
  private double[] dewDensityArray = new double[0];
  private double[] dewEntropyArray = new double[0];
  private double[] bubTempArray = new double[0];
  private double[] bubPresArray = new double[0];
  private double[] bubEnthalpyArray = new double[0];
  private double[] bubDensityArray = new double[0];
  private double[] bubEntropyArray = new double[0];

  /** Structured per-branch view of the envelope (built after run). */
  private List<EnvelopeSegment> segments = Collections.emptyList();

  // --- Critical point and characteristic points ---
  private double[] cricondenTherm = new double[3];
  private double[] cricondenBar = new double[3];
  private double[] cricondenThermX = new double[100];
  private double[] cricondenThermY = new double[100];
  private double[] cricondenBarX = new double[100];
  private double[] cricondenBarY = new double[100];

  // --- Saved first-pass data for merging after restart ---
  private double[] cricondenThermFirst = new double[3];
  private double[] cricondenBarFirst = new double[3];
  private double[] cricondenThermXFirst = new double[100];
  private double[] cricondenThermYFirst = new double[100];
  private double[] cricondenBarXFirst = new double[100];
  private double[] cricondenBarYFirst = new double[100];

  // --- Critical points (supports multiple CPs for complex mixtures) ---
  private ArrayList<double[]> criticalPoints = new ArrayList<double[]>();

  // --- Quality line data (keyed by "qualityT_0.5", "qualityP_0.5", etc.) ---
  private Map<String, double[]> qualityLineData = new HashMap<String, double[]>();
  /** Traced quality line beta values. */
  private double[] qualityBetaValues = new double[0];

  // --- Three-phase stability analysis results ---
  /** Temperatures of envelope points where 3+ phases detected. */
  private double[] threePhaseRegionT = new double[0];
  /** Pressures of envelope points where 3+ phases detected. */
  private double[] threePhaseRegionP = new double[0];
  /** Non-destructive secondary-phase stability diagnostics along the traced boundary. */
  private List<SecondaryStabilitySample> secondaryStabilitySamples = Collections.emptyList();
  /** Adjacent samples that bracket entry to or exit from a three-phase region. */
  private List<SecondaryStabilityBracket> secondaryStabilityBrackets = Collections.emptyList();
  /** Refined two-phase/three-phase transition candidates. */
  private List<ThreePhasePointCandidate> threePhasePointCandidates = Collections.emptyList();

  /**
   * Default constructor.
   */
  public PTPhaseEnvelopeMichelsen() {
  }

  /**
   * Constructor for PTPhaseEnvelopeMichelsen.
   *
   * @param system the thermodynamic system
   * @param name output file name (unused, kept for API compatibility)
   * @param phaseFraction initial phase fraction (near 0 = bubble, near 1 = dew)
   * @param lowPres starting low pressure in bara
   * @param bubfirst if true, trace bubble point curve first
   */
  public PTPhaseEnvelopeMichelsen(SystemInterface system, String name, double phaseFraction, double lowPres,
      boolean bubfirst) {
    this.system = system;
    this.phaseFraction = phaseFraction;
    this.lowPres = lowPres;
    this.bubblePointFirst = bubfirst;
  }

  /** {@inheritDoc} */
  @Override
  public SystemInterface getThermoSystem() {
    return system;
  }

  /**
   * Remove components with negligible overall mole fractions from a private clone of the working system. These
   * components make the continuation Jacobian singular because its equations contain terms divided by the overall mole
   * fraction.
   *
   * @return {@code true} if one or more components were removed
   */
  private boolean filterZeroFractionComponents() {
    List<String> componentNames = new ArrayList<String>();
    for (int i = 0; i < system.getNumberOfComponents(); i++) {
      if (system.getComponent(i).getz() < ZERO_FRACTION_THRESHOLD) {
        componentNames.add(system.getComponent(i).getComponentName());
      }
    }
    if (componentNames.isEmpty()) {
      return false;
    }

    logger.warn("Phase envelope: removing {} component(s) with z < {}: {}", componentNames.size(),
        ZERO_FRACTION_THRESHOLD, componentNames);
    SystemInterface filteredSystem = system.clone();
    for (String componentName : componentNames) {
      filteredSystem.removeComponent(componentName);
    }
    system = filteredSystem;
    return true;
  }

  /**
   * {@inheritDoc}
   *
   * <p>
   * Traces the phase envelope using a two-pass approach. The primary pass traces from the starting side (dew or bubble)
   * through the full envelope. If it crashes, a second pass traces from the opposite side to fill in the gap. This
   * avoids the recursive restart used in legacy implementations.
   * </p>
   */
  /**
   * Reports whether the two phases currently differ, so the trivial {@code K = 1} root is rejected.
   *
   * <p>
   * A saturation flash on a water-dominated feed readily "converges" with both phases equal to the feed. That state
   * satisfies the equations but is not a phase boundary, and the continuation Jacobian is singular there, so the
   * envelope then fails to advance at all. Measured on the real 23-component fluid: the trivial root reports an
   * {@code ln K} spread of 0 with an incipient water fraction identical to the feed, while a genuine dew point reports
   * a spread near 80.
   * </p>
   *
   * @return {@code true} when the phase compositions differ enough to be a real split
   */
  private boolean hasNonTrivialPhaseSplit() {
    if (system.getNumberOfPhases() < 2) {
      return false;
    }
    double maximumAbsoluteLogK = 0.0;
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      double vapourFraction = system.getPhase(0).getComponent(i).getx();
      double liquidFraction = system.getPhase(1).getComponent(i).getx();
      if (!(vapourFraction > 0.0) || !(liquidFraction > 0.0)) {
        continue;
      }
      maximumAbsoluteLogK = Math.max(maximumAbsoluteLogK, Math.abs(Math.log(vapourFraction / liquidFraction)));
    }
    return maximumAbsoluteLogK > MINIMUM_NON_TRIVIAL_LOG_K;
  }

  /**
   * Searches a wide temperature range for a first saturation point when the Wilson-anchored guess is too far off.
   *
   * <p>
   * Used only as a fallback. The scan sweeps outward from the current estimate in both directions so a true saturation
   * temperature several hundred kelvin below the guess is still reachable, and rejects a flash that returned its own
   * starting value untouched, which is how a non-converging saturation flash reports success.
   * </p>
   *
   * @param testOps operations bound to the working system
   * @param beta specified phase fraction; below 0.5 selects the bubble-point flash
   * @return converged saturation temperature in K, or NaN when the scan found none
   */
  private double scanForFirstSaturationPoint(ThermodynamicOperations testOps, double beta) {
    double centre = system.getTemperature();
    if (!(centre > FIRST_POINT_SCAN_MINIMUM_K) || !(centre < FIRST_POINT_SCAN_MAXIMUM_K)) {
      centre = 0.5 * (FIRST_POINT_SCAN_MINIMUM_K + FIRST_POINT_SCAN_MAXIMUM_K);
    }
    int steps = (int) Math.ceil((FIRST_POINT_SCAN_MAXIMUM_K - FIRST_POINT_SCAN_MINIMUM_K) / FIRST_POINT_SCAN_STEP_K);
    for (int offset = 0; offset <= steps; offset++) {
      for (int direction = 0; direction < 2; direction++) {
        if (offset == 0 && direction == 1) {
          continue;
        }
        double start = centre + (direction == 0 ? offset : -offset) * FIRST_POINT_SCAN_STEP_K;
        if (start < FIRST_POINT_SCAN_MINIMUM_K || start > FIRST_POINT_SCAN_MAXIMUM_K) {
          continue;
        }
        try {
          // The failed attempts that led here leave the phases in a broken state (NaN compressibility), which
          // makes every later flash fail no matter the temperature. Re-seed before each probe.
          system.setTemperature(start);
          system.setPressure(lowPres);
          system.init(0);
          resetKValuesWithWilson();
          if (beta < 0.5) {
            testOps.bubblePointTemperatureFlash();
          } else {
            testOps.dewPointTemperatureFlash();
          }
        } catch (Exception ignored) {
          continue;
        }
        double converged = system.getTemperature();
        if (!Double.isNaN(converged) && converged > FIRST_POINT_SCAN_MINIMUM_K && converged < FIRST_POINT_SCAN_MAXIMUM_K
            && Math.abs(converged - start) > 1.0e-6 && hasNonTrivialPhaseSplit()) {
          return converged;
        }
      }
    }
    return Double.NaN;
  }

  @Override
  public void run() {
    iterationLimitReached = false;
    terminationReasons.clear();
    boundaryStateSeeds.clear();
    double initialTemp = system.getTemperature();
    double initialPres = system.getPressure();

    // isDewPhase determines which list receives each traced point.
    // When starting from dew side (bubblePointFirst=false, phaseFraction~1),
    // points go to dew lists. At CP, they switch to bubble lists.
    isDewPhase = true;

    boolean needRestart = false;
    double restartTmin = 0.0;

    // === Two-pass loop: primary trace + optional restart ===
    for (int pass = 0; pass < 2; pass++) {
      if (pass == 1) {
        if (!needRestart) {
          break;
        }

        // Save first-pass cricondentherm/bar for later merge
        cricondenThermFirst = cricondenTherm.clone();
        cricondenBarFirst = cricondenBar.clone();
        cricondenThermXFirst = cricondenThermX.clone();
        cricondenThermYFirst = cricondenThermY.clone();
        cricondenBarXFirst = cricondenBarX.clone();
        cricondenBarYFirst = cricondenBarY.clone();

        // Reset tracking for second pass
        cricondenTherm = new double[3];
        cricondenBar = new double[3];
        cricondenThermX = new double[100];
        cricondenThermY = new double[100];
        cricondenBarX = new double[100];
        cricondenBarY = new double[100];

        // Flip conditions for second pass
        phaseFraction = 1.0 - phaseFraction;
        bubblePointFirst = !bubblePointFirst;
        isDewPhase = false;
        kValuesDivergedAfterCP = false;

        // Insert NaN "break" markers so plotters do not draw a straight line
        // from the last first-pass point to the first second-pass point
        // (the two passes cover disjoint parts of the envelope and must be
        // rendered as separate polylines).
        addBranchBreak();

        // Reset K-values using Wilson correlation
        resetKValuesWithWilson();
      }

      // === Standard initialization ===
      speceq = 0;
      system.init(0);
      if (pass == 0 && filterZeroFractionComponents()) {
        system.init(0);
      }

      for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
        if (system.getComponent(i).getz() < 1e-10) {
          continue;
        }
        if (system.getPhase(0).getComponent(i).getIonicCharge() == 0) {
          if (bubblePointFirst
              && system.getPhase(0).getComponent(speceq).getTC() > system.getPhase(0).getComponent(i).getTC()) {
            speceq = system.getPhase(0).getComponent(i).getComponentNumber();
          }
          if (!bubblePointFirst
              && system.getPhase(0).getComponent(speceq).getTC() < system.getPhase(0).getComponent(i).getTC()) {
            speceq = system.getPhase(0).getComponent(i).getComponentNumber();
          }
        }
      }

      // Estimate initial temperature using Wilson correlation
      double temp = tempKWilson(phaseFraction, lowPres);
      // Diagnostic override: lets a probe seed the first point on a known branch to separate "the continuation
      // cannot march" from "the continuation starts on the wrong saturation branch". Never set in production.
      String seedOverride = System.getProperty("neqsim.probe.seedTemperatureK");
      if (seedOverride != null) {
        temp = Double.parseDouble(seedOverride);
      }
      if (Double.isNaN(temp)) {
        temp = system.getPhase(0).getComponent(speceq).getTC() - 20.0;
      }
      system.setTemperature(temp);
      system.setPressure(lowPres);

      // Converge first point using saturation flash
      ThermodynamicOperations testOps = new ThermodynamicOperations(system);
      boolean firstPointConverged = false;
      for (int attempt = 0; attempt < FIRST_POINT_ATTEMPTS; attempt++) {
        try {
          if (phaseFraction < 0.5) {
            temp += attempt * FIRST_POINT_STEP;
            system.setTemperature(temp);
            testOps.bubblePointTemperatureFlash();
          } else {
            temp += attempt * FIRST_POINT_STEP;
            system.setTemperature(temp);
            testOps.dewPointTemperatureFlash();
          }
        } catch (Exception ex) {
          continue;
        }
        double tempNy = system.getTemperature();
        if (!Double.isNaN(tempNy)) {
          temp = tempNy;
          firstPointConverged = true;
          break;
        }
      }
      if (!firstPointConverged) {
        // The Wilson estimate anchors on the extreme-Tc component, which is meaningless once a nearly immiscible
        // component such as water dominates the composition: for a 65% water heavy oil it overshoots the true dew
        // temperature by roughly 300 K, while the march above only searches 20 K upward. Fall back to a wide
        // bidirectional scan. This runs only after the original attempts have already failed, so fluids that
        // converged before keep converging on exactly the same first point.
        double scanned = scanForFirstSaturationPoint(testOps, phaseFraction);
        if (System.getProperty("neqsim.probe.firstPoint") != null) {
          System.err.printf("[probe] pass=%d beta=%.3g lowPres=%.4g wilsonTemp=%.2f scan=%s%n", pass, phaseFraction,
              lowPres, temp, Double.isNaN(scanned) ? "FAILED" : String.format("%.2f K", scanned));
        }
        if (!Double.isNaN(scanned)) {
          temp = scanned;
          firstPointConverged = true;
        }
      }
      if (!firstPointConverged) {
        logger.warn("Could not converge first envelope point for pass={}, beta={}", pass, phaseFraction);
        continue;
      }

      // Set up for continuation
      system.setBeta(phaseFraction);
      system.setPressure(lowPres);
      system.setTemperature(temp);

      SysNewtonRhapsonPhaseEnvelope nonLinSolver = new SysNewtonRhapsonPhaseEnvelope(system, 2,
          system.getPhase(0).getNumberOfComponents());
      nonLinSolver.dTmax = this.dTmax;
      nonLinSolver.dPmax = this.dPmax;
      nonLinSolver.setu();
      boolean passedCricoT = false;

      // === Main continuation loop ===
      int np;
      for (np = 1; np < maximumEnvelopeIterations; np++) {
        try {
          nonLinSolver.calcInc(np);
          nonLinSolver.solve(np);
        } catch (Exception e0) {
          if (System.getProperty("neqsim.probe.firstPoint") != null) {
            System.err.printf("[probe] continuation failed at np=%d T=%.2f P=%.4f : %s: %s%n", np,
                system.getTemperature(), system.getPressure(), e0.getClass().getSimpleName(), e0.getMessage());
            StackTraceElement[] frames = e0.getStackTrace();
            for (int f = 0; f < Math.min(4, frames.length); f++) {
              System.err.println("[probe]     at " + frames[f]);
            }
          }
          if (pass == 0) {
            // Primary trace crashed: schedule restart from opposite side
            needRestart = true;
            if (np > 2) {
              // Use recent stored temperature as Tmin for second pass
              ArrayList<Double> tempList = isDewPhase ? dewPointTemperatures : bubblePointTemperatures;
              if (!tempList.isEmpty()) {
                restartTmin = tempList.get(tempList.size() - 1);
              } else {
                restartTmin = system.getTemperature();
              }
            } else {
              restartTmin = system.getTemperature();
            }
          }
          np = np - 1;
          break;
        }

        double currentT = system.getTemperature();
        double currentP = system.getPressure();

        // === Critical point detection via K-value convergence ===
        double Kvallc = system.getPhase(0).getComponent(nonLinSolver.lc).getx()
            / system.getPhase(1).getComponent(nonLinSolver.lc).getx();
        double Kvalhc = system.getPhase(0).getComponent(nonLinSolver.hc).getx()
            / system.getPhase(1).getComponent(nonLinSolver.hc).getx();

        if (!nonLinSolver.etterCP) {
          if (Kvallc < 1.05 && Kvalhc > 0.95) {
            nonLinSolver.npCrit = np;
            system.invertPhaseTypes();
            nonLinSolver.etterCP = true;
            isDewPhase = !isDewPhase;
            nonLinSolver.calcCrit();
            criticalPoints.add(new double[] { system.getTC(), system.getPC() });
            kValuesDivergedAfterCP = false;
            addBranchBreak();
          }
        } else {
          // Track K-value divergence after CP: K-values must move significantly
          // away from 1 before we allow re-detection. This prevents false triggers
          // when K-values linger near 1 just past the first CP, or naturally
          // approach 1 at low pressures at the tail end of the envelope.
          if (!kValuesDivergedAfterCP) {
            if (Kvallc < 0.5 || Kvalhc > 2.0) {
              kValuesDivergedAfterCP = true;
            }
          } else {
            // K-values have diverged strongly and are now re-converging — second
            // CP.
            // Additional safety guards: must be far enough from first CP (>50
            // steps)
            // and at pressure above 5 bar (not at the tail end of the envelope).
            if (Kvallc < 1.05 && Kvalhc > 0.95 && (np - nonLinSolver.npCrit) > 50 && currentP > 5.0) {
              nonLinSolver.npCrit = np;
              system.invertPhaseTypes();
              isDewPhase = !isDewPhase;
              nonLinSolver.calcCrit();
              criticalPoints.add(new double[] { system.getTC(), system.getPC() });
              kValuesDivergedAfterCP = false;
              addBranchBreak();
            }
          }
        }

        // === Cricondentherm tracking (maximum temperature) ===
        if (currentT > cricondenTherm[0]) {
          cricondenTherm[1] = currentP;
          cricondenTherm[0] = currentT;
          for (int ii = 0; ii < nonLinSolver.numberOfComponents; ii++) {
            cricondenThermX[ii] = system.getPhase(1).getComponent(ii).getx();
            cricondenThermY[ii] = system.getPhase(0).getComponent(ii).getx();
          }
        } else {
          nonLinSolver.ettercricoT = true;
          passedCricoT = true;
        }

        // === Cricondenbar tracking (maximum pressure) ===
        if (currentP > cricondenBar[1]) {
          cricondenBar[0] = currentT;
          cricondenBar[1] = currentP;
          for (int ii = 0; ii < nonLinSolver.numberOfComponents; ii++) {
            cricondenBarX[ii] = system.getPhase(1).getComponent(ii).getx();
            cricondenBarY[ii] = system.getPhase(0).getComponent(ii).getx();
          }
        }

        // === Exit criteria ===
        if (currentP < minPressure && passedCricoT) {
          terminationReasons.add(TerminationReason.PRESSURE_FLOOR);
          break;
        }
        if (currentP > maxPressure) {
          terminationReasons.add(TerminationReason.PRESSURE_CEILING);
          break;
        }
        if (pass == 1 && restartTmin > 0 && currentT > restartTmin) {
          terminationReasons.add(TerminationReason.RESTART_OVERLAP);
          break;
        }

        // === Store the point ===
        if (currentT > 1e-6 && currentP > 1e-6 && !Double.isNaN(currentT) && !Double.isNaN(currentP)) {
          double enthalpy = system.getPhase(1).getEnthalpy() / system.getPhase(1).getNumberOfMolesInPhase()
              / system.getPhase(1).getMolarMass() / 1e3;
          double density = system.getPhase(1).getDensity();
          double entropy = system.getPhase(1).getEntropy() / system.getPhase(1).getNumberOfMolesInPhase()
              / system.getPhase(1).getMolarMass() / 1e3;

          if (isDewPhase) {
            dewPointTemperatures.add(currentT);
            dewPointPressures.add(currentP);
            dewPointEnthalpies.add(enthalpy);
            dewPointDensities.add(density);
            dewPointEntropies.add(entropy);
            boundaryStateSeeds.add(BoundaryStateSeed.capture(EnvelopeSegment.PhaseType.DEW, system));
          } else {
            bubblePointTemperatures.add(currentT);
            bubblePointPressures.add(currentP);
            bubblePointEnthalpies.add(enthalpy);
            bubblePointDensities.add(density);
            bubblePointEntropies.add(entropy);
            boundaryStateSeeds.add(BoundaryStateSeed.capture(EnvelopeSegment.PhaseType.BUBBLE, system));
          }
        }
      }
      if (np >= maximumEnvelopeIterations) {
        iterationLimitReached = true;
        terminationReasons.add(TerminationReason.ITERATION_LIMIT);
        logger.warn("Phase-envelope continuation stopped at configured iteration limit {} on pass {}",
            maximumEnvelopeIterations, pass);
      } else if (terminationReasons.size() <= pass) {
        terminationReasons.add(TerminationReason.CONTINUATION_END);
      }

      // Set critical point on the system
      system.setTemperature(system.getTC());
      system.setPressure(system.getPC());
    }

    // === Merge cricondentherm/bar from first and second passes ===
    if (needRestart) {
      if (cricondenThermFirst[0] > cricondenTherm[0]) {
        cricondenTherm = cricondenThermFirst;
        cricondenThermX = cricondenThermXFirst;
        cricondenThermY = cricondenThermYFirst;
      }
      if (cricondenBarFirst[1] > cricondenBar[1]) {
        cricondenBar = cricondenBarFirst;
        cricondenBarX = cricondenBarXFirst;
        cricondenBarY = cricondenBarYFirst;
      }
    }

    // Validate final cricondenbar and cricondentherm
    if (!Double.isFinite(cricondenBar[0]) || !Double.isFinite(cricondenBar[1])
        || (cricondenBar[0] == 0.0 && cricondenBar[1] == 0.0)) {
      cricondenBar[0] = initialTemp;
      cricondenBar[1] = initialPres;
    }
    if (!Double.isFinite(cricondenTherm[0]) || !Double.isFinite(cricondenTherm[1])
        || (cricondenTherm[0] == 0.0 && cricondenTherm[1] == 0.0)) {
      cricondenTherm[0] = initialTemp;
      cricondenTherm[1] = initialPres;
    }

    // Convert ArrayLists to output arrays
    buildOutputArrays();
  }

  /**
   * Check whether the phase envelope is closed.
   *
   * <p>
   * Both branches must carry at least three points, the two branches must be joined by at least one critical point, and
   * no pass may have been cut short by a solver or domain limit. Having points on both sides is not on its own evidence
   * of closure: a branch that ran into the maximum pressure while still climbing produces exactly that signature while
   * leaving most of the boundary untraced.
   * </p>
   *
   * @return true if both branches were traced to a physical end and joined at a critical point
   */
  public boolean isEnvelopeClosed() {
    return dewTempArray.length >= 3 && bubTempArray.length >= 3 && !criticalPoints.isEmpty() && !isTruncated();
  }

  /**
   * Check whether any continuation pass was cut short by a solver or domain limit.
   *
   * <p>
   * A restart pass that runs into the temperature range the first pass already traced is not truncated: the two passes
   * met, which is how a closed envelope finishes.
   * </p>
   *
   * @return true if a pass stopped at the pressure ceiling or ran out of iterations
   */
  public boolean isTruncated() {
    return terminationReasons.contains(TerminationReason.PRESSURE_CEILING)
        || terminationReasons.contains(TerminationReason.ITERATION_LIMIT);
  }

  /**
   * Getter for the reason each continuation pass stopped.
   *
   * @return one entry per executed pass, in pass order
   */
  public List<TerminationReason> getTerminationReasons() {
    return Collections.unmodifiableList(terminationReasons);
  }

  /**
   * Returns the number of critical points detected on the phase envelope. Typical mixtures have exactly one critical
   * point. Highly asymmetric mixtures (e.g., CH4 + nC16+, CO2 + heavy hydrocarbons) may exhibit two or more critical
   * points where the K-values converge to unity.
   *
   * @return number of detected critical points (0 if none detected)
   */
  public int getNumberOfCriticalPoints() {
    return criticalPoints.size();
  }

  /**
   * Trace quality lines (constant vapor fraction curves) inside the two-phase region. Each quality line is traced using
   * the same continuation method as the phase boundary, but with a fixed molar vapor fraction beta. At each point, the
   * volume fraction and mass fraction are also computed.
   *
   * <p>
   * Results are stored internally and accessible via {@link #get(String)} with keys:
   * </p>
   * <ul>
   * <li>{@code "qualityT_X"} - temperatures (K) along quality line at molar fraction X</li>
   * <li>{@code "qualityP_X"} - pressures (bara) along quality line at molar fraction X</li>
   * <li>{@code "qualityVolFrac_X"} - volume vapor fractions along the line</li>
   * <li>{@code "qualityMassFrac_X"} - mass vapor fractions along the line</li>
   * </ul>
   *
   * @param betaValues array of molar vapor fractions to trace (between 0 and 1 exclusive)
   */
  public void calcQualityLines(double[] betaValues) {
    this.qualityBetaValues = betaValues.clone();
    double maxEnvelopeP = Math.max(cricondenBar[1], 1.0);

    for (double beta : betaValues) {
      if (beta <= 0.0 || beta >= 1.0) {
        continue;
      }

      SystemInterface clonedSystem = system.clone();

      // Estimate initial temperature at low pressure for this beta
      double temp = tempKWilsonForSystem(clonedSystem, beta, lowPres);
      if (Double.isNaN(temp)) {
        continue;
      }
      clonedSystem.setTemperature(temp);
      clonedSystem.setPressure(lowPres);

      // Converge first point using saturation flash
      ThermodynamicOperations flashOps = new ThermodynamicOperations(clonedSystem);
      boolean converged = false;
      for (int attempt = 0; attempt < FIRST_POINT_ATTEMPTS; attempt++) {
        try {
          double tempAttempt = temp + attempt * FIRST_POINT_STEP;
          clonedSystem.setTemperature(tempAttempt);
          if (beta < 0.5) {
            flashOps.bubblePointTemperatureFlash();
          } else {
            flashOps.dewPointTemperatureFlash();
          }
        } catch (Exception ex) {
          continue;
        }
        double tempNy = clonedSystem.getTemperature();
        if (!Double.isNaN(tempNy)) {
          temp = tempNy;
          converged = true;
          break;
        }
      }
      if (!converged) {
        logger.debug("Could not converge first quality line point for beta={}", beta);
        continue;
      }

      clonedSystem.setBeta(beta);
      clonedSystem.setPressure(lowPres);
      clonedSystem.setTemperature(temp);

      SysNewtonRhapsonPhaseEnvelope solver = new SysNewtonRhapsonPhaseEnvelope(clonedSystem, 2,
          clonedSystem.getPhase(0).getNumberOfComponents());
      solver.dTmax = this.dTmax;
      solver.dPmax = this.dPmax;
      solver.setu();

      ArrayList<Double> qTemps = new ArrayList<Double>();
      ArrayList<Double> qPress = new ArrayList<Double>();
      ArrayList<Double> qVolFracs = new ArrayList<Double>();
      ArrayList<Double> qMassFracs = new ArrayList<Double>();

      boolean pastPeak = false;
      double peakP = 0.0;

      for (int np = 1; np < MAX_QUALITY_LINE_POINTS; np++) {
        try {
          solver.calcInc(np);
          solver.solve(np);
        } catch (Exception e) {
          break;
        }

        double currentT = clonedSystem.getTemperature();
        double currentP = clonedSystem.getPressure();

        if (Double.isNaN(currentT) || Double.isNaN(currentP) || currentT < 1e-6 || currentP < 1e-6) {
          break;
        }

        // Track pressure peak for exit criteria
        if (currentP > peakP) {
          peakP = currentP;
        } else {
          pastPeak = true;
        }

        // Exit if pressure exceeds envelope maximum
        if (currentP > maxEnvelopeP * 1.05) {
          break;
        }
        // Exit if pressure drops below minimum after passing peak
        if (pastPeak && currentP < minPressure) {
          break;
        }

        qTemps.add(currentT);
        qPress.add(currentP);

        // Compute volume fraction: betaV = beta * Vm_vap / (beta * Vm_vap + (1-beta) * Vm_liq)
        double mwVap = clonedSystem.getPhase(0).getMolarMass();
        double mwLiq = clonedSystem.getPhase(1).getMolarMass();
        double densVap = clonedSystem.getPhase(0).getDensity();
        double densLiq = clonedSystem.getPhase(1).getDensity();

        double volFrac = 0.5;
        if (densVap > 1e-10 && densLiq > 1e-10) {
          double vmVap = mwVap / densVap;
          double vmLiq = mwLiq / densLiq;
          volFrac = beta * vmVap / (beta * vmVap + (1.0 - beta) * vmLiq);
        }
        qVolFracs.add(volFrac);

        // Compute mass fraction: betaW = beta * Mw_vap / (beta * Mw_vap + (1-beta) * Mw_liq)
        double massFrac = beta * mwVap / (beta * mwVap + (1.0 - beta) * mwLiq);
        qMassFracs.add(massFrac);
      }

      if (!qTemps.isEmpty()) {
        String betaKey = formatBetaKey(beta);
        qualityLineData.put("qualityT_" + betaKey, toDoubleArray(qTemps));
        qualityLineData.put("qualityP_" + betaKey, toDoubleArray(qPress));
        qualityLineData.put("qualityVolFrac_" + betaKey, toDoubleArray(qVolFracs));
        qualityLineData.put("qualityMassFrac_" + betaKey, toDoubleArray(qMassFracs));
      }
    }
  }

  /**
   * Format a beta value as a clean string key (e.g. 0.5 becomes "0.5", 0.25 becomes "0.25").
   *
   * @param beta the molar vapor fraction
   * @return formatted string key
   */
  private String formatBetaKey(double beta) {
    if (beta == (int) beta) {
      return String.valueOf((int) beta);
    }
    String s = String.valueOf(beta);
    // Remove trailing zeros but keep at least one decimal
    while (s.endsWith("0") && !s.endsWith(".0")) {
      s = s.substring(0, s.length() - 1);
    }
    return s;
  }

  /**
   * Estimate the initial temperature using Wilson correlation for a given system.
   *
   * @param sys the thermodynamic system to use
   * @param beta overall vapor fraction
   * @param pressure pressure in bara
   * @return estimated temperature in Kelvin
   */
  private double tempKWilsonForSystem(SystemInterface sys, double beta, double pressure) {
    int numberOfComponents = sys.getPhase(0).getNumberOfComponents();
    int lcIdx = 0;
    int hcIdx = 0;
    double minTc = 1e10;
    double maxTc = 0;

    for (int i = 0; i < numberOfComponents; i++) {
      if (sys.getPhase(0).getComponent(i).getTC() > maxTc) {
        maxTc = sys.getPhase(0).getComponent(i).getTC();
        hcIdx = i;
      }
      if (sys.getPhase(0).getComponent(i).getTC() < minTc) {
        minTc = sys.getPhase(0).getComponent(i).getTC();
        lcIdx = i;
      }
    }

    int refIdx = beta <= 0.5 ? lcIdx : hcIdx;
    double refTc = sys.getPhase(0).getComponent(refIdx).getTC();
    double refPc = sys.getPhase(0).getComponent(refIdx).getPC();
    double refAc = sys.getPhase(0).getComponent(refIdx).getAcentricFactor();

    double lnPr = Math.log(pressure / refPc);
    double tEst = refTc * WILSON_CONST * (1 + refAc) / (WILSON_CONST * (1 + refAc) - lnPr);
    double tOld = 0;

    try {
      double[] kw = new double[numberOfComponents];
      for (int iter = 0; iter < MAX_WILSON_ITERATIONS; iter++) {
        double f = 0;
        double df = 0;
        for (int j = 0; j < numberOfComponents; j++) {
          kw[j] = sys.getPhase(0).getComponent(j).getPC() / pressure
              * Math.exp(WILSON_CONST * (1.0 + sys.getPhase(0).getComponent(j).getAcentricFactor())
                  * (1.0 - sys.getPhase(0).getComponent(j).getTC() / tEst));
        }
        for (int j = 0; j < numberOfComponents; j++) {
          double zj = sys.getPhase(0).getComponent(j).getz();
          double tc = sys.getPhase(0).getComponent(j).getTC();
          double ac = sys.getPhase(0).getComponent(j).getAcentricFactor();
          if (beta < 0.5) {
            f += zj * kw[j];
            df += zj * kw[j] * WILSON_CONST * (1 + ac) * tc / (tEst * tEst);
          } else {
            f += zj / kw[j];
            df -= zj / kw[j] * WILSON_CONST * (1 + ac) * tc / (tEst * tEst);
          }
        }
        f -= 1.0;
        if (Math.abs(f / df) > 0.1 * tEst) {
          tEst -= 0.001 * f / df;
        } else {
          tEst -= f / df;
        }
        if (Math.abs(tEst - tOld) < 1e-5) {
          return tEst;
        }
        tOld = tEst;
      }
    } catch (Exception ex) {
      lnPr = Math.log(pressure / refPc);
      tEst = refTc * WILSON_CONST * (1 + refAc) / (WILSON_CONST * (1 + refAc) - lnPr);
    }

    if (Double.isNaN(tEst) || Double.isInfinite(tEst)) {
      tEst = refTc * WILSON_CONST * (1 + refAc) / (WILSON_CONST * (1 + refAc) - Math.log(pressure / refPc));
    }
    return tEst;
  }

  /**
   * Get the list of quality line beta values that were traced.
   *
   * @return array of molar vapor fraction values
   */
  public double[] getQualityBetaValues() {
    return qualityBetaValues;
  }

  /**
   * Get data for a specific quality line.
   *
   * @param beta the molar vapor fraction of the quality line
   * @return array of [temperatures(K), pressures(bara), volumeFractions, massFractions], or null if not traced
   */
  public double[][] getQualityLine(double beta) {
    String key = formatBetaKey(beta);
    double[] qT = qualityLineData.get("qualityT_" + key);
    double[] qP = qualityLineData.get("qualityP_" + key);
    double[] qV = qualityLineData.get("qualityVolFrac_" + key);
    double[] qM = qualityLineData.get("qualityMassFrac_" + key);
    if (qT == null) {
      return null;
    }
    return new double[][] { qT, qP, qV, qM };
  }

  /**
   * Perform stability analysis along the traced phase envelope to detect potential three-phase (VLLE) regions. At each
   * sampled point on the envelope, both a multi-phase TP flash and an independent tangent-plane minimization are
   * performed. The latter is required because a conventional TP flash can miss an incipient water or hydrocarbon phase
   * before its amount becomes numerically visible.
   *
   * <p>
   * This is a post-processing step that should be called after the envelope has been traced (after {@link #run()}). It
   * detects three-phase regions such as:
   * </p>
   * <ul>
   * <li>Vapor-liquid-liquid equilibrium (VLLE) — common in water + hydrocarbon or CO2 systems</li>
   * <li>Systems approaching liquid-liquid boundaries embedded within the VLE region</li>
   * </ul>
   *
   * <p>
   * Results are accessible via {@link #get(String)} with keys {@code "threePhaseT"} and {@code "threePhaseP"}, or via
   * {@link #getThreePhaseRegionT()} and {@link #getThreePhaseRegionP()}.
   * </p>
   *
   * <p>
   * Reference: Cismondi &amp; Michelsen, "Global calculation of phase equilibrium", Fluid Phase Equilibria, 259,
   * 228-234 (2007).
   * </p>
   */
  public void checkStabilityAlongEnvelope() {
    ArrayList<Double> threePhaseTemps = new ArrayList<Double>();
    ArrayList<Double> threePhasePress = new ArrayList<Double>();
    ArrayList<SecondaryStabilitySample> samples = new ArrayList<SecondaryStabilitySample>();

    int sampleStep = Math.max(1, boundaryStateSeeds.size() / 100);
    for (int seedIndex = 0; seedIndex < boundaryStateSeeds.size(); seedIndex += sampleStep) {
      BoundaryStateSeed seed = boundaryStateSeeds.get(seedIndex);
      if (seed == null) {
        continue;
      }
      SecondaryStabilitySample sample = checkSecondaryStability(seed);
      samples.add(sample);
      if (sample.indicatesThreePhaseRegion()) {
        threePhaseTemps.add(seed.temperature);
        threePhasePress.add(seed.pressure);
      }
    }

    threePhaseRegionT = toDoubleArray(threePhaseTemps);
    threePhaseRegionP = toDoubleArray(threePhasePress);
    secondaryStabilitySamples = Collections.unmodifiableList(samples);
    secondaryStabilityBrackets = buildSecondaryStabilityBrackets(samples);
    threePhasePointCandidates = Collections.emptyList();
  }

  private List<SecondaryStabilityBracket> buildSecondaryStabilityBrackets(List<SecondaryStabilitySample> samples) {
    ArrayList<SecondaryStabilityBracket> brackets = new ArrayList<SecondaryStabilityBracket>();
    for (int sampleIndex = 1; sampleIndex < samples.size(); sampleIndex++) {
      SecondaryStabilitySample first = samples.get(sampleIndex - 1);
      SecondaryStabilitySample second = samples.get(sampleIndex);
      boolean contiguous = first.getBranch() == second.getBranch()
          && Math.abs(first.getTemperature() - second.getTemperature()) <= JUMP_BREAK_FACTOR * dTmax
          && Math.abs(first.getPressure() - second.getPressure()) <= JUMP_BREAK_FACTOR * dPmax;
      if (contiguous && first.getFailureMessage() == null && second.getFailureMessage() == null
          && first.isPrimaryBoundaryValid() && second.isPrimaryBoundaryValid()
          && first.getIncipientPhase() == IncipientPhaseStabilityAnalyzer.CandidatePhase.AQUEOUS
          && second.getIncipientPhase() == IncipientPhaseStabilityAnalyzer.CandidatePhase.AQUEOUS
          && Double.isFinite(first.getStabilityFunction()) && Double.isFinite(second.getStabilityFunction())
          && first.indicatesThreePhaseRegion() != second.indicatesThreePhaseRegion()) {
        brackets.add(new SecondaryStabilityBracket(first, second));
      }
    }
    return Collections.unmodifiableList(brackets);
  }

  private SecondaryStabilitySample checkSecondaryStability(BoundaryStateSeed seed) {
    int equilibriumPhaseCount = checkPhaseCount(seed.temperature, seed.pressure);
    try {
      if (seed.motherPhase == IncipientPhaseStabilityAnalyzer.CandidatePhase.AQUEOUS
          || seed.primaryIncipientPhase == IncipientPhaseStabilityAnalyzer.CandidatePhase.AQUEOUS
          || seed.motherPhase == seed.primaryIncipientPhase) {
        return SecondaryStabilitySample.failure(seed.branch, seed.temperature, seed.pressure, equilibriumPhaseCount,
            "conventional hydrocarbon boundary does not contain distinct gas and oil phase families");
      }
      SystemInterface mother = singleMotherPhase(seed);
      IncipientPhaseStabilityAnalyzer analyzer = new IncipientPhaseStabilityAnalyzer(mother).setMaximumIterations(300)
          .setDampingFactor(0.35);
      IncipientPhaseStabilityAnalyzer.Candidate primary = analyzer.analyzeCandidate(seed.primaryIncipientPhase,
          seed.primaryComposition);
      IncipientPhaseStabilityAnalyzer.Candidate aqueous = analyzer
          .analyzeCandidate(IncipientPhaseStabilityAnalyzer.CandidatePhase.AQUEOUS, aqueousSeed(mother));
      return SecondaryStabilitySample.success(seed.branch, seed.temperature, seed.pressure, equilibriumPhaseCount,
          seed.motherPhase, seed.primaryIncipientPhase, primary.getStabilityFunction(), primary.isConverged(),
          aqueous.getPhase(), aqueous.getStabilityFunction(), aqueous.getTangentPlaneDistance(),
          aqueous.getStationarityResidual(), seed.beta, seed.logK);
    } catch (RuntimeException error) {
      return SecondaryStabilitySample.failure(seed.branch, seed.temperature, seed.pressure, equilibriumPhaseCount,
          error.getMessage());
    }
  }

  private SystemInterface singleMotherPhase(BoundaryStateSeed seed) {
    SystemInterface mother = system.clone();
    mother.setNumberOfPhases(1);
    mother.setTemperature(seed.temperature);
    mother.setPressure(seed.pressure);
    mother.setPhaseType(0,
        seed.motherPhase == IncipientPhaseStabilityAnalyzer.CandidatePhase.GAS ? neqsim.thermo.phase.PhaseType.GAS
            : neqsim.thermo.phase.PhaseType.OIL);
    for (int componentIndex = 0; componentIndex < mother.getPhase(0).getNumberOfComponents(); componentIndex++) {
      mother.getPhase(0).getComponent(componentIndex)
          .setx(Math.max(mother.getPhase(0).getComponent(componentIndex).getz(), 1.0e-100));
    }
    mother.getPhase(0).normalize();
    mother.init(1);
    return mother;
  }

  private static double[] aqueousSeed(SystemInterface mother) {
    double[] composition = new double[mother.getPhase(0).getNumberOfComponents()];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      neqsim.thermo.component.ComponentInterface component = mother.getPhase(0).getComponent(componentIndex);
      if (component.getComponentName().equalsIgnoreCase("water")) {
        composition[componentIndex] = Math.max(component.getz(), 0.99);
      } else if (component.isHydrocarbon()) {
        composition[componentIndex] = Math.max(component.getz() * 1.0e-8, 1.0e-100);
      } else {
        composition[componentIndex] = Math.max(component.getz(), 1.0e-100);
      }
      total += composition[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }

  private static double[] logKVector(SystemInterface boundaryState) {
    double[] logK = new double[boundaryState.getPhase(0).getNumberOfComponents()];
    for (int componentIndex = 0; componentIndex < logK.length; componentIndex++) {
      double k = boundaryState.getPhase(0).getComponent(componentIndex).getK();
      logK[componentIndex] = k > 0.0 && Double.isFinite(k) ? Math.log(k) : Double.NaN;
    }
    return logK;
  }

  private SecondaryStabilitySample checkSecondaryStabilityAtFixedTemperature(EnvelopeSegment.PhaseType branch,
      double temperature, double pressureGuess, double beta, double[] initialLogK) {
    try {
      SystemInterface restricted = system.clone();
      restricted.setMultiPhaseCheck(false);
      restricted.setMaxNumberOfPhases(2);
      if (restricted.getNumberOfPhases() != 2) {
        restricted.setNumberOfPhases(2);
      }
      restricted.setTemperature(temperature);
      restricted.setPressure(pressureGuess);
      restricted.setBeta(Math.max(1.0e-10, Math.min(1.0 - 1.0e-10, beta)));
      if (initialLogK.length != restricted.getPhase(0).getNumberOfComponents()) {
        throw new IllegalArgumentException("initial log(K) vector has wrong component count");
      }
      for (int componentIndex = 0; componentIndex < initialLogK.length; componentIndex++) {
        double k = Math.exp(initialLogK[componentIndex]);
        restricted.getPhase(0).getComponent(componentIndex).setK(k);
        restricted.getPhase(1).getComponent(componentIndex).setK(k);
      }
      SysNewtonRhapsonPhaseEnvelope corrector = new SysNewtonRhapsonPhaseEnvelope(restricted, 2,
          restricted.getPhase(0).getNumberOfComponents());
      if (!corrector.solveAtFixedTemperature(temperature, 40, 1.0e-8)) {
        return SecondaryStabilitySample.failure(branch, temperature, pressureGuess, -1,
            "fixed-temperature two-phase corrector did not converge");
      }
      restricted.init(1);
      double correctedPressure = restricted.getPressure();
      double[] correctedLogK = logKVector(restricted);
      double maximumAbsoluteLogK = 0.0;
      for (double value : correctedLogK) {
        maximumAbsoluteLogK = Math.max(maximumAbsoluteLogK, Math.abs(value));
      }
      double maximumLocalPressureChange = Math.max(2.0 * dPmax, 0.35 * pressureGuess);
      if (maximumAbsoluteLogK < 1.0e-4) {
        return SecondaryStabilitySample.failure(branch, temperature, pressureGuess, -1,
            "fixed-temperature corrector reached the trivial K=1 solution");
      }
      if (!Double.isFinite(correctedPressure) || correctedPressure <= 0.0
          || Math.abs(correctedPressure - pressureGuess) > maximumLocalPressureChange) {
        return SecondaryStabilitySample.failure(branch, temperature, pressureGuess, -1,
            "fixed-temperature corrector left the local pressure bracket");
      }
      return checkSecondaryStability(BoundaryStateSeed.capture(branch, restricted));
    } catch (RuntimeException error) {
      return SecondaryStabilitySample.failure(branch, temperature, pressureGuess, -1, error.getMessage());
    }
  }

  /**
   * Check how many equilibrium phases exist at a given (T, P) for the system's overall composition. Uses a multi-phase
   * TP flash with stability analysis enabled.
   *
   * @param temperature temperature in Kelvin
   * @param pressure pressure in bara
   * @return number of equilibrium phases (1, 2, or 3+), or -1 if the flash failed
   */
  public int checkPhaseCount(double temperature, double pressure) {
    try {
      SystemInterface testSys = system.clone();
      testSys.setTemperature(temperature);
      testSys.setPressure(pressure);
      testSys.setMultiPhaseCheck(true);
      ThermodynamicOperations testOps = new ThermodynamicOperations(testSys);
      testOps.TPflash();
      return testSys.getNumberOfPhases();
    } catch (Exception e) {
      return -1;
    }
  }

  /**
   * Get temperatures of envelope points where three or more phases were detected. Requires
   * {@link #checkStabilityAlongEnvelope()} to have been called first.
   *
   * @return array of temperatures (K) at three-phase points, empty if none found
   */
  public double[] getThreePhaseRegionT() {
    return threePhaseRegionT;
  }

  /**
   * Get pressures of envelope points where three or more phases were detected. Requires
   * {@link #checkStabilityAlongEnvelope()} to have been called first.
   *
   * @return array of pressures (bara) at three-phase points, empty if none found
   */
  public double[] getThreePhaseRegionP() {
    return threePhaseRegionP;
  }

  /**
   * Returns tangent-plane stability samples collected along the conventional two-phase boundary.
   *
   * @return immutable samples, empty until {@link #checkStabilityAlongEnvelope()} is called
   */
  public List<SecondaryStabilitySample> getSecondaryStabilitySamples() {
    return secondaryStabilitySamples;
  }

  /** @return adjacent stability samples that bracket a two-phase/three-phase topology transition */
  public List<SecondaryStabilityBracket> getSecondaryStabilityBrackets() {
    return secondaryStabilityBrackets;
  }

  /**
   * Refines every coarse topology bracket while projecting each trial back onto the two-phase boundary.
   *
   * @param maximumBisections maximum fixed-temperature bisections per bracket
   * @param temperatureToleranceK terminal temperature width in kelvin
   * @param pressureToleranceBara terminal pressure width in bara
   * @return immutable refined transition candidates
   */
  public List<ThreePhasePointCandidate> refineSecondaryStabilityBrackets(int maximumBisections,
      double temperatureToleranceK, double pressureToleranceBara) {
    if (maximumBisections < 1 || !Double.isFinite(temperatureToleranceK) || temperatureToleranceK <= 0.0
        || !Double.isFinite(pressureToleranceBara) || pressureToleranceBara <= 0.0) {
      throw new IllegalArgumentException("invalid three-phase bracket refinement configuration");
    }
    ArrayList<ThreePhasePointCandidate> refined = new ArrayList<ThreePhasePointCandidate>();
    for (SecondaryStabilityBracket bracket : secondaryStabilityBrackets) {
      SecondaryStabilitySample twoPhase = bracket.getTwoPhaseSide();
      SecondaryStabilitySample threePhase = bracket.getThreePhaseSide();
      String failureMessage = null;
      int iterations = 0;
      while (iterations < maximumBisections
          && (Math.abs(twoPhase.getTemperature() - threePhase.getTemperature()) > temperatureToleranceK
              || Math.abs(twoPhase.getPressure() - threePhase.getPressure()) > pressureToleranceBara)) {
        iterations++;
        double temperature = 0.5 * (twoPhase.getTemperature() + threePhase.getTemperature());
        double pressureGuess = 0.5 * (twoPhase.getPressure() + threePhase.getPressure());
        SecondaryStabilitySample midpoint = refineBracketMidpoint(twoPhase, threePhase, temperature, pressureGuess);
        if (midpoint.getFailureMessage() != null) {
          failureMessage = midpoint.getFailureMessage();
          break;
        }
        if (midpoint.indicatesThreePhaseRegion()) {
          threePhase = midpoint;
        } else {
          twoPhase = midpoint;
        }
      }
      boolean toleranceSatisfied = Math
          .abs(twoPhase.getTemperature() - threePhase.getTemperature()) <= temperatureToleranceK
          && Math.abs(twoPhase.getPressure() - threePhase.getPressure()) <= pressureToleranceBara;
      refined.add(new ThreePhasePointCandidate(twoPhase, threePhase, iterations, toleranceSatisfied, failureMessage));
    }
    threePhasePointCandidates = Collections.unmodifiableList(refined);
    return threePhasePointCandidates;
  }

  private SecondaryStabilitySample refineBracketMidpoint(SecondaryStabilitySample twoPhase,
      SecondaryStabilitySample threePhase, double temperature, double pressureGuess) {
    double branchBeta = twoPhase.getBranch() == EnvelopeSegment.PhaseType.DEW ? 1.0 - 1.0e-10 : 1.0e-10;
    double[][] logKSeeds = new double[][] { interpolateLogK(twoPhase.getLogK(), threePhase.getLogK()),
        twoPhase.getLogK(), threePhase.getLogK() };
    double[] betaSeeds = new double[] { 0.5 * (twoPhase.getBeta() + threePhase.getBeta()), branchBeta,
        1.0 - branchBeta };
    SecondaryStabilitySample lastFailure = null;
    for (double[] logKSeed : logKSeeds) {
      for (double betaSeed : betaSeeds) {
        SecondaryStabilitySample attempt = checkSecondaryStabilityAtFixedTemperature(twoPhase.getBranch(), temperature,
            pressureGuess, betaSeed, logKSeed);
        if (attempt.getFailureMessage() == null) {
          return attempt;
        }
        lastFailure = attempt;
      }
    }
    return lastFailure == null ? SecondaryStabilitySample.failure(twoPhase.getBranch(), temperature, pressureGuess, -1,
        "no fixed-temperature corrector seed was available") : lastFailure;
  }

  private static double[] interpolateLogK(double[] first, double[] second) {
    if (first.length == 0 || first.length != second.length) {
      throw new IllegalArgumentException("cannot interpolate incompatible log(K) vectors");
    }
    double[] midpoint = new double[first.length];
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      midpoint[componentIndex] = 0.5 * (first[componentIndex] + second[componentIndex]);
    }
    return midpoint;
  }

  /** @return most recently refined transition candidates */
  public List<ThreePhasePointCandidate> getThreePhasePointCandidates() {
    return threePhasePointCandidates;
  }

  /** A coarse bracket around one entry to or exit from a three-phase region. */
  public static final class SecondaryStabilityBracket implements java.io.Serializable {
    private static final long serialVersionUID = 1L;
    private final SecondaryStabilitySample first;
    private final SecondaryStabilitySample second;

    private SecondaryStabilityBracket(SecondaryStabilitySample first, SecondaryStabilitySample second) {
      this.first = first;
      this.second = second;
    }

    /** @return first adjacent sample */
    public SecondaryStabilitySample getFirst() {
      return first;
    }

    /** @return second adjacent sample */
    public SecondaryStabilitySample getSecond() {
      return second;
    }

    /** @return sample on the conventional two-phase side */
    public SecondaryStabilitySample getTwoPhaseSide() {
      return first.indicatesThreePhaseRegion() ? second : first;
    }

    /** @return sample on the detected three-phase side */
    public SecondaryStabilitySample getThreePhaseSide() {
      return first.indicatesThreePhaseRegion() ? first : second;
    }
  }

  /** One refined candidate where a two-phase boundary enters or exits a three-phase region. */
  public static final class ThreePhasePointCandidate implements java.io.Serializable {
    private static final long serialVersionUID = 1L;
    private final SecondaryStabilitySample twoPhaseSide;
    private final SecondaryStabilitySample threePhaseSide;
    private final int iterations;
    private final boolean toleranceSatisfied;
    private final String failureMessage;

    private ThreePhasePointCandidate(SecondaryStabilitySample twoPhaseSide, SecondaryStabilitySample threePhaseSide,
        int iterations, boolean toleranceSatisfied, String failureMessage) {
      this.twoPhaseSide = twoPhaseSide;
      this.threePhaseSide = threePhaseSide;
      this.iterations = iterations;
      this.toleranceSatisfied = toleranceSatisfied;
      this.failureMessage = failureMessage;
    }

    /** @return midpoint temperature estimate in kelvin */
    public double getTemperature() {
      return 0.5 * (twoPhaseSide.getTemperature() + threePhaseSide.getTemperature());
    }

    /** @return midpoint pressure estimate in bara */
    public double getPressure() {
      return 0.5 * (twoPhaseSide.getPressure() + threePhaseSide.getPressure());
    }

    /** @return final temperature bracket width in kelvin */
    public double getTemperatureWidth() {
      return Math.abs(twoPhaseSide.getTemperature() - threePhaseSide.getTemperature());
    }

    /** @return final pressure bracket width in bara */
    public double getPressureWidth() {
      return Math.abs(twoPhaseSide.getPressure() - threePhaseSide.getPressure());
    }

    /** @return incipient phase on the detected three-phase side */
    public IncipientPhaseStabilityAnalyzer.CandidatePhase getIncipientPhase() {
      return threePhaseSide.getIncipientPhase();
    }

    /** @return number of bisection/correction attempts */
    public int getIterations() {
      return iterations;
    }

    /** @return true when no correction failed and both requested bracket tolerances were reached */
    public boolean isConverged() {
      return failureMessage == null && toleranceSatisfied;
    }

    /** @return true when both requested temperature and pressure widths were reached */
    public boolean isToleranceSatisfied() {
      return toleranceSatisfied;
    }

    /** @return failure diagnostic, or {@code null} */
    public String getFailureMessage() {
      return failureMessage;
    }

    /** @return last point classified on the conventional two-phase side */
    public SecondaryStabilitySample getTwoPhaseSide() {
      return twoPhaseSide;
    }

    /** @return last point classified on the three-phase side */
    public SecondaryStabilitySample getThreePhaseSide() {
      return threePhaseSide;
    }
  }

  /** One non-destructive stability diagnostic at a traced PT boundary point. */
  public static final class SecondaryStabilitySample implements java.io.Serializable {
    private static final long serialVersionUID = 1L;
    private static final double INSTABILITY_TOLERANCE = -1.0e-8;

    private final EnvelopeSegment.PhaseType branch;
    private final double temperature;
    private final double pressure;
    private final int equilibriumPhaseCount;
    private final IncipientPhaseStabilityAnalyzer.CandidatePhase motherPhase;
    private final IncipientPhaseStabilityAnalyzer.CandidatePhase primaryIncipientPhase;
    private final double primaryBoundaryStabilityFunction;
    private final boolean primaryBoundaryConverged;
    private final IncipientPhaseStabilityAnalyzer.CandidatePhase incipientPhase;
    private final double stabilityFunction;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;
    private final double beta;
    private final double[] logK;
    private final String failureMessage;

    private SecondaryStabilitySample(EnvelopeSegment.PhaseType branch, double temperature, double pressure,
        int equilibriumPhaseCount, IncipientPhaseStabilityAnalyzer.CandidatePhase motherPhase,
        IncipientPhaseStabilityAnalyzer.CandidatePhase primaryIncipientPhase, double primaryBoundaryStabilityFunction,
        boolean primaryBoundaryConverged, IncipientPhaseStabilityAnalyzer.CandidatePhase incipientPhase,
        double stabilityFunction, double tangentPlaneDistance, double stationarityResidual, double beta, double[] logK,
        String failureMessage) {
      this.branch = branch;
      this.temperature = temperature;
      this.pressure = pressure;
      this.equilibriumPhaseCount = equilibriumPhaseCount;
      this.motherPhase = motherPhase;
      this.primaryIncipientPhase = primaryIncipientPhase;
      this.primaryBoundaryStabilityFunction = primaryBoundaryStabilityFunction;
      this.primaryBoundaryConverged = primaryBoundaryConverged;
      this.incipientPhase = incipientPhase;
      this.stabilityFunction = stabilityFunction;
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.stationarityResidual = stationarityResidual;
      this.beta = beta;
      this.logK = logK.clone();
      this.failureMessage = failureMessage;
    }

    private static SecondaryStabilitySample success(EnvelopeSegment.PhaseType branch, double temperature,
        double pressure, int equilibriumPhaseCount, IncipientPhaseStabilityAnalyzer.CandidatePhase motherPhase,
        IncipientPhaseStabilityAnalyzer.CandidatePhase primaryIncipientPhase, double primaryBoundaryStabilityFunction,
        boolean primaryBoundaryConverged, IncipientPhaseStabilityAnalyzer.CandidatePhase incipientPhase,
        double stabilityFunction, double tangentPlaneDistance, double stationarityResidual, double beta,
        double[] logK) {
      return new SecondaryStabilitySample(branch, temperature, pressure, equilibriumPhaseCount, motherPhase,
          primaryIncipientPhase, primaryBoundaryStabilityFunction, primaryBoundaryConverged, incipientPhase,
          stabilityFunction, tangentPlaneDistance, stationarityResidual, beta, logK, null);
    }

    private static SecondaryStabilitySample failure(EnvelopeSegment.PhaseType branch, double temperature,
        double pressure, int equilibriumPhaseCount, String failureMessage) {
      return new SecondaryStabilitySample(branch, temperature, pressure, equilibriumPhaseCount, null, null, Double.NaN,
          false, null, Double.NaN, Double.NaN, Double.NaN, Double.NaN, new double[0], failureMessage);
    }

    /** @return dew or bubble branch */
    public EnvelopeSegment.PhaseType getBranch() {
      return branch;
    }

    /** @return temperature in K */
    public double getTemperature() {
      return temperature;
    }

    /** @return pressure in bara */
    public double getPressure() {
      return pressure;
    }

    /** @return stable multiphase TP-flash phase count, or {@code -1} on failure */
    public int getEquilibriumPhaseCount() {
      return equilibriumPhaseCount;
    }

    /** @return fixed-composition mother phase used for the TPD reference */
    public IncipientPhaseStabilityAnalyzer.CandidatePhase getMotherPhase() {
      return motherPhase;
    }

    /** @return hydrocarbon phase already known to be incipient on the conventional boundary */
    public IncipientPhaseStabilityAnalyzer.CandidatePhase getPrimaryIncipientPhase() {
      return primaryIncipientPhase;
    }

    /** @return stability function of the primary hydrocarbon incipient phase */
    public double getPrimaryBoundaryStabilityFunction() {
      return primaryBoundaryStabilityFunction;
    }

    /** @return true when the primary hydrocarbon stationary iteration converged */
    public boolean isPrimaryBoundaryConverged() {
      return primaryBoundaryConverged;
    }

    /** @return true when the saved continuation point still satisfies the primary zero-TPD boundary */
    public boolean isPrimaryBoundaryValid() {
      return primaryBoundaryConverged && Double.isFinite(primaryBoundaryStabilityFunction)
          && Math.abs(primaryBoundaryStabilityFunction) <= 1.0e-4;
    }

    /** @return lowest-TPD non-trivial candidate phase, or {@code null} */
    public IncipientPhaseStabilityAnalyzer.CandidatePhase getIncipientPhase() {
      return incipientPhase;
    }

    /** @return Michelsen stability function {@code 1 - sum(W)} */
    public double getStabilityFunction() {
      return stabilityFunction;
    }

    /** @return tangent-plane distance at the stationary trial */
    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    /** @return aqueous stationary-point logarithmic residual */
    public double getStationarityResidual() {
      return stationarityResidual;
    }

    /** @return phase-0 mole fraction used by the reconstructed two-phase boundary state */
    public double getBeta() {
      return beta;
    }

    /** @return logarithmic K-values of the reconstructed two-phase boundary state */
    public double[] getLogK() {
      return logK.clone();
    }

    /** @return failure diagnostic, or {@code null} */
    public String getFailureMessage() {
      return failureMessage;
    }

    /** @return true when a valid primary boundary has negative aqueous TPD */
    public boolean indicatesThreePhaseRegion() {
      return failureMessage == null && isPrimaryBoundaryValid()
          && incipientPhase == IncipientPhaseStabilityAnalyzer.CandidatePhase.AQUEOUS
          && Double.isFinite(stabilityFunction) && stabilityFunction < INSTABILITY_TOLERANCE;
    }
  }

  /** Exact two-phase continuation state used to reproduce a primary incipient boundary without another TP flash. */
  private static final class BoundaryStateSeed {
    private final EnvelopeSegment.PhaseType branch;
    private final double temperature;
    private final double pressure;
    private final double beta;
    private final IncipientPhaseStabilityAnalyzer.CandidatePhase motherPhase;
    private final IncipientPhaseStabilityAnalyzer.CandidatePhase primaryIncipientPhase;
    private final double[] primaryComposition;
    private final double[] logK;

    private BoundaryStateSeed(EnvelopeSegment.PhaseType branch, double temperature, double pressure, double beta,
        IncipientPhaseStabilityAnalyzer.CandidatePhase motherPhase,
        IncipientPhaseStabilityAnalyzer.CandidatePhase primaryIncipientPhase, double[] primaryComposition,
        double[] logK) {
      this.branch = branch;
      this.temperature = temperature;
      this.pressure = pressure;
      this.beta = beta;
      this.motherPhase = motherPhase;
      this.primaryIncipientPhase = primaryIncipientPhase;
      this.primaryComposition = primaryComposition;
      this.logK = logK;
    }

    private static BoundaryStateSeed capture(EnvelopeSegment.PhaseType branch, SystemInterface boundaryState) {
      int componentCount = boundaryState.getPhase(0).getNumberOfComponents();
      double[] overall = new double[componentCount];
      double[] first = new double[componentCount];
      double[] second = new double[componentCount];
      double firstDistance = 0.0;
      double secondDistance = 0.0;
      for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
        overall[componentIndex] = boundaryState.getPhase(0).getComponent(componentIndex).getz();
        first[componentIndex] = boundaryState.getPhase(0).getComponent(componentIndex).getx();
        second[componentIndex] = boundaryState.getPhase(1).getComponent(componentIndex).getx();
        firstDistance += Math.abs(first[componentIndex] - overall[componentIndex]);
        secondDistance += Math.abs(second[componentIndex] - overall[componentIndex]);
      }
      boolean firstIsMother = firstDistance <= secondDistance;
      int motherIndex = firstIsMother ? 0 : 1;
      int incipientIndex = firstIsMother ? 1 : 0;
      return new BoundaryStateSeed(branch, boundaryState.getTemperature(), boundaryState.getPressure(),
          boundaryState.getBeta(), phaseFamily(boundaryState.getPhase(motherIndex).getType()),
          phaseFamily(boundaryState.getPhase(incipientIndex).getType()), firstIsMother ? second : first,
          logKVector(boundaryState));
    }

    private static IncipientPhaseStabilityAnalyzer.CandidatePhase phaseFamily(neqsim.thermo.phase.PhaseType phaseType) {
      if (phaseType == neqsim.thermo.phase.PhaseType.GAS) {
        return IncipientPhaseStabilityAnalyzer.CandidatePhase.GAS;
      }
      if (phaseType == neqsim.thermo.phase.PhaseType.AQUEOUS) {
        return IncipientPhaseStabilityAnalyzer.CandidatePhase.AQUEOUS;
      }
      return IncipientPhaseStabilityAnalyzer.CandidatePhase.OIL;
    }
  }

  /**
   * Reset K-values on the system using the Wilson correlation at the configured low pressure. Called before the restart
   * pass to ensure fresh initial estimates after a crash.
   */
  private void resetKValuesWithWilson() {
    double restartTemp = tempKWilson(phaseFraction, lowPres);
    if (Double.isNaN(restartTemp)) {
      restartTemp = system.getPhase(0).getComponent(speceq).getTC() - 20.0;
    }
    for (int ic = 0; ic < system.getPhase(0).getNumberOfComponents(); ic++) {
      double Kwil = system.getPhase(0).getComponent(ic).getPC() / lowPres
          * Math.exp(WILSON_CONST * (1.0 + system.getPhase(0).getComponent(ic).getAcentricFactor())
              * (1.0 - system.getPhase(0).getComponent(ic).getTC() / restartTemp));
      system.getPhase(0).getComponent(ic).setK(Kwil);
      system.getPhase(1).getComponent(ic).setK(Kwil);
    }
  }

  /**
   * Insert a NaN sentinel into every per-point list so that plotting libraries render the dew and bubble curves as
   * disjoint polylines across branch transitions (primary-to-restart pass and critical-point crossings). Without this
   * break plotters draw a spurious straight line across the two-phase region.
   */
  private void addBranchBreak() {
    dewPointTemperatures.add(Double.NaN);
    dewPointPressures.add(Double.NaN);
    dewPointEnthalpies.add(Double.NaN);
    dewPointDensities.add(Double.NaN);
    dewPointEntropies.add(Double.NaN);
    bubblePointTemperatures.add(Double.NaN);
    bubblePointPressures.add(Double.NaN);
    bubblePointEnthalpies.add(Double.NaN);
    bubblePointDensities.add(Double.NaN);
    bubblePointEntropies.add(Double.NaN);
    boundaryStateSeeds.add(null);
  }

  /**
   * Convert all ArrayList results to primitive arrays for efficient access via get(), and build the structured
   * per-branch segment list.
   */
  private void buildOutputArrays() {
    dewTempArray = toDoubleArray(dewPointTemperatures);
    dewPresArray = toDoubleArray(dewPointPressures);
    dewEnthalpyArray = toDoubleArray(dewPointEnthalpies);
    dewDensityArray = toDoubleArray(dewPointDensities);
    dewEntropyArray = toDoubleArray(dewPointEntropies);

    bubTempArray = toDoubleArray(bubblePointTemperatures);
    bubPresArray = toDoubleArray(bubblePointPressures);
    bubEnthalpyArray = toDoubleArray(bubblePointEnthalpies);
    bubDensityArray = toDoubleArray(bubblePointDensities);
    bubEntropyArray = toDoubleArray(bubblePointEntropies);

    // Re-break each branch wherever consecutive points jump much further than a single
    // continuation step can travel. The tracer already inserts NaN sentinels at
    // critical-point crossings and restart passes, but a disjoint sub-branch (e.g. a
    // near-critical segment separated from the low-pressure tail) can still leave two far
    // apart points adjacent in the same list; without a break a plotter draws a straight
    // chord across the gap. Inserting NaN here makes get()/getSegments() robust regardless
    // of how the caller renders the curve.
    insertJumpBreaks();

    segments = buildSegments();
  }

  /**
   * Insert NaN break sentinels into both branches wherever consecutive non-NaN points are separated by more than
   * {@link #JUMP_BREAK_FACTOR} times the maximum continuation step. This guarantees that disjoint sub-branches are
   * rendered as separate polylines even when the tracer stored them adjacently in the same list.
   */
  private void insertJumpBreaks() {
    double dTbreak = JUMP_BREAK_FACTOR * dTmax;
    double dPbreak = JUMP_BREAK_FACTOR * dPmax;

    double[][] dew = insertJumpBreaks(dTbreak, dPbreak, dewTempArray, dewPresArray, dewEnthalpyArray, dewDensityArray,
        dewEntropyArray);
    dewTempArray = dew[0];
    dewPresArray = dew[1];
    dewEnthalpyArray = dew[2];
    dewDensityArray = dew[3];
    dewEntropyArray = dew[4];

    double[][] bub = insertJumpBreaks(dTbreak, dPbreak, bubTempArray, bubPresArray, bubEnthalpyArray, bubDensityArray,
        bubEntropyArray);
    bubTempArray = bub[0];
    bubPresArray = bub[1];
    bubEnthalpyArray = bub[2];
    bubDensityArray = bub[3];
    bubEntropyArray = bub[4];
  }

  /**
   * Copy a single branch's parallel arrays, inserting a NaN sentinel into every array at the same index wherever the
   * temperature or pressure step between two consecutive non-NaN points exceeds the supplied thresholds.
   *
   * @param dTbreak temperature jump (K) above which a break is inserted
   * @param dPbreak pressure jump (bara) above which a break is inserted
   * @param T temperatures (may already contain NaN sentinels)
   * @param P pressures, same length as T
   * @param H mass enthalpies, same length as T
   * @param D mass densities, same length as T
   * @param S mass entropies, same length as T
   * @return five new arrays {T, P, H, D, S} with NaN breaks inserted at large discontinuities
   */
  private double[][] insertJumpBreaks(double dTbreak, double dPbreak, double[] T, double[] P, double[] H, double[] D,
      double[] S) {
    int n = T.length;
    ArrayList<Double> oT = new ArrayList<Double>(n);
    ArrayList<Double> oP = new ArrayList<Double>(n);
    ArrayList<Double> oH = new ArrayList<Double>(n);
    ArrayList<Double> oD = new ArrayList<Double>(n);
    ArrayList<Double> oS = new ArrayList<Double>(n);
    for (int i = 0; i < n; i++) {
      if (i > 0 && !Double.isNaN(T[i]) && !Double.isNaN(T[i - 1])
          && (Math.abs(T[i] - T[i - 1]) > dTbreak || Math.abs(P[i] - P[i - 1]) > dPbreak)) {
        oT.add(Double.NaN);
        oP.add(Double.NaN);
        oH.add(Double.NaN);
        oD.add(Double.NaN);
        oS.add(Double.NaN);
      }
      oT.add(T[i]);
      oP.add(P[i]);
      oH.add(H[i]);
      oD.add(D[i]);
      oS.add(S[i]);
    }
    return new double[][] { toDoubleArray(oT), toDoubleArray(oP), toDoubleArray(oH), toDoubleArray(oD),
        toDoubleArray(oS) };
  }

  /**
   * Split the NaN-delimited per-point arrays into contiguous branch segments.
   *
   * @return unmodifiable list of dew and bubble segments in the order they were traced
   */
  private List<EnvelopeSegment> buildSegments() {
    ArrayList<EnvelopeSegment> out = new ArrayList<EnvelopeSegment>();
    extractSegmentsInto(out, EnvelopeSegment.PhaseType.DEW, dewTempArray, dewPresArray, dewEnthalpyArray,
        dewDensityArray, dewEntropyArray);
    extractSegmentsInto(out, EnvelopeSegment.PhaseType.BUBBLE, bubTempArray, bubPresArray, bubEnthalpyArray,
        bubDensityArray, bubEntropyArray);
    return Collections.unmodifiableList(out);
  }

  /**
   * Append every contiguous non-NaN run of a flat branch array as a segment.
   *
   * @param out list to append to
   * @param type phase type for the segments being extracted
   * @param T temperatures (possibly containing NaN break sentinels)
   * @param P pressures, same length as T
   * @param H mass enthalpies, same length as T
   * @param D mass densities, same length as T
   * @param S mass entropies, same length as T
   */
  private static void extractSegmentsInto(ArrayList<EnvelopeSegment> out, EnvelopeSegment.PhaseType type, double[] T,
      double[] P, double[] H, double[] D, double[] S) {
    int i = 0;
    int n = T.length;
    while (i < n) {
      // skip leading NaN sentinels
      while (i < n && Double.isNaN(T[i])) {
        i++;
      }
      int start = i;
      while (i < n && !Double.isNaN(T[i])) {
        i++;
      }
      int len = i - start;
      if (len <= 0) {
        continue;
      }
      double[] segT = new double[len];
      double[] segP = new double[len];
      double[] segH = new double[len];
      double[] segD = new double[len];
      double[] segS = new double[len];
      System.arraycopy(T, start, segT, 0, len);
      System.arraycopy(P, start, segP, 0, len);
      System.arraycopy(H, start, segH, 0, len);
      System.arraycopy(D, start, segD, 0, len);
      System.arraycopy(S, start, segS, 0, len);
      out.add(new EnvelopeSegment(type, segT, segP, segH, segD, segS));
    }
  }

  /**
   * Convert an ArrayList of Double to a primitive double array.
   *
   * @param list the list to convert
   * @return primitive double array
   */
  private double[] toDoubleArray(ArrayList<Double> list) {
    double[] result = new double[list.size()];
    for (int i = 0; i < list.size(); i++) {
      result[i] = list.get(i).doubleValue();
    }
    return result;
  }

  /**
   * Estimate the initial temperature using the Wilson correlation.
   *
   * <p>
   * Iteratively solves for the temperature at which the Wilson K-value correlation satisfies the Rachford-Rice equation
   * for the given phase fraction and pressure.
   * </p>
   *
   * @param beta overall vapor fraction
   * @param P pressure in bara
   * @return estimated temperature in Kelvin
   */
  private double tempKWilson(double beta, double P) {
    int numberOfComponents = system.getPhase(0).getNumberOfComponents();
    int lc = 0;
    int hc = 0;
    double min = 1e10;
    double max = 0;

    for (int i = 0; i < numberOfComponents; i++) {
      if (system.getPhase(0).getComponent(i).getTC() > max) {
        max = system.getPhase(0).getComponent(i).getTC();
        hc = i;
      }
      if (system.getPhase(0).getComponent(i).getTC() < min) {
        min = system.getPhase(0).getComponent(i).getTC();
        lc = i;
      }
    }

    double initTc;
    double initPc;
    double initAc;
    if (beta <= 0.5) {
      initTc = system.getPhase(0).getComponent(lc).getTC();
      initPc = system.getPhase(0).getComponent(lc).getPC();
      initAc = system.getPhase(0).getComponent(lc).getAcentricFactor();
    } else {
      initTc = system.getPhase(0).getComponent(hc).getTC();
      initPc = system.getPhase(0).getComponent(hc).getPC();
      initAc = system.getPhase(0).getComponent(hc).getAcentricFactor();
    }

    double lnPratio = Math.log(P / initPc);
    double Tstart = initTc * WILSON_CONST * (1 + initAc) / (WILSON_CONST * (1 + initAc) - lnPratio);
    double Tstartold = 0;

    try {
      double[] Kwil = new double[numberOfComponents];
      for (int i = 0; i < MAX_WILSON_ITERATIONS; i++) {
        double initT = 0;
        double dinitT = 0;
        for (int j = 0; j < numberOfComponents; j++) {
          Kwil[j] = system.getPhase(0).getComponent(j).getPC() / P
              * Math.exp(WILSON_CONST * (1.0 + system.getPhase(0).getComponent(j).getAcentricFactor())
                  * (1.0 - system.getPhase(0).getComponent(j).getTC() / Tstart));
        }

        for (int j = 0; j < numberOfComponents; j++) {
          if (beta < 0.5) {
            initT += system.getPhase(0).getComponent(j).getz() * Kwil[j];
            dinitT += system.getPhase(0).getComponent(j).getz() * Kwil[j] * WILSON_CONST
                * (1 + system.getPhase(0).getComponent(j).getAcentricFactor())
                * system.getPhase(0).getComponent(j).getTC() / (Tstart * Tstart);
          } else {
            initT += system.getPhase(0).getComponent(j).getz() / Kwil[j];
            dinitT -= system.getPhase(0).getComponent(j).getz() / Kwil[j] * WILSON_CONST
                * (1 + system.getPhase(0).getComponent(j).getAcentricFactor())
                * system.getPhase(0).getComponent(j).getTC() / (Tstart * Tstart);
          }
        }

        initT -= 1.0;
        if (Math.abs(initT / dinitT) > 0.1 * Tstart) {
          Tstart -= 0.001 * initT / dinitT;
        } else {
          Tstart -= initT / dinitT;
        }
        if (Math.abs(Tstart - Tstartold) < 1e-5) {
          return Tstart;
        }
        Tstartold = Tstart;
      }
    } catch (Exception ex) {
      lnPratio = Math.log(P / initPc);
      Tstart = initTc * WILSON_CONST * (1 + initAc) / (WILSON_CONST * (1 + initAc) - lnPratio);
    }

    if (Double.isNaN(Tstart) || Double.isInfinite(Tstart)) {
      lnPratio = Math.log(P / initPc);
      Tstart = initTc * WILSON_CONST * (1 + initAc) / (WILSON_CONST * (1 + initAc) - lnPratio);
    }
    return Tstart;
  }

  // ==================== Configuration setters ====================

  /**
   * Set the maximum pressure limit for the phase envelope.
   *
   * @param maxPressure maximum pressure in bara
   */
  public void setMaxPressure(double maxPressure) {
    this.maxPressure = maxPressure;
  }

  /**
   * Set the minimum pressure limit for the phase envelope.
   *
   * @param minPressure minimum pressure in bara
   */
  public void setMinPressure(double minPressure) {
    this.minPressure = minPressure;
  }

  /**
   * Set the maximum temperature step per iteration.
   *
   * @param dTmax max temperature change in K per step
   */
  public void setDTmax(double dTmax) {
    this.dTmax = dTmax;
  }

  /**
   * Set the maximum pressure step per iteration.
   *
   * @param dPmax max pressure change in bar per step
   */
  public void setDPmax(double dPmax) {
    this.dPmax = dPmax;
  }

  /**
   * Set a hard continuation-iteration limit for each pass.
   *
   * @param maximumEnvelopeIterations maximum iterations per continuation pass
   */
  public void setMaximumEnvelopeIterations(int maximumEnvelopeIterations) {
    if (maximumEnvelopeIterations < 10) {
      throw new IllegalArgumentException("maximumEnvelopeIterations must be at least 10");
    }
    this.maximumEnvelopeIterations = maximumEnvelopeIterations;
  }

  /** @return true when a continuation pass stopped at the configured hard iteration limit */
  public boolean isIterationLimitReached() {
    return iterationLimitReached;
  }

  // ==================== Result getters ====================

  /**
   * Get the dew point temperatures in Kelvin.
   *
   * @return array of dew point temperatures
   */
  public double[] getDewPointTemperatures() {
    return dewTempArray;
  }

  /**
   * Get the dew point pressures in bara.
   *
   * @return array of dew point pressures
   */
  public double[] getDewPointPressures() {
    return dewPresArray;
  }

  /**
   * Get the bubble point temperatures in Kelvin.
   *
   * @return array of bubble point temperatures
   */
  public double[] getBubblePointTemperatures() {
    return bubTempArray;
  }

  /**
   * Get the bubble point pressures in bara.
   *
   * @return array of bubble point pressures
   */
  public double[] getBubblePointPressures() {
    return bubPresArray;
  }

  /**
   * Get the cricondentherm values.
   *
   * @return array [T(K), P(bara), 0] at max temperature on envelope
   */
  public double[] getCricondenTherm() {
    return cricondenTherm;
  }

  /**
   * Get the cricondenbar values.
   *
   * @return array [T(K), P(bara), 0] at max pressure on envelope
   */
  public double[] getCricondenBar() {
    return cricondenBar;
  }

  /**
   * Get the critical point temperature.
   *
   * @return critical temperature in Kelvin
   */
  public double getCriticalTemperature() {
    return system.getTC();
  }

  /**
   * Get the critical point pressure.
   *
   * @return critical pressure in bara
   */
  public double getCriticalPressure() {
    return system.getPC();
  }

  // ==================== OperationInterface implementation ====================

  /** {@inheritDoc} */
  @Override
  public void displayResult() {
    // No GUI display - use get() methods to retrieve data
  }

  /** {@inheritDoc} */
  @Override
  public double[][] getPoints(int i) {
    return new double[][] { dewTempArray, dewPresArray, bubTempArray, bubPresArray };
  }

  /**
   * Return the envelope as a list of contiguous branch segments.
   *
   * <p>
   * Preferred over the flat {@link #get(String)} arrays for plotting, machine-readable export, and any code that needs
   * to reason about individual dew / bubble branches. Each segment is a polyline with uniform phase type (dew or
   * bubble) and never contains NaN.
   * </p>
   *
   * <p>
   * The flat arrays returned by {@code get("dewT")}, {@code get("dewP")}, etc. remain available for backward
   * compatibility and use NaN sentinels as branch-break markers.
   * </p>
   *
   * @return unmodifiable list of segments (possibly empty if the envelope trace produced no points)
   */
  public List<EnvelopeSegment> getSegments() {
    return segments;
  }

  /** {@inheritDoc} */
  @Override
  public double[] get(String name) {
    if (name.equals("dewT")) {
      return dewTempArray;
    }
    if (name.equals("dewP")) {
      return dewPresArray;
    }
    if (name.equals("bubT")) {
      return bubTempArray;
    }
    if (name.equals("bubP")) {
      return bubPresArray;
    }
    if (name.equals("dewH")) {
      return dewEnthalpyArray;
    }
    if (name.equals("dewDens")) {
      return dewDensityArray;
    }
    if (name.equals("dewS")) {
      return dewEntropyArray;
    }
    if (name.equals("bubH")) {
      return bubEnthalpyArray;
    }
    if (name.equals("bubDens")) {
      return bubDensityArray;
    }
    if (name.equals("bubS")) {
      return bubEntropyArray;
    }
    if (name.equals("cricondentherm")) {
      return cricondenTherm;
    }
    if (name.equals("cricondenthermX")) {
      return cricondenThermX;
    }
    if (name.equals("cricondenthermY")) {
      return cricondenThermY;
    }
    if (name.equals("cricondenbar")) {
      return cricondenBar;
    }
    if (name.equals("cricondenbarX")) {
      return cricondenBarX;
    }
    if (name.equals("cricondenbarY")) {
      return cricondenBarY;
    }
    if (name.equals("dewT2") || name.equals("dewP2") || name.equals("bubT2") || name.equals("bubP2")) {
      // Return null to match legacy PTphaseEnvelope behavior when no second pass exists.
      // The Michelsen method merges all points into dewT/dewP/bubT/bubP, so separate
      // "pass 2" arrays are not applicable. Returning null ensures downstream consumers
      // (e.g., NeqSimAPI's SPhaseopt.calculate_cricondenbar) correctly fall back to dewP.
      return null;
    }
    if (name.equals("criticalPoint1")) {
      if (!criticalPoints.isEmpty()) {
        return criticalPoints.get(0);
      }
      return new double[] { system.getTC(), system.getPC() };
    }
    if (name.equals("criticalPoint2")) {
      if (criticalPoints.size() >= 2) {
        return criticalPoints.get(1);
      }
      return new double[] { 0, 0 };
    }
    if (name.equals("criticalPoint3")) {
      if (criticalPoints.size() >= 3) {
        return criticalPoints.get(2);
      }
      return new double[] { 0, 0 };
    }
    // Quality line keys: qualityT_X, qualityP_X, qualityVolFrac_X, qualityMassFrac_X
    if (name.startsWith("quality")) {
      return qualityLineData.get(name);
    }
    // Three-phase stability analysis results
    if (name.equals("threePhaseT")) {
      return threePhaseRegionT;
    }
    if (name.equals("threePhaseP")) {
      return threePhaseRegionP;
    }
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public String[][] getResultTable() {
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public void printToFile(String name) {
  }
}
