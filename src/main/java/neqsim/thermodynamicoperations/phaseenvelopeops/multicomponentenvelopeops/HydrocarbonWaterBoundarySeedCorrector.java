package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBranchSeedScanner.Seed;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Routes a coarse three-phase grid seed into one rigorously specified two-to-three-phase boundary family. */
public final class HydrocarbonWaterBoundarySeedCorrector {
  private final SystemInterface template;
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase incipientPhase;
  private int scanIntervals = 80;
  private int maximumBisections = 80;
  private double temperatureToleranceK = 1.0e-5;
  private double tangentPlaneTolerance = 1.0e-8;

  /** Creates a corrector for one explicit topology transition, such as {@code GO -> GOW}. */
  public HydrocarbonWaterBoundarySeedCorrector(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase) {
    if (template == null || retainedPhaseZero == null || retainedPhaseOne == null || incipientPhase == null) {
      throw new IllegalArgumentException("template and all three phase families are required");
    }
    if (retainedPhaseZero == retainedPhaseOne || retainedPhaseZero == incipientPhase
        || retainedPhaseOne == incipientPhase) {
      throw new IllegalArgumentException("retained and incipient phase families must be distinct");
    }
    this.template = template.clone();
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.incipientPhase = incipientPhase;
  }

  /** Sets the strict fixed-pressure scan and root-refinement controls. */
  public HydrocarbonWaterBoundarySeedCorrector setNumericalControls(int scanIntervals, int maximumBisections,
      double temperatureToleranceK, double tangentPlaneTolerance) {
    if (scanIntervals < 2 || maximumBisections < 1 || !Double.isFinite(temperatureToleranceK)
        || temperatureToleranceK <= 0.0 || !Double.isFinite(tangentPlaneTolerance) || tangentPlaneTolerance <= 0.0) {
      throw new IllegalArgumentException("invalid hydrocarbon-water seed-correction controls");
    }
    this.scanIntervals = scanIntervals;
    this.maximumBisections = maximumBisections;
    this.temperatureToleranceK = temperatureToleranceK;
    this.tangentPlaneTolerance = tangentPlaneTolerance;
    return this;
  }

  /**
   * Corrects a stable three-phase grid seed to every strict zero-TPD root of the requested boundary family.
   *
   * <p>
   * The seed contributes only same-state phase compositions and a retained-phase fraction guess. The returned points
   * are recomputed from the specified two-phase equilibrium equations and an independent incipient-phase stability
   * condition. Consequently, an uncorrected TP-grid point can never enter an engineering envelope through this API.
   * </p>
   *
   * @param seed stable grid state containing all requested phases
   * @param minimumTemperatureK lower fixed-pressure search bound
   * @param maximumTemperatureK upper fixed-pressure search bound
   * @return all strict corrected roots and their diagnostics
   */
  public TwoToThreePhaseBoundaryPointSolver.RootSet correctAtSeedPressure(Seed seed, double minimumTemperatureK,
      double maximumTemperatureK) {
    if (seed == null) {
      throw new IllegalArgumentException("a three-phase grid seed is required");
    }
    return correctAtSeedPressure(Collections.singletonList(seed), minimumTemperatureK, maximumTemperatureK);
  }

  /**
   * Corrects all stationary branches represented by same-pressure grid seeds in one multistart solve.
   *
   * <p>
   * One seed with the most balanced retained-phase inventory initializes the specified two-phase flash. Every distinct
   * incipient composition is passed to the stationary-point multistart search, so disconnected third-phase roots are
   * preserved rather than overwritten by the last grid state.
   * </p>
   */
  public TwoToThreePhaseBoundaryPointSolver.RootSet correctAtSeedPressure(List<Seed> seeds, double minimumTemperatureK,
      double maximumTemperatureK) {
    if (seeds == null || seeds.isEmpty()) {
      throw new IllegalArgumentException("at least one three-phase grid seed is required");
    }
    Seed retainedSeed = null;
    double bestRetainedBalance = Double.NEGATIVE_INFINITY;
    double pressureBara = seeds.get(0).getPressureBara();
    java.util.ArrayList<double[]> incipientSeeds = new java.util.ArrayList<double[]>();
    for (Seed seed : seeds) {
      if (seed == null || !seed.hasPhase(retainedPhaseZero) || !seed.hasPhase(retainedPhaseOne)
          || !seed.hasPhase(incipientPhase)) {
        throw new IllegalArgumentException("every grid seed must contain the requested retained and incipient phases");
      }
      if (Math.abs(seed.getPressureBara() - pressureBara) > 1.0e-10 * Math.max(1.0, pressureBara)) {
        throw new IllegalArgumentException("all grid seeds must have the same pressure");
      }
      double retainedFraction = seed.getPhaseFraction(retainedPhaseZero) + seed.getPhaseFraction(retainedPhaseOne);
      double retainedBalance = retainedFraction > 0.0
          ? Math.min(seed.getPhaseFraction(retainedPhaseZero), seed.getPhaseFraction(retainedPhaseOne))
              / retainedFraction
          : Double.NEGATIVE_INFINITY;
      if (retainedBalance > bestRetainedBalance) {
        retainedSeed = seed;
        bestRetainedBalance = retainedBalance;
      }
      incipientSeeds.add(seed.getPhaseComposition(incipientPhase));
    }
    if (retainedSeed == null) {
      throw new IllegalArgumentException("grid seeds have no finite retained-phase inventory");
    }
    double retainedFraction = retainedSeed.getPhaseFraction(retainedPhaseZero)
        + retainedSeed.getPhaseFraction(retainedPhaseOne);
    if (!(retainedFraction > 0.0) || !Double.isFinite(retainedFraction)) {
      throw new IllegalArgumentException("grid seed has no finite retained-phase inventory");
    }
    double betaSeed = retainedSeed.getPhaseFraction(retainedPhaseZero) / retainedFraction;
    if (!(betaSeed > 0.0 && betaSeed < 1.0)) {
      throw new IllegalArgumentException("grid seed does not contain both retained phases");
    }
    return new TwoToThreePhaseBoundaryPointSolver(template, retainedPhaseZero, retainedPhaseOne, incipientPhase)
        .setNumericalControls(scanIntervals, maximumBisections, temperatureToleranceK, tangentPlaneTolerance)
        .solveAll(pressureBara, minimumTemperatureK, maximumTemperatureK, betaSeed,
            retainedSeed.getPhaseComposition(retainedPhaseZero), retainedSeed.getPhaseComposition(retainedPhaseOne),
            incipientSeeds);
  }
}
