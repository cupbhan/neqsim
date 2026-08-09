package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Coarse TP-grid scanner that discovers composition seeds without converting grid states into envelope curves. */
public final class HydrocarbonWaterBranchSeedScanner {
  private final SystemInterface template;
  private final CandidatePhase requestedPhase;
  private double compositionDeduplicationTolerance = 1.0e-4;

  /** Creates a non-destructive seed scanner for one physical phase family. */
  public HydrocarbonWaterBranchSeedScanner(SystemInterface template, CandidatePhase requestedPhase) {
    if (template == null || requestedPhase == null) {
      throw new IllegalArgumentException("template and requested phase are required");
    }
    this.template = template.clone();
    this.requestedPhase = requestedPhase;
  }

  /** Sets the L1 composition threshold used only to remove equivalent seed candidates. */
  public HydrocarbonWaterBranchSeedScanner setCompositionDeduplicationTolerance(double tolerance) {
    if (!Double.isFinite(tolerance) || tolerance <= 0.0) {
      throw new IllegalArgumentException("composition deduplication tolerance must be positive");
    }
    this.compositionDeduplicationTolerance = tolerance;
    return this;
  }

  /**
   * Scans the Cartesian TP grid and retains every distinct requested-phase composition returned by a stable TP flash.
   *
   * <p>
   * These points are seed evidence only. They contain no adjacency and must never be rendered as a phase boundary. Each
   * final boundary point must subsequently be corrected by the appropriate zero-TPD continuation equations.
   * </p>
   */
  public Result scan(double[] temperaturesK, double[] pressuresBara) {
    validateGrid(temperaturesK, pressuresBara);
    List<Seed> seeds = new ArrayList<Seed>();
    List<GridFailure> failures = new ArrayList<GridFailure>();
    int evaluatedPoints = 0;
    for (double pressureBara : pressuresBara) {
      for (double temperatureK : temperaturesK) {
        evaluatedPoints++;
        SystemInterface working = template.clone();
        working.setTemperature(temperatureK);
        working.setPressure(pressureBara);
        working.setMultiPhaseCheck(true);
        working.setMaxNumberOfPhases(3);
        try {
          new ThermodynamicOperations(working).TPflash();
          working.init(1);
          double[][] phaseCompositions = new double[CandidatePhase.values().length][];
          double[] phaseFractions = new double[CandidatePhase.values().length];
          boolean requestedPhaseFound = false;
          for (int phaseIndex = 0; phaseIndex < working.getNumberOfPhases(); phaseIndex++) {
            CandidatePhase physicalPhase = candidatePhase(working.getPhase(phaseIndex).getType());
            phaseCompositions[physicalPhase.ordinal()] = composition(working, phaseIndex);
            phaseFractions[physicalPhase.ordinal()] = working.getBeta(phaseIndex);
            if (physicalPhase != requestedPhase) {
              continue;
            }
            requestedPhaseFound = true;
          }
          if (requestedPhaseFound) {
            addDistinctSeed(seeds,
                new Seed(temperatureK, pressureBara, requestedPhase, phaseCompositions[requestedPhase.ordinal()],
                    working.getNumberOfPhases(), phaseCompositions, phaseFractions));
          } else {
            failures.add(new GridFailure(temperatureK, pressureBara, "REQUESTED_PHASE_ABSENT"));
          }
        } catch (RuntimeException error) {
          failures.add(new GridFailure(temperatureK, pressureBara,
              error.getMessage() == null ? error.getClass().getSimpleName() : error.getMessage()));
        }
      }
    }
    return new Result(seeds, failures, evaluatedPoints);
  }

  private void addDistinctSeed(List<Seed> seeds, Seed candidate) {
    for (Seed seed : seeds) {
      if (compositionDistance(seed.composition, candidate.composition) <= compositionDeduplicationTolerance) {
        return;
      }
    }
    seeds.add(candidate);
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

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

  private static void validateGrid(double[] temperaturesK, double[] pressuresBara) {
    if (temperaturesK == null || temperaturesK.length == 0 || pressuresBara == null || pressuresBara.length == 0) {
      throw new IllegalArgumentException("non-empty temperature and pressure grids are required");
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

  /** One requested-phase composition found at a stable TP-grid state. */
  public static final class Seed {
    private final double temperatureK;
    private final double pressureBara;
    private final CandidatePhase phase;
    private final double[] composition;
    private final int stablePhaseCount;
    private final double[][] phaseCompositions;
    private final double[] phaseFractions;

    private Seed(double temperatureK, double pressureBara, CandidatePhase phase, double[] composition,
        int stablePhaseCount, double[][] phaseCompositions, double[] phaseFractions) {
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.phase = phase;
      this.composition = composition.clone();
      this.stablePhaseCount = stablePhaseCount;
      this.phaseCompositions = new double[phaseCompositions.length][];
      for (int phaseIndex = 0; phaseIndex < phaseCompositions.length; phaseIndex++) {
        this.phaseCompositions[phaseIndex] = phaseCompositions[phaseIndex] == null ? new double[0]
            : phaseCompositions[phaseIndex].clone();
      }
      this.phaseFractions = phaseFractions.clone();
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public CandidatePhase getPhase() {
      return phase;
    }

    public double[] getComposition() {
      return composition.clone();
    }

    public int getStablePhaseCount() {
      return stablePhaseCount;
    }

    public boolean hasPhase(CandidatePhase candidatePhase) {
      return candidatePhase != null && phaseCompositions[candidatePhase.ordinal()].length > 0;
    }

    public double[] getPhaseComposition(CandidatePhase candidatePhase) {
      if (candidatePhase == null) {
        throw new IllegalArgumentException("candidate phase is required");
      }
      return phaseCompositions[candidatePhase.ordinal()].clone();
    }

    public double getPhaseFraction(CandidatePhase candidatePhase) {
      if (candidatePhase == null) {
        throw new IllegalArgumentException("candidate phase is required");
      }
      return phaseFractions[candidatePhase.ordinal()];
    }
  }

  /** Explicit failure or absence diagnostic for one grid point. */
  public static final class GridFailure {
    private final double temperatureK;
    private final double pressureBara;
    private final String reason;

    private GridFailure(double temperatureK, double pressureBara, String reason) {
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

  /** Immutable seed scan preserving all grid failures for topology auditing. */
  public static final class Result {
    private final List<Seed> seeds;
    private final List<GridFailure> failures;
    private final int evaluatedPointCount;

    private Result(List<Seed> seeds, List<GridFailure> failures, int evaluatedPointCount) {
      this.seeds = Collections.unmodifiableList(new ArrayList<Seed>(seeds));
      this.failures = Collections.unmodifiableList(new ArrayList<GridFailure>(failures));
      this.evaluatedPointCount = evaluatedPointCount;
    }

    public List<Seed> getSeeds() {
      return seeds;
    }

    public List<GridFailure> getFailures() {
      return failures;
    }

    public int getEvaluatedPointCount() {
      return evaluatedPointCount;
    }
  }
}
