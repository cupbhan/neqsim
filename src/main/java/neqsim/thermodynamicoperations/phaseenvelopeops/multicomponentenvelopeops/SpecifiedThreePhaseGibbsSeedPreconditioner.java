package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.flashops.TPmultiflash;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Produces a Gibbs-descent seed for a specified three-phase Newton flash.
 *
 * <p>
 * This class deliberately invokes only {@link TPmultiflash#solveBeta()} on a private, caller-seeded three-slot system.
 * It does not run the ordinary stability analysis, add phases, remove phases, or classify duplicate phase families.
 * Consequently, two independent oil slots can move away from the coincident solution before the rigorous specified
 * flash corrects fugacity equality and material balance.
 * </p>
 */
public final class SpecifiedThreePhaseGibbsSeedPreconditioner {
  private final SystemInterface template;
  private final CandidatePhase phaseZero;
  private final CandidatePhase phaseOne;
  private final CandidatePhase phaseTwo;
  private int maximumCalls = 4;
  private double errorTolerance = 1.0e-10;

  /** Creates a preconditioner for three specified phase slots; repeated families are permitted. */
  public SpecifiedThreePhaseGibbsSeedPreconditioner(SystemInterface template, CandidatePhase phaseZero,
      CandidatePhase phaseOne, CandidatePhase phaseTwo) {
    if (template == null || phaseZero == null || phaseOne == null || phaseTwo == null) {
      throw new IllegalArgumentException("template and all three phase slots are required");
    }
    this.template = template.clone();
    this.phaseZero = phaseZero;
    this.phaseOne = phaseOne;
    this.phaseTwo = phaseTwo;
  }

  /** Sets the number of bounded Gibbs calls and the reported convergence tolerance. */
  public SpecifiedThreePhaseGibbsSeedPreconditioner setNumericalControls(int maximumCalls, double errorTolerance) {
    if (maximumCalls < 1 || !Double.isFinite(errorTolerance) || errorTolerance <= 0.0) {
      throw new IllegalArgumentException("invalid Gibbs seed preconditioner controls");
    }
    this.maximumCalls = maximumCalls;
    this.errorTolerance = errorTolerance;
    return this;
  }

  /** Runs bounded Gibbs phase-fraction/composition updates without mutating the template. */
  public Result precondition(double temperatureK, double pressureBara, double[] phaseFractions,
      double[] phaseZeroComposition, double[] phaseOneComposition, double[] phaseTwoComposition) {
    validateInputs(temperatureK, pressureBara, phaseFractions, phaseZeroComposition, phaseOneComposition,
        phaseTwoComposition);
    SystemInterface working;
    try {
      working = createWorkingSystem(temperatureK, pressureBara, phaseFractions,
          new double[][] {phaseZeroComposition, phaseOneComposition, phaseTwoComposition});
    } catch (RuntimeException error) {
      return Result.failure(phaseZero, phaseOne, phaseTwo, temperatureK, pressureBara,
          template.getPhase(0).getNumberOfComponents(), error.getMessage());
    }
    TPmultiflash operation = new TPmultiflash(working, false);
    operation.setDoubleArrays();
    double error = Double.POSITIVE_INFINITY;
    int calls = 0;
    String failureMessage = null;
    while (calls < maximumCalls && error > errorTolerance) {
      calls++;
      try {
        error = operation.solveBeta();
      } catch (RuntimeException exception) {
        failureMessage = "Gibbs beta preconditioning failed: " + exception.getMessage();
        break;
      }
      if (!Double.isFinite(error)) {
        failureMessage = "Gibbs beta preconditioning returned a non-finite error";
        break;
      }
      if (working.getNumberOfPhases() != 3) {
        failureMessage = "Gibbs beta preconditioning changed the specified phase count";
        break;
      }
    }
    double[] fractions = new double[3];
    double[][] compositions = new double[3][template.getPhase(0).getNumberOfComponents()];
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      fractions[phaseIndex] = working.getBeta(phaseIndex);
      for (int componentIndex = 0; componentIndex < compositions[phaseIndex].length; componentIndex++) {
        compositions[phaseIndex][componentIndex] = working.getPhase(phaseIndex).getComponent(componentIndex).getx();
      }
    }
    return new Result(phaseZero, phaseOne, phaseTwo, temperatureK, pressureBara, fractions, compositions, calls, error,
        error <= errorTolerance, failureMessage);
  }

  private SystemInterface createWorkingSystem(double temperatureK, double pressureBara, double[] fractions,
      double[][] compositions) {
    SystemInterface working = template.clone();
    working.setMultiPhaseCheck(false);
    working.setMaxNumberOfPhases(3);
    working.setNumberOfPhases(3);
    PhaseInterface phaseTemplate = template.getPhase(0);
    CandidatePhase[] phaseSlots = new CandidatePhase[] {phaseZero, phaseOne, phaseTwo};
    double fractionTotal = phaseFractionsTotal(fractions);
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      working.setPhase(phaseTemplate.clone(), phaseIndex);
      working.setPhaseType(phaseIndex, toPhaseType(phaseSlots[phaseIndex]));
      working.setBeta(phaseIndex, fractions[phaseIndex] / fractionTotal);
      for (int componentIndex = 0; componentIndex < compositions[phaseIndex].length; componentIndex++) {
        working.getPhase(phaseIndex).getComponent(componentIndex).setx(compositions[phaseIndex][componentIndex]);
      }
      working.getPhase(phaseIndex).normalize();
      working.setPhaseType(phaseIndex, toPhaseType(phaseSlots[phaseIndex]));
      working.init(1, phaseIndex);
    }
    return working;
  }

  private void validateInputs(double temperatureK, double pressureBara, double[] phaseFractions,
      double[] phaseZeroComposition, double[] phaseOneComposition, double[] phaseTwoComposition) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    if (!Double.isFinite(temperatureK) || temperatureK < 50.0 || !Double.isFinite(pressureBara) || pressureBara <= 0.0
        || phaseFractions == null || phaseFractions.length != 3 || phaseZeroComposition == null
        || phaseOneComposition == null || phaseTwoComposition == null || phaseZeroComposition.length != componentCount
        || phaseOneComposition.length != componentCount || phaseTwoComposition.length != componentCount) {
      throw new IllegalArgumentException("invalid specified three-phase Gibbs seed input");
    }
    phaseFractionsTotal(phaseFractions);
    validateComposition(phaseZeroComposition);
    validateComposition(phaseOneComposition);
    validateComposition(phaseTwoComposition);
  }

  private static double phaseFractionsTotal(double[] fractions) {
    double total = 0.0;
    for (double fraction : fractions) {
      if (!Double.isFinite(fraction) || fraction <= 0.0) {
        throw new IllegalArgumentException("all phase fractions must be positive");
      }
      total += fraction;
    }
    return total;
  }

  private static void validateComposition(double[] composition) {
    double total = 0.0;
    for (double value : composition) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException("phase composition contains an invalid value");
      }
      total += value;
    }
    if (!(total > 0.0)) {
      throw new IllegalArgumentException("phase composition cannot be normalized");
    }
  }

  private static PhaseType toPhaseType(CandidatePhase phase) {
    switch (phase) {
    case GAS:
      return PhaseType.GAS;
    case OIL:
      return PhaseType.OIL;
    case AQUEOUS:
      return PhaseType.AQUEOUS;
    default:
      throw new IllegalArgumentException("unsupported phase family " + phase);
    }
  }

  /** Immutable Gibbs seed and convergence diagnostics. */
  public static final class Result {
    private final CandidatePhase phaseZero;
    private final CandidatePhase phaseOne;
    private final CandidatePhase phaseTwo;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] phaseFractions;
    private final double[][] phaseCompositions;
    private final int calls;
    private final double finalError;
    private final boolean converged;
    private final String failureMessage;

    private Result(CandidatePhase phaseZero, CandidatePhase phaseOne, CandidatePhase phaseTwo, double temperatureK,
        double pressureBara, double[] phaseFractions, double[][] phaseCompositions, int calls, double finalError,
        boolean converged, String failureMessage) {
      this.phaseZero = phaseZero;
      this.phaseOne = phaseOne;
      this.phaseTwo = phaseTwo;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.phaseFractions = phaseFractions.clone();
      this.phaseCompositions = copy(phaseCompositions);
      this.calls = calls;
      this.finalError = finalError;
      this.converged = converged;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase phaseZero, CandidatePhase phaseOne, CandidatePhase phaseTwo,
        double temperatureK, double pressureBara, int componentCount, String failureMessage) {
      return new Result(phaseZero, phaseOne, phaseTwo, temperatureK, pressureBara, new double[3],
          new double[3][componentCount], 0, Double.NaN, false, failureMessage);
    }

    public boolean isConverged() {
      return converged && failureMessage == null;
    }

    public CandidatePhase getPhaseZero() {
      return phaseZero;
    }

    public CandidatePhase getPhaseOne() {
      return phaseOne;
    }

    public CandidatePhase getPhaseTwo() {
      return phaseTwo;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double[] getPhaseFractions() {
      return phaseFractions.clone();
    }

    public double[] getPhaseComposition(int phaseIndex) {
      return phaseCompositions[phaseIndex].clone();
    }

    public double getPhaseCompositionDistance(int firstPhaseIndex, int secondPhaseIndex) {
      double distance = 0.0;
      for (int componentIndex = 0; componentIndex < phaseCompositions[firstPhaseIndex].length; componentIndex++) {
        distance += Math.abs(
            phaseCompositions[firstPhaseIndex][componentIndex] - phaseCompositions[secondPhaseIndex][componentIndex]);
      }
      return distance;
    }

    public int getCalls() {
      return calls;
    }

    public double getFinalError() {
      return finalError;
    }

    public String getFailureMessage() {
      return failureMessage;
    }

    private static double[][] copy(double[][] matrix) {
      double[][] result = new double[matrix.length][];
      for (int row = 0; row < matrix.length; row++) {
        result[row] = matrix[row].clone();
      }
      return result;
    }
  }
}
