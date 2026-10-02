package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import Jama.Matrix;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Damped Newton three-phase TP flash with caller-specified phase slots.
 *
 * <p>
 * Unlike the ordinary NeqSim multiphase flash, phase families are allowed to repeat. In particular, a gas phase and two
 * compositionally distinct oil phases can be represented as {@code GAS/OIL/OIL}. The two oil slots retain independent
 * compositions and fractions even though both use the oil EOS root.
 * </p>
 *
 * <p>
 * The unknowns are two sets of component log K-values relative to phase zero and two phase-fraction logits. Component
 * material balances are eliminated analytically. Two remaining equations enforce normalization of all three phase
 * compositions, while the other equations enforce component fugacity equality.
 * </p>
 */
public final class SpecifiedThreePhaseFlashSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;
  private static final double MINIMUM_PHASE_FRACTION = 1.0e-14;

  private final SystemInterface template;
  private final CandidatePhase phaseZero;
  private final CandidatePhase phaseOne;
  private final CandidatePhase phaseTwo;
  private int maximumIterations = 100;
  private double residualTolerance = 1.0e-9;
  private double finiteDifferenceStep = 2.0e-5;
  private int deflatedPhaseOne = -1;
  private int deflatedPhaseTwo = -1;
  private double deflationPower = 1.0;
  private double deflationShift = 1.0;

  /** Creates a specified three-phase flash; repeated physical phase families are permitted. */
  public SpecifiedThreePhaseFlashSolver(SystemInterface template, CandidatePhase phaseZero, CandidatePhase phaseOne,
      CandidatePhase phaseTwo) {
    if (template == null || phaseZero == null || phaseOne == null || phaseTwo == null) {
      throw new IllegalArgumentException("template and all three phase slots must be specified");
    }
    this.template = template.clone();
    this.phaseZero = phaseZero;
    this.phaseOne = phaseOne;
    this.phaseTwo = phaseTwo;
  }

  /** Sets Newton and numerical-Jacobian controls. */
  public SpecifiedThreePhaseFlashSolver setNumericalControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    if (maximumIterations < 1 || !positive(residualTolerance) || !positive(finiteDifferenceStep)) {
      throw new IllegalArgumentException("invalid specified three-phase flash numerical controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /**
   * Deflates the known solution manifold where phase slots one and two have identical compositions.
   *
   * <p>
   * The physical equations are unchanged. Their residual vector is multiplied by the shifted deflation factor
   * {@code distance^-power + shift}, using the distance between the two log K-value blocks. This removes the
   * coincident-phase solution as a Newton attractor while retaining every distinct equilibrium root.
   * </p>
   */
  public SpecifiedThreePhaseFlashSolver setCoincidentPhaseDeflation(int firstPhaseIndex, int secondPhaseIndex,
      double power, double shift) {
    if (firstPhaseIndex != 1 || secondPhaseIndex != 2 || !positive(power) || !Double.isFinite(shift) || shift < 0.0) {
      throw new IllegalArgumentException("current log-K deflation requires phase slots one and two");
    }
    this.deflatedPhaseOne = firstPhaseIndex;
    this.deflatedPhaseTwo = secondPhaseIndex;
    this.deflationPower = power;
    this.deflationShift = shift;
    return this;
  }

  /**
   * Solves one fixed-temperature, fixed-pressure three-phase equilibrium state.
   *
   * @param temperatureK temperature in kelvin
   * @param pressureBara pressure in bara
   * @param phaseFractions three positive fractions that sum to one after normalization
   * @param phaseZeroComposition phase-zero composition seed
   * @param phaseOneComposition phase-one composition seed
   * @param phaseTwoComposition phase-two composition seed
   * @return immutable flash result and numerical diagnostics
   */
  public Result solve(double temperatureK, double pressureBara, double[] phaseFractions, double[] phaseZeroComposition,
      double[] phaseOneComposition, double[] phaseTwoComposition) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    validateInputs(temperatureK, pressureBara, phaseFractions, phaseZeroComposition, phaseOneComposition,
        phaseTwoComposition, componentCount);
    double[] normalizedFractions = normalizedFractions(phaseFractions);
    double[] variables = new double[2 * componentCount + 2];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double reference = Math.max(phaseZeroComposition[componentIndex], MINIMUM_COMPOSITION);
      variables[componentIndex] = Math
          .log(Math.max(phaseOneComposition[componentIndex], MINIMUM_COMPOSITION) / reference);
      variables[componentCount + componentIndex] = Math
          .log(Math.max(phaseTwoComposition[componentIndex], MINIMUM_COMPOSITION) / reference);
    }
    variables[2 * componentCount] = Math.log(normalizedFractions[1] / normalizedFractions[0]);
    variables[2 * componentCount + 1] = Math.log(normalizedFractions[2] / normalizedFractions[0]);

    Evaluation evaluation;
    try {
      evaluation = evaluate(temperatureK, pressureBara, variables);
    } catch (RuntimeException error) {
      return Result.failure(phaseZero, phaseOne, phaseTwo, temperatureK, pressureBara, componentCount,
          error.getMessage());
    }
    double initialResidual = evaluation.maximumResidual;
    String failureMessage = null;
    int iterations = 0;
    double conditionNumber = Double.NaN;
    while (iterations < maximumIterations && evaluation.maximumResidual > residualTolerance) {
      iterations++;
      Matrix jacobian;
      try {
        jacobian = numericalJacobian(temperatureK, pressureBara, variables, evaluation.residual);
        conditionNumber = jacobian.cond();
      } catch (RuntimeException error) {
        failureMessage = "Jacobian evaluation failed: " + error.getMessage();
        break;
      }
      Matrix correction;
      try {
        correction = solveCorrection(jacobian, evaluation.residual, conditionNumber);
      } catch (RuntimeException error) {
        failureMessage = "Jacobian solve failed: " + error.getMessage();
        break;
      }
      limitCorrection(correction, componentCount);
      boolean accepted = false;
      double damping = 1.0;
      for (int lineSearch = 0; lineSearch < 24; lineSearch++) {
        double[] trialVariables = variables.clone();
        for (int variableIndex = 0; variableIndex < trialVariables.length; variableIndex++) {
          double lower = variableIndex < 2 * componentCount ? -700.0 : -32.0;
          double upper = variableIndex < 2 * componentCount ? 700.0 : 32.0;
          trialVariables[variableIndex] = Math.max(lower,
              Math.min(upper, variables[variableIndex] - damping * correction.get(variableIndex, 0)));
        }
        try {
          Evaluation trial = evaluate(temperatureK, pressureBara, trialVariables);
          if (trial.maximumResidual < evaluation.maximumResidual) {
            variables = trialVariables;
            evaluation = trial;
            accepted = true;
            break;
          }
        } catch (RuntimeException error) {
          // Try a shorter Newton step.
        }
        damping *= 0.5;
      }
      if (!accepted) {
        failureMessage = "line search could not reduce the specified three-phase flash residual";
        break;
      }
    }
    boolean converged = evaluation.maximumResidual <= residualTolerance;
    if (failureMessage == null && !converged) {
      failureMessage = "maximum specified three-phase flash iteration count reached";
    }
    return new Result(phaseZero, phaseOne, phaseTwo, temperatureK, pressureBara, evaluation.phaseFractions,
        evaluation.phaseCompositions, initialResidual, evaluation.maximumResidual,
        evaluation.thermodynamicMaximumResidual,
        materialBalanceResidual(evaluation.phaseFractions, evaluation.phaseCompositions), iterations, conditionNumber,
        converged, failureMessage);
  }

  /** Reconstructs a private three-phase thermodynamic system from a converged result. */
  public SystemInterface toThermodynamicSystem(Result result) {
    if (result == null || result.getPhaseZero() != phaseZero || result.getPhaseOne() != phaseOne
        || result.getPhaseTwo() != phaseTwo || !result.isConverged()) {
      throw new IllegalArgumentException("a converged result with matching phase slots is required");
    }
    return createWorkingSystem(result.getTemperatureK(), result.getPressureBara(), result.getPhaseFractions(),
        result.phaseCompositions);
  }

  private Matrix numericalJacobian(double temperatureK, double pressureBara, double[] variables,
      double[] baseResidual) {
    int variableCount = variables.length;
    Matrix jacobian = new Matrix(variableCount, variableCount);
    for (int column = 0; column < variableCount; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(variables[column]));
      double[] plus = variables.clone();
      double[] minus = variables.clone();
      plus[column] += step;
      minus[column] -= step;
      double[] plusResidual = null;
      double[] minusResidual = null;
      try {
        plusResidual = evaluate(temperatureK, pressureBara, plus).residual;
      } catch (RuntimeException error) {
        // A one-sided column can still be used.
      }
      try {
        minusResidual = evaluate(temperatureK, pressureBara, minus).residual;
      } catch (RuntimeException error) {
        // A one-sided column can still be used.
      }
      if (plusResidual == null && minusResidual == null) {
        throw new IllegalStateException("both finite-difference evaluations failed");
      }
      double denominator;
      if (plusResidual == null) {
        plusResidual = baseResidual;
        denominator = variables[column] - minus[column];
      } else if (minusResidual == null) {
        minusResidual = baseResidual;
        denominator = plus[column] - variables[column];
      } else {
        denominator = plus[column] - minus[column];
      }
      for (int row = 0; row < variableCount; row++) {
        jacobian.set(row, column, (plusResidual[row] - minusResidual[row]) / denominator);
      }
    }
    return jacobian;
  }

  private static Matrix solveCorrection(Matrix jacobian, double[] residual, double conditionNumber) {
    Matrix rightHandSide = new Matrix(residual, residual.length);
    if (Double.isFinite(conditionNumber) && conditionNumber < 1.0e11) {
      try {
        return jacobian.solve(rightHandSide);
      } catch (RuntimeException error) {
        // Fall through to a regularized least-squares correction.
      }
    }
    Matrix transpose = jacobian.transpose();
    Matrix normal = transpose.times(jacobian);
    Matrix normalRightHandSide = transpose.times(rightHandSide);
    double lambda = 1.0e-10;
    RuntimeException lastFailure = null;
    for (int attempt = 0; attempt < 8; attempt++) {
      try {
        return normal.plus(Matrix.identity(normal.getRowDimension(), normal.getColumnDimension()).times(lambda))
            .solve(normalRightHandSide);
      } catch (RuntimeException error) {
        lastFailure = error;
        lambda *= 100.0;
      }
    }
    throw lastFailure == null ? new IllegalStateException("regularized Jacobian solve failed") : lastFailure;
  }

  /** Numerically safe {@code log(sum(exp(values)))}, used so the inventory never has to be formed directly. */
  private static double logSumOfExponentials(double[] values) {
    double maximum = Double.NEGATIVE_INFINITY;
    for (int index = 0; index < values.length; index++) {
      maximum = Math.max(maximum, values[index]);
    }
    if (!Double.isFinite(maximum)) {
      return maximum;
    }
    double total = 0.0;
    for (int index = 0; index < values.length; index++) {
      total += Math.exp(values[index] - maximum);
    }
    return maximum + Math.log(total);
  }

  private Evaluation evaluate(double temperatureK, double pressureBara, double[] variables) {
    int componentCount = (variables.length - 2) / 2;
    double[] fractions = softmaxFractions(variables[2 * componentCount], variables[2 * componentCount + 1]);
    double[] overall = overallComposition();
    double[][] compositions = new double[3][componentCount];
    double[] sums = new double[3];
    // Evaluated in log space so the composition floor never enters the generalized Rachford-Rice residual. A heavy
    // pseudo-component against water reaches K ~ 1e99 and a true reference-phase fraction near 1e-102, below the
    // floor; flooring it propagates into the other phases as a large spurious inventory and pins the residual to
    // the floor rather than to the equations. The floor belongs only at the EOS input.
    double[][] logCompositions = new double[3][componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double firstLogK = Math.max(-700.0, Math.min(700.0, variables[componentIndex]));
      double secondLogK = Math.max(-700.0, Math.min(700.0, variables[componentCount + componentIndex]));
      double logDenominator = logSumOfExponentials(new double[] {Math.log(fractions[0]),
          Math.log(fractions[1]) + firstLogK, Math.log(fractions[2]) + secondLogK});
      if (!Double.isFinite(logDenominator)) {
        throw new IllegalStateException("invalid generalized Rachford-Rice denominator");
      }
      logCompositions[0][componentIndex] = Math.log(overall[componentIndex]) - logDenominator;
      logCompositions[1][componentIndex] = firstLogK + logCompositions[0][componentIndex];
      logCompositions[2][componentIndex] = secondLogK + logCompositions[0][componentIndex];
      for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
        compositions[phaseIndex][componentIndex] = Math.exp(logCompositions[phaseIndex][componentIndex]);
      }
    }
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      sums[phaseIndex] = Math.exp(logSumOfExponentials(logCompositions[phaseIndex]));
      if (!Double.isFinite(sums[phaseIndex]) || !(sums[phaseIndex] > 0.0)) {
        throw new IllegalStateException("specified three-phase flash inventory left the physical domain");
      }
    }
    SystemInterface working = createWorkingSystem(temperatureK, pressureBara, fractions, compositions);
    double[] residual = new double[2 * componentCount + 2];
    double thermodynamicMaximumResidual = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      residual[componentIndex] = variables[componentIndex]
          + working.getPhase(1).getComponent(componentIndex).getLogFugacityCoefficient()
          - working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
      residual[componentCount + componentIndex] = variables[componentCount + componentIndex]
          + working.getPhase(2).getComponent(componentIndex).getLogFugacityCoefficient()
          - working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
      thermodynamicMaximumResidual = Math.max(thermodynamicMaximumResidual, Math.abs(residual[componentIndex]));
      thermodynamicMaximumResidual = Math.max(thermodynamicMaximumResidual,
          Math.abs(residual[componentCount + componentIndex]));
    }
    residual[2 * componentCount] = sums[1] - sums[0];
    residual[2 * componentCount + 1] = sums[2] - sums[0];
    thermodynamicMaximumResidual = Math.max(thermodynamicMaximumResidual, Math.abs(residual[2 * componentCount]));
    thermodynamicMaximumResidual = Math.max(thermodynamicMaximumResidual, Math.abs(residual[2 * componentCount + 1]));
    double multiplier = deflationMultiplier(variables, componentCount);
    for (int residualIndex = 0; residualIndex < residual.length; residualIndex++) {
      residual[residualIndex] *= multiplier;
    }
    double maximumResidual = thermodynamicMaximumResidual * multiplier;
    if (!Double.isFinite(maximumResidual) || !Double.isFinite(thermodynamicMaximumResidual)) {
      throw new IllegalStateException("specified three-phase flash residual is non-finite");
    }
    return new Evaluation(residual, maximumResidual, thermodynamicMaximumResidual, fractions, compositions);
  }

  private double deflationMultiplier(double[] variables, int componentCount) {
    if (deflatedPhaseOne < 0 || deflatedPhaseTwo < 0) {
      return 1.0;
    }
    double squaredDistance = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double difference = variables[componentIndex] - variables[componentCount + componentIndex];
      squaredDistance += difference * difference;
    }
    double distance = Math.sqrt(squaredDistance);
    if (!(distance > 1.0e-12) || !Double.isFinite(distance)) {
      throw new IllegalStateException("deflated phase slots coalesced");
    }
    return Math.pow(distance, -deflationPower) + deflationShift;
  }

  private SystemInterface createWorkingSystem(double temperatureK, double pressureBara, double[] fractions,
      double[][] compositions) {
    SystemInterface working = template.clone();
    working.setMultiPhaseCheck(false);
    working.setMaxNumberOfPhases(3);
    working.setNumberOfPhases(3);
    PhaseInterface phaseTemplate = template.getPhase(0);
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      working.setPhase(phaseTemplate.clone(), phaseIndex);
    }
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    CandidatePhase[] phaseSlots = new CandidatePhase[] {phaseZero, phaseOne, phaseTwo};
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      working.setPhaseType(phaseIndex, toPhaseType(phaseSlots[phaseIndex]));
      working.setBeta(phaseIndex, fractions[phaseIndex]);
      for (int componentIndex = 0; componentIndex < compositions[phaseIndex].length; componentIndex++) {
        // The EOS requires a strictly positive mole fraction; this is the only place the floor belongs.
        working.getPhase(phaseIndex).getComponent(componentIndex)
            .setx(Math.max(compositions[phaseIndex][componentIndex], MINIMUM_COMPOSITION));
      }
      working.getPhase(phaseIndex).normalize();
      working.setPhaseType(phaseIndex, toPhaseType(phaseSlots[phaseIndex]));
      working.init(1, phaseIndex);
    }
    return working;
  }

  private double[] overallComposition() {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    double[] composition = new double[componentCount];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(),
          MINIMUM_COMPOSITION);
      total += composition[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }

  private double materialBalanceResidual(double[] fractions, double[][] compositions) {
    double[] overall = overallComposition();
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < overall.length; componentIndex++) {
      double reconstructed = 0.0;
      for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
        reconstructed += fractions[phaseIndex] * compositions[phaseIndex][componentIndex];
      }
      maximum = Math.max(maximum, Math.abs(reconstructed - overall[componentIndex]));
    }
    return maximum;
  }

  private static void limitCorrection(Matrix correction, int componentCount) {
    double scale = 1.0;
    for (int variableIndex = 0; variableIndex < correction.getRowDimension(); variableIndex++) {
      double maximumStep = variableIndex < 2 * componentCount ? 2.0 : 1.0;
      double magnitude = Math.abs(correction.get(variableIndex, 0));
      if (magnitude > maximumStep) {
        scale = Math.min(scale, maximumStep / magnitude);
      }
    }
    if (scale < 1.0) {
      correction.timesEquals(scale);
    }
  }

  private static double[] softmaxFractions(double firstLogit, double secondLogit) {
    double maximum = Math.max(0.0, Math.max(firstLogit, secondLogit));
    double zero = Math.exp(-maximum);
    double one = Math.exp(firstLogit - maximum);
    double two = Math.exp(secondLogit - maximum);
    double total = zero + one + two;
    double[] fractions = new double[] {zero / total, one / total, two / total};
    for (int phaseIndex = 0; phaseIndex < fractions.length; phaseIndex++) {
      fractions[phaseIndex] = Math.max(fractions[phaseIndex], MINIMUM_PHASE_FRACTION);
    }
    return normalizedFractions(fractions);
  }

  private static double[] normalizedFractions(double[] phaseFractions) {
    double[] result = phaseFractions.clone();
    double total = 0.0;
    for (double fraction : result) {
      if (!positive(fraction)) {
        throw new IllegalArgumentException("all specified phase fractions must be positive");
      }
      total += fraction;
    }
    for (int phaseIndex = 0; phaseIndex < result.length; phaseIndex++) {
      result[phaseIndex] /= total;
    }
    return result;
  }

  private static void validateInputs(double temperatureK, double pressureBara, double[] phaseFractions,
      double[] phaseZeroComposition, double[] phaseOneComposition, double[] phaseTwoComposition, int componentCount) {
    if (!Double.isFinite(temperatureK) || temperatureK < 50.0 || !positive(pressureBara) || phaseFractions == null
        || phaseFractions.length != 3 || phaseZeroComposition == null || phaseOneComposition == null
        || phaseTwoComposition == null || phaseZeroComposition.length != componentCount
        || phaseOneComposition.length != componentCount || phaseTwoComposition.length != componentCount) {
      throw new IllegalArgumentException("invalid specified three-phase flash input");
    }
    normalizedFractions(phaseFractions);
    validateComposition(phaseZeroComposition);
    validateComposition(phaseOneComposition);
    validateComposition(phaseTwoComposition);
  }

  private static void validateComposition(double[] composition) {
    double total = 0.0;
    for (double value : composition) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException("specified phase composition contains an invalid value");
      }
      total += value;
    }
    if (!(total > 0.0)) {
      throw new IllegalArgumentException("specified phase composition cannot be normalized");
    }
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
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

  private static final class Evaluation {
    private final double[] residual;
    private final double maximumResidual;
    private final double thermodynamicMaximumResidual;
    private final double[] phaseFractions;
    private final double[][] phaseCompositions;

    private Evaluation(double[] residual, double maximumResidual, double thermodynamicMaximumResidual,
        double[] phaseFractions, double[][] phaseCompositions) {
      this.residual = residual;
      this.maximumResidual = maximumResidual;
      this.thermodynamicMaximumResidual = thermodynamicMaximumResidual;
      this.phaseFractions = phaseFractions;
      this.phaseCompositions = phaseCompositions;
    }
  }

  /** Immutable specified three-phase TP flash result. */
  public static final class Result {
    private final CandidatePhase phaseZero;
    private final CandidatePhase phaseOne;
    private final CandidatePhase phaseTwo;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] phaseFractions;
    private final double[][] phaseCompositions;
    private final double initialMaximumResidual;
    private final double maximumResidual;
    private final double thermodynamicMaximumResidual;
    private final double materialBalanceResidual;
    private final int iterations;
    private final double jacobianConditionNumber;
    private final boolean converged;
    private final String failureMessage;

    private Result(CandidatePhase phaseZero, CandidatePhase phaseOne, CandidatePhase phaseTwo, double temperatureK,
        double pressureBara, double[] phaseFractions, double[][] phaseCompositions, double initialMaximumResidual,
        double maximumResidual, double thermodynamicMaximumResidual, double materialBalanceResidual, int iterations,
        double jacobianConditionNumber, boolean converged, String failureMessage) {
      this.phaseZero = phaseZero;
      this.phaseOne = phaseOne;
      this.phaseTwo = phaseTwo;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.phaseFractions = phaseFractions.clone();
      this.phaseCompositions = copy(phaseCompositions);
      this.initialMaximumResidual = initialMaximumResidual;
      this.maximumResidual = maximumResidual;
      this.thermodynamicMaximumResidual = thermodynamicMaximumResidual;
      this.materialBalanceResidual = materialBalanceResidual;
      this.iterations = iterations;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.converged = converged;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase phaseZero, CandidatePhase phaseOne, CandidatePhase phaseTwo,
        double temperatureK, double pressureBara, int componentCount, String failureMessage) {
      return new Result(phaseZero, phaseOne, phaseTwo, temperatureK, pressureBara, new double[3],
          new double[3][componentCount], Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0, Double.NaN, false,
          failureMessage);
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

    public double getPhaseFraction(int phaseIndex) {
      return phaseFractions[phaseIndex];
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

    public double getInitialMaximumResidual() {
      return initialMaximumResidual;
    }

    public double getMaximumResidual() {
      return maximumResidual;
    }

    /** @return unscaled fugacity/normalization residual before optional branch deflation */
    public double getThermodynamicMaximumResidual() {
      return thermodynamicMaximumResidual;
    }

    public double getMaterialBalanceResidual() {
      return materialBalanceResidual;
    }

    public int getIterations() {
      return iterations;
    }

    public double getJacobianConditionNumber() {
      return jacobianConditionNumber;
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
