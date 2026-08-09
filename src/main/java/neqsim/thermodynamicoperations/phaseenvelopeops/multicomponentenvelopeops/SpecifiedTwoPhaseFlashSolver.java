package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import Jama.Matrix;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Damped Newton two-phase TP flash that preserves caller-specified physical phase families. */
public final class SpecifiedTwoPhaseFlashSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;
  private static final double MINIMUM_BETA = 1.0e-10;

  private final SystemInterface template;
  private final SystemInterface working;
  private final CandidatePhase phaseZero;
  private final CandidatePhase phaseOne;
  private int maximumIterations = 60;
  private double residualTolerance = 1.0e-9;
  private double finiteDifferenceStep = 2.0e-5;

  /**
   * Creates a specified-phase flash solver.
   *
   * @param template configured EOS/mixing-rule system
   * @param phaseZero first physical phase family; its mole fraction is beta
   * @param phaseOne second physical phase family
   */
  public SpecifiedTwoPhaseFlashSolver(SystemInterface template, CandidatePhase phaseZero, CandidatePhase phaseOne) {
    if (template == null || phaseZero == null || phaseOne == null) {
      throw new IllegalArgumentException("template and both phase families must be specified");
    }
    if (phaseZero == phaseOne) {
      throw new IllegalArgumentException("specified phase families must be distinct");
    }
    this.template = template.clone();
    this.working = this.template.clone();
    this.working.setNumberOfPhases(2);
    this.phaseZero = phaseZero;
    this.phaseOne = phaseOne;
  }

  /** Sets Newton and numerical-Jacobian controls. */
  public SpecifiedTwoPhaseFlashSolver setNumericalControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    if (maximumIterations < 1 || !Double.isFinite(residualTolerance) || residualTolerance <= 0.0
        || !Double.isFinite(finiteDifferenceStep) || finiteDifferenceStep <= 0.0) {
      throw new IllegalArgumentException("invalid specified two-phase flash numerical controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /**
   * Solves at fixed T/P from beta and two composition seeds.
   *
   * @param temperatureK temperature in kelvin
   * @param pressureBara pressure in bara
   * @param betaSeed phase-zero mole fraction in (0,1)
   * @param phaseZeroComposition phase-zero seed
   * @param phaseOneComposition phase-one seed
   * @return immutable specified-phase flash result
   */
  public Result solve(double temperatureK, double pressureBara, double betaSeed, double[] phaseZeroComposition,
      double[] phaseOneComposition) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    validateInputs(temperatureK, pressureBara, betaSeed, phaseZeroComposition, phaseOneComposition, componentCount);
    double[] variables = new double[componentCount + 1];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      variables[componentIndex] = Math.log(Math.max(phaseZeroComposition[componentIndex], MINIMUM_COMPOSITION)
          / Math.max(phaseOneComposition[componentIndex], MINIMUM_COMPOSITION));
    }
    variables[componentCount] = betaSeed;

    Evaluation evaluation;
    try {
      evaluation = evaluate(temperatureK, pressureBara, variables);
    } catch (RuntimeException error) {
      return Result.failure(phaseZero, phaseOne, temperatureK, pressureBara, componentCount, error.getMessage());
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
        correction = jacobian.solve(new Matrix(evaluation.residual, evaluation.residual.length));
      } catch (RuntimeException error) {
        failureMessage = "Jacobian solve failed: " + error.getMessage();
        break;
      }
      limitCorrection(correction, componentCount);
      boolean accepted = false;
      double damping = 1.0;
      for (int lineSearch = 0; lineSearch < 16; lineSearch++) {
        double[] trialVariables = variables.clone();
        for (int variableIndex = 0; variableIndex < componentCount; variableIndex++) {
          trialVariables[variableIndex] = Math.max(-700.0,
              Math.min(700.0, variables[variableIndex] - damping * correction.get(variableIndex, 0)));
        }
        trialVariables[componentCount] = Math.max(MINIMUM_BETA,
            Math.min(1.0 - MINIMUM_BETA, variables[componentCount] - damping * correction.get(componentCount, 0)));
        try {
          Evaluation trial = evaluate(temperatureK, pressureBara, trialVariables);
          if (trial.maximumResidual < evaluation.maximumResidual) {
            variables = trialVariables;
            evaluation = trial;
            accepted = true;
            break;
          }
        } catch (RuntimeException error) {
          // Try a shorter step.
        }
        damping *= 0.5;
      }
      if (!accepted) {
        failureMessage = "line search could not reduce the specified-phase flash residual";
        break;
      }
    }
    boolean converged = evaluation.maximumResidual <= residualTolerance;
    if (failureMessage == null && !converged) {
      failureMessage = "maximum specified-phase flash iteration count reached";
    }
    return new Result(phaseZero, phaseOne, temperatureK, pressureBara, evaluation.beta, evaluation.phaseZeroComposition,
        evaluation.phaseOneComposition, initialResidual, evaluation.maximumResidual, iterations, conditionNumber,
        converged, failureMessage);
  }

  /**
   * Reconstructs a converged result as a private thermodynamic system for third-phase stability analysis.
   *
   * @param result result produced by this solver with matching phase families
   * @return independent two-phase thermodynamic system
   */
  public SystemInterface toThermodynamicSystem(Result result) {
    if (result == null || result.getPhaseZero() != phaseZero || result.getPhaseOne() != phaseOne
        || !result.isConverged()) {
      throw new IllegalArgumentException("a converged result with matching phase families is required");
    }
    return createWorkingSystem(result.getTemperatureK(), result.getPressureBara(), result.getPhaseZeroComposition(),
        result.getPhaseOneComposition(), result.getBeta()).clone();
  }

  private Matrix numericalJacobian(double temperatureK, double pressureBara, double[] variables,
      double[] baseResidual) {
    int variableCount = variables.length;
    Matrix jacobian = new Matrix(variableCount, variableCount);
    for (int column = 0; column < variableCount; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(variables[column]));
      double[] plus = variables.clone();
      double[] minus = variables.clone();
      if (column == variableCount - 1) {
        plus[column] = Math.min(1.0 - MINIMUM_BETA, plus[column] + step);
        minus[column] = Math.max(MINIMUM_BETA, minus[column] - step);
      } else {
        plus[column] += step;
        minus[column] -= step;
      }
      double denominator = plus[column] - minus[column];
      double[] plusResidual = evaluate(temperatureK, pressureBara, plus).residual;
      double[] minusResidual;
      try {
        minusResidual = evaluate(temperatureK, pressureBara, minus).residual;
      } catch (RuntimeException error) {
        minusResidual = baseResidual;
        denominator = plus[column] - variables[column];
      }
      for (int row = 0; row < variableCount; row++) {
        jacobian.set(row, column, (plusResidual[row] - minusResidual[row]) / denominator);
      }
    }
    return jacobian;
  }

  private Evaluation evaluate(double temperatureK, double pressureBara, double[] variables) {
    int componentCount = variables.length - 1;
    double beta = variables[componentCount];
    if (!(beta > 0.0 && beta < 1.0)) {
      throw new IllegalStateException("phase fraction left the physical interval");
    }
    // Evaluated in log space so the composition floor never enters the Rachford-Rice residual. A heavy pseudo
    // component against water reaches K ~ 1e99, which puts its true liquid mole fraction near 1e-102 - below the
    // floor. Flooring it there propagates into first[] = K * second[] as a large spurious inventory and pins the
    // residual to the floor rather than to the equations, so the line search cannot reduce it no matter how far
    // it backtracks. The floor is applied only where the EOS needs a strictly positive mole fraction.
    double[] z = overallComposition();
    double[] first = new double[componentCount];
    double[] second = new double[componentCount];
    double[] logFirst = new double[componentCount];
    double[] logSecond = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double logK = Math.max(-700.0, Math.min(700.0, variables[componentIndex]));
      double logDenominator = rachfordRiceLogDenominator(beta, logK);
      if (!Double.isFinite(logDenominator)) {
        throw new IllegalStateException("invalid Rachford-Rice denominator");
      }
      logSecond[componentIndex] = Math.log(z[componentIndex]) - logDenominator;
      logFirst[componentIndex] = logK + logSecond[componentIndex];
      second[componentIndex] = Math.exp(logSecond[componentIndex]);
      first[componentIndex] = Math.exp(logFirst[componentIndex]);
    }
    double sumFirst = Math.exp(logSumExp(logFirst));
    double sumSecond = Math.exp(logSumExp(logSecond));
    if (!Double.isFinite(sumFirst) || !Double.isFinite(sumSecond) || !(sumFirst > 0.0) || !(sumSecond > 0.0)) {
      throw new IllegalStateException("specified-phase flash inventory left the physical domain");
    }
    SystemInterface working = createWorkingSystem(temperatureK, pressureBara, first, second, beta);
    double[] residual = new double[componentCount + 1];
    double maximumResidual = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      residual[componentIndex] = variables[componentIndex]
          + working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient()
          - working.getPhase(1).getComponent(componentIndex).getLogFugacityCoefficient();
      maximumResidual = Math.max(maximumResidual, Math.abs(residual[componentIndex]));
    }
    residual[componentCount] = sumFirst - sumSecond;
    maximumResidual = Math.max(maximumResidual, Math.abs(residual[componentCount]));
    if (!Double.isFinite(maximumResidual)) {
      throw new IllegalStateException("specified-phase flash residual is non-finite");
    }
    return new Evaluation(residual, maximumResidual, beta, normalizedComposition(working, 0, componentCount),
        normalizedComposition(working, 1, componentCount));
  }

  /** Evaluates {@code log(1 - beta + beta * exp(logK))} without forming {@code exp(logK)} when it would overflow. */
  private static double rachfordRiceLogDenominator(double beta, double logK) {
    double liquidFraction = 1.0 - beta;
    if (logK > 0.0) {
      return logK + Math.log(beta + liquidFraction * Math.exp(-logK));
    }
    return Math.log(liquidFraction + beta * Math.exp(logK));
  }

  private static double logSumExp(double[] values) {
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

  private static double[] normalizedComposition(SystemInterface working, int phaseIndex, int componentCount) {
    double[] composition = new double[componentCount];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = Math.max(working.getPhase(phaseIndex).getComponent(componentIndex).getx(),
          MINIMUM_COMPOSITION);
      total += composition[componentIndex];
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalStateException("specified-phase flash composition cannot be normalized");
    }
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }

  private SystemInterface createWorkingSystem(double temperatureK, double pressureBara, double[] first, double[] second,
      double beta) {
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    working.setPhaseType(0, toPhaseType(phaseZero));
    working.setPhaseType(1, toPhaseType(phaseOne));
    working.setBeta(0, beta);
    working.setBeta(1, 1.0 - beta);
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      // The EOS requires a strictly positive mole fraction; this is the only place the floor belongs. A component
      // this dilute is indistinguishable from zero for mixture properties, whereas flooring the inventory itself
      // would corrupt the Rachford-Rice residual being driven to zero.
      working.getPhase(0).getComponent(componentIndex).setx(Math.max(first[componentIndex], MINIMUM_COMPOSITION));
      working.getPhase(1).getComponent(componentIndex).setx(Math.max(second[componentIndex], MINIMUM_COMPOSITION));
    }
    working.getPhase(0).normalize();
    working.getPhase(1).normalize();
    working.setPhaseType(0, toPhaseType(phaseZero));
    working.setPhaseType(1, toPhaseType(phaseOne));
    working.init(1, 0);
    working.init(1, 1);
    return working;
  }

  private double[] overallComposition() {
    double[] z = new double[template.getPhase(0).getNumberOfComponents()];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < z.length; componentIndex++) {
      z[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(), MINIMUM_COMPOSITION);
      total += z[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < z.length; componentIndex++) {
      z[componentIndex] /= total;
    }
    return z;
  }

  private static void limitCorrection(Matrix correction, int componentCount) {
    double scale = 1.0;
    for (int variableIndex = 0; variableIndex < correction.getRowDimension(); variableIndex++) {
      double maximumStep = variableIndex < componentCount ? 2.0 : 0.1;
      double magnitude = Math.abs(correction.get(variableIndex, 0));
      if (magnitude > maximumStep) {
        scale = Math.min(scale, maximumStep / magnitude);
      }
    }
    if (scale < 1.0) {
      correction.timesEquals(scale);
    }
  }

  private static void validateInputs(double temperatureK, double pressureBara, double betaSeed, double[] first,
      double[] second, int componentCount) {
    if (!Double.isFinite(temperatureK) || temperatureK < 50.0 || !Double.isFinite(pressureBara) || pressureBara <= 0.0
        || !Double.isFinite(betaSeed) || betaSeed <= 0.0 || betaSeed >= 1.0 || first == null || second == null
        || first.length != componentCount || second.length != componentCount) {
      throw new IllegalArgumentException("invalid specified two-phase flash input");
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

  private static final class Evaluation {
    private final double[] residual;
    private final double maximumResidual;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;

    private Evaluation(double[] residual, double maximumResidual, double beta, double[] phaseZeroComposition,
        double[] phaseOneComposition) {
      this.residual = residual;
      this.maximumResidual = maximumResidual;
      this.beta = beta;
      this.phaseZeroComposition = phaseZeroComposition;
      this.phaseOneComposition = phaseOneComposition;
    }
  }

  /** Immutable specified-phase TP flash result. */
  public static final class Result {
    private final CandidatePhase phaseZero;
    private final CandidatePhase phaseOne;
    private final double temperatureK;
    private final double pressureBara;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double initialMaximumResidual;
    private final double maximumResidual;
    private final int iterations;
    private final double jacobianConditionNumber;
    private final boolean converged;
    private final String failureMessage;

    private Result(CandidatePhase phaseZero, CandidatePhase phaseOne, double temperatureK, double pressureBara,
        double beta, double[] phaseZeroComposition, double[] phaseOneComposition, double initialMaximumResidual,
        double maximumResidual, int iterations, double jacobianConditionNumber, boolean converged,
        String failureMessage) {
      this.phaseZero = phaseZero;
      this.phaseOne = phaseOne;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.beta = beta;
      this.phaseZeroComposition = phaseZeroComposition.clone();
      this.phaseOneComposition = phaseOneComposition.clone();
      this.initialMaximumResidual = initialMaximumResidual;
      this.maximumResidual = maximumResidual;
      this.iterations = iterations;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.converged = converged;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase phaseZero, CandidatePhase phaseOne, double temperatureK,
        double pressureBara, int componentCount, String failureMessage) {
      return new Result(phaseZero, phaseOne, temperatureK, pressureBara, Double.NaN, new double[componentCount],
          new double[componentCount], Double.NaN, Double.NaN, 0, Double.NaN, false, failureMessage);
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

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double getBeta() {
      return beta;
    }

    public double[] getPhaseZeroComposition() {
      return phaseZeroComposition.clone();
    }

    public double[] getPhaseOneComposition() {
      return phaseOneComposition.clone();
    }

    public double getInitialMaximumResidual() {
      return initialMaximumResidual;
    }

    public double getMaximumResidual() {
      return maximumResidual;
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
  }
}
