package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.Arrays;
import Jama.Matrix;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Coupled corrector for a point where two different phases become incipient from one mother phase.
 *
 * <p>
 * The unknown vector is {@code [ln(W_b), ln(W_c), ln(T), ln(P)]}. For each trial phase, {@code W} is the unnormalised
 * Michelsen stability variable and the EOS is evaluated at {@code x = W/sum(W)}. The equations are the component
 * fugacity equalities to the fixed-composition mother phase and {@code ln(sum(W)) = 0} for both trial phases.
 * Consequently, a converged solution is not merely a TP flash with three visible phases: both distinct incipient phases
 * independently have zero tangent-plane distance at the same temperature and pressure.
 * </p>
 */
public final class ThreePhasePointSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;
  private static final double MINIMUM_LOG_COMPOSITION = -700.0;
  private static final double MAXIMUM_LOG_COMPOSITION = 100.0;
  private static final double MINIMUM_TEMPERATURE_K = 50.0;
  private static final double MAXIMUM_TEMPERATURE_K = 2500.0;
  private static final double MINIMUM_PRESSURE_BARA = 1.0e-6;
  private static final double MAXIMUM_PRESSURE_BARA = 1.0e6;
  private static final double DISTINCT_PHASE_TOLERANCE = 1.0e-5;

  private final SystemInterface template;
  private final CandidatePhase motherPhase;
  private final CandidatePhase firstIncipientPhase;
  private final CandidatePhase secondIncipientPhase;
  private int maximumIterations = 40;
  private double residualTolerance = 1.0e-8;
  private double finiteDifferenceStep = 2.0e-5;

  /**
   * Creates a coupled three-phase point corrector.
   *
   * @param template fully configured EOS and mixing-rule system
   * @param motherPhase phase whose composition equals the fixed overall composition
   * @param firstIncipientPhase first distinct incipient phase
   * @param secondIncipientPhase second distinct incipient phase
   */
  public ThreePhasePointSolver(SystemInterface template, CandidatePhase motherPhase, CandidatePhase firstIncipientPhase,
      CandidatePhase secondIncipientPhase) {
    if (template == null || motherPhase == null || firstIncipientPhase == null || secondIncipientPhase == null) {
      throw new IllegalArgumentException("template and all phase families must be specified");
    }
    if (motherPhase == firstIncipientPhase || motherPhase == secondIncipientPhase
        || firstIncipientPhase == secondIncipientPhase) {
      throw new IllegalArgumentException("mother and incipient phase families must be distinct");
    }
    this.template = template.clone();
    this.motherPhase = motherPhase;
    this.firstIncipientPhase = firstIncipientPhase;
    this.secondIncipientPhase = secondIncipientPhase;
  }

  /**
   * Sets numerical controls.
   *
   * @param maximumIterations maximum damped Newton steps
   * @param residualTolerance infinity-norm residual tolerance
   * @param finiteDifferenceStep central-difference step in transformed variables
   * @return this solver
   */
  public ThreePhasePointSolver setNumericalControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    if (maximumIterations < 1 || !Double.isFinite(residualTolerance) || residualTolerance <= 0.0
        || !Double.isFinite(finiteDifferenceStep) || finiteDifferenceStep <= 0.0) {
      throw new IllegalArgumentException("invalid three-phase point numerical controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /**
   * Solves from physical T/P and automatically generated gas/oil/water trial compositions.
   *
   * @param temperatureK initial temperature in kelvin
   * @param pressureBara initial pressure in bara
   * @return immutable coupled solution and diagnostics
   */
  public Result solve(double temperatureK, double pressureBara) {
    validateTemperaturePressure(temperatureK, pressureBara);
    double[] firstComposition = initialComposition(firstIncipientPhase, temperatureK, pressureBara);
    double[] secondComposition = initialComposition(secondIncipientPhase, temperatureK, pressureBara);
    return solve(temperatureK, pressureBara, firstComposition, secondComposition);
  }

  /**
   * Solves from caller-supplied incipient-phase composition seeds.
   *
   * @param temperatureK initial temperature in kelvin
   * @param pressureBara initial pressure in bara
   * @param firstComposition first incipient phase mole fractions
   * @param secondComposition second incipient phase mole fractions
   * @return immutable coupled solution and diagnostics
   */
  public Result solve(double temperatureK, double pressureBara, double[] firstComposition, double[] secondComposition) {
    validateTemperaturePressure(temperatureK, pressureBara);
    int componentCount = template.getPhase(0).getNumberOfComponents();
    validateComposition(firstComposition, componentCount, "firstComposition");
    validateComposition(secondComposition, componentCount, "secondComposition");
    double[] variables = initialVariables(firstComposition, secondComposition, temperatureK, pressureBara);
    double[][] variableBounds = localVariableBounds(componentCount, temperatureK, pressureBara);
    Evaluation evaluation;
    try {
      evaluation = evaluate(variables);
    } catch (RuntimeException error) {
      return Result.failure(motherPhase, firstIncipientPhase, secondIncipientPhase, temperatureK, pressureBara,
          componentCount, 0, error.getMessage());
    }

    double initialResidual = evaluation.maximumResidual;
    String failureMessage = null;
    int iterations = 0;
    double finalJacobianConditionNumber = Double.NaN;
    while (iterations < maximumIterations && evaluation.maximumResidual > residualTolerance) {
      iterations++;
      Matrix jacobian;
      try {
        jacobian = numericalJacobian(variables, evaluation.residual, variableBounds[0], variableBounds[1]);
        finalJacobianConditionNumber = jacobian.cond();
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
      if (!Double.isFinite(correction.normInf())) {
        failureMessage = "Newton correction is non-finite";
        break;
      }
      limitCorrection(correction, componentCount);

      boolean accepted = false;
      double damping = 1.0;
      for (int lineSearch = 0; lineSearch < 16; lineSearch++) {
        double[] trialVariables = variables.clone();
        for (int variableIndex = 0; variableIndex < trialVariables.length; variableIndex++) {
          trialVariables[variableIndex] = clampToBounds(
              trialVariables[variableIndex] - damping * correction.get(variableIndex, 0),
              variableBounds[0][variableIndex], variableBounds[1][variableIndex]);
        }
        try {
          Evaluation trialEvaluation = evaluate(trialVariables);
          if (trialEvaluation.maximumResidual < evaluation.maximumResidual) {
            variables = trialVariables;
            evaluation = trialEvaluation;
            accepted = true;
            break;
          }
        } catch (RuntimeException error) {
          // Continue with a shorter line-search step.
        }
        damping *= 0.5;
      }
      if (!accepted) {
        failureMessage = "line search could not reduce the coupled residual";
        break;
      }
    }

    boolean residualConverged = evaluation.maximumResidual <= residualTolerance;
    double firstDistance = compositionDistance(evaluation.motherComposition, evaluation.firstComposition);
    double secondDistance = compositionDistance(evaluation.motherComposition, evaluation.secondComposition);
    double incipientDistance = compositionDistance(evaluation.firstComposition, evaluation.secondComposition);
    boolean distinctPhases = firstDistance > DISTINCT_PHASE_TOLERANCE && secondDistance > DISTINCT_PHASE_TOLERANCE
        && incipientDistance > DISTINCT_PHASE_TOLERANCE;
    if (failureMessage == null && !residualConverged) {
      failureMessage = "maximum coupled iteration count reached";
    }
    if (failureMessage == null && !distinctPhases) {
      failureMessage = "coupled solution is a trivial or duplicate phase solution";
    }
    return new Result(motherPhase, firstIncipientPhase, secondIncipientPhase, evaluation.temperatureK,
        evaluation.pressureBara, evaluation.motherComposition, evaluation.firstComposition,
        evaluation.secondComposition, initialResidual, evaluation.maximumResidual, iterations,
        finalJacobianConditionNumber, residualConverged, distinctPhases, failureMessage);
  }

  private Matrix numericalJacobian(double[] variables, double[] baseResidual, double[] lowerBounds,
      double[] upperBounds) {
    int variableCount = variables.length;
    Matrix jacobian = new Matrix(variableCount, variableCount);
    for (int column = 0; column < variableCount; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(variables[column]));
      double[] plus = variables.clone();
      double[] minus = variables.clone();
      plus[column] = clampToBounds(plus[column] + step, lowerBounds[column], upperBounds[column]);
      minus[column] = clampToBounds(minus[column] - step, lowerBounds[column], upperBounds[column]);
      double denominator = plus[column] - minus[column];
      if (!(denominator > 0.0)) {
        throw new IllegalStateException("finite-difference variable is pinned at a bound");
      }
      double[] plusResidual = evaluate(plus).residual;
      double[] minusResidual;
      try {
        minusResidual = evaluate(minus).residual;
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

  private Evaluation evaluate(double[] variables) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    double temperatureK = Math.exp(variables[2 * componentCount]);
    double pressureBara = Math.exp(variables[2 * componentCount + 1]);
    validateTemperaturePressure(temperatureK, pressureBara);

    double firstLogSum = logSumExp(variables, 0, componentCount);
    double secondLogSum = logSumExp(variables, componentCount, 2 * componentCount);
    double[] firstComposition = normalizedExponentials(variables, 0, componentCount, firstLogSum);
    double[] secondComposition = normalizedExponentials(variables, componentCount, 2 * componentCount, secondLogSum);
    double[] motherComposition = overallComposition();
    SystemInterface working = createThreePhaseSystem(temperatureK, pressureBara, motherComposition, firstComposition,
        secondComposition);

    double[] residual = new double[2 * componentCount + 2];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double reference = Math.log(Math.max(motherComposition[componentIndex], MINIMUM_COMPOSITION))
          + working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
      residual[componentIndex] = variables[componentIndex]
          + working.getPhase(1).getComponent(componentIndex).getLogFugacityCoefficient() - reference;
      residual[componentCount + componentIndex] = variables[componentCount + componentIndex]
          + working.getPhase(2).getComponent(componentIndex).getLogFugacityCoefficient() - reference;
    }
    residual[2 * componentCount] = firstLogSum;
    residual[2 * componentCount + 1] = secondLogSum;
    double maximumResidual = 0.0;
    for (double value : residual) {
      if (!Double.isFinite(value)) {
        throw new IllegalStateException("non-finite coupled residual");
      }
      maximumResidual = Math.max(maximumResidual, Math.abs(value));
    }
    return new Evaluation(residual, maximumResidual, temperatureK, pressureBara, motherComposition, firstComposition,
        secondComposition);
  }

  private SystemInterface createThreePhaseSystem(double temperatureK, double pressureBara, double[] motherComposition,
      double[] firstComposition, double[] secondComposition) {
    SystemInterface working = template.clone();
    working.setMultiPhaseCheck(false);
    working.setMaxNumberOfPhases(3);
    working.setNumberOfPhases(3);
    PhaseInterface phaseTemplate = template.getPhase(0);
    working.setPhase(phaseTemplate.clone(), 0);
    working.setPhase(phaseTemplate.clone(), 1);
    working.setPhase(phaseTemplate.clone(), 2);
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    setPhase(working, 0, motherPhase, motherComposition);
    setPhase(working, 1, firstIncipientPhase, firstComposition);
    setPhase(working, 2, secondIncipientPhase, secondComposition);
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      working.init(1, phaseIndex);
    }
    return working;
  }

  private static void setPhase(SystemInterface system, int phaseIndex, CandidatePhase phase, double[] composition) {
    system.setPhaseType(phaseIndex, toNeqSimPhaseType(phase));
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      system.getPhase(phaseIndex).getComponent(componentIndex).setx(composition[componentIndex]);
    }
    system.getPhase(phaseIndex).normalize();
    system.setPhaseType(phaseIndex, toNeqSimPhaseType(phase));
  }

  private double[] initialComposition(CandidatePhase phase, double temperatureK, double pressureBara) {
    if (phase != CandidatePhase.AQUEOUS) {
      double[] flashComposition = twoPhaseFlashComposition(phase, temperatureK, pressureBara);
      if (flashComposition != null) {
        return flashComposition;
      }
    }
    try {
      SystemInterface mother = template.clone();
      mother.setNumberOfPhases(1);
      mother.setTemperature(temperatureK);
      mother.setPressure(pressureBara);
      setPhase(mother, 0, motherPhase, overallComposition());
      mother.init(1);
      Candidate candidate = new IncipientPhaseStabilityAnalyzer(mother).setMaximumIterations(300).setDampingFactor(0.35)
          .analyze().getCandidate(phase);
      if (candidate != null && candidate.isConverged() && !candidate.isTrivial()) {
        return candidate.getComposition();
      }
    } catch (RuntimeException error) {
      // Fall through to a physical Wilson/water seed.
    }
    return physicalSeed(phase, temperatureK, pressureBara);
  }

  private double[] twoPhaseFlashComposition(CandidatePhase requestedPhase, double temperatureK, double pressureBara) {
    try {
      SystemInterface twoPhase = template.clone();
      twoPhase.setTemperature(temperatureK);
      twoPhase.setPressure(pressureBara);
      twoPhase.setMultiPhaseCheck(false);
      twoPhase.setMaxNumberOfPhases(2);
      if (twoPhase.getNumberOfPhases() > 2) {
        twoPhase.setNumberOfPhases(2);
      }
      new ThermodynamicOperations(twoPhase).TPflash();
      for (int phaseIndex = 0; phaseIndex < twoPhase.getNumberOfPhases(); phaseIndex++) {
        PhaseType phaseType = twoPhase.getPhase(phaseIndex).getType();
        boolean matches = requestedPhase == CandidatePhase.GAS && phaseType == PhaseType.GAS
            || requestedPhase == CandidatePhase.OIL && (phaseType == PhaseType.OIL || phaseType == PhaseType.LIQUID);
        if (!matches) {
          continue;
        }
        double[] composition = new double[twoPhase.getPhase(phaseIndex).getNumberOfComponents()];
        for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
          composition[componentIndex] = Math.max(twoPhase.getPhase(phaseIndex).getComponent(componentIndex).getx(),
              MINIMUM_COMPOSITION);
        }
        normalize(composition);
        if (compositionDistance(composition, overallComposition()) > DISTINCT_PHASE_TOLERANCE) {
          return composition;
        }
      }
    } catch (RuntimeException error) {
      return null;
    }
    return null;
  }

  private double[] physicalSeed(CandidatePhase phase, double temperatureK, double pressureBara) {
    double[] composition = overallComposition();
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      ComponentInterface component = template.getPhase(0).getComponent(componentIndex);
      double wilsonK = component.getPC() / pressureBara
          * Math.exp(5.373 * (1.0 + component.getAcentricFactor()) * (1.0 - component.getTC() / temperatureK));
      wilsonK = Math.max(1.0e-20, Math.min(1.0e20, wilsonK));
      if (phase == CandidatePhase.GAS) {
        composition[componentIndex] *= wilsonK;
      } else if (phase == CandidatePhase.OIL) {
        composition[componentIndex] /= wilsonK;
      } else if (component.getComponentName().equalsIgnoreCase("water")) {
        composition[componentIndex] = Math.max(composition[componentIndex], 0.99);
      } else if (component.isHydrocarbon()) {
        composition[componentIndex] *= 1.0e-8;
      }
      composition[componentIndex] = Math.max(composition[componentIndex], MINIMUM_COMPOSITION);
    }
    normalize(composition);
    return composition;
  }

  private double[] overallComposition() {
    double[] composition = new double[template.getPhase(0).getNumberOfComponents()];
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(),
          MINIMUM_COMPOSITION);
    }
    normalize(composition);
    return composition;
  }

  private static double[] initialVariables(double[] firstComposition, double[] secondComposition, double temperatureK,
      double pressureBara) {
    double[] first = firstComposition.clone();
    double[] second = secondComposition.clone();
    normalize(first);
    normalize(second);
    double[] variables = new double[2 * first.length + 2];
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      variables[componentIndex] = Math.log(Math.max(first[componentIndex], MINIMUM_COMPOSITION));
      variables[first.length + componentIndex] = Math.log(Math.max(second[componentIndex], MINIMUM_COMPOSITION));
    }
    variables[2 * first.length] = Math.log(temperatureK);
    variables[2 * first.length + 1] = Math.log(pressureBara);
    return variables;
  }

  private static double[][] localVariableBounds(int componentCount, double temperatureK, double pressureBara) {
    int variableCount = 2 * componentCount + 2;
    double[] lower = new double[variableCount];
    double[] upper = new double[variableCount];
    Arrays.fill(lower, MINIMUM_LOG_COMPOSITION);
    Arrays.fill(upper, MAXIMUM_LOG_COMPOSITION);
    lower[2 * componentCount] = Math.log(Math.max(MINIMUM_TEMPERATURE_K, 0.8 * temperatureK));
    upper[2 * componentCount] = Math.log(Math.min(MAXIMUM_TEMPERATURE_K, 1.2 * temperatureK));
    lower[2 * componentCount + 1] = Math.log(Math.max(MINIMUM_PRESSURE_BARA, 0.5 * pressureBara));
    upper[2 * componentCount + 1] = Math.log(Math.min(MAXIMUM_PRESSURE_BARA, 1.5 * pressureBara));
    return new double[][] {lower, upper};
  }

  private static double clampToBounds(double value, double lower, double upper) {
    return Math.max(lower, Math.min(upper, value));
  }

  private static void limitCorrection(Matrix correction, int componentCount) {
    double scale = 1.0;
    for (int variableIndex = 0; variableIndex < correction.getRowDimension(); variableIndex++) {
      double maximumStep = variableIndex < 2 * componentCount ? 2.0 : 0.08;
      double magnitude = Math.abs(correction.get(variableIndex, 0));
      if (magnitude > maximumStep) {
        scale = Math.min(scale, maximumStep / magnitude);
      }
    }
    if (scale < 1.0) {
      correction.timesEquals(scale);
    }
  }

  private static double logSumExp(double[] values, int start, int end) {
    double maximum = Double.NEGATIVE_INFINITY;
    for (int index = start; index < end; index++) {
      maximum = Math.max(maximum, values[index]);
    }
    double total = 0.0;
    for (int index = start; index < end; index++) {
      total += Math.exp(values[index] - maximum);
    }
    return maximum + Math.log(total);
  }

  private static double[] normalizedExponentials(double[] values, int start, int end, double logSum) {
    double[] composition = new double[end - start];
    for (int index = start; index < end; index++) {
      composition[index - start] = Math.exp(values[index] - logSum);
    }
    return composition;
  }

  private static void normalize(double[] values) {
    double total = 0.0;
    for (double value : values) {
      total += value;
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalArgumentException("composition cannot be normalized");
    }
    for (int index = 0; index < values.length; index++) {
      values[index] /= total;
    }
  }

  private static void validateComposition(double[] composition, int componentCount, String name) {
    if (composition == null || composition.length != componentCount) {
      throw new IllegalArgumentException(name + " must contain one entry per component");
    }
    for (double value : composition) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException(name + " contains an invalid mole fraction");
      }
    }
  }

  private static void validateTemperaturePressure(double temperatureK, double pressureBara) {
    if (!Double.isFinite(temperatureK) || temperatureK < MINIMUM_TEMPERATURE_K || temperatureK > MAXIMUM_TEMPERATURE_K
        || !Double.isFinite(pressureBara) || pressureBara < MINIMUM_PRESSURE_BARA
        || pressureBara > MAXIMUM_PRESSURE_BARA) {
      throw new IllegalArgumentException("temperature or pressure is outside the solver domain");
    }
  }

  private static PhaseType toNeqSimPhaseType(CandidatePhase phase) {
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

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }

  private static final class Evaluation {
    private final double[] residual;
    private final double maximumResidual;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] motherComposition;
    private final double[] firstComposition;
    private final double[] secondComposition;

    private Evaluation(double[] residual, double maximumResidual, double temperatureK, double pressureBara,
        double[] motherComposition, double[] firstComposition, double[] secondComposition) {
      this.residual = residual;
      this.maximumResidual = maximumResidual;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.motherComposition = motherComposition;
      this.firstComposition = firstComposition;
      this.secondComposition = secondComposition;
    }
  }

  /** Immutable coupled three-phase point solution and numerical diagnostics. */
  public static final class Result {
    private final CandidatePhase motherPhase;
    private final CandidatePhase firstIncipientPhase;
    private final CandidatePhase secondIncipientPhase;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] motherComposition;
    private final double[] firstComposition;
    private final double[] secondComposition;
    private final double initialMaximumResidual;
    private final double maximumResidual;
    private final int iterations;
    private final double jacobianConditionNumber;
    private final boolean residualConverged;
    private final boolean distinctPhases;
    private final String failureMessage;

    private Result(CandidatePhase motherPhase, CandidatePhase firstIncipientPhase, CandidatePhase secondIncipientPhase,
        double temperatureK, double pressureBara, double[] motherComposition, double[] firstComposition,
        double[] secondComposition, double initialMaximumResidual, double maximumResidual, int iterations,
        double jacobianConditionNumber, boolean residualConverged, boolean distinctPhases, String failureMessage) {
      this.motherPhase = motherPhase;
      this.firstIncipientPhase = firstIncipientPhase;
      this.secondIncipientPhase = secondIncipientPhase;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.motherComposition = motherComposition.clone();
      this.firstComposition = firstComposition.clone();
      this.secondComposition = secondComposition.clone();
      this.initialMaximumResidual = initialMaximumResidual;
      this.maximumResidual = maximumResidual;
      this.iterations = iterations;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.residualConverged = residualConverged;
      this.distinctPhases = distinctPhases;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase motherPhase, CandidatePhase firstIncipientPhase,
        CandidatePhase secondIncipientPhase, double temperatureK, double pressureBara, int componentCount,
        int iterations, String failureMessage) {
      return new Result(motherPhase, firstIncipientPhase, secondIncipientPhase, temperatureK, pressureBara,
          new double[componentCount], new double[componentCount], new double[componentCount], Double.NaN, Double.NaN,
          iterations, Double.NaN, false, false, failureMessage);
    }

    /** @return true only when residual and distinct-phase gates both pass */
    public boolean isConverged() {
      return residualConverged && distinctPhases && failureMessage == null;
    }

    /** @return mother phase family */
    public CandidatePhase getMotherPhase() {
      return motherPhase;
    }

    /** @return first incipient phase family */
    public CandidatePhase getFirstIncipientPhase() {
      return firstIncipientPhase;
    }

    /** @return second incipient phase family */
    public CandidatePhase getSecondIncipientPhase() {
      return secondIncipientPhase;
    }

    /** @return temperature in kelvin */
    public double getTemperatureK() {
      return temperatureK;
    }

    /** @return pressure in bara */
    public double getPressureBara() {
      return pressureBara;
    }

    /** @return fixed mother composition */
    public double[] getMotherComposition() {
      return motherComposition.clone();
    }

    /** @return first incipient composition */
    public double[] getFirstComposition() {
      return firstComposition.clone();
    }

    /** @return second incipient composition */
    public double[] getSecondComposition() {
      return secondComposition.clone();
    }

    /** @return initial infinity-norm residual */
    public double getInitialMaximumResidual() {
      return initialMaximumResidual;
    }

    /** @return final infinity-norm residual */
    public double getMaximumResidual() {
      return maximumResidual;
    }

    /** @return damped Newton iteration count */
    public int getIterations() {
      return iterations;
    }

    /** @return final numerical Jacobian condition number */
    public double getJacobianConditionNumber() {
      return jacobianConditionNumber;
    }

    /** @return true when the residual tolerance was met */
    public boolean isResidualConverged() {
      return residualConverged;
    }

    /** @return true when all three compositions are physically distinct */
    public boolean hasDistinctPhases() {
      return distinctPhases;
    }

    /** @return failure diagnostic, or {@code null} */
    public String getFailureMessage() {
      return failureMessage;
    }

    @Override
    public String toString() {
      return "ThreePhasePointResult{T=" + temperatureK + ", P=" + pressureBara + ", residual=" + maximumResidual
          + ", initialResidual=" + initialMaximumResidual + ", iterations=" + iterations + ", phases=" + motherPhase
          + "/" + firstIncipientPhase + "/" + secondIncipientPhase + ", converged=" + isConverged() + ", failure="
          + failureMessage + ", xMother=" + Arrays.toString(motherComposition) + "}";
    }
  }
}
