package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import Jama.Matrix;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Bordered Newton corrector that grows a finite phase from a retained-phase spinodal.
 *
 * <p>
 * A prescribed phase fraction does not remove the coincident-phase solution manifold. This corrector instead fixes the
 * signed separation of the retained and bifurcating compositions along the spinodal eigenvector, while solving the new
 * phase fraction as an additional unknown. At fixed pressure, temperature, all compositions, both retained fractions,
 * and the new phase fraction are corrected simultaneously.
 *
 *
 * @author NeqSim contributors
 * @version 1.0
 */
public final class RetainedPhaseBifurcationCorrector {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;
  private static final double MINIMUM_BETA = 1.0e-10;

  private final SystemInterface template;
  private final Object resultOwner = new Object();
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase bifurcatingPhase;
  private int maximumIterations = 120;
  private double residualTolerance = 1.0e-9;
  private double finiteDifferenceStep = 2.0e-5;

  /**
   * Creates a bordered corrector; the bifurcating family must match one retained phase.
   *
   * @param template configured fluid copied privately by this solver
   * @param retainedPhaseZero retained phase zero
   * @param retainedPhaseOne retained phase one
   * @param bifurcatingPhase bifurcating phase
   */
  public RetainedPhaseBifurcationCorrector(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase bifurcatingPhase) {
    if (template == null || template.getNumberOfComponents() < 2 || retainedPhaseZero == null
        || retainedPhaseOne == null || bifurcatingPhase == null) {
      throw new IllegalArgumentException("template and all phase identities are required");
    }
    if (retainedPhaseZero == retainedPhaseOne
        || bifurcatingPhase != retainedPhaseZero && bifurcatingPhase != retainedPhaseOne) {
      throw new IllegalArgumentException("bifurcating phase must match exactly one retained phase family");
    }
    this.template = template.clone();
    this.template.init(0);
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.bifurcatingPhase = bifurcatingPhase;
  }

  /**
   * Sets damped-Newton and numerical-Jacobian controls.
   *
   * @param maximumIterations positive nonlinear iteration limit
   * @param residualTolerance positive maximum equilibrium residual
   * @param finiteDifferenceStep positive relative Jacobian perturbation
   * @return computed set numerical controls result
   */
  public RetainedPhaseBifurcationCorrector setNumericalControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    if (maximumIterations < 1 || !positive(residualTolerance) || !positive(finiteDifferenceStep)) {
      throw new IllegalArgumentException("invalid retained-phase bifurcation controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /**
   * Corrects a finite separation from one retained two-phase state at fixed pressure.
   *
   * @param retainedState converged retained two-phase state with a bifurcating composition seed in its third slot
   * @param mode spinodal stationarity-Jacobian mode
   * @param pressureBara fixed pressure
   * @param targetSeparation signed log-composition separation along the mode
   * @param bifurcatingFractionSeed strictly positive seed for the new phase fraction
   *
   * @return computed correct result
   */
  public Result correct(TwoToThreePhaseArcLengthCorrector.State retainedState,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode, double pressureBara, double targetSeparation,
      double bifurcatingFractionSeed) {
    validateInputs(retainedState, mode, pressureBara, targetSeparation, bifurcatingFractionSeed);
    double[] variables = variables(retainedState, bifurcatingFractionSeed);
    variables[2 * componentCount() + 2] = Math.log(pressureBara);
    clampVariables(variables);
    Evaluation evaluation;
    SystemInterface working;
    try {
      working = createWorkingSystem(retainedState.getTemperatureK(), pressureBara, retainedState.getBeta(),
          bifurcatingFractionSeed, retainedState.getPhaseZeroComposition(), retainedState.getPhaseOneComposition(),
          retainedState.getIncipientComposition());
      evaluation = evaluate(working, variables, mode, pressureBara, targetSeparation);
    } catch (RuntimeException error) {
      return Result.failure(resultOwner, retainedPhaseZero, retainedPhaseOne, bifurcatingPhase,
          retainedState.getTemperatureK(), pressureBara, componentCount(), error.getMessage());
    }
    double initialResidual = evaluation.maximumResidual;
    double conditionNumber = Double.NaN;
    int iterations = 0;
    String failureMessage = null;
    while (iterations < maximumIterations && evaluation.maximumResidual > residualTolerance) {
      iterations++;
      Matrix jacobian;
      try {
        jacobian = numericalJacobian(working, variables, evaluation.residual, mode, pressureBara, targetSeparation);
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
      limitCorrection(correction);
      boolean accepted = false;
      double damping = 1.0;
      for (int lineSearch = 0; lineSearch < 24; lineSearch++) {
        double[] trialVariables = variables.clone();
        for (int variableIndex = 0; variableIndex < variables.length; variableIndex++) {
          trialVariables[variableIndex] -= damping * correction.get(variableIndex, 0);
        }
        clampVariables(trialVariables);
        try {
          Evaluation trial = evaluate(working, trialVariables, mode, pressureBara, targetSeparation);
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
        failureMessage = "line search could not reduce the bordered bifurcation residual";
        break;
      }
    }
    boolean residualConverged = evaluation.maximumResidual <= residualTolerance;
    if (failureMessage == null && !residualConverged) {
      failureMessage = "maximum bordered bifurcation iteration count reached";
    }
    boolean distinct = compositionDistance(evaluation.duplicatedRetainedComposition,
        evaluation.bifurcatingComposition) > 1.0e-5;
    if (failureMessage == null && !distinct) {
      failureMessage = "bordered bifurcation correction returned coincident phase compositions";
    }
    boolean positiveFraction = evaluation.bifurcatingFraction > 1.0e-12
        && evaluation.bifurcatingFraction < 1.0 - 1.0e-12;
    if (failureMessage == null && !positiveFraction) {
      failureMessage = "bordered bifurcation correction returned a vanishing phase fraction";
    }
    if (failureMessage == null && (!Double.isFinite(evaluation.materialBalanceResidual)
        || evaluation.materialBalanceResidual > Math.max(residualTolerance, 1.0e-10))) {
      failureMessage = "normalized result failed component material balance";
    }
    return new Result(resultOwner, retainedPhaseZero, retainedPhaseOne, bifurcatingPhase, evaluation.temperatureK,
        evaluation.pressureBara, evaluation.retainedBeta, evaluation.bifurcatingFraction,
        evaluation.phaseZeroComposition, evaluation.phaseOneComposition, evaluation.bifurcatingComposition,
        initialResidual, evaluation.maximumResidual, evaluation.equilibriumMaximumResidual, evaluation.pressureResidual,
        evaluation.separationResidual, evaluation.materialBalanceResidual, iterations, conditionNumber,
        residualConverged, distinct, positiveFraction, failureMessage);
  }

  /**
   * Reconstructs a private numerical root for diagnostics; this does not certify phase stability.
   *
   * @param result converged result produced by this corrector instance
   * @return private system carrying the three phase compositions and fractions
   * @throws IllegalArgumentException when the result did not converge or has another owner
   */
  public SystemInterface toThermodynamicSystem(Result result) {
    if (result == null || result.owner != resultOwner || !result.isConverged()) {
      throw new IllegalArgumentException("a converged result from this corrector is required");
    }
    return createWorkingSystem(result.temperatureK, result.pressureBara, result.retainedBeta,
        result.bifurcatingFraction, result.phaseZeroComposition, result.phaseOneComposition,
        result.bifurcatingComposition);
  }

  /**
   * Independently checks a finite-phase root, including phase identities and stability searches.
   *
   * @param result converged result produced by this corrector instance
   * @return scoped physical acceptance diagnostics
   * @throws IllegalArgumentException when the result cannot be reconstructed
   */
  public SpecifiedPhaseEquilibriumValidator.Result validateEquilibrium(Result result) {
    return SpecifiedPhaseEquilibriumValidator.validate(template, toThermodynamicSystem(result),
        new CandidatePhase[] {retainedPhaseZero, retainedPhaseOne, bifurcatingPhase});
  }

  /**
   * Exports a finite-phase root only after replay and stability checks pass.
   *
   * @param result converged result produced by this corrector instance
   * @return private accepted thermodynamic system
   * @throws IllegalArgumentException when the result cannot be reconstructed
   * @throws IllegalStateException when physical acceptance fails
   */
  public SystemInterface toValidatedThermodynamicSystem(Result result) {
    SpecifiedPhaseEquilibriumValidator.Result validation = validateEquilibrium(result);
    if (!validation.isAccepted()) {
      throw new IllegalStateException("bifurcation equilibrium rejected: " + validation.getViolations());
    }
    return toThermodynamicSystem(result);
  }

  /**
   * Evaluates the residual Jacobian using central or one-sided finite differences.
   *
   * @param variables variables
   * @param baseResidual base residual
   * @param mode mode
   * @param pressureBara absolute pressure in bara
   * @param targetSeparation signed log-composition projection along the selected mode
   * @return computed numerical jacobian result
   *
   * @param working working
   */
  private Matrix numericalJacobian(SystemInterface working, double[] variables, double[] baseResidual,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode, double pressureBara, double targetSeparation) {
    Matrix jacobian = new Matrix(variables.length, variables.length);
    for (int column = 0; column < variables.length; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(variables[column]));
      double[] plus = variables.clone();
      double[] minus = variables.clone();
      plus[column] += step;
      minus[column] -= step;
      clampVariables(plus);
      clampVariables(minus);
      double[] plusResidual = null;
      double[] minusResidual = null;
      try {
        plusResidual = evaluate(working, plus, mode, pressureBara, targetSeparation).residual;
      } catch (RuntimeException error) {
        // One-sided difference remains available.
      }
      try {
        minusResidual = evaluate(working, minus, mode, pressureBara, targetSeparation).residual;
      } catch (RuntimeException error) {
        // One-sided difference remains available.
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
      for (int row = 0; row < variables.length; row++) {
        jacobian.set(row, column, (plusResidual[row] - minusResidual[row]) / denominator);
      }
    }
    return jacobian;
  }

  /**
   * Evaluates the nonlinear equilibrium equations at the supplied solver coordinates.
   *
   * @param variables variables
   * @param mode mode
   * @param specifiedPressureBara fixed absolute pressure in bara
   * @param targetSeparation signed log-composition projection along the selected mode
   * @return computed evaluate result
   *
   * @param working working
   */
  private Evaluation evaluate(SystemInterface working, double[] variables,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode, double specifiedPressureBara, double targetSeparation) {
    int componentCount = componentCount();
    int betaIndex = componentCount;
    int logWStart = componentCount + 1;
    int temperatureIndex = 2 * componentCount + 1;
    int pressureIndex = 2 * componentCount + 2;
    int fractionIndex = 2 * componentCount + 3;
    double retainedBeta = variables[betaIndex];
    double temperatureK = Math.exp(variables[temperatureIndex]);
    double pressureBara = Math.exp(variables[pressureIndex]);
    double bifurcatingFraction = logistic(variables[fractionIndex]);
    if (!(retainedBeta > 0.0 && retainedBeta < 1.0) || !positive(temperatureK) || temperatureK < 50.0
        || temperatureK > 2500.0 || !positive(pressureBara) || pressureBara > 1.0e6
        || !(bifurcatingFraction > 0.0 && bifurcatingFraction < 1.0)) {
      throw new IllegalStateException("bordered bifurcation variables left the physical domain");
    }
    double logSumW = logSumExp(variables, logWStart, logWStart + componentCount);
    double[] third = normalizedExponentials(variables, logWStart, logWStart + componentCount, logSumW);
    double[] overall = overallComposition();
    // Evaluated in log space so the composition floor never enters the material balance. A heavy pseudo-component
    // against water reaches K ~ 1e99 and a true retained-liquid fraction near 1e-102, below the floor; flooring it
    // propagates into first[] = K * second[] as a large spurious inventory and makes the converged point a root of
    // a perturbed system. The floor belongs only where the EOS needs a strictly positive mole fraction.
    double[] first = new double[componentCount];
    double[] second = new double[componentCount];
    double[] logFirst = new double[componentCount];
    double[] logSecond = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double retainedOverall = (overall[componentIndex] - bifurcatingFraction * third[componentIndex])
          / (1.0 - bifurcatingFraction);
      if (!(retainedOverall > 0.0) || !Double.isFinite(retainedOverall)) {
        throw new IllegalStateException("bifurcating phase exceeds an overall component inventory");
      }
      double logK = Math.max(-700.0, Math.min(100.0, variables[componentIndex]));
      double logDenominator = retainedPhaseLogDenominator(retainedBeta, logK);
      if (!Double.isFinite(logDenominator)) {
        throw new IllegalStateException("invalid retained-phase material-balance denominator");
      }
      logSecond[componentIndex] = Math.log(retainedOverall) - logDenominator;
      logFirst[componentIndex] = logK + logSecond[componentIndex];
      second[componentIndex] = Math.exp(logSecond[componentIndex]);
      first[componentIndex] = Math.exp(logFirst[componentIndex]);
    }
    double sumFirst = Math.exp(logSumExp(logFirst, 0, componentCount));
    double sumSecond = Math.exp(logSumExp(logSecond, 0, componentCount));
    if (!Double.isFinite(sumFirst) || !Double.isFinite(sumSecond) || !(sumFirst > 0.0) || !(sumSecond > 0.0)) {
      throw new IllegalStateException("retained-phase inventory left the physical domain");
    }
    updateWorkingSystem(working, temperatureK, pressureBara, retainedBeta, bifurcatingFraction, first, second, third);
    double[] residual = new double[2 * componentCount + 4];
    double equilibriumMaximumResidual = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      residual[componentIndex] = variables[componentIndex]
          + working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient()
          - working.getPhase(1).getComponent(componentIndex).getLogFugacityCoefficient();
      equilibriumMaximumResidual = Math.max(equilibriumMaximumResidual, Math.abs(residual[componentIndex]));
    }
    residual[componentCount] = sumFirst - sumSecond;
    equilibriumMaximumResidual = Math.max(equilibriumMaximumResidual, Math.abs(residual[componentCount]));
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double reference = Math
          .log(Math.max(working.getPhase(0).getComponent(componentIndex).getx(), MINIMUM_COMPOSITION))
          + working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
      residual[componentCount + 1 + componentIndex] = variables[logWStart + componentIndex]
          + working.getPhase(2).getComponent(componentIndex).getLogFugacityCoefficient() - reference;
      equilibriumMaximumResidual = Math.max(equilibriumMaximumResidual,
          Math.abs(residual[componentCount + 1 + componentIndex]));
    }
    residual[2 * componentCount + 1] = logSumW;
    equilibriumMaximumResidual = Math.max(equilibriumMaximumResidual, Math.abs(logSumW));
    double pressureResidual = Math.log(pressureBara / specifiedPressureBara);
    residual[2 * componentCount + 2] = pressureResidual;
    double[] normalizedFirst = normalized(first);
    double[] normalizedSecond = normalized(second);
    double[] duplicatedRetained = bifurcatingPhase == retainedPhaseZero ? normalizedFirst : normalizedSecond;
    double separation = modeSeparation(duplicatedRetained, third, mode);
    double separationResidual = separation - targetSeparation;
    residual[2 * componentCount + 3] = separationResidual;
    double maximumResidual = Math.max(equilibriumMaximumResidual,
        Math.max(Math.abs(pressureResidual), Math.abs(separationResidual)));
    if (!Double.isFinite(maximumResidual)) {
      throw new IllegalStateException("bordered bifurcation residual is non-finite");
    }
    return new Evaluation(residual, maximumResidual, equilibriumMaximumResidual, pressureResidual, separationResidual,
        temperatureK, pressureBara, retainedBeta, bifurcatingFraction, normalizedFirst, normalizedSecond, third,
        duplicatedRetained,
        materialBalanceResidual(retainedBeta, bifurcatingFraction, normalizedFirst, normalizedSecond, third));
  }

  /**
   * Builds a private EOS state from the requested phase compositions and fractions.
   *
   * @param temperatureK temperature in kelvin
   * @param pressureBara absolute pressure in bara
   * @param retainedBeta retained beta
   * @param bifurcatingFraction bifurcating fraction
   * @param first first
   * @param second second
   * @param third third
   * @return computed create working system result
   */
  private SystemInterface createWorkingSystem(double temperatureK, double pressureBara, double retainedBeta,
      double bifurcatingFraction, double[] first, double[] second, double[] third) {
    SystemInterface working = template.clone();
    working.setMultiPhaseCheck(false);
    working.setMaxNumberOfPhases(3);
    working.setNumberOfPhases(3);
    PhaseInterface phaseTemplate = template.getPhase(0);
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      working.getPhases()[phaseIndex] = phaseTemplate.clone();
    }
    updateWorkingSystem(working, temperatureK, pressureBara, retainedBeta, bifurcatingFraction, first, second, third);
    return working;
  }

  /**
   * Replaces all trial-dependent values in the private workspace for this correction.
   *
   * @param working solve-local EOS workspace
   * @param temperatureK temperature in kelvin
   * @param pressureBara absolute pressure in bara
   * @param retainedBeta first retained phase fraction within the retained inventory
   * @param bifurcatingFraction new phase mole fraction of the complete inventory
   * @param first first retained composition
   * @param second second retained composition
   * @param third bifurcating composition
   */
  private void updateWorkingSystem(SystemInterface working, double temperatureK, double pressureBara,
      double retainedBeta, double bifurcatingFraction, double[] first, double[] second, double[] third) {
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    working.setBeta(0, retainedBeta * (1.0 - bifurcatingFraction));
    working.setBeta(1, (1.0 - retainedBeta) * (1.0 - bifurcatingFraction));
    working.setBeta(2, bifurcatingFraction);
    setPhase(working, 0, retainedPhaseZero, first);
    setPhase(working, 1, retainedPhaseOne, second);
    setPhase(working, 2, bifurcatingPhase, third);
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      working.init(1, phaseIndex);
    }
  }

  /**
   * Solves a Newton step, using regularized least squares for ill-conditioned systems.
   *
   * @param jacobian jacobian
   * @param residual residual
   * @param conditionNumber condition number
   * @return computed solve correction result
   */
  private static Matrix solveCorrection(Matrix jacobian, double[] residual, double conditionNumber) {
    Matrix rightHandSide = new Matrix(residual, residual.length);
    if (Double.isFinite(conditionNumber) && conditionNumber < 1.0e8) {
      try {
        return jacobian.solve(rightHandSide);
      } catch (RuntimeException error) {
        // Fall through to regularized least squares.
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

  /**
   * Limits the Newton step without changing its direction.
   *
   * @param correction correction
   */
  private void limitCorrection(Matrix correction) {
    int componentCount = componentCount();
    double scale = 1.0;
    for (int variableIndex = 0; variableIndex < correction.getRowDimension(); variableIndex++) {
      double maximumStep;
      if (variableIndex == componentCount) {
        maximumStep = 0.05;
      } else if (variableIndex == 2 * componentCount + 1 || variableIndex == 2 * componentCount + 2) {
        maximumStep = 0.04;
      } else if (variableIndex == 2 * componentCount + 3) {
        maximumStep = 1.0;
      } else {
        maximumStep = 1.0;
      }
      double magnitude = Math.abs(correction.get(variableIndex, 0));
      if (magnitude > maximumStep) {
        scale = Math.min(scale, maximumStep / magnitude);
      }
    }
    if (scale < 1.0) {
      correction.timesEquals(scale);
    }
  }

  /**
   * Constructs logarithmic solver coordinates from the phase seeds.
   *
   * @param state state
   * @param bifurcatingFraction bifurcating fraction
   * @return computed variables result
   */
  private double[] variables(TwoToThreePhaseArcLengthCorrector.State state, double bifurcatingFraction) {
    int componentCount = componentCount();
    double[] variables = new double[2 * componentCount + 4];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      variables[componentIndex] = Math
          .log(Math.max(state.getPhaseZeroComposition()[componentIndex], MINIMUM_COMPOSITION)
              / Math.max(state.getPhaseOneComposition()[componentIndex], MINIMUM_COMPOSITION));
      variables[componentCount + 1 + componentIndex] = Math
          .log(Math.max(state.getIncipientComposition()[componentIndex], MINIMUM_COMPOSITION));
    }
    variables[componentCount] = state.getBeta();
    variables[2 * componentCount + 1] = Math.log(state.getTemperatureK());
    variables[2 * componentCount + 2] = Math.log(state.getPressureBara());
    variables[2 * componentCount + 3] = Math.log(bifurcatingFraction / (1.0 - bifurcatingFraction));
    return variables;
  }

  /**
   * Restricts solver coordinates to their numerical and physical bounds.
   *
   * @param variables variables
   */
  private void clampVariables(double[] variables) {
    int componentCount = componentCount();
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      variables[componentIndex] = clamp(variables[componentIndex], -700.0, 100.0);
      variables[componentCount + 1 + componentIndex] = clamp(variables[componentCount + 1 + componentIndex], -700.0,
          100.0);
    }
    variables[componentCount] = clamp(variables[componentCount], MINIMUM_BETA, 1.0 - MINIMUM_BETA);
    variables[2 * componentCount + 1] = clamp(variables[2 * componentCount + 1], Math.log(50.0), Math.log(2500.0));
    variables[2 * componentCount + 2] = clamp(variables[2 * componentCount + 2], Math.log(1.0e-6), Math.log(1.0e6));
    variables[2 * componentCount + 3] = clamp(variables[2 * componentCount + 3], -32.0, 20.0);
  }

  /**
   * Checks dimensions and physical domains before starting a numerical solve.
   *
   * @param state state
   * @param mode mode
   * @param pressureBara absolute pressure in bara
   * @param targetSeparation signed log-composition projection along the selected mode
   * @param bifurcatingFractionSeed initial mole fraction of the newly formed phase
   */
  private void validateInputs(TwoToThreePhaseArcLengthCorrector.State state,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode, double pressureBara, double targetSeparation,
      double bifurcatingFractionSeed) {
    int componentCount = componentCount();
    if (state == null || state.getRetainedPhaseZero() != retainedPhaseZero
        || state.getRetainedPhaseOne() != retainedPhaseOne || state.getIncipientPhase() != bifurcatingPhase
        || state.getPhaseZeroComposition().length != componentCount
        || state.getPhaseOneComposition().length != componentCount
        || state.getIncipientComposition().length != componentCount || mode == null || !mode.hasRealBifurcationMode()
        || mode.getReferenceComponentIndex() < 0 || mode.getReferenceComponentIndex() >= componentCount
        || mode.getBifurcationEigenvector().length != componentCount - 1 || !positive(pressureBara)
        || pressureBara < 1.0e-6 || pressureBara > 1.0e6 || !Double.isFinite(state.getTemperatureK())
        || state.getTemperatureK() < 50.0 || state.getTemperatureK() > 2500.0 || !positive(state.getBeta())
        || state.getBeta() >= 1.0 || !Double.isFinite(targetSeparation) || Math.abs(targetSeparation) < 1.0e-8
        || !positive(bifurcatingFractionSeed) || bifurcatingFractionSeed >= 1.0) {
      throw new IllegalArgumentException("invalid bordered retained-phase bifurcation input");
    }
    validateComposition(state.getPhaseZeroComposition());
    validateComposition(state.getPhaseOneComposition());
    validateComposition(state.getIncipientComposition());
  }

  /**
   * Requires finite, normalized phase seeds before constructing logarithmic variables.
   *
   * @param composition phase mole fractions
   * @throws IllegalArgumentException when any value is invalid or the sum is not unity
   */
  private static void validateComposition(double[] composition) {
    double total = 0.0;
    for (double value : composition) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException("bifurcation phase seeds must be finite and nonnegative");
      }
      total += value;
    }
    if (!Double.isFinite(total) || Math.abs(total - 1.0) > 1.0e-8) {
      throw new IllegalArgumentException("bifurcation phase seeds must sum to one");
    }
  }

  /**
   * Measures the largest deviation from the original overall component inventory.
   *
   * @param retainedBeta retained beta
   * @param bifurcatingFraction bifurcating fraction
   * @param first first
   * @param second second
   * @param third third
   * @return computed material balance residual result
   */
  private double materialBalanceResidual(double retainedBeta, double bifurcatingFraction, double[] first,
      double[] second, double[] third) {
    double[] overall = overallComposition();
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < overall.length; componentIndex++) {
      double reconstructed = (1.0 - bifurcatingFraction)
          * (retainedBeta * first[componentIndex] + (1.0 - retainedBeta) * second[componentIndex])
          + bifurcatingFraction * third[componentIndex];
      maximum = Math.max(maximum, Math.abs(reconstructed - overall[componentIndex]));
    }
    return maximum;
  }

  /**
   * Returns the normalized overall component inventory of the private template.
   *
   * @return computed overall composition result
   */
  private double[] overallComposition() {
    int componentCount = componentCount();
    double[] overall = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      overall[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(), MINIMUM_COMPOSITION);
    }
    return normalized(overall);
  }

  /**
   * Returns the number of components in the private template.
   *
   * @return computed component count result
   */
  private int componentCount() {
    return template.getPhase(0).getNumberOfComponents();
  }

  /**
   * Projects the signed log-composition separation onto the bifurcation eigenvector.
   *
   * @param retained retained
   * @param bifurcating bifurcating
   * @param mode mode
   * @return computed mode separation result
   */
  private static double modeSeparation(double[] retained, double[] bifurcating,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode) {
    int referenceIndex = mode.getReferenceComponentIndex();
    double retainedReference = Math.max(retained[referenceIndex], MINIMUM_COMPOSITION);
    double bifurcatingReference = Math.max(bifurcating[referenceIndex], MINIMUM_COMPOSITION);
    double[] eigenvector = mode.getBifurcationEigenvector();
    double separation = 0.0;
    int coordinateIndex = 0;
    for (int componentIndex = 0; componentIndex < retained.length; componentIndex++) {
      if (componentIndex == referenceIndex) {
        continue;
      }
      double difference = Math.log(Math.max(bifurcating[componentIndex], MINIMUM_COMPOSITION) / bifurcatingReference)
          - Math.log(Math.max(retained[componentIndex], MINIMUM_COMPOSITION) / retainedReference);
      separation += eigenvector[coordinateIndex++] * difference;
    }
    return separation;
  }

  /**
   * Initializes a private EOS phase with normalized composition and the requested type.
   *
   * @param system system
   * @param phaseIndex zero-based phase index
   * @param phase phase
   * @param composition composition
   */
  private static void setPhase(SystemInterface system, int phaseIndex, CandidatePhase phase, double[] composition) {
    system.setPhaseType(phaseIndex, toPhaseType(phase));
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      // The EOS requires a strictly positive mole fraction; this is the only place the floor belongs.
      system.getPhase(phaseIndex).getComponent(componentIndex)
          .setx(Math.max(composition[componentIndex], MINIMUM_COMPOSITION));
    }
    system.getPhase(phaseIndex).normalize();
    system.setPhaseType(phaseIndex, toPhaseType(phase));
  }

  /**
   * Evaluates {@code log(1 - beta + beta * exp(logK))} without forming {@code exp(logK)} when it would overflow.
   *
   * @param beta beta
   * @param logK log k
   * @return computed retained phase log denominator result
   */
  private static double retainedPhaseLogDenominator(double beta, double logK) {
    double retainedFraction = 1.0 - beta;
    if (logK > 0.0) {
      return logK + Math.log(beta + retainedFraction * Math.exp(-logK));
    }
    return Math.log(retainedFraction + beta * Math.exp(logK));
  }

  /**
   * Computes a stable log-sum-exp over the requested interval.
   *
   * @param values values
   * @param start start
   * @param end end
   * @return computed log sum exp result
   */
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

  /**
   * Recovers mole fractions from logarithms and their log-sum.
   *
   * @param values values
   * @param start start
   * @param end end
   * @param logSum log sum
   * @return computed normalized exponentials result
   */
  private static double[] normalizedExponentials(double[] values, int start, int end, double logSum) {
    double[] composition = new double[end - start];
    for (int index = start; index < end; index++) {
      composition[index - start] = Math.exp(values[index] - logSum);
    }
    return composition;
  }

  /**
   * Returns a normalized defensive copy of a composition.
   *
   * @param values values
   * @return computed normalized result
   */
  private static double[] normalized(double[] values) {
    double[] result = values.clone();
    double total = 0.0;
    for (double value : result) {
      total += value;
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalStateException("composition cannot be normalized");
    }
    for (int index = 0; index < result.length; index++) {
      result[index] /= total;
    }
    return result;
  }

  /**
   * Checks or computes logistic.
   *
   * @param value value
   * @return computed logistic result
   */
  private static double logistic(double value) {
    if (value >= 0.0) {
      double exponential = Math.exp(-value);
      return 1.0 / (1.0 + exponential);
    }
    double exponential = Math.exp(value);
    return exponential / (1.0 + exponential);
  }

  /**
   * Computes the L1 distance between phase compositions.
   *
   * @param first first
   * @param second second
   * @return computed composition distance result
   */
  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

  /**
   * Restricts a scalar to an inclusive interval.
   *
   * @param value value
   * @param lower lower
   * @param upper upper
   * @return computed clamp result
   */
  private static double clamp(double value, double lower, double upper) {
    return Math.max(lower, Math.min(upper, value));
  }

  /**
   * Checks that a scalar is finite and strictly positive.
   *
   * @param value value
   * @return true when the documented condition holds
   */
  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /**
   * Maps the requested physical family to its EOS phase type.
   *
   * @param phase phase
   * @return computed to phase type result
   */
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
    private final double equilibriumMaximumResidual;
    private final double pressureResidual;
    private final double separationResidual;
    private final double temperatureK;
    private final double pressureBara;
    private final double retainedBeta;
    private final double bifurcatingFraction;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double[] bifurcatingComposition;
    private final double[] duplicatedRetainedComposition;
    private final double materialBalanceResidual;

    /**
     * Creates an immutable evaluation.
     *
     * @param residual residual
     * @param maximumResidual maximum residual
     * @param equilibriumMaximumResidual equilibrium maximum residual
     * @param pressureResidual pressure residual
     * @param separationResidual separation residual
     * @param temperatureK temperature in kelvin
     * @param pressureBara absolute pressure in bara
     * @param retainedBeta retained beta
     * @param bifurcatingFraction bifurcating fraction
     * @param phaseZeroComposition phase zero composition
     * @param phaseOneComposition phase one composition
     * @param bifurcatingComposition bifurcating composition
     * @param duplicatedRetainedComposition duplicated retained composition
     * @param materialBalanceResidual material balance residual
     */
    private Evaluation(double[] residual, double maximumResidual, double equilibriumMaximumResidual,
        double pressureResidual, double separationResidual, double temperatureK, double pressureBara,
        double retainedBeta, double bifurcatingFraction, double[] phaseZeroComposition, double[] phaseOneComposition,
        double[] bifurcatingComposition, double[] duplicatedRetainedComposition, double materialBalanceResidual) {
      this.residual = residual;
      this.maximumResidual = maximumResidual;
      this.equilibriumMaximumResidual = equilibriumMaximumResidual;
      this.pressureResidual = pressureResidual;
      this.separationResidual = separationResidual;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.retainedBeta = retainedBeta;
      this.bifurcatingFraction = bifurcatingFraction;
      this.phaseZeroComposition = phaseZeroComposition;
      this.phaseOneComposition = phaseOneComposition;
      this.bifurcatingComposition = bifurcatingComposition;
      this.duplicatedRetainedComposition = duplicatedRetainedComposition;
      this.materialBalanceResidual = materialBalanceResidual;
    }
  }

  /**
   * Immutable bordered bifurcation correction and physical diagnostics.
   *
   * @author NeqSim contributors
   * @version 1.0
   */
  public static final class Result {
    private final Object owner;
    private final CandidatePhase retainedPhaseZero;
    private final CandidatePhase retainedPhaseOne;
    private final CandidatePhase bifurcatingPhase;
    private final double temperatureK;
    private final double pressureBara;
    private final double retainedBeta;
    private final double bifurcatingFraction;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double[] bifurcatingComposition;
    private final double initialMaximumResidual;
    private final double maximumResidual;
    private final double equilibriumMaximumResidual;
    private final double pressureResidual;
    private final double separationResidual;
    private final double materialBalanceResidual;
    private final int iterations;
    private final double jacobianConditionNumber;
    private final boolean residualConverged;
    private final boolean distinct;
    private final boolean positiveFraction;
    private final String failureMessage;

    /**
     * Creates an immutable result.
     *
     * @param owner private identity of the solver that produced this result
     * @param retainedPhaseZero retained phase zero
     * @param retainedPhaseOne retained phase one
     * @param bifurcatingPhase bifurcating phase
     * @param temperatureK temperature in kelvin
     * @param pressureBara absolute pressure in bara
     * @param retainedBeta retained beta
     * @param bifurcatingFraction bifurcating fraction
     * @param phaseZeroComposition phase zero composition
     * @param phaseOneComposition phase one composition
     * @param bifurcatingComposition bifurcating composition
     * @param initialMaximumResidual initial maximum residual
     * @param maximumResidual maximum residual
     * @param equilibriumMaximumResidual equilibrium maximum residual
     * @param pressureResidual pressure residual
     * @param separationResidual separation residual
     * @param materialBalanceResidual material balance residual
     * @param iterations iterations
     * @param jacobianConditionNumber jacobian condition number
     * @param residualConverged residual converged
     * @param distinct distinct
     * @param positiveFraction positive fraction
     * @param failureMessage failure diagnostic, or null for success
     */
    private Result(Object owner, CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase bifurcatingPhase, double temperatureK, double pressureBara, double retainedBeta,
        double bifurcatingFraction, double[] phaseZeroComposition, double[] phaseOneComposition,
        double[] bifurcatingComposition, double initialMaximumResidual, double maximumResidual,
        double equilibriumMaximumResidual, double pressureResidual, double separationResidual,
        double materialBalanceResidual, int iterations, double jacobianConditionNumber, boolean residualConverged,
        boolean distinct, boolean positiveFraction, String failureMessage) {
      this.owner = owner;
      this.retainedPhaseZero = retainedPhaseZero;
      this.retainedPhaseOne = retainedPhaseOne;
      this.bifurcatingPhase = bifurcatingPhase;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.retainedBeta = retainedBeta;
      this.bifurcatingFraction = bifurcatingFraction;
      this.phaseZeroComposition = phaseZeroComposition.clone();
      this.phaseOneComposition = phaseOneComposition.clone();
      this.bifurcatingComposition = bifurcatingComposition.clone();
      this.initialMaximumResidual = initialMaximumResidual;
      this.maximumResidual = maximumResidual;
      this.equilibriumMaximumResidual = equilibriumMaximumResidual;
      this.pressureResidual = pressureResidual;
      this.separationResidual = separationResidual;
      this.materialBalanceResidual = materialBalanceResidual;
      this.iterations = iterations;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.residualConverged = residualConverged;
      this.distinct = distinct;
      this.positiveFraction = positiveFraction;
      this.failureMessage = failureMessage;
    }

    /**
     * Constructs an immutable unsuccessful result with no accepted equilibrium.
     *
     * @param owner private identity of the solver that produced this result
     * @param retainedPhaseZero retained phase zero
     * @param retainedPhaseOne retained phase one
     * @param bifurcatingPhase bifurcating phase
     * @param temperatureK temperature in kelvin
     * @param pressureBara absolute pressure in bara
     * @param componentCount number of components
     * @param failureMessage failure diagnostic, or null for success
     * @return computed failure result
     */
    private static Result failure(Object owner, CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase bifurcatingPhase, double temperatureK, double pressureBara, int componentCount,
        String failureMessage) {
      return new Result(owner, retainedPhaseZero, retainedPhaseOne, bifurcatingPhase, temperatureK, pressureBara,
          Double.NaN, Double.NaN, new double[componentCount], new double[componentCount], new double[componentCount],
          Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0, Double.NaN, false, false, false,
          failureMessage);
    }

    /**
     * Reports numerical convergence only; use validateEquilibrium before accepting a physical state.
     *
     * @return true when the documented condition holds
     */
    public boolean isConverged() {
      return residualConverged && distinct && positiveFraction && failureMessage == null;
    }

    /**
     * Returns the retained phase zero.
     *
     * @return retained phase zero
     */
    public CandidatePhase getRetainedPhaseZero() {
      return retainedPhaseZero;
    }

    /**
     * Returns the retained phase one.
     *
     * @return retained phase one
     */
    public CandidatePhase getRetainedPhaseOne() {
      return retainedPhaseOne;
    }

    /**
     * Returns the bifurcating phase.
     *
     * @return bifurcating phase
     */
    public CandidatePhase getBifurcatingPhase() {
      return bifurcatingPhase;
    }

    /**
     * Returns the temperature k.
     *
     * @return temperature k
     */
    public double getTemperatureK() {
      return temperatureK;
    }

    /**
     * Returns the pressure bara.
     *
     * @return pressure bara
     */
    public double getPressureBara() {
      return pressureBara;
    }

    /**
     * Returns the retained beta.
     *
     * @return retained beta
     */
    public double getRetainedBeta() {
      return retainedBeta;
    }

    /**
     * Returns the bifurcating fraction.
     *
     * @return bifurcating fraction
     */
    public double getBifurcatingFraction() {
      return bifurcatingFraction;
    }

    /**
     * Returns the phase zero composition.
     *
     * @return phase zero composition
     */
    public double[] getPhaseZeroComposition() {
      return phaseZeroComposition.clone();
    }

    /**
     * Returns the phase one composition.
     *
     * @return phase one composition
     */
    public double[] getPhaseOneComposition() {
      return phaseOneComposition.clone();
    }

    /**
     * Returns the bifurcating composition.
     *
     * @return bifurcating composition
     */
    public double[] getBifurcatingComposition() {
      return bifurcatingComposition.clone();
    }

    /**
     * Returns the initial maximum residual.
     *
     * @return initial maximum residual
     */
    public double getInitialMaximumResidual() {
      return initialMaximumResidual;
    }

    /**
     * Returns the maximum residual.
     *
     * @return maximum residual
     */
    public double getMaximumResidual() {
      return maximumResidual;
    }

    /**
     * Returns the equilibrium maximum residual.
     *
     * @return equilibrium maximum residual
     */
    public double getEquilibriumMaximumResidual() {
      return equilibriumMaximumResidual;
    }

    /**
     * Returns the pressure residual.
     *
     * @return pressure residual
     */
    public double getPressureResidual() {
      return pressureResidual;
    }

    /**
     * Returns the separation residual.
     *
     * @return separation residual
     */
    public double getSeparationResidual() {
      return separationResidual;
    }

    /**
     * Returns the material balance residual.
     *
     * @return material balance residual
     */
    public double getMaterialBalanceResidual() {
      return materialBalanceResidual;
    }

    /**
     * Returns the iterations.
     *
     * @return iterations
     */
    public int getIterations() {
      return iterations;
    }

    /**
     * Returns the jacobian condition number.
     *
     * @return jacobian condition number
     */
    public double getJacobianConditionNumber() {
      return jacobianConditionNumber;
    }

    /**
     * Returns the bifurcating phase distance.
     *
     * @return bifurcating phase distance
     */
    public double getBifurcatingPhaseDistance() {
      double[] duplicated = bifurcatingPhase == retainedPhaseZero ? phaseZeroComposition : phaseOneComposition;
      return compositionDistance(duplicated, bifurcatingComposition);
    }

    /**
     * Returns the failure message.
     *
     * @return failure message
     */
    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
