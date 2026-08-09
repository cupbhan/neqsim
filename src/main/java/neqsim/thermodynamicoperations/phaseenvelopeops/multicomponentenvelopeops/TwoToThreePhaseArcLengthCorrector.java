package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import Jama.Matrix;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Pseudo-arclength predictor-corrector for a specified two-to-three-phase PT boundary.
 *
 * <p>
 * The unknown vector is {@code [ln(K), beta, ln(W_3), ln(T), ln(P)]}. The thermodynamic equations are the retained
 * two-phase fugacity equalities and material-balance closure, followed by the third-phase fugacity equalities and
 * {@code ln(sum(W_3)) = 0}. These leave a one-dimensional PT manifold. A secant-tangent hyperplane equation closes the
 * Newton system without fixing either temperature or pressure, allowing continuation through pressure or temperature
 * turning points.
 * </p>
 */
public final class TwoToThreePhaseArcLengthCorrector {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;
  private static final double MINIMUM_BETA = 1.0e-10;
  private static final double MINIMUM_LOG_VALUE = -700.0;
  private static final double MAXIMUM_LOG_K_VALUE = 700.0;
  private static final double MAXIMUM_LOG_W_VALUE = 100.0;

  private final SystemInterface template;
  private final SystemInterface working;
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase incipientPhase;
  private final boolean retainedPhaseBifurcation;
  private int maximumIterations = 50;
  private double residualTolerance = 1.0e-8;
  private double finiteDifferenceStep = 2.0e-5;
  private double deflationPower = 0.0;
  private double deflationShift = 1.0;

  /** Creates a corrector for three distinct physical phase families. */
  public TwoToThreePhaseArcLengthCorrector(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase) {
    this(template, retainedPhaseZero, retainedPhaseOne, incipientPhase, false);
  }

  /**
   * Creates a finite-fraction homotopy corrector for a new phase that bifurcates from a retained phase family.
   *
   * <p>
   * The returned corrector uses the same equations as the ordinary two-to-three-phase boundary, but permits the third
   * phase slot to share the GAS or OIL EOS family of one retained phase. The two slots still carry independent
   * compositions and fractions, and duplicate compositions remain rejected by the result quality check.
   * </p>
   */
  public static TwoToThreePhaseArcLengthCorrector forRetainedPhaseBifurcation(SystemInterface template,
      CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne, CandidatePhase bifurcatingPhase) {
    return new TwoToThreePhaseArcLengthCorrector(template, retainedPhaseZero, retainedPhaseOne, bifurcatingPhase, true);
  }

  private TwoToThreePhaseArcLengthCorrector(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase, boolean retainedPhaseBifurcation) {
    if (template == null || retainedPhaseZero == null || retainedPhaseOne == null || incipientPhase == null) {
      throw new IllegalArgumentException("template and all three phase families must be specified");
    }
    if (retainedPhaseZero == retainedPhaseOne
        || !retainedPhaseBifurcation && (retainedPhaseZero == incipientPhase || retainedPhaseOne == incipientPhase)) {
      throw new IllegalArgumentException("the three physical phase families must be distinct");
    }
    if (retainedPhaseBifurcation && incipientPhase != retainedPhaseZero && incipientPhase != retainedPhaseOne) {
      throw new IllegalArgumentException("the bifurcating phase must match one retained phase family");
    }
    this.template = template.clone();
    this.working = this.template.clone();
    this.working.setMultiPhaseCheck(false);
    this.working.setMaxNumberOfPhases(3);
    this.working.setNumberOfPhases(3);
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.incipientPhase = incipientPhase;
    this.retainedPhaseBifurcation = retainedPhaseBifurcation;
  }

  /** Sets damped-Newton and numerical-Jacobian controls. */
  public TwoToThreePhaseArcLengthCorrector setNumericalControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    if (maximumIterations < 1 || !Double.isFinite(residualTolerance) || residualTolerance <= 0.0
        || !Double.isFinite(finiteDifferenceStep) || finiteDifferenceStep <= 0.0) {
      throw new IllegalArgumentException("invalid pseudo-arclength corrector controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /** Sets coincident-phase deflation used only by the retained-phase bifurcation factory. */
  public TwoToThreePhaseArcLengthCorrector setCoincidentPhaseDeflation(double power, double shift) {
    if (!retainedPhaseBifurcation) {
      throw new IllegalStateException("coincident-phase deflation is only available for retained-phase bifurcation");
    }
    if (!Double.isFinite(power) || power <= 0.0 || !Double.isFinite(shift) || shift < 0.0) {
      throw new IllegalArgumentException("invalid coincident-phase deflation controls");
    }
    this.deflationPower = power;
    this.deflationShift = shift;
    return this;
  }

  /**
   * Predicts from two accepted points and corrects on the pseudo-arclength hyperplane.
   *
   * @param previous previous accepted point
   * @param current current accepted point
   * @param stepFactor multiplier on the transformed secant distance between the two points
   * @return corrected point and diagnostics
   */
  public Result correct(State previous, State current, double stepFactor) {
    return correctInternal(previous, current, stepFactor, true);
  }

  /**
   * Corrects with an absolute transformed-space pseudo-arclength step.
   *
   * @param previous previous accepted point
   * @param current current accepted point
   * @param arcStep positive absolute step in the transformed unknown vector
   * @return corrected point and diagnostics
   */
  public Result correctWithArcStep(State previous, State current, double arcStep) {
    return correctInternal(previous, current, arcStep, false);
  }

  /**
   * Corrects the full two-to-three-phase equations while fixing pressure instead of pseudo-arclength.
   *
   * <p>
   * The existing final hyperplane equation is reused with a unit vector in {@code ln(P)}, making its residual exactly
   * {@code ln(P/P_target)}. Temperature and every equilibrium/composition variable remain simultaneous Newton unknowns.
   * This avoids nested TPD root finding when a stationary-point branch folds or switches.
   * </p>
   */
  public Result correctAtPressure(State initial, double pressureBara) {
    return correctAtPressureWithIncipientFraction(initial, pressureBara, 0.0);
  }

  /**
   * Corrects the full two-to-three-phase equations while fixing temperature instead of pseudo-arclength.
   *
   * <p>
   * This is the complementary local parameterization to {@link #correctAtPressure(State, double)}. It is required to
   * construct a second strict continuation seed when a boundary is locally vertical in the PT plane or folds in
   * pressure.
   * </p>
   *
   * @param initial full corrected boundary seed
   * @param temperatureK fixed absolute temperature
   * @return corrected zero-third-phase boundary state
   */
  public Result correctAtTemperature(State initial, double temperatureK) {
    validateState(initial);
    if (!Double.isFinite(temperatureK) || temperatureK < 50.0 || temperatureK > 2500.0) {
      throw new IllegalArgumentException("fixed temperature must be finite and physical");
    }
    double[] predictor = transformedVariables(initial);
    int temperatureIndex = predictor.length - 2;
    predictor[temperatureIndex] = Math.log(temperatureK);
    double[] temperatureSpecification = new double[predictor.length];
    temperatureSpecification[temperatureIndex] = 1.0;
    clampPhysicalVariables(predictor);
    return correctFromPredictor(predictor, temperatureSpecification, 0.0, 0.0);
  }

  /**
   * Builds a one-point pseudo-arclength predictor from the null space of the boundary Jacobian and corrects it.
   *
   * <p>
   * The thermodynamic boundary has one fewer independent equation than transformed unknowns. Supplying a zero final
   * hyperplane row exposes that one-dimensional null space; the right singular vector associated with the smallest
   * singular value is the local continuation tangent. This initializer remains valid at pressure and temperature folds,
   * where fixed-P and fixed-T parameterizations can both fail.
   * </p>
   *
   * @param initial one converged zero-third-phase boundary state
   * @param arcStep positive transformed-space predictor distance
   * @param orientation either {@code +1} or {@code -1} for the two tangent directions
   * @return corrected neighboring boundary point
   */
  public Result correctFromLocalTangent(State initial, double arcStep, int orientation) {
    validateState(initial);
    if (!Double.isFinite(arcStep) || arcStep <= 0.0 || orientation != -1 && orientation != 1) {
      throw new IllegalArgumentException("a positive local arc step and orientation of +1 or -1 are required");
    }
    double[] baseVariables = transformedVariables(initial);
    double[] zeroHyperplane = new double[baseVariables.length];
    Evaluation base;
    Matrix jacobian;
    try {
      base = evaluate(baseVariables, baseVariables, zeroHyperplane, 0.0);
      jacobian = numericalJacobian(baseVariables, baseVariables, zeroHyperplane, base.residual, 0.0);
    } catch (RuntimeException error) {
      return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase, arcStep,
          "local tangent Jacobian evaluation failed: " + error.getMessage());
    }
    Matrix rightSingularVectors;
    try {
      rightSingularVectors = jacobian.svd().getV();
    } catch (RuntimeException error) {
      return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase, arcStep,
          "local tangent singular-value decomposition failed: " + error.getMessage());
    }
    int tangentColumn = rightSingularVectors.getColumnDimension() - 1;
    double[] tangent = new double[baseVariables.length];
    double norm = 0.0;
    for (int index = 0; index < tangent.length; index++) {
      tangent[index] = orientation * rightSingularVectors.get(index, tangentColumn);
      norm += tangent[index] * tangent[index];
    }
    norm = Math.sqrt(norm);
    if (!(norm > 1.0e-12) || !Double.isFinite(norm)) {
      return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase, arcStep,
          "local boundary Jacobian did not expose a finite tangent");
    }
    double[] predictor = baseVariables.clone();
    for (int index = 0; index < tangent.length; index++) {
      tangent[index] /= norm;
      predictor[index] += arcStep * tangent[index];
    }
    clampPhysicalVariables(predictor);
    return correctFromPredictor(predictor, tangent, arcStep, 0.0);
  }

  /**
   * Corrects at fixed pressure with a prescribed finite incipient-phase fraction for phase-fraction homotopy.
   *
   * @param initial full three-phase seed
   * @param pressureBara fixed pressure
   * @param incipientFraction prescribed total fraction of phase three in {@code [0, 1)}
   * @return corrected finite-fraction state; zero fraction is the physical 2P-to-3P boundary
   */
  public Result correctAtPressureWithIncipientFraction(State initial, double pressureBara, double incipientFraction) {
    validateState(initial);
    if (!Double.isFinite(pressureBara) || pressureBara <= 1.0e-6 || pressureBara > 1.0e6) {
      throw new IllegalArgumentException("fixed pressure must be finite and positive");
    }
    if (!Double.isFinite(incipientFraction) || incipientFraction < 0.0 || incipientFraction >= 1.0) {
      throw new IllegalArgumentException("incipient phase fraction must be in [0, 1)");
    }
    double[] predictor = transformedVariables(initial);
    int pressureIndex = predictor.length - 1;
    predictor[pressureIndex] = Math.log(pressureBara);
    double[] pressureSpecification = new double[predictor.length];
    pressureSpecification[pressureIndex] = 1.0;
    clampPhysicalVariables(predictor);
    return correctFromPredictor(predictor, pressureSpecification, 0.0, incipientFraction);
  }

  private Result correctInternal(State previous, State current, double requestedStep, boolean relativeToSecant) {
    validateState(previous);
    validateState(current);
    if (!Double.isFinite(requestedStep) || requestedStep <= 0.0) {
      throw new IllegalArgumentException("pseudo-arclength step must be positive");
    }
    double[] previousVariables = transformedVariables(previous);
    double[] currentVariables = transformedVariables(current);
    double[] direction = new double[currentVariables.length];
    double[] metricWeights = continuationMetricWeights(current);
    double secantNorm = 0.0;
    for (int variableIndex = 0; variableIndex < direction.length; variableIndex++) {
      direction[variableIndex] = currentVariables[variableIndex] - previousVariables[variableIndex];
      double scaledDifference = metricWeights[variableIndex] * direction[variableIndex];
      secantNorm += scaledDifference * scaledDifference;
    }
    secantNorm = Math.sqrt(secantNorm);
    if (!(secantNorm > 1.0e-12) || !Double.isFinite(secantNorm)) {
      throw new IllegalArgumentException("accepted points do not define a finite continuation tangent");
    }
    double[] hyperplaneNormal = new double[direction.length];
    double normalNorm = 0.0;
    for (int variableIndex = 0; variableIndex < direction.length; variableIndex++) {
      direction[variableIndex] /= secantNorm;
      hyperplaneNormal[variableIndex] = metricWeights[variableIndex] * metricWeights[variableIndex]
          * direction[variableIndex];
      normalNorm += hyperplaneNormal[variableIndex] * hyperplaneNormal[variableIndex];
    }
    normalNorm = Math.sqrt(normalNorm);
    if (!(normalNorm > 1.0e-12) || !Double.isFinite(normalNorm)) {
      throw new IllegalArgumentException("accepted points do not define a finite scaled continuation hyperplane");
    }
    for (int variableIndex = 0; variableIndex < hyperplaneNormal.length; variableIndex++) {
      hyperplaneNormal[variableIndex] /= normalNorm;
    }
    double arcStep = relativeToSecant ? requestedStep * secantNorm : requestedStep;
    double[] predictor = currentVariables.clone();
    for (int variableIndex = 0; variableIndex < predictor.length; variableIndex++) {
      predictor[variableIndex] += arcStep * direction[variableIndex];
    }
    clampPhysicalVariables(predictor);
    return correctFromPredictor(predictor, hyperplaneNormal, arcStep);
  }

  private Result correctFromPredictor(double[] predictor, double[] tangent, double arcStep) {
    return correctFromPredictor(predictor, tangent, arcStep, 0.0);
  }

  private Result correctFromPredictor(double[] predictor, double[] tangent, double arcStep, double incipientFraction) {
    double[] variables = predictor.clone();
    Evaluation evaluation;
    try {
      evaluation = evaluate(variables, predictor, tangent, incipientFraction);
    } catch (RuntimeException error) {
      return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase, arcStep, error.getMessage());
    }
    double initialResidual = evaluation.maximumResidual;
    int iterations = 0;
    double conditionNumber = Double.NaN;
    String failureMessage = null;
    while (iterations < maximumIterations && evaluation.maximumResidual > residualTolerance) {
      iterations++;
      Matrix jacobian;
      try {
        jacobian = numericalJacobian(variables, predictor, tangent, evaluation.residual, incipientFraction);
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
      if (!Double.isFinite(correction.normInf())) {
        failureMessage = "pseudo-arclength Newton correction is non-finite";
        break;
      }
      limitCorrection(correction);
      boolean accepted = false;
      double damping = 1.0;
      for (int lineSearch = 0; lineSearch < 24; lineSearch++) {
        double[] trialVariables = variables.clone();
        for (int variableIndex = 0; variableIndex < trialVariables.length; variableIndex++) {
          trialVariables[variableIndex] -= damping * correction.get(variableIndex, 0);
        }
        clampPhysicalVariables(trialVariables);
        try {
          Evaluation trial = evaluate(trialVariables, predictor, tangent, incipientFraction);
          if (trial.maximumResidual < evaluation.maximumResidual) {
            variables = trialVariables;
            evaluation = trial;
            accepted = true;
            break;
          }
        } catch (RuntimeException error) {
          // Retry with a shorter line-search step.
        }
        damping *= 0.5;
      }
      if (!accepted) {
        failureMessage = "line search could not reduce the pseudo-arclength residual";
        break;
      }
    }
    boolean residualConverged = evaluation.maximumResidual <= residualTolerance;
    if (failureMessage == null && !residualConverged) {
      failureMessage = "maximum pseudo-arclength iteration count reached";
    }
    boolean distinct = distinctPhases(evaluation);
    if (failureMessage == null && !distinct) {
      failureMessage = "corrected state contains a trivial or duplicate physical phase";
    }
    State state = new State(retainedPhaseZero, retainedPhaseOne, incipientPhase, evaluation.temperatureK,
        evaluation.pressureBara, evaluation.beta, evaluation.phaseZeroComposition, evaluation.phaseOneComposition,
        evaluation.incipientComposition, variables);
    EquationResiduals independentReplay = null;
    EquationResiduals returnedStateReplay = null;
    StateFingerprint returnedStateReplayFingerprint = null;
    StateFingerprint independentReplayFingerprint = null;
    if (failureMessage == null && residualConverged && distinct) {
      try {
        // Replaying on this same instance isolates the loss caused purely by rebuilding the variables from the
        // returned state; the independent instance additionally exposes template-clone and phase-slot leakage.
        Evaluation returnedStateReplayEvaluation = replayEvaluation(state);
        returnedStateReplay = returnedStateReplayEvaluation.equationResiduals;
        returnedStateReplayFingerprint = returnedStateReplayEvaluation.fingerprint;
        TwoToThreePhaseArcLengthCorrector independent = new TwoToThreePhaseArcLengthCorrector(template,
            retainedPhaseZero, retainedPhaseOne, incipientPhase, retainedPhaseBifurcation);
        Evaluation independentReplayEvaluation = independent.replayEvaluation(state);
        independentReplay = independentReplayEvaluation.equationResiduals;
        independentReplayFingerprint = independentReplayEvaluation.fingerprint;
        if (independentReplay.getThermodynamicMaximumResidual() > residualTolerance) {
          failureMessage = "corrected state failed independent thermodynamic replay: residual="
              + independentReplay.getThermodynamicMaximumResidual() + " tolerance=" + residualTolerance;
        }
      } catch (RuntimeException error) {
        failureMessage = "corrected state independent thermodynamic replay failed: " + error.getMessage();
      }
    }
    return new Result(state, initialResidual, evaluation.maximumResidual, evaluation.thermodynamicMaximumResidual,
        evaluation.arcLengthResidual, evaluation.equationResiduals, independentReplay, returnedStateReplay, arcStep,
        iterations, conditionNumber, residualConverged, distinct, failureMessage, evaluation.fingerprint,
        returnedStateReplayFingerprint, independentReplayFingerprint);
  }

  private static Matrix solveCorrection(Matrix jacobian, double[] residual, double conditionNumber) {
    Matrix rightHandSide = new Matrix(residual, residual.length);
    if (Double.isFinite(conditionNumber) && conditionNumber < 1.0e8) {
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

  private Matrix numericalJacobian(double[] variables, double[] predictor, double[] tangent, double[] baseResidual,
      double incipientFraction) {
    int variableCount = variables.length;
    Matrix jacobian = new Matrix(variableCount, variableCount);
    for (int column = 0; column < variableCount; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(variables[column]));
      double[] plus = variables.clone();
      double[] minus = variables.clone();
      plus[column] += step;
      minus[column] -= step;
      clampPhysicalVariables(plus);
      clampPhysicalVariables(minus);
      double denominator = plus[column] - minus[column];
      if (!(Math.abs(denominator) > 1.0e-20)) {
        throw new IllegalStateException("finite-difference variable is pinned at a physical bound");
      }
      double[] plusResidual = evaluate(plus, predictor, tangent, incipientFraction).residual;
      double[] minusResidual;
      try {
        minusResidual = evaluate(minus, predictor, tangent, incipientFraction).residual;
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

  private Evaluation evaluate(double[] variables, double[] predictor, double[] tangent, double incipientFraction) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    int betaIndex = componentCount;
    int logWStart = componentCount + 1;
    int temperatureIndex = 2 * componentCount + 1;
    int pressureIndex = 2 * componentCount + 2;
    double beta = variables[betaIndex];
    double temperatureK = Math.exp(variables[temperatureIndex]);
    double pressureBara = Math.exp(variables[pressureIndex]);
    if (!(beta > 0.0 && beta < 1.0) || !Double.isFinite(temperatureK) || temperatureK < 50.0 || temperatureK > 2500.0
        || !Double.isFinite(pressureBara) || pressureBara <= 1.0e-6 || pressureBara > 1.0e6) {
      throw new IllegalStateException("pseudo-arclength variables left the physical domain");
    }

    double[] z = overallComposition();
    double logSumW = logSumExp(variables, logWStart, logWStart + componentCount);
    double[] third = normalizedExponentials(variables, logWStart, logWStart + componentCount, logSumW);
    // The retained-phase inventory is evaluated in log space. A heavy pseudo-component in the aqueous phase
    // legitimately sits far below any fixed composition floor (C50-C80 at K ~ 1e99 lands near 1e-102), and flooring
    // it here would propagate into first[] = K * second[] as a large spurious inventory, break the Rachford-Rice
    // identity sum(first) = sum(second) = 1 and make the converged point a root of a perturbed system. The floor is
    // applied only where the EOS needs a strictly positive mole fraction, in setPhase.
    double[] first = new double[componentCount];
    double[] second = new double[componentCount];
    double[] logFirst = new double[componentCount];
    double[] logSecond = new double[componentCount];
    double[] equilibriumRatios = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double logK = clamp(variables[componentIndex], MINIMUM_LOG_VALUE, MAXIMUM_LOG_K_VALUE);
      equilibriumRatios[componentIndex] = Math.exp(logK);
      double logDenominator = retainedPhaseLogDenominator(beta, logK);
      if (!Double.isFinite(logDenominator)) {
        throw new IllegalStateException("invalid retained-phase material-balance denominator");
      }
      double retainedOverall = (z[componentIndex] - incipientFraction * third[componentIndex])
          / (1.0 - incipientFraction);
      if (!(retainedOverall > 0.0) || !Double.isFinite(retainedOverall)) {
        throw new IllegalStateException("prescribed incipient fraction exceeds the overall component inventory");
      }
      logSecond[componentIndex] = Math.log(retainedOverall) - logDenominator;
      logFirst[componentIndex] = logK + logSecond[componentIndex];
      second[componentIndex] = Math.exp(logSecond[componentIndex]);
      first[componentIndex] = Math.exp(logFirst[componentIndex]);
    }
    double logSumFirst = logSumExp(logFirst, 0, componentCount);
    double logSumSecond = logSumExp(logSecond, 0, componentCount);
    double sumFirst = Math.exp(logSumFirst);
    double sumSecond = Math.exp(logSumSecond);
    if (!Double.isFinite(sumFirst) || !Double.isFinite(sumSecond) || !(sumFirst > 0.0) || !(sumSecond > 0.0)) {
      throw new IllegalStateException("retained-phase inventory left the physical domain");
    }
    SystemInterface working = createWorkingSystem(temperatureK, pressureBara, beta, first, second, third,
        incipientFraction);

    double[] residual = new double[2 * componentCount + 3];
    double thermodynamicMaximumResidual = 0.0;
    double retainedFugacityMaximumResidual = 0.0;
    int retainedFugacityWorstComponentIndex = -1;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      residual[componentIndex] = variables[componentIndex]
          + working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient()
          - working.getPhase(1).getComponent(componentIndex).getLogFugacityCoefficient();
      thermodynamicMaximumResidual = Math.max(thermodynamicMaximumResidual, Math.abs(residual[componentIndex]));
      if (Math.abs(residual[componentIndex]) > retainedFugacityMaximumResidual) {
        retainedFugacityMaximumResidual = Math.abs(residual[componentIndex]);
        retainedFugacityWorstComponentIndex = componentIndex;
      }
    }
    residual[componentCount] = sumFirst - sumSecond;
    double retainedMaterialBalanceResidual = Math.abs(residual[componentCount]);
    thermodynamicMaximumResidual = Math.max(thermodynamicMaximumResidual, retainedMaterialBalanceResidual);
    double incipientFugacityMaximumResidual = 0.0;
    int incipientFugacityWorstComponentIndex = -1;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double reference = Math
          .log(Math.max(working.getPhase(0).getComponent(componentIndex).getx(), MINIMUM_COMPOSITION))
          + working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
      residual[componentCount + 1 + componentIndex] = variables[logWStart + componentIndex]
          + working.getPhase(2).getComponent(componentIndex).getLogFugacityCoefficient() - reference;
      thermodynamicMaximumResidual = Math.max(thermodynamicMaximumResidual,
          Math.abs(residual[componentCount + 1 + componentIndex]));
      if (Math.abs(residual[componentCount + 1 + componentIndex]) > incipientFugacityMaximumResidual) {
        incipientFugacityMaximumResidual = Math.abs(residual[componentCount + 1 + componentIndex]);
        incipientFugacityWorstComponentIndex = componentIndex;
      }
    }
    residual[2 * componentCount + 1] = logSumW;
    double incipientNormalizationResidual = Math.abs(logSumW);
    thermodynamicMaximumResidual = Math.max(thermodynamicMaximumResidual, incipientNormalizationResidual);
    EquationResiduals equationResiduals = new EquationResiduals(retainedFugacityMaximumResidual,
        retainedMaterialBalanceResidual, incipientFugacityMaximumResidual, incipientNormalizationResidual,
        thermodynamicMaximumResidual, componentName(retainedFugacityWorstComponentIndex),
        componentName(incipientFugacityWorstComponentIndex));
    double deflationMultiplier = retainedPhaseBifurcation && deflationPower > 0.0
        ? retainedPhaseDeflationMultiplier(first, second, third)
        : 1.0;
    if (deflationMultiplier != 1.0) {
      for (int residualIndex = 0; residualIndex <= 2 * componentCount + 1; residualIndex++) {
        residual[residualIndex] *= deflationMultiplier;
      }
    }
    double arcResidual = 0.0;
    for (int variableIndex = 0; variableIndex < variables.length; variableIndex++) {
      arcResidual += (variables[variableIndex] - predictor[variableIndex]) * tangent[variableIndex];
    }
    residual[2 * componentCount + 2] = arcResidual;
    double maximumResidual = Math.max(thermodynamicMaximumResidual * deflationMultiplier, Math.abs(arcResidual));
    if (!Double.isFinite(maximumResidual)) {
      throw new IllegalStateException("pseudo-arclength residual is non-finite");
    }
    StateFingerprint fingerprint = new StateFingerprint(componentCount, componentNames(componentCount), z, beta,
        temperatureK, pressureBara, equilibriumRatios, first, second, third, residual[componentCount],
        phaseModelIdentities(working));
    return new Evaluation(residual, maximumResidual, thermodynamicMaximumResidual, arcResidual, equationResiduals,
        temperatureK, pressureBara, beta, normalizedExponentials(logFirst, 0, componentCount, logSumFirst),
        normalizedExponentials(logSecond, 0, componentCount, logSumSecond), third, fingerprint);
  }

  /**
   * Evaluates {@code log(1 - beta + beta * exp(logK))} without forming {@code exp(logK)} when it would overflow.
   *
   * <p>
   * Heavy pseudo-components reach equilibrium ratios near {@code 1e99} against water, so the linear form overflows long
   * before the physics does.
   * </p>
   */
  private static double retainedPhaseLogDenominator(double beta, double logK) {
    double retainedFraction = 1.0 - beta;
    if (logK > 0.0) {
      return logK + Math.log(beta + retainedFraction * Math.exp(-logK));
    }
    return Math.log(retainedFraction + beta * Math.exp(logK));
  }

  private String[] componentNames(int componentCount) {
    String[] names = new String[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      names[componentIndex] = template.getPhase(0).getComponent(componentIndex).getComponentName();
    }
    return names;
  }

  private static String[] phaseModelIdentities(SystemInterface system) {
    String[] identities = new String[3];
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      PhaseInterface phase = system.getPhase(phaseIndex);
      StringBuilder identity = new StringBuilder(phase.getClass().getSimpleName());
      identity.append('/').append(phase.getType());
      if (phase instanceof PhaseEos) {
        identity.append('/').append(((PhaseEos) phase).getMixingRuleName());
      }
      identities[phaseIndex] = identity.toString();
    }
    return identities;
  }

  private String componentName(int componentIndex) {
    return componentIndex < 0 ? null : template.getPhase(0).getComponent(componentIndex).getComponentName();
  }

  private double retainedPhaseDeflationMultiplier(double[] first, double[] second, double[] third) {
    double[] retained = incipientPhase == retainedPhaseZero ? first : second;
    double squaredDistance = 0.0;
    for (int componentIndex = 0; componentIndex < retained.length; componentIndex++) {
      double difference = retained[componentIndex] - third[componentIndex];
      squaredDistance += difference * difference;
    }
    double distance = Math.sqrt(squaredDistance);
    if (!(distance > 1.0e-12) || !Double.isFinite(distance)) {
      throw new IllegalStateException("deflated retained and bifurcating phase slots coalesced");
    }
    return Math.pow(distance, -deflationPower) + deflationShift;
  }

  private SystemInterface createWorkingSystem(double temperatureK, double pressureBara, double beta, double[] first,
      double[] second, double[] third, double incipientFraction) {
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    double numericalThirdPhaseFraction = Math.max(incipientFraction, MINIMUM_BETA);
    working.setBeta(0, beta * (1.0 - numericalThirdPhaseFraction));
    working.setBeta(1, (1.0 - beta) * (1.0 - numericalThirdPhaseFraction));
    working.setBeta(2, numericalThirdPhaseFraction);
    setPhase(working, 0, retainedPhaseZero, first);
    setPhase(working, 1, retainedPhaseOne, second);
    setPhase(working, 2, incipientPhase, third);
    for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
      working.init(1, phaseIndex);
    }
    return working;
  }

  private static void setPhase(SystemInterface system, int phaseIndex, CandidatePhase phase, double[] composition) {
    system.setPhaseType(phaseIndex, toPhaseType(phase));
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      // The EOS requires a strictly positive mole fraction. This floor is deliberately confined to the EOS input:
      // a component this dilute is indistinguishable from zero for mixture properties, whereas flooring the
      // material-balance inventory itself would corrupt the equations being solved.
      system.getPhase(phaseIndex).getComponent(componentIndex)
          .setx(Math.max(composition[componentIndex], MINIMUM_COMPOSITION));
    }
    system.getPhase(phaseIndex).normalize();
    system.setPhaseType(phaseIndex, toPhaseType(phase));
  }

  /** Returns the value itself when it is a usable positive number, and only then falls back to the floor. */
  private static double positiveOrFloor(double composition) {
    return composition > 0.0 && Double.isFinite(composition) ? composition : MINIMUM_COMPOSITION;
  }

  static double[] transformedVariables(State state) {
    if (state.solverVariables != null) {
      return state.solverVariables.clone();
    }
    int componentCount = state.phaseZeroComposition.length;
    double[] variables = new double[2 * componentCount + 3];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      // The floor guards against a composition of exactly zero and nothing else. Applying it to a value that is
      // merely tiny truncates it: a heavy pseudo-component against water sits near 1e-102, which double precision
      // represents without difficulty, and clamping it to 1e-100 shifts ln(K) by ln(1e-100 / x). That shift is
      // what made an otherwise clean boundary point replay with a residual of 2.7686 instead of ~1e-13.
      variables[componentIndex] = Math.log(positiveOrFloor(state.phaseZeroComposition[componentIndex])
          / positiveOrFloor(state.phaseOneComposition[componentIndex]));
      variables[componentCount + 1 + componentIndex] = Math
          .log(positiveOrFloor(state.incipientComposition[componentIndex]));
    }
    variables[componentCount] = state.beta;
    variables[2 * componentCount + 1] = Math.log(state.temperatureK);
    variables[2 * componentCount + 2] = Math.log(state.pressureBara);
    return variables;
  }

  /** Re-evaluates the thermodynamic equations of a continuation state without a predictor contribution. */
  double replayThermodynamicResidual(State state) {
    return replayEquationResiduals(state).getThermodynamicMaximumResidual();
  }

  /** Re-evaluates each physical equation group of a continuation state. */
  EquationResiduals replayEquationResiduals(State state) {
    return replayEvaluation(state).equationResiduals;
  }

  /**
   * Rebuilds the transformed variables from a continuation state and re-evaluates every equation group, keeping the
   * full inventory fingerprint so the reconstruction can be compared slot by slot against the solving path.
   */
  Evaluation replayEvaluation(State state) {
    if (state == null || state.getRetainedPhaseZero() != retainedPhaseZero
        || state.getRetainedPhaseOne() != retainedPhaseOne || state.getIncipientPhase() != incipientPhase) {
      throw new IllegalArgumentException("continuation state does not match this corrector phase topology");
    }
    double[] variables = transformedVariables(state);
    return evaluate(variables, variables, new double[variables.length], 0.0);
  }

  static double[] continuationMetricWeights(State state) {
    int componentCount = state.phaseZeroComposition.length;
    double[] weights = new double[2 * componentCount + 3];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double retainedRelevance = 0.5
          * (state.phaseZeroComposition[componentIndex] + state.phaseOneComposition[componentIndex]);
      weights[componentIndex] = Math.sqrt(Math.max(retainedRelevance, 1.0e-12));
      weights[componentCount + 1 + componentIndex] = Math
          .sqrt(Math.max(state.incipientComposition[componentIndex], 1.0e-12));
    }
    weights[componentCount] = 1.0;
    weights[2 * componentCount + 1] = 1.0;
    weights[2 * componentCount + 2] = 1.0;
    return weights;
  }

  private void validateState(State state) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    if (state == null || state.retainedPhaseZero != retainedPhaseZero || state.retainedPhaseOne != retainedPhaseOne
        || state.incipientPhase != incipientPhase || state.phaseZeroComposition.length != componentCount
        || state.phaseOneComposition.length != componentCount || state.incipientComposition.length != componentCount
        || !Double.isFinite(state.temperatureK) || state.temperatureK < 50.0 || !Double.isFinite(state.pressureBara)
        || state.pressureBara <= 0.0 || !Double.isFinite(state.beta) || state.beta <= 0.0 || state.beta >= 1.0) {
      throw new IllegalArgumentException("state does not match the pseudo-arclength corrector");
    }
  }

  private void clampPhysicalVariables(double[] variables) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      variables[componentIndex] = clamp(variables[componentIndex], MINIMUM_LOG_VALUE, MAXIMUM_LOG_K_VALUE);
      variables[componentCount + 1 + componentIndex] = clamp(variables[componentCount + 1 + componentIndex],
          MINIMUM_LOG_VALUE, MAXIMUM_LOG_W_VALUE);
    }
    variables[componentCount] = clamp(variables[componentCount], MINIMUM_BETA, 1.0 - MINIMUM_BETA);
    variables[2 * componentCount + 1] = clamp(variables[2 * componentCount + 1], Math.log(50.0), Math.log(2500.0));
    variables[2 * componentCount + 2] = clamp(variables[2 * componentCount + 2], Math.log(1.0e-6), Math.log(1.0e6));
  }

  private void limitCorrection(Matrix correction) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    double scale = 1.0;
    for (int variableIndex = 0; variableIndex < correction.getRowDimension(); variableIndex++) {
      double maximumStep;
      if (variableIndex == componentCount) {
        maximumStep = 0.05;
      } else if (variableIndex >= 2 * componentCount + 1) {
        maximumStep = 0.04;
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

  private double[] overallComposition() {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    double[] composition = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(),
          MINIMUM_COMPOSITION);
    }
    return normalized(composition);
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

  private static double[] normalized(double[] values) {
    double[] normalized = values.clone();
    double total = 0.0;
    for (double value : normalized) {
      total += value;
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalStateException("composition cannot be normalized");
    }
    for (int index = 0; index < normalized.length; index++) {
      normalized[index] /= total;
    }
    return normalized;
  }

  private static boolean distinctPhases(Evaluation evaluation) {
    return compositionDistance(evaluation.phaseZeroComposition, evaluation.phaseOneComposition) > 1.0e-5
        && compositionDistance(evaluation.phaseZeroComposition, evaluation.incipientComposition) > 1.0e-5
        && compositionDistance(evaluation.phaseOneComposition, evaluation.incipientComposition) > 1.0e-5;
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }

  private static double clamp(double value, double lower, double upper) {
    return Math.max(lower, Math.min(upper, value));
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
    private final double arcLengthResidual;
    private final EquationResiduals equationResiduals;
    private final double temperatureK;
    private final double pressureBara;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double[] incipientComposition;
    private final StateFingerprint fingerprint;

    private Evaluation(double[] residual, double maximumResidual, double thermodynamicMaximumResidual,
        double arcLengthResidual, EquationResiduals equationResiduals, double temperatureK, double pressureBara,
        double beta, double[] phaseZeroComposition, double[] phaseOneComposition, double[] incipientComposition,
        StateFingerprint fingerprint) {
      this.fingerprint = fingerprint;
      this.residual = residual;
      this.maximumResidual = maximumResidual;
      this.thermodynamicMaximumResidual = thermodynamicMaximumResidual;
      this.arcLengthResidual = arcLengthResidual;
      this.equationResiduals = equationResiduals;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.beta = beta;
      this.phaseZeroComposition = phaseZeroComposition;
      this.phaseOneComposition = phaseOneComposition;
      this.incipientComposition = incipientComposition;
    }
  }

  /**
   * Immutable record of the exact inventory and variables one residual evaluation actually used.
   *
   * <p>
   * Every array is stored in the internal component order of the corrector template, so two fingerprints taken on
   * different construction paths can be compared slot by slot. This is what distinguishes a genuine thermodynamic
   * disagreement from a component-order, overall-composition or state-reconstruction inconsistency.
   * </p>
   */
  public static final class StateFingerprint {
    private static final String EQUILIBRIUM_RATIO_DEFINITION = "K[i] = exp(clamp(variables[i])) = unnormalizedPhaseZeroComposition[i] / unnormalizedPhaseOneComposition[i]";

    private final int componentCount;
    private final String[] componentNames;
    private final double[] overallComposition;
    private final double overallCompositionSum;
    private final double beta;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] equilibriumRatios;
    private final double[] unnormalizedPhaseZeroComposition;
    private final double[] unnormalizedPhaseOneComposition;
    private final double[] incipientComposition;
    private final double[] materialBalanceByComponent;
    private final double rachfordRiceResidual;
    private final String[] phaseModelIdentities;

    private StateFingerprint(int componentCount, String[] componentNames, double[] overallComposition, double beta,
        double temperatureK, double pressureBara, double[] equilibriumRatios, double[] unnormalizedPhaseZeroComposition,
        double[] unnormalizedPhaseOneComposition, double[] incipientComposition, double rachfordRiceResidual,
        String[] phaseModelIdentities) {
      this.componentCount = componentCount;
      this.componentNames = componentNames.clone();
      this.overallComposition = overallComposition.clone();
      double sum = 0.0;
      for (int componentIndex = 0; componentIndex < overallComposition.length; componentIndex++) {
        sum += overallComposition[componentIndex];
      }
      this.overallCompositionSum = sum;
      this.beta = beta;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.equilibriumRatios = equilibriumRatios.clone();
      this.unnormalizedPhaseZeroComposition = unnormalizedPhaseZeroComposition.clone();
      this.unnormalizedPhaseOneComposition = unnormalizedPhaseOneComposition.clone();
      this.incipientComposition = incipientComposition.clone();
      this.materialBalanceByComponent = new double[componentCount];
      for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
        this.materialBalanceByComponent[componentIndex] = unnormalizedPhaseZeroComposition[componentIndex]
            - unnormalizedPhaseOneComposition[componentIndex];
      }
      this.rachfordRiceResidual = rachfordRiceResidual;
      this.phaseModelIdentities = phaseModelIdentities.clone();
    }

    public int getComponentCount() {
      return componentCount;
    }

    public String[] getComponentNames() {
      return componentNames.clone();
    }

    public double[] getOverallComposition() {
      return overallComposition.clone();
    }

    public double getOverallCompositionSum() {
      return overallCompositionSum;
    }

    public double getBeta() {
      return beta;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double[] getEquilibriumRatios() {
      return equilibriumRatios.clone();
    }

    public String getEquilibriumRatioDefinition() {
      return EQUILIBRIUM_RATIO_DEFINITION;
    }

    public double[] getUnnormalizedPhaseZeroComposition() {
      return unnormalizedPhaseZeroComposition.clone();
    }

    public double[] getUnnormalizedPhaseOneComposition() {
      return unnormalizedPhaseOneComposition.clone();
    }

    public double[] getIncipientComposition() {
      return incipientComposition.clone();
    }

    public double[] getMaterialBalanceByComponent() {
      return materialBalanceByComponent.clone();
    }

    public double getRachfordRiceResidual() {
      return rachfordRiceResidual;
    }

    public String[] getPhaseModelIdentities() {
      return phaseModelIdentities.clone();
    }
  }

  /** Immutable unscaled residuals for the four thermodynamic equation groups. */
  public static final class EquationResiduals {
    private final double retainedFugacityMaximumResidual;
    private final double retainedMaterialBalanceResidual;
    private final double incipientFugacityMaximumResidual;
    private final double incipientNormalizationResidual;
    private final double thermodynamicMaximumResidual;
    private final String retainedFugacityWorstComponent;
    private final String incipientFugacityWorstComponent;

    private EquationResiduals(double retainedFugacityMaximumResidual, double retainedMaterialBalanceResidual,
        double incipientFugacityMaximumResidual, double incipientNormalizationResidual,
        double thermodynamicMaximumResidual, String retainedFugacityWorstComponent,
        String incipientFugacityWorstComponent) {
      this.retainedFugacityMaximumResidual = retainedFugacityMaximumResidual;
      this.retainedMaterialBalanceResidual = retainedMaterialBalanceResidual;
      this.incipientFugacityMaximumResidual = incipientFugacityMaximumResidual;
      this.incipientNormalizationResidual = incipientNormalizationResidual;
      this.thermodynamicMaximumResidual = thermodynamicMaximumResidual;
      this.retainedFugacityWorstComponent = retainedFugacityWorstComponent;
      this.incipientFugacityWorstComponent = incipientFugacityWorstComponent;
    }

    public double getRetainedFugacityMaximumResidual() {
      return retainedFugacityMaximumResidual;
    }

    public double getRetainedMaterialBalanceResidual() {
      return retainedMaterialBalanceResidual;
    }

    public double getIncipientFugacityMaximumResidual() {
      return incipientFugacityMaximumResidual;
    }

    public double getIncipientNormalizationResidual() {
      return incipientNormalizationResidual;
    }

    public double getThermodynamicMaximumResidual() {
      return thermodynamicMaximumResidual;
    }

    public String getRetainedFugacityWorstComponent() {
      return retainedFugacityWorstComponent;
    }

    public String getIncipientFugacityWorstComponent() {
      return incipientFugacityWorstComponent;
    }
  }

  /** Immutable accepted boundary state used as a pseudo-arclength seed. */
  public static final class State {
    private final CandidatePhase retainedPhaseZero;
    private final CandidatePhase retainedPhaseOne;
    private final CandidatePhase incipientPhase;
    private final double temperatureK;
    private final double pressureBara;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double[] incipientComposition;
    private final double[] solverVariables;

    /** Builds a seed state from a converged fixed-pressure boundary point. */
    public static State from(TwoToThreePhaseBoundaryPointSolver.Result point) {
      if (point == null || !point.isConverged()) {
        throw new IllegalArgumentException("a converged fixed-pressure boundary point is required");
      }
      return new State(point.getRetainedPhaseZero(), point.getRetainedPhaseOne(), point.getIncipientPhase(),
          point.getTemperatureK(), point.getPressureBara(), point.getBeta(), point.getPhaseZeroComposition(),
          point.getPhaseOneComposition(), point.getIncipientComposition());
    }

    /** Builds an explicit seed from retained two-phase and candidate incipient compositions. */
    public static State create(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase incipientPhase, double temperatureK, double pressureBara, double beta,
        double[] phaseZeroComposition, double[] phaseOneComposition, double[] incipientComposition) {
      return create(retainedPhaseZero, retainedPhaseOne, incipientPhase, temperatureK, pressureBara, beta,
          phaseZeroComposition, phaseOneComposition, incipientComposition, null);
    }

    /**
     * Builds a seed that also carries the transformed solver variables it was produced from.
     *
     * <p>
     * Compositions alone cannot reproduce {@code ln(K)} once any mole fraction falls below the composition floor used
     * when writing the EOS input: a heavy pseudo-component in the aqueous phase reaches {@code 1e-102}, so
     * reconstructing {@code K = x0 / x1} from stored compositions silently truncates it. Carrying the variables makes a
     * replay an exact re-evaluation of the same equations instead of a lossy round trip.
     * </p>
     */
    public static State create(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase incipientPhase, double temperatureK, double pressureBara, double beta,
        double[] phaseZeroComposition, double[] phaseOneComposition, double[] incipientComposition,
        double[] solverVariables) {
      if (retainedPhaseZero == null || retainedPhaseOne == null || incipientPhase == null
          || phaseZeroComposition == null || phaseOneComposition == null || incipientComposition == null) {
        throw new IllegalArgumentException("all phase identities and compositions are required");
      }
      if (solverVariables != null && solverVariables.length != 2 * phaseZeroComposition.length + 3) {
        throw new IllegalArgumentException("solver variables must hold 2 * componentCount + 3 entries");
      }
      return new State(retainedPhaseZero, retainedPhaseOne, incipientPhase, temperatureK, pressureBara, beta,
          phaseZeroComposition, phaseOneComposition, incipientComposition, solverVariables);
    }

    private State(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase,
        double temperatureK, double pressureBara, double beta, double[] phaseZeroComposition,
        double[] phaseOneComposition, double[] incipientComposition) {
      this(retainedPhaseZero, retainedPhaseOne, incipientPhase, temperatureK, pressureBara, beta, phaseZeroComposition,
          phaseOneComposition, incipientComposition, null);
    }

    private State(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase,
        double temperatureK, double pressureBara, double beta, double[] phaseZeroComposition,
        double[] phaseOneComposition, double[] incipientComposition, double[] solverVariables) {
      this.solverVariables = solverVariables == null ? null : solverVariables.clone();
      this.retainedPhaseZero = retainedPhaseZero;
      this.retainedPhaseOne = retainedPhaseOne;
      this.incipientPhase = incipientPhase;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.beta = beta;
      this.phaseZeroComposition = phaseZeroComposition.clone();
      this.phaseOneComposition = phaseOneComposition.clone();
      this.incipientComposition = incipientComposition.clone();
    }

    /**
     * Transformed solver variables {@code [ln(K), beta, ln(W_3), ln(T), ln(P)]} this state was produced from, or
     * {@code null} for seeds that only carry compositions.
     */
    public double[] getSolverVariables() {
      return solverVariables == null ? null : solverVariables.clone();
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public CandidatePhase getRetainedPhaseZero() {
      return retainedPhaseZero;
    }

    public CandidatePhase getRetainedPhaseOne() {
      return retainedPhaseOne;
    }

    public CandidatePhase getIncipientPhase() {
      return incipientPhase;
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

    public double[] getIncipientComposition() {
      return incipientComposition.clone();
    }
  }

  /** Immutable pseudo-arclength correction result and residual diagnostics. */
  public static final class Result {
    private final State state;
    private final double initialMaximumResidual;
    private final double maximumResidual;
    private final double thermodynamicMaximumResidual;
    private final double arcLengthResidual;
    private final EquationResiduals equationResiduals;
    private final EquationResiduals independentReplayEquationResiduals;
    private final EquationResiduals returnedStateReplayEquationResiduals;
    private final double arcStep;
    private final int iterations;
    private final double jacobianConditionNumber;
    private final boolean residualConverged;
    private final boolean distinctPhases;
    private final String failureMessage;
    private final StateFingerprint solveFingerprint;
    private final StateFingerprint returnedStateReplayFingerprint;
    private final StateFingerprint independentReplayFingerprint;

    private Result(State state, double initialMaximumResidual, double maximumResidual,
        double thermodynamicMaximumResidual, double arcLengthResidual, EquationResiduals equationResiduals,
        EquationResiduals independentReplayEquationResiduals, EquationResiduals returnedStateReplayEquationResiduals,
        double arcStep, int iterations, double jacobianConditionNumber, boolean residualConverged,
        boolean distinctPhases, String failureMessage, StateFingerprint solveFingerprint,
        StateFingerprint returnedStateReplayFingerprint, StateFingerprint independentReplayFingerprint) {
      this.returnedStateReplayEquationResiduals = returnedStateReplayEquationResiduals;
      this.solveFingerprint = solveFingerprint;
      this.returnedStateReplayFingerprint = returnedStateReplayFingerprint;
      this.independentReplayFingerprint = independentReplayFingerprint;
      this.state = state;
      this.initialMaximumResidual = initialMaximumResidual;
      this.maximumResidual = maximumResidual;
      this.thermodynamicMaximumResidual = thermodynamicMaximumResidual;
      this.arcLengthResidual = arcLengthResidual;
      this.equationResiduals = equationResiduals;
      this.independentReplayEquationResiduals = independentReplayEquationResiduals;
      this.arcStep = arcStep;
      this.iterations = iterations;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.residualConverged = residualConverged;
      this.distinctPhases = distinctPhases;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase incipientPhase, double arcStep, String failureMessage) {
      State failureState = new State(retainedPhaseZero, retainedPhaseOne, incipientPhase, Double.NaN, Double.NaN,
          Double.NaN, new double[0], new double[0], new double[0]);
      return new Result(failureState, Double.NaN, Double.NaN, Double.NaN, Double.NaN, null, null, null, arcStep, 0,
          Double.NaN, false, false, failureMessage, null, null, null);
    }

    public boolean isConverged() {
      return residualConverged && distinctPhases && failureMessage == null;
    }

    public State getState() {
      return state;
    }

    public double getInitialMaximumResidual() {
      return initialMaximumResidual;
    }

    public double getMaximumResidual() {
      return maximumResidual;
    }

    public double getThermodynamicMaximumResidual() {
      return thermodynamicMaximumResidual;
    }

    public double getArcLengthResidual() {
      return arcLengthResidual;
    }

    public EquationResiduals getEquationResiduals() {
      return equationResiduals;
    }

    public EquationResiduals getIndependentReplayEquationResiduals() {
      return independentReplayEquationResiduals;
    }

    /** Residuals obtained by rebuilding the variables from the returned state on the solving corrector instance. */
    public EquationResiduals getReturnedStateReplayEquationResiduals() {
      return returnedStateReplayEquationResiduals;
    }

    /** Inventory actually used by the converged solving evaluation. */
    public StateFingerprint getSolveFingerprint() {
      return solveFingerprint;
    }

    /** Inventory rebuilt from the returned state on the solving corrector instance. */
    public StateFingerprint getReturnedStateReplayFingerprint() {
      return returnedStateReplayFingerprint;
    }

    /** Inventory rebuilt from the returned state on a freshly constructed corrector instance. */
    public StateFingerprint getIndependentReplayFingerprint() {
      return independentReplayFingerprint;
    }

    public boolean isIndependentlyReproducible() {
      return independentReplayEquationResiduals != null && failureMessage == null;
    }

    public double getArcStep() {
      return arcStep;
    }

    public int getIterations() {
      return iterations;
    }

    public double getJacobianConditionNumber() {
      return jacobianConditionNumber;
    }

    public boolean hasDistinctPhases() {
      return distinctPhases;
    }

    public boolean isResidualConverged() {
      return residualConverged;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
