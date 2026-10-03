package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import Jama.Matrix;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Damped-Newton solver for non-trivial TPD stationary points in unnormalised log-mole variables. */
public final class IncipientPhaseStationaryPointSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;
  private static final double MINIMUM_LOG_VALUE = -700.0;
  private static final double MAXIMUM_LOG_VALUE = 100.0;

  private final SystemInterface reference;
  private final CandidatePhase candidatePhase;
  private int maximumIterations = 60;
  private double residualTolerance = 1.0e-9;
  private double finiteDifferenceStep = 2.0e-5;

  /** Creates a stationary-point solver for one requested physical trial phase. */
  public IncipientPhaseStationaryPointSolver(SystemInterface reference, CandidatePhase candidatePhase) {
    if (reference == null || candidatePhase == null) {
      throw new IllegalArgumentException("reference system and candidate phase are required");
    }
    this.reference = reference.clone();
    this.candidatePhase = candidatePhase;
  }

  /** Sets Newton iteration, residual, and numerical-Jacobian controls. */
  public IncipientPhaseStationaryPointSolver setNumericalControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    if (maximumIterations < 1 || !Double.isFinite(residualTolerance) || residualTolerance <= 0.0
        || !Double.isFinite(finiteDifferenceStep) || finiteDifferenceStep <= 0.0) {
      throw new IllegalArgumentException("invalid stationary-point solver controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /** Solves the Michelsen stationary equations from an explicit composition seed. */
  public Result solve(double[] compositionSeed) {
    int componentCount = reference.getPhase(0).getNumberOfComponents();
    if (compositionSeed == null || compositionSeed.length != componentCount) {
      throw new IllegalArgumentException("composition seed must contain one value per component");
    }
    double[] normalizedSeed = compositionSeed.clone();
    normalize(normalizedSeed);
    double[] logW = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      logW[componentIndex] = Math.log(Math.max(normalizedSeed[componentIndex], MINIMUM_COMPOSITION));
    }

    EvaluationContext context;
    Evaluation evaluation;
    try {
      context = createContext();
      evaluation = evaluate(context, logW);
    } catch (RuntimeException error) {
      return Result.failure(candidatePhase, componentCount, error.getMessage());
    }
    double initialResidual = evaluation.maximumResidual;
    double conditionNumber = Double.NaN;
    int iterations = 0;
    String failureMessage = null;
    while (iterations < maximumIterations && evaluation.maximumResidual > residualTolerance) {
      iterations++;
      Matrix jacobian;
      try {
        jacobian = numericalJacobian(context, logW, evaluation.residual);
        conditionNumber = jacobian.cond();
      } catch (RuntimeException error) {
        failureMessage = "stationary-point Jacobian evaluation failed: " + error.getMessage();
        break;
      }
      Matrix correction;
      try {
        correction = jacobian.solve(new Matrix(evaluation.residual, componentCount));
      } catch (RuntimeException error) {
        correction = regularizedCorrection(jacobian, evaluation.residual);
      }
      if (correction == null || !Double.isFinite(correction.normInf())) {
        failureMessage = "stationary-point Newton correction is non-finite";
        break;
      }
      limitCorrection(correction, 2.0);
      boolean accepted = false;
      double damping = 1.0;
      for (int lineSearch = 0; lineSearch < 20; lineSearch++) {
        double[] trialLogW = logW.clone();
        for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
          trialLogW[componentIndex] = clamp(trialLogW[componentIndex] - damping * correction.get(componentIndex, 0));
        }
        try {
          Evaluation trial = evaluate(context, trialLogW);
          if (trial.maximumResidual < evaluation.maximumResidual) {
            logW = trialLogW;
            evaluation = trial;
            accepted = true;
            break;
          }
        } catch (RuntimeException error) {
          // Retry with a shorter step.
        }
        damping *= 0.5;
      }
      if (!accepted) {
        failureMessage = "line search could not reduce the stationary-point residual";
        break;
      }
    }
    boolean converged = evaluation.maximumResidual <= residualTolerance;
    if (!converged && failureMessage == null) {
      failureMessage = "maximum stationary-point Newton iteration count reached";
    }
    boolean trivial = isTrivial(evaluation.composition);
    CandidatePhase physicalPhase = evaluation.physicalPhase;
    return new Result(candidatePhase, physicalPhase, evaluation.composition, evaluation.tangentPlaneDistance,
        evaluation.maximumResidual, initialResidual, iterations, conditionNumber, converged, trivial, failureMessage);
  }

  private EvaluationContext createContext() {
    SystemInterface initializedReference = reference.clone();
    initializedReference.init(1);
    int componentCount = initializedReference.getPhase(0).getNumberOfComponents();
    double[] logReferenceFugacity = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double largestComposition = -1.0;
      double selectedLogFugacity = Double.NaN;
      for (int phaseIndex = 0; phaseIndex < initializedReference.getNumberOfPhases(); phaseIndex++) {
        double composition = initializedReference.getPhase(phaseIndex).getComponent(componentIndex).getx();
        if (composition > largestComposition && composition > MINIMUM_COMPOSITION) {
          double value = Math.log(composition)
              + initializedReference.getPhase(phaseIndex).getComponent(componentIndex).getLogFugacityCoefficient();
          if (Double.isFinite(value)) {
            largestComposition = composition;
            selectedLogFugacity = value;
          }
        }
      }
      if (!Double.isFinite(selectedLogFugacity)) {
        throw new IllegalStateException("no finite reference fugacity for component " + componentIndex);
      }
      logReferenceFugacity[componentIndex] = selectedLogFugacity;
    }
    SystemInterface trial = initializedReference.clone();
    if (trial.getNumberOfPhases() < 2) {
      trial.setNumberOfPhases(2);
    }
    trial.setPhaseType(1, toPhaseType(candidatePhase));
    return new EvaluationContext(initializedReference, trial, logReferenceFugacity);
  }

  private Evaluation evaluate(EvaluationContext context, double[] logW) {
    double logSumW = logSumExp(logW);
    double[] composition = new double[logW.length];
    for (int componentIndex = 0; componentIndex < logW.length; componentIndex++) {
      composition[componentIndex] = Math.exp(Math.max(MINIMUM_LOG_VALUE, logW[componentIndex] - logSumW));
      context.trial.getPhase(1).getComponent(componentIndex).setx(composition[componentIndex]);
    }
    context.trial.getPhase(1).normalize();
    context.trial.setPhaseType(1, toPhaseType(candidatePhase));
    context.trial.init(1, 1);
    double[] residual = new double[logW.length];
    double maximumResidual = 0.0;
    for (int componentIndex = 0; componentIndex < logW.length; componentIndex++) {
      residual[componentIndex] = logW[componentIndex]
          + context.trial.getPhase(1).getComponent(componentIndex).getLogFugacityCoefficient()
          - context.logReferenceFugacity[componentIndex];
      maximumResidual = Math.max(maximumResidual, Math.abs(residual[componentIndex]));
    }
    if (!Double.isFinite(maximumResidual)) {
      throw new IllegalStateException("stationary-point residual is non-finite");
    }
    return new Evaluation(residual, maximumResidual, -logSumW, composition,
        IncipientPhaseStabilityAnalyzer.classifyPhysicalPhase(context.trial.getPhase(1)));
  }

  private Matrix numericalJacobian(EvaluationContext context, double[] logW, double[] baseResidual) {
    int componentCount = logW.length;
    Matrix jacobian = new Matrix(componentCount, componentCount);
    for (int column = 0; column < componentCount; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(logW[column]));
      double[] plus = logW.clone();
      double[] minus = logW.clone();
      plus[column] = clamp(plus[column] + step);
      minus[column] = clamp(minus[column] - step);
      double denominator = plus[column] - minus[column];
      double[] plusResidual = evaluate(context, plus).residual;
      double[] minusResidual;
      try {
        minusResidual = evaluate(context, minus).residual;
      } catch (RuntimeException error) {
        minusResidual = baseResidual;
        denominator = plus[column] - logW[column];
      }
      for (int row = 0; row < componentCount; row++) {
        jacobian.set(row, column, (plusResidual[row] - minusResidual[row]) / denominator);
      }
    }
    return jacobian;
  }

  private static Matrix regularizedCorrection(Matrix jacobian, double[] residual) {
    try {
      Matrix transpose = jacobian.transpose();
      Matrix normal = transpose.times(jacobian);
      double regularization = 1.0e-8;
      for (int diagonal = 0; diagonal < normal.getRowDimension(); diagonal++) {
        normal.set(diagonal, diagonal, normal.get(diagonal, diagonal) + regularization);
      }
      return normal.solve(transpose.times(new Matrix(residual, residual.length)));
    } catch (RuntimeException error) {
      return null;
    }
  }

  private boolean isTrivial(double[] composition) {
    for (int phaseIndex = 0; phaseIndex < reference.getNumberOfPhases(); phaseIndex++) {
      double distance = 0.0;
      for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
        distance += Math
            .abs(composition[componentIndex] - reference.getPhase(phaseIndex).getComponent(componentIndex).getx());
      }
      if (distance <= 1.0e-7) {
        return true;
      }
    }
    return false;
  }

  private static double logSumExp(double[] values) {
    double maximum = Double.NEGATIVE_INFINITY;
    for (double value : values) {
      maximum = Math.max(maximum, value);
    }
    double sum = 0.0;
    for (double value : values) {
      sum += Math.exp(value - maximum);
    }
    return maximum + Math.log(sum);
  }

  private static void normalize(double[] composition) {
    double total = 0.0;
    for (double value : composition) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException("composition seed contains an invalid value");
      }
      total += value;
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalArgumentException("composition seed cannot be normalized");
    }
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] /= total;
    }
  }

  private static void limitCorrection(Matrix correction, double maximumMagnitude) {
    double scale = 1.0;
    for (int row = 0; row < correction.getRowDimension(); row++) {
      double magnitude = Math.abs(correction.get(row, 0));
      if (magnitude > maximumMagnitude) {
        scale = Math.min(scale, maximumMagnitude / magnitude);
      }
    }
    if (scale < 1.0) {
      correction.timesEquals(scale);
    }
  }

  private static double clamp(double value) {
    return Math.max(MINIMUM_LOG_VALUE, Math.min(MAXIMUM_LOG_VALUE, value));
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

  private static final class EvaluationContext {
    private final SystemInterface reference;
    private final SystemInterface trial;
    private final double[] logReferenceFugacity;

    private EvaluationContext(SystemInterface reference, SystemInterface trial, double[] logReferenceFugacity) {
      this.reference = reference;
      this.trial = trial;
      this.logReferenceFugacity = logReferenceFugacity;
    }
  }

  private static final class Evaluation {
    private final double[] residual;
    private final double maximumResidual;
    private final double tangentPlaneDistance;
    private final double[] composition;
    private final CandidatePhase physicalPhase;

    /**
     * Captures diagnostics and root identity before subsequent Jacobian trials mutate the workspace.
     *
     * @param residual stationarity equations
     * @param maximumResidual largest absolute residual
     * @param tangentPlaneDistance trial tangent-plane distance
     * @param composition normalized trial composition
     * @param physicalPhase evaluated EOS root family
     */
    private Evaluation(double[] residual, double maximumResidual, double tangentPlaneDistance, double[] composition,
        CandidatePhase physicalPhase) {
      this.physicalPhase = physicalPhase;
      this.residual = residual;
      this.maximumResidual = maximumResidual;
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.composition = composition.clone();
    }
  }

  /** Immutable stationary-point solution and numerical diagnostics. */
  public static final class Result {
    private final CandidatePhase seedPhase;
    private final CandidatePhase physicalPhase;
    private final double[] composition;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;
    private final double initialResidual;
    private final int iterations;
    private final double jacobianConditionNumber;
    private final boolean converged;
    private final boolean trivial;
    private final String failureMessage;

    private Result(CandidatePhase seedPhase, CandidatePhase physicalPhase, double[] composition,
        double tangentPlaneDistance, double stationarityResidual, double initialResidual, int iterations,
        double jacobianConditionNumber, boolean converged, boolean trivial, String failureMessage) {
      this.seedPhase = seedPhase;
      this.physicalPhase = physicalPhase;
      this.composition = composition.clone();
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.stationarityResidual = stationarityResidual;
      this.initialResidual = initialResidual;
      this.iterations = iterations;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.converged = converged;
      this.trivial = trivial;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase phase, int componentCount, String failureMessage) {
      return new Result(phase, phase, new double[componentCount], Double.NaN, Double.POSITIVE_INFINITY,
          Double.POSITIVE_INFINITY, 0, Double.NaN, false, false, failureMessage);
    }

    public boolean isConverged() {
      return converged && failureMessage == null;
    }

    public boolean isTrivial() {
      return trivial;
    }

    public CandidatePhase getSeedPhase() {
      return seedPhase;
    }

    public CandidatePhase getPhysicalPhase() {
      return physicalPhase;
    }

    public double[] getComposition() {
      return composition.clone();
    }

    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    public double getStationarityResidual() {
      return stationarityResidual;
    }

    public double getInitialResidual() {
      return initialResidual;
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
