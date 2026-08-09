package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

/**
 * Solves a nonlinear least-squares problem with variable scaling and an SVD trust region.
 *
 * <p>
 * The local step minimizes {@code ||r + J D q||} subject to {@code ||q|| <= delta}, where {@code D} contains caller
 * supplied variable scales. A Levenberg parameter is selected in singular-value coordinates, so the implementation does
 * not form {@code J'J}. This is important for phase-boundary equations close to critical and multiphase bifurcation
 * points, where a normal-equation solve squares the Jacobian condition number.
 * </p>
 */
public final class ScaledTrustRegionLeastSquaresSolver {
  /** Residual callback for one trial vector. */
  public interface ResidualFunction {
    /**
     * Evaluates the unweighted residual vector.
     *
     * @param variables trial variables
     * @return finite residual vector
     */
    double[] evaluate(double[] variables);
  }

  private final ResidualFunction function;
  private int maximumIterations = 100;
  private double residualTolerance = 1.0e-9;
  private double finiteDifferenceStep = 2.0e-5;
  private double initialTrustRadius = 1.0;
  private double maximumTrustRadius = 100.0;
  private double[] variableScales;
  private double[] lowerBounds;
  private double[] upperBounds;

  /** Creates a trust-region least-squares solver. */
  public ScaledTrustRegionLeastSquaresSolver(ResidualFunction function) {
    if (function == null) {
      throw new IllegalArgumentException("residual function is required");
    }
    this.function = function;
  }

  /** Sets iteration, convergence, finite-difference, and trust-radius controls. */
  public ScaledTrustRegionLeastSquaresSolver setNumericalControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep, double initialTrustRadius, double maximumTrustRadius) {
    if (maximumIterations < 1 || !positive(residualTolerance) || !positive(finiteDifferenceStep)
        || !positive(initialTrustRadius) || !positive(maximumTrustRadius) || initialTrustRadius > maximumTrustRadius) {
      throw new IllegalArgumentException("invalid trust-region numerical controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    this.initialTrustRadius = initialTrustRadius;
    this.maximumTrustRadius = maximumTrustRadius;
    return this;
  }

  /** Sets positive scales used to measure the trust-region step. */
  public ScaledTrustRegionLeastSquaresSolver setVariableScales(double[] variableScales) {
    validatePositiveVector(variableScales, "variable scales");
    this.variableScales = variableScales.clone();
    return this;
  }

  /** Sets finite lower and upper variable bounds. */
  public ScaledTrustRegionLeastSquaresSolver setBounds(double[] lowerBounds, double[] upperBounds) {
    if (lowerBounds == null || upperBounds == null || lowerBounds.length == 0
        || lowerBounds.length != upperBounds.length) {
      throw new IllegalArgumentException("matching variable bounds are required");
    }
    for (int index = 0; index < lowerBounds.length; index++) {
      if (!Double.isFinite(lowerBounds[index]) || !Double.isFinite(upperBounds[index])
          || lowerBounds[index] > upperBounds[index]) {
        throw new IllegalArgumentException("invalid variable bound at index " + index);
      }
    }
    this.lowerBounds = lowerBounds.clone();
    this.upperBounds = upperBounds.clone();
    return this;
  }

  /** Solves from one caller-supplied initial vector. */
  public Result solve(double[] initialVariables) {
    if (initialVariables == null || initialVariables.length == 0) {
      throw new IllegalArgumentException("initial variables are required");
    }
    int variableCount = initialVariables.length;
    double[] scales = scales(variableCount);
    validateConfiguration(variableCount);
    double[] variables = initialVariables.clone();
    clamp(variables);
    double[] residual;
    try {
      residual = checkedResidual(variables, -1);
    } catch (RuntimeException error) {
      return Result.failure(variables, "initial residual evaluation failed: " + error.getMessage());
    }
    double initialMaximumResidual = maximumAbsolute(residual);
    double trustRadius = initialTrustRadius;
    double conditionNumber = Double.NaN;
    int acceptedSteps = 0;
    int rejectedSteps = 0;
    String failureMessage = null;
    int iteration = 0;
    for (; iteration < maximumIterations && maximumAbsolute(residual) > residualTolerance; iteration++) {
      Matrix jacobian;
      try {
        jacobian = numericalJacobian(variables, residual);
      } catch (RuntimeException error) {
        failureMessage = "Jacobian evaluation failed: " + error.getMessage();
        break;
      }
      Matrix scaledJacobian = scaleColumns(jacobian, scales);
      SingularValueDecomposition decomposition;
      try {
        decomposition = scaledJacobian.svd();
        conditionNumber = conditionNumber(decomposition.getSingularValues());
      } catch (RuntimeException error) {
        failureMessage = "SVD failed: " + error.getMessage();
        break;
      }
      double[] scaledStep = trustRegionStep(decomposition, residual, trustRadius);
      if (!(norm(scaledStep) > 1.0e-14)) {
        failureMessage = "trust-region model returned a zero step before residual convergence";
        break;
      }
      double[] trial = variables.clone();
      for (int index = 0; index < variableCount; index++) {
        trial[index] += scales[index] * scaledStep[index];
      }
      clamp(trial);
      double[] actualScaledStep = new double[variableCount];
      for (int index = 0; index < variableCount; index++) {
        actualScaledStep[index] = (trial[index] - variables[index]) / scales[index];
      }
      double stepNorm = norm(actualScaledStep);
      if (!(stepNorm > 1.0e-14)) {
        failureMessage = "active bounds prevented a nonzero trust-region step";
        break;
      }
      double[] trialResidual;
      try {
        trialResidual = checkedResidual(trial, residual.length);
      } catch (RuntimeException error) {
        trialResidual = null;
      }
      double predictedReduction = predictedReduction(residual, scaledJacobian, actualScaledStep);
      double actualReduction = trialResidual == null ? Double.NEGATIVE_INFINITY
          : 0.5 * (squaredNorm(residual) - squaredNorm(trialResidual));
      double ratio = predictedReduction > 0.0 && Double.isFinite(actualReduction) ? actualReduction / predictedReduction
          : Double.NEGATIVE_INFINITY;
      if (ratio < 0.25) {
        trustRadius = Math.max(1.0e-10, 0.25 * trustRadius);
      } else if (ratio > 0.75 && stepNorm >= 0.8 * trustRadius) {
        trustRadius = Math.min(maximumTrustRadius, 2.0 * trustRadius);
      }
      if (ratio > 1.0e-4 && actualReduction > 0.0 && trialResidual != null) {
        variables = trial;
        residual = trialResidual;
        acceptedSteps++;
      } else {
        rejectedSteps++;
      }
      if (trustRadius <= 1.0e-10 && maximumAbsolute(residual) > residualTolerance) {
        failureMessage = "trust-region radius collapsed before residual convergence";
        break;
      }
    }
    boolean converged = maximumAbsolute(residual) <= residualTolerance;
    if (!converged && failureMessage == null) {
      failureMessage = "maximum trust-region iteration count reached";
    }
    return new Result(variables, residual, initialMaximumResidual, maximumAbsolute(residual), iteration, acceptedSteps,
        rejectedSteps, trustRadius, conditionNumber, converged, failureMessage);
  }

  private Matrix numericalJacobian(double[] variables, double[] baseResidual) {
    Matrix jacobian = new Matrix(baseResidual.length, variables.length);
    for (int column = 0; column < variables.length; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(variables[column]));
      double[] plus = variables.clone();
      double[] minus = variables.clone();
      plus[column] += step;
      minus[column] -= step;
      clamp(plus);
      clamp(minus);
      double[] plusResidual = tryResidual(plus, baseResidual.length);
      double[] minusResidual = tryResidual(minus, baseResidual.length);
      double denominator;
      if (plusResidual == null && minusResidual == null) {
        throw new IllegalStateException("both finite-difference evaluations failed at column " + column);
      } else if (plusResidual == null) {
        plusResidual = baseResidual;
        denominator = variables[column] - minus[column];
      } else if (minusResidual == null) {
        minusResidual = baseResidual;
        denominator = plus[column] - variables[column];
      } else {
        denominator = plus[column] - minus[column];
      }
      if (!(Math.abs(denominator) > 0.0)) {
        throw new IllegalStateException("finite-difference bound collapsed at column " + column);
      }
      for (int row = 0; row < baseResidual.length; row++) {
        jacobian.set(row, column, (plusResidual[row] - minusResidual[row]) / denominator);
      }
    }
    return jacobian;
  }

  private double[] tryResidual(double[] variables, int expectedLength) {
    try {
      return checkedResidual(variables, expectedLength);
    } catch (RuntimeException error) {
      return null;
    }
  }

  private double[] checkedResidual(double[] variables, int expectedLength) {
    double[] residual = function.evaluate(variables.clone());
    if (residual == null || residual.length == 0 || expectedLength >= 0 && residual.length != expectedLength) {
      throw new IllegalStateException("residual function changed vector size");
    }
    for (double value : residual) {
      if (!Double.isFinite(value)) {
        throw new IllegalStateException("residual function returned a non-finite value");
      }
    }
    return residual;
  }

  private static Matrix scaleColumns(Matrix jacobian, double[] scales) {
    Matrix result = jacobian.copy();
    for (int column = 0; column < result.getColumnDimension(); column++) {
      for (int row = 0; row < result.getRowDimension(); row++) {
        result.set(row, column, result.get(row, column) * scales[column]);
      }
    }
    return result;
  }

  private static double[] trustRegionStep(SingularValueDecomposition decomposition, double[] residual,
      double trustRadius) {
    Matrix u = decomposition.getU();
    Matrix v = decomposition.getV();
    double[] singularValues = decomposition.getSingularValues();
    int rankDimension = Math.min(singularValues.length, Math.min(u.getColumnDimension(), v.getColumnDimension()));
    double[] projected = new double[rankDimension];
    for (int column = 0; column < rankDimension; column++) {
      for (int row = 0; row < residual.length; row++) {
        projected[column] += u.get(row, column) * residual[row];
      }
    }
    double[] gaussNewton = spectralStep(v, singularValues, projected, rankDimension, 0.0);
    if (norm(gaussNewton) <= trustRadius) {
      return gaussNewton;
    }
    double lower = 0.0;
    double upper = 1.0;
    while (norm(spectralStep(v, singularValues, projected, rankDimension, upper)) > trustRadius && upper < 1.0e30) {
      upper *= 10.0;
    }
    for (int iteration = 0; iteration < 80; iteration++) {
      double middle = 0.5 * (lower + upper);
      double middleNorm = norm(spectralStep(v, singularValues, projected, rankDimension, middle));
      if (middleNorm > trustRadius) {
        lower = middle;
      } else {
        upper = middle;
      }
    }
    return spectralStep(v, singularValues, projected, rankDimension, upper);
  }

  private static double[] spectralStep(Matrix v, double[] singularValues, double[] projected, int rankDimension,
      double damping) {
    double[] step = new double[v.getRowDimension()];
    double largest = rankDimension == 0 ? 0.0 : singularValues[0];
    double threshold = Math.max(1.0e-14, largest * 1.0e-13);
    for (int singularIndex = 0; singularIndex < rankDimension; singularIndex++) {
      double singular = singularValues[singularIndex];
      if (singular <= threshold) {
        continue;
      }
      double coefficient = -singular * projected[singularIndex] / (singular * singular + damping);
      for (int row = 0; row < step.length; row++) {
        step[row] += v.get(row, singularIndex) * coefficient;
      }
    }
    return step;
  }

  private static double predictedReduction(double[] residual, Matrix scaledJacobian, double[] scaledStep) {
    double[] linearResidual = residual.clone();
    for (int row = 0; row < scaledJacobian.getRowDimension(); row++) {
      for (int column = 0; column < scaledJacobian.getColumnDimension(); column++) {
        linearResidual[row] += scaledJacobian.get(row, column) * scaledStep[column];
      }
    }
    return 0.5 * (squaredNorm(residual) - squaredNorm(linearResidual));
  }

  private double[] scales(int variableCount) {
    if (variableScales == null) {
      double[] result = new double[variableCount];
      for (int index = 0; index < result.length; index++) {
        result[index] = 1.0;
      }
      return result;
    }
    return variableScales.clone();
  }

  private void validateConfiguration(int variableCount) {
    if (variableScales != null && variableScales.length != variableCount) {
      throw new IllegalArgumentException("variable scale count does not match initial vector");
    }
    if (lowerBounds != null && lowerBounds.length != variableCount) {
      throw new IllegalArgumentException("variable bound count does not match initial vector");
    }
  }

  private void clamp(double[] variables) {
    if (lowerBounds == null) {
      return;
    }
    for (int index = 0; index < variables.length; index++) {
      variables[index] = Math.max(lowerBounds[index], Math.min(upperBounds[index], variables[index]));
    }
  }

  private static double conditionNumber(double[] singularValues) {
    if (singularValues.length == 0 || !(singularValues[0] > 0.0)) {
      return Double.POSITIVE_INFINITY;
    }
    double smallest = Double.POSITIVE_INFINITY;
    for (double singular : singularValues) {
      if (singular > singularValues[0] * 1.0e-15) {
        smallest = Math.min(smallest, singular);
      }
    }
    return Double.isFinite(smallest) ? singularValues[0] / smallest : Double.POSITIVE_INFINITY;
  }

  private static double maximumAbsolute(double[] vector) {
    double maximum = 0.0;
    for (double value : vector) {
      maximum = Math.max(maximum, Math.abs(value));
    }
    return maximum;
  }

  private static double squaredNorm(double[] vector) {
    double result = 0.0;
    for (double value : vector) {
      result += value * value;
    }
    return result;
  }

  private static double norm(double[] vector) {
    return Math.sqrt(squaredNorm(vector));
  }

  private static void validatePositiveVector(double[] vector, String name) {
    if (vector == null || vector.length == 0) {
      throw new IllegalArgumentException(name + " are required");
    }
    for (double value : vector) {
      if (!positive(value)) {
        throw new IllegalArgumentException(name + " must contain positive finite values");
      }
    }
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Immutable trust-region result and diagnostics. */
  public static final class Result {
    private final double[] variables;
    private final double[] residual;
    private final double initialMaximumResidual;
    private final double maximumResidual;
    private final int iterations;
    private final int acceptedSteps;
    private final int rejectedSteps;
    private final double finalTrustRadius;
    private final double jacobianConditionNumber;
    private final boolean converged;
    private final String failureMessage;

    private Result(double[] variables, double[] residual, double initialMaximumResidual, double maximumResidual,
        int iterations, int acceptedSteps, int rejectedSteps, double finalTrustRadius, double jacobianConditionNumber,
        boolean converged, String failureMessage) {
      this.variables = variables.clone();
      this.residual = residual.clone();
      this.initialMaximumResidual = initialMaximumResidual;
      this.maximumResidual = maximumResidual;
      this.iterations = iterations;
      this.acceptedSteps = acceptedSteps;
      this.rejectedSteps = rejectedSteps;
      this.finalTrustRadius = finalTrustRadius;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.converged = converged;
      this.failureMessage = failureMessage;
    }

    private static Result failure(double[] variables, String failureMessage) {
      return new Result(variables, new double[0], Double.NaN, Double.NaN, 0, 0, 0, Double.NaN, Double.NaN, false,
          failureMessage);
    }

    public boolean isConverged() {
      return converged && failureMessage == null;
    }

    public double[] getVariables() {
      return variables.clone();
    }

    public double[] getResidual() {
      return residual.clone();
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

    public int getAcceptedSteps() {
      return acceptedSteps;
    }

    public int getRejectedSteps() {
      return rejectedSteps;
    }

    public double getFinalTrustRadius() {
      return finalTrustRadius;
    }

    public double getJacobianConditionNumber() {
      return jacobianConditionNumber;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
