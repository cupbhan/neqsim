package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.Arrays;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Local bifurcation analysis of stationary-phase chemical-potential-difference equations.
 *
 * <p>
 * The reduced residual subtracts one component equation from every other stationary equation in log-ratio composition
 * coordinates. A phase compared with itself therefore has an exactly zero residual even when a mixing-rule fugacity
 * implementation is not a numerically exact derivative of a scalar TPD function. A zero eigenvalue of the symmetric
 * reduced Jacobian identifies loss of local composition stability and birth of a secondary stationary branch. The
 * antisymmetry norm is retained as a thermodynamic-consistency diagnostic.
 * </p>
 */
public final class IncipientPhaseStationarityJacobianAnalyzer {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;

  private final SystemInterface template;
  private final CandidatePhase referencePhase;
  private final CandidatePhase trialPhase;
  private double finiteDifferenceStep = 2.0e-4;

  /** Creates a reduced stationary-equation Jacobian analyzer. */
  public IncipientPhaseStationarityJacobianAnalyzer(SystemInterface template, CandidatePhase referencePhase,
      CandidatePhase trialPhase) {
    if (template == null || referencePhase == null || trialPhase == null) {
      throw new IllegalArgumentException("template and reference/trial phase families are required");
    }
    this.template = template.clone();
    this.referencePhase = referencePhase;
    this.trialPhase = trialPhase;
  }

  /** Sets the central finite-difference step in log-ratio composition coordinates. */
  public IncipientPhaseStationarityJacobianAnalyzer setFiniteDifferenceStep(double finiteDifferenceStep) {
    if (!Double.isFinite(finiteDifferenceStep) || finiteDifferenceStep <= 0.0) {
      throw new IllegalArgumentException("finiteDifferenceStep must be positive");
    }
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /** Evaluates the reduced residual and symmetric stationary Jacobian at one boundary state. */
  public Result analyze(TwoToThreePhaseArcLengthCorrector.State state) {
    if (state == null || !Double.isFinite(state.getTemperatureK()) || !Double.isFinite(state.getPressureBara())) {
      throw new IllegalArgumentException("a finite corrected boundary state is required");
    }
    double[] referenceComposition = referenceComposition(state);
    double[] trialComposition = state.getIncipientComposition();
    if (referenceComposition.length != trialComposition.length || trialComposition.length < 2) {
      throw new IllegalArgumentException("state compositions do not match the configured fluid");
    }
    int referenceComponentIndex = largestComponentIndex(trialComposition);
    double[] coordinates = logRatioCoordinates(trialComposition, referenceComponentIndex);
    EvaluationContext context = createContext(state.getTemperatureK(), state.getPressureBara(), referenceComposition);
    double[] residual = stationaryResidual(context, coordinates, referenceComponentIndex);
    int dimension = coordinates.length;
    Matrix jacobian = new Matrix(dimension, dimension);
    for (int column = 0; column < dimension; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(coordinates[column]));
      double[] plus = coordinates.clone();
      double[] minus = coordinates.clone();
      plus[column] += step;
      minus[column] -= step;
      double[] plusResidual = stationaryResidual(context, plus, referenceComponentIndex);
      double[] minusResidual = stationaryResidual(context, minus, referenceComponentIndex);
      for (int row = 0; row < dimension; row++) {
        jacobian.set(row, column, (plusResidual[row] - minusResidual[row]) / (2.0 * step));
      }
    }
    Matrix symmetric = jacobian.plus(jacobian.transpose()).times(0.5);
    EigenvalueDecomposition decomposition = symmetric.eig();
    double[] rawEigenvalues = decomposition.getRealEigenvalues();
    Integer[] order = new Integer[rawEigenvalues.length];
    for (int index = 0; index < order.length; index++) {
      order[index] = index;
    }
    Arrays.sort(order, (first, second) -> Double.compare(rawEigenvalues[first], rawEigenvalues[second]));
    double[] eigenvalues = new double[rawEigenvalues.length];
    double[] minimumEigenvector = new double[dimension];
    Matrix eigenvectors = decomposition.getV();
    for (int index = 0; index < order.length; index++) {
      eigenvalues[index] = rawEigenvalues[order[index]];
    }
    if (dimension > 0) {
      double norm = 0.0;
      for (int row = 0; row < dimension; row++) {
        minimumEigenvector[row] = eigenvectors.get(row, order[0]);
        norm += minimumEigenvector[row] * minimumEigenvector[row];
      }
      norm = Math.sqrt(norm);
      for (int row = 0; row < dimension; row++) {
        minimumEigenvector[row] /= norm;
      }
    }
    double maximumResidual = maximumAbsolute(residual);
    double maximumAntisymmetry = 0.0;
    for (int row = 0; row < dimension; row++) {
      for (int column = 0; column < row; column++) {
        maximumAntisymmetry = Math.max(maximumAntisymmetry,
            Math.abs(jacobian.get(row, column) - jacobian.get(column, row)));
      }
    }
    EigenvalueDecomposition rawDecomposition = jacobian.eig();
    double[] rawRealEigenvalues = rawDecomposition.getRealEigenvalues();
    double[] rawImaginaryEigenvalues = rawDecomposition.getImagEigenvalues();
    int closestEigenvalueIndex = closestRealEigenvalueIndex(rawRealEigenvalues, rawImaginaryEigenvalues);
    double closestMagnitude = Double.POSITIVE_INFINITY;
    if (closestEigenvalueIndex >= 0) {
      closestMagnitude = Math.abs(rawRealEigenvalues[closestEigenvalueIndex]);
    }
    double bifurcationEigenvalue = closestEigenvalueIndex < 0 ? Double.NaN : rawRealEigenvalues[closestEigenvalueIndex];
    double[] bifurcationEigenvector = new double[dimension];
    if (Double.isFinite(bifurcationEigenvalue) && dimension > 0) {
      Matrix rawEigenvectors = rawDecomposition.getV();
      double norm = 0.0;
      for (int row = 0; row < dimension; row++) {
        bifurcationEigenvector[row] = rawEigenvectors.get(row, closestEigenvalueIndex);
        norm += bifurcationEigenvector[row] * bifurcationEigenvector[row];
      }
      norm = Math.sqrt(norm);
      if (norm > 0.0 && Double.isFinite(norm)) {
        for (int row = 0; row < dimension; row++) {
          bifurcationEigenvector[row] /= norm;
        }
      }
    }
    double[] singularValues = jacobian.svd().getSingularValues();
    double minimumSingularValue = singularValues.length == 0 ? Double.NaN : singularValues[singularValues.length - 1];
    return new Result(residual, jacobian.getArrayCopy(), symmetric.getArrayCopy(), eigenvalues, minimumEigenvector,
        rawRealEigenvalues, rawImaginaryEigenvalues, bifurcationEigenvalue, bifurcationEigenvector,
        minimumSingularValue, jacobian.det(), maximumResidual, maximumAntisymmetry, referenceComponentIndex,
        finiteDifferenceStep, allFinite(residual) && allFinite(jacobian.getArray()) && allFinite(rawRealEigenvalues)
            && allFinite(rawImaginaryEigenvalues),
        closestMagnitude);
  }

  /**
   * Selects the real stationary mode nearest zero.
   *
   * <p>
   * A real, non-symmetric stationarity Jacobian may have a complex-conjugate pair whose magnitude is smaller than every
   * real eigenvalue. Such a pair is not a static composition bifurcation and must not hide an available real
   * zero-crossing mode. The previous implementation selected the smallest complex magnitude first and then returned
   * {@code NaN}; for an odd-dimensional 22-component reduction this incorrectly rejected an otherwise finite
   * retained-oil spectrum.
   * </p>
   */
  static int closestRealEigenvalueIndex(double[] realEigenvalues, double[] imaginaryEigenvalues) {
    if (realEigenvalues == null || imaginaryEigenvalues == null
        || realEigenvalues.length != imaginaryEigenvalues.length) {
      throw new IllegalArgumentException("matching real and imaginary eigenvalue arrays are required");
    }
    int selected = -1;
    double selectedMagnitude = Double.POSITIVE_INFINITY;
    for (int index = 0; index < realEigenvalues.length; index++) {
      double real = realEigenvalues[index];
      double imaginary = imaginaryEigenvalues[index];
      if (!Double.isFinite(real) || !Double.isFinite(imaginary)
          || Math.abs(imaginary) > 1.0e-7 * Math.max(1.0, Math.abs(real))) {
        continue;
      }
      double magnitude = Math.abs(real);
      if (magnitude < selectedMagnitude) {
        selected = index;
        selectedMagnitude = magnitude;
      }
    }
    return selected;
  }

  private EvaluationContext createContext(double temperatureK, double pressureBara, double[] referenceComposition) {
    PhaseInterface phaseTemplate = template.getPhase(0);
    SystemInterface reference = template.clone();
    reference.setMultiPhaseCheck(false);
    reference.setNumberOfPhases(1);
    reference.setPhase(phaseTemplate.clone(), 0);
    reference.setTemperature(temperatureK);
    reference.setPressure(pressureBara);
    reference.setBeta(0, 1.0);
    setPhase(reference, 0, referencePhase, referenceComposition);
    reference.init(1, 0);
    double[] logReferenceFugacity = new double[referenceComposition.length];
    for (int componentIndex = 0; componentIndex < referenceComposition.length; componentIndex++) {
      logReferenceFugacity[componentIndex] = Math
          .log(Math.max(reference.getPhase(0).getComponent(componentIndex).getx(), MINIMUM_COMPOSITION))
          + reference.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
    }

    SystemInterface trial = template.clone();
    trial.setMultiPhaseCheck(false);
    trial.setNumberOfPhases(1);
    trial.setPhase(phaseTemplate.clone(), 0);
    trial.setTemperature(temperatureK);
    trial.setPressure(pressureBara);
    trial.setBeta(0, 1.0);
    setPhase(trial, 0, trialPhase, referenceComposition);
    trial.init(1, 0);
    return new EvaluationContext(trial, logReferenceFugacity, toPhaseType(trialPhase));
  }

  private static double[] stationaryResidual(EvaluationContext context, double[] coordinates,
      int referenceComponentIndex) {
    double[] composition = compositionFromLogRatios(coordinates, referenceComponentIndex,
        context.system.getPhase(0).getNumberOfComponents());
    setPhase(context.system, 0, null, composition);
    context.system.setPhaseType(0, context.trialPhaseType);
    context.system.init(1, 0);
    double referenceEquation = Math.log(Math.max(composition[referenceComponentIndex], MINIMUM_COMPOSITION))
        + context.system.getPhase(0).getComponent(referenceComponentIndex).getLogFugacityCoefficient()
        - context.logReferenceFugacity[referenceComponentIndex];
    double[] residual = new double[composition.length - 1];
    int residualIndex = 0;
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      if (componentIndex == referenceComponentIndex) {
        continue;
      }
      residual[residualIndex++] = Math.log(Math.max(composition[componentIndex], MINIMUM_COMPOSITION))
          + context.system.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient()
          - context.logReferenceFugacity[componentIndex] - referenceEquation;
    }
    return residual;
  }

  private double[] referenceComposition(TwoToThreePhaseArcLengthCorrector.State state) {
    if (referencePhase == state.getRetainedPhaseZero()) {
      return state.getPhaseZeroComposition();
    }
    if (referencePhase == state.getRetainedPhaseOne()) {
      return state.getPhaseOneComposition();
    }
    throw new IllegalArgumentException("reference phase is not one of the supplied state's retained phases");
  }

  private static void setPhase(SystemInterface system, int phaseIndex, CandidatePhase phase, double[] composition) {
    if (phase != null) {
      system.setPhaseType(phaseIndex, toPhaseType(phase));
    }
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      system.getPhase(phaseIndex).getComponent(componentIndex).setx(composition[componentIndex]);
    }
    system.getPhase(phaseIndex).normalize();
    if (phase != null) {
      system.setPhaseType(phaseIndex, toPhaseType(phase));
    }
  }

  private static int largestComponentIndex(double[] composition) {
    int index = 0;
    for (int componentIndex = 1; componentIndex < composition.length; componentIndex++) {
      if (composition[componentIndex] > composition[index]) {
        index = componentIndex;
      }
    }
    return index;
  }

  private static double[] logRatioCoordinates(double[] composition, int referenceComponentIndex) {
    double[] coordinates = new double[composition.length - 1];
    double reference = Math.max(composition[referenceComponentIndex], MINIMUM_COMPOSITION);
    int coordinateIndex = 0;
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      if (componentIndex != referenceComponentIndex) {
        coordinates[coordinateIndex++] = Math
            .log(Math.max(composition[componentIndex], MINIMUM_COMPOSITION) / reference);
      }
    }
    return coordinates;
  }

  private static double[] compositionFromLogRatios(double[] coordinates, int referenceComponentIndex,
      int componentCount) {
    double[] logarithms = new double[componentCount];
    int coordinateIndex = 0;
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      logarithms[componentIndex] = componentIndex == referenceComponentIndex ? 0.0 : coordinates[coordinateIndex++];
      maximum = Math.max(maximum, logarithms[componentIndex]);
    }
    double total = 0.0;
    double[] composition = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = Math.exp(Math.max(-700.0, logarithms[componentIndex] - maximum));
      total += composition[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }

  private static double maximumAbsolute(double[] values) {
    double maximum = 0.0;
    for (double value : values) {
      maximum = Math.max(maximum, Math.abs(value));
    }
    return maximum;
  }

  private static boolean allFinite(double[] values) {
    for (double value : values) {
      if (!Double.isFinite(value)) {
        return false;
      }
    }
    return true;
  }

  private static boolean allFinite(double[][] values) {
    for (double[] row : values) {
      if (!allFinite(row)) {
        return false;
      }
    }
    return true;
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
    private final SystemInterface system;
    private final double[] logReferenceFugacity;
    private final PhaseType trialPhaseType;

    private EvaluationContext(SystemInterface system, double[] logReferenceFugacity, PhaseType trialPhaseType) {
      this.system = system;
      this.logReferenceFugacity = logReferenceFugacity;
      this.trialPhaseType = trialPhaseType;
    }
  }

  /** Immutable stationary residual, Jacobian, symmetric spectrum, and consistency diagnostics. */
  public static final class Result {
    private final double[] residual;
    private final double[][] jacobian;
    private final double[][] symmetricJacobian;
    private final double[] eigenvalues;
    private final double[] minimumEigenvector;
    private final double[] rawRealEigenvalues;
    private final double[] rawImaginaryEigenvalues;
    private final double bifurcationEigenvalue;
    private final double[] bifurcationEigenvector;
    private final double minimumSingularValue;
    private final double determinant;
    private final double maximumResidual;
    private final double maximumAntisymmetry;
    private final int referenceComponentIndex;
    private final double finiteDifferenceStep;
    private final boolean finiteEvaluation;
    private final double closestRealEigenvalueMagnitude;

    private Result(double[] residual, double[][] jacobian, double[][] symmetricJacobian, double[] eigenvalues,
        double[] minimumEigenvector, double[] rawRealEigenvalues, double[] rawImaginaryEigenvalues,
        double bifurcationEigenvalue, double[] bifurcationEigenvector, double minimumSingularValue, double determinant,
        double maximumResidual, double maximumAntisymmetry, int referenceComponentIndex, double finiteDifferenceStep,
        boolean finiteEvaluation, double closestRealEigenvalueMagnitude) {
      this.residual = residual.clone();
      this.jacobian = copy(jacobian);
      this.symmetricJacobian = copy(symmetricJacobian);
      this.eigenvalues = eigenvalues.clone();
      this.minimumEigenvector = minimumEigenvector.clone();
      this.rawRealEigenvalues = rawRealEigenvalues.clone();
      this.rawImaginaryEigenvalues = rawImaginaryEigenvalues.clone();
      this.bifurcationEigenvalue = bifurcationEigenvalue;
      this.bifurcationEigenvector = bifurcationEigenvector.clone();
      this.minimumSingularValue = minimumSingularValue;
      this.determinant = determinant;
      this.maximumResidual = maximumResidual;
      this.maximumAntisymmetry = maximumAntisymmetry;
      this.referenceComponentIndex = referenceComponentIndex;
      this.finiteDifferenceStep = finiteDifferenceStep;
      this.finiteEvaluation = finiteEvaluation;
      this.closestRealEigenvalueMagnitude = closestRealEigenvalueMagnitude;
    }

    public double[] getResidual() {
      return residual.clone();
    }

    public double[][] getJacobian() {
      return copy(jacobian);
    }

    public double[][] getSymmetricJacobian() {
      return copy(symmetricJacobian);
    }

    public double[] getEigenvalues() {
      return eigenvalues.clone();
    }

    public double getMinimumEigenvalue() {
      return eigenvalues.length == 0 ? Double.NaN : eigenvalues[0];
    }

    public double[] getMinimumEigenvector() {
      return minimumEigenvector.clone();
    }

    public double[] getRawRealEigenvalues() {
      return rawRealEigenvalues.clone();
    }

    public double[] getRawImaginaryEigenvalues() {
      return rawImaginaryEigenvalues.clone();
    }

    /** @return real raw-Jacobian eigenvalue closest to zero, or NaN when the closest pair is complex */
    public double getBifurcationEigenvalue() {
      return bifurcationEigenvalue;
    }

    /**
     * Returns the normalized right eigenvector associated with {@link #getBifurcationEigenvalue()}.
     *
     * <p>
     * Its coordinates are log composition ratios in component order with the reference component omitted. Both signs
     * describe the same local bifurcation direction and must be tried when constructing a new phase branch.
     * </p>
     */
    public double[] getBifurcationEigenvector() {
      return bifurcationEigenvector.clone();
    }

    public double getMinimumSingularValue() {
      return minimumSingularValue;
    }

    public double getDeterminant() {
      return determinant;
    }

    public double getMaximumResidual() {
      return maximumResidual;
    }

    public double getMaximumAntisymmetry() {
      return maximumAntisymmetry;
    }

    public int getReferenceComponentIndex() {
      return referenceComponentIndex;
    }

    public double getFiniteDifferenceStep() {
      return finiteDifferenceStep;
    }

    /** @return true when residual, finite-difference Jacobian, and raw spectrum are all finite */
    public boolean isFiniteEvaluation() {
      return finiteEvaluation;
    }

    /** @return true when the finite raw spectrum contains a real static-bifurcation mode */
    public boolean hasRealBifurcationMode() {
      return finiteEvaluation && Double.isFinite(bifurcationEigenvalue);
    }

    /** @return magnitude of the selected real eigenvalue, or positive infinity when none exists */
    public double getClosestRealEigenvalueMagnitude() {
      return closestRealEigenvalueMagnitude;
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
