package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.Arrays;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Numerical TPD curvature analysis for an incipient phase on a corrected two-to-three-phase boundary. */
public final class IncipientPhaseCurvatureAnalyzer {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;

  private final SystemInterface template;
  private final CandidatePhase referencePhase;
  private final CandidatePhase incipientPhase;
  private double finiteDifferenceStep = 2.0e-4;

  /**
   * Creates a curvature analyzer using one retained phase as the equilibrium chemical-potential reference.
   *
   * <p>
   * The trial phase may equal the reference phase when evaluating homogeneous-mixture criticality. A distinct trial
   * phase is used for ordinary third-phase stationary-point diagnostics.
   * </p>
   */
  public IncipientPhaseCurvatureAnalyzer(SystemInterface template, CandidatePhase referencePhase,
      CandidatePhase incipientPhase) {
    if (template == null || referencePhase == null || incipientPhase == null) {
      throw new IllegalArgumentException("template and reference/trial phase families are required");
    }
    this.template = template.clone();
    this.referencePhase = referencePhase;
    this.incipientPhase = incipientPhase;
  }

  /** Sets the central finite-difference step in log-ratio composition coordinates. */
  public IncipientPhaseCurvatureAnalyzer setFiniteDifferenceStep(double finiteDifferenceStep) {
    if (!Double.isFinite(finiteDifferenceStep) || finiteDifferenceStep <= 0.0) {
      throw new IllegalArgumentException("finiteDifferenceStep must be positive");
    }
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /** Calculates the TPD gradient and Hessian eigenvalues at one corrected boundary state. */
  public Result analyze(TwoToThreePhaseArcLengthCorrector.State state) {
    if (state == null || !Double.isFinite(state.getTemperatureK()) || !Double.isFinite(state.getPressureBara())) {
      throw new IllegalArgumentException("a finite corrected boundary state is required");
    }
    double[] referenceComposition = referenceComposition(state);
    double[] incipientComposition = state.getIncipientComposition();
    if (referenceComposition.length != incipientComposition.length || incipientComposition.length < 2) {
      throw new IllegalArgumentException("state compositions do not match the configured fluid");
    }
    int referenceComponentIndex = largestComponentIndex(incipientComposition);
    double[] coordinates = logRatioCoordinates(incipientComposition, referenceComponentIndex);
    EvaluationContext context = createContext(state.getTemperatureK(), state.getPressureBara(), referenceComposition,
        incipientComposition);
    double baseTpd = tangentPlaneDistance(context, coordinates, referenceComponentIndex);
    int dimension = coordinates.length;
    double[] gradient = new double[dimension];
    double[][] hessian = new double[dimension][dimension];
    for (int row = 0; row < dimension; row++) {
      double rowStep = finiteDifferenceStep * Math.max(1.0, Math.abs(coordinates[row]));
      double[] plus = coordinates.clone();
      double[] minus = coordinates.clone();
      plus[row] += rowStep;
      minus[row] -= rowStep;
      double plusTpd = tangentPlaneDistance(context, plus, referenceComponentIndex);
      double minusTpd = tangentPlaneDistance(context, minus, referenceComponentIndex);
      gradient[row] = (plusTpd - minusTpd) / (2.0 * rowStep);
      hessian[row][row] = (plusTpd - 2.0 * baseTpd + minusTpd) / (rowStep * rowStep);
      for (int column = 0; column < row; column++) {
        double columnStep = finiteDifferenceStep * Math.max(1.0, Math.abs(coordinates[column]));
        double[] plusPlus = coordinates.clone();
        double[] plusMinus = coordinates.clone();
        double[] minusPlus = coordinates.clone();
        double[] minusMinus = coordinates.clone();
        plusPlus[row] += rowStep;
        plusPlus[column] += columnStep;
        plusMinus[row] += rowStep;
        plusMinus[column] -= columnStep;
        minusPlus[row] -= rowStep;
        minusPlus[column] += columnStep;
        minusMinus[row] -= rowStep;
        minusMinus[column] -= columnStep;
        double value = (tangentPlaneDistance(context, plusPlus, referenceComponentIndex)
            - tangentPlaneDistance(context, plusMinus, referenceComponentIndex)
            - tangentPlaneDistance(context, minusPlus, referenceComponentIndex)
            + tangentPlaneDistance(context, minusMinus, referenceComponentIndex)) / (4.0 * rowStep * columnStep);
        hessian[row][column] = value;
        hessian[column][row] = value;
      }
    }
    EigenvalueDecomposition eigen = new Matrix(hessian).eig();
    double[] rawEigenvalues = eigen.getRealEigenvalues();
    Integer[] eigenvalueOrder = new Integer[rawEigenvalues.length];
    for (int index = 0; index < eigenvalueOrder.length; index++) {
      eigenvalueOrder[index] = index;
    }
    Arrays.sort(eigenvalueOrder, (first, second) -> Double.compare(rawEigenvalues[first], rawEigenvalues[second]));
    double[] eigenvalues = new double[rawEigenvalues.length];
    for (int index = 0; index < eigenvalues.length; index++) {
      eigenvalues[index] = rawEigenvalues[eigenvalueOrder[index]];
    }
    double[] minimumEigenvector = new double[dimension];
    Matrix eigenvectors = eigen.getV();
    double eigenvectorNorm = 0.0;
    if (dimension > 0) {
      int minimumColumn = eigenvalueOrder[0];
      for (int row = 0; row < dimension; row++) {
        minimumEigenvector[row] = eigenvectors.get(row, minimumColumn);
        eigenvectorNorm += minimumEigenvector[row] * minimumEigenvector[row];
      }
      eigenvectorNorm = Math.sqrt(eigenvectorNorm);
      for (int row = 0; row < dimension; row++) {
        minimumEigenvector[row] /= eigenvectorNorm;
      }
    }
    double thirdDirectionalDerivative = thirdDirectionalDerivative(context, coordinates, referenceComponentIndex,
        minimumEigenvector, Math.max(1.0e-3, 5.0 * finiteDifferenceStep));
    double maximumGradient = 0.0;
    for (double value : gradient) {
      maximumGradient = Math.max(maximumGradient, Math.abs(value));
    }
    return new Result(baseTpd, gradient, maximumGradient, eigenvalues, minimumEigenvector, thirdDirectionalDerivative,
        referenceComponentIndex, finiteDifferenceStep);
  }

  private static double thirdDirectionalDerivative(EvaluationContext context, double[] coordinates,
      int referenceComponentIndex, double[] direction, double step) {
    if (direction.length == 0) {
      return Double.NaN;
    }
    double[] plusOne = displaced(coordinates, direction, step);
    double[] minusOne = displaced(coordinates, direction, -step);
    double[] plusTwo = displaced(coordinates, direction, 2.0 * step);
    double[] minusTwo = displaced(coordinates, direction, -2.0 * step);
    return (tangentPlaneDistance(context, plusTwo, referenceComponentIndex)
        - 2.0 * tangentPlaneDistance(context, plusOne, referenceComponentIndex)
        + 2.0 * tangentPlaneDistance(context, minusOne, referenceComponentIndex)
        - tangentPlaneDistance(context, minusTwo, referenceComponentIndex)) / (2.0 * step * step * step);
  }

  private static double[] displaced(double[] coordinates, double[] direction, double distance) {
    double[] displaced = coordinates.clone();
    for (int index = 0; index < displaced.length; index++) {
      displaced[index] += distance * direction[index];
    }
    return displaced;
  }

  private EvaluationContext createContext(double temperatureK, double pressureBara, double[] referenceComposition,
      double[] incipientComposition) {
    SystemInterface referenceSystem = template.clone();
    referenceSystem.setMultiPhaseCheck(false);
    referenceSystem.setNumberOfPhases(1);
    PhaseInterface phaseTemplate = template.getPhase(0);
    referenceSystem.setPhase(phaseTemplate.clone(), 0);
    referenceSystem.setTemperature(temperatureK);
    referenceSystem.setPressure(pressureBara);
    referenceSystem.setBeta(0, 1.0);
    setPhase(referenceSystem, 0, referencePhase, referenceComposition);
    referenceSystem.init(1, 0);
    double[] logReferenceFugacity = new double[referenceComposition.length];
    for (int componentIndex = 0; componentIndex < referenceComposition.length; componentIndex++) {
      logReferenceFugacity[componentIndex] = Math
          .log(Math.max(referenceSystem.getPhase(0).getComponent(componentIndex).getx(), MINIMUM_COMPOSITION))
          + referenceSystem.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
    }
    SystemInterface trialSystem = template.clone();
    trialSystem.setMultiPhaseCheck(false);
    trialSystem.setNumberOfPhases(1);
    trialSystem.setPhase(phaseTemplate.clone(), 0);
    trialSystem.setTemperature(temperatureK);
    trialSystem.setPressure(pressureBara);
    trialSystem.setBeta(0, 1.0);
    setPhase(trialSystem, 0, incipientPhase, incipientComposition);
    trialSystem.init(1, 0);
    return new EvaluationContext(trialSystem, logReferenceFugacity);
  }

  private static double tangentPlaneDistance(EvaluationContext context, double[] coordinates,
      int referenceComponentIndex) {
    double[] composition = compositionFromLogRatios(coordinates, referenceComponentIndex,
        context.system.getPhase(0).getNumberOfComponents());
    setPhase(context.system, 0, null, composition);
    context.system.setPhaseType(0, context.incipientPhaseType);
    context.system.init(1, 0);
    double tpd = 0.0;
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      double x = Math.max(composition[componentIndex], MINIMUM_COMPOSITION);
      double logTrialFugacity = Math.log(x)
          + context.system.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
      tpd += x * (logTrialFugacity - context.logReferenceFugacity[componentIndex]);
    }
    return tpd;
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

  private double[] referenceComposition(TwoToThreePhaseArcLengthCorrector.State state) {
    if (referencePhase == state.getRetainedPhaseZero()) {
      return state.getPhaseZeroComposition();
    }
    if (referencePhase == state.getRetainedPhaseOne()) {
      return state.getPhaseOneComposition();
    }
    throw new IllegalArgumentException("reference phase is not one of the supplied state's retained phases");
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
      if (componentIndex == referenceComponentIndex) {
        logarithms[componentIndex] = 0.0;
      } else {
        logarithms[componentIndex] = coordinates[coordinateIndex++];
      }
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
    private final PhaseType incipientPhaseType;

    private EvaluationContext(SystemInterface system, double[] logReferenceFugacity) {
      this.system = system;
      this.logReferenceFugacity = logReferenceFugacity;
      this.incipientPhaseType = system.getPhase(0).getType();
    }
  }

  /** Immutable TPD value, stationarity residual, and sorted Hessian spectrum. */
  public static final class Result {
    private final double tangentPlaneDistance;
    private final double[] gradient;
    private final double maximumGradient;
    private final double[] eigenvalues;
    private final double[] minimumEigenvector;
    private final double thirdDirectionalDerivative;
    private final int referenceComponentIndex;
    private final double finiteDifferenceStep;

    private Result(double tangentPlaneDistance, double[] gradient, double maximumGradient, double[] eigenvalues,
        double[] minimumEigenvector, double thirdDirectionalDerivative, int referenceComponentIndex,
        double finiteDifferenceStep) {
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.gradient = gradient.clone();
      this.maximumGradient = maximumGradient;
      this.eigenvalues = eigenvalues.clone();
      this.minimumEigenvector = minimumEigenvector.clone();
      this.thirdDirectionalDerivative = thirdDirectionalDerivative;
      this.referenceComponentIndex = referenceComponentIndex;
      this.finiteDifferenceStep = finiteDifferenceStep;
    }

    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    public double getMaximumGradient() {
      return maximumGradient;
    }

    public double[] getGradient() {
      return gradient.clone();
    }

    public double[] getEigenvalues() {
      return eigenvalues.clone();
    }

    public double getMinimumEigenvalue() {
      return eigenvalues.length == 0 ? Double.NaN : eigenvalues[0];
    }

    public double getSecondEigenvalue() {
      return eigenvalues.length < 2 ? Double.NaN : eigenvalues[1];
    }

    public double[] getMinimumEigenvector() {
      return minimumEigenvector.clone();
    }

    public double getThirdDirectionalDerivative() {
      return thirdDirectionalDerivative;
    }

    public int getReferenceComponentIndex() {
      return referenceComponentIndex;
    }

    public double getFiniteDifferenceStep() {
      return finiteDifferenceStep;
    }
  }
}
