package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import Jama.Matrix;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Solves a critical endpoint where the two retained phases coalesce while the third phase remains distinct and
 * incipient.
 *
 * <p>
 * This topology is different from an incipient-with-retained critical endpoint. At retained-pair coalescence the
 * specified two-phase flash is singular and its phase fraction is indeterminate, so the criticality equations are
 * evaluated on one homogeneous phase at the fixed overall composition. Temperature and pressure satisfy the minimum
 * homogeneous TPD-Hessian eigenvalue and third directional derivative equations. The requested third phase is then
 * independently required to be a distinct zero-TPD stationary point at the same state.
 * </p>
 */
public final class HydrocarbonWaterRetainedPairCriticalEndpointSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;

  private final SystemInterface template;
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase incipientPhase;
  private final CandidatePhase criticalPhase;
  private int maximumIterations = 24;
  private double criticalityTolerance = 1.0e-7;
  private double finiteDifferenceStep = 1.0e-3;
  private double curvatureFiniteDifferenceStep = 2.0e-4;
  private double stationarityTolerance = 1.0e-7;
  private double minimumDistinctPhaseDistance = 1.0e-4;

  /** Creates a homogeneous retained-pair critical endpoint solver for one explicit topology. */
  public HydrocarbonWaterRetainedPairCriticalEndpointSolver(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase, CandidatePhase criticalPhase) {
    if (template == null || retainedPhaseZero == null || retainedPhaseOne == null || incipientPhase == null
        || criticalPhase == null) {
      throw new IllegalArgumentException("template and all phase identities are required");
    }
    if (retainedPhaseZero == retainedPhaseOne || retainedPhaseZero == incipientPhase
        || retainedPhaseOne == incipientPhase) {
      throw new IllegalArgumentException("the retained and incipient phase families must be distinct");
    }
    if (criticalPhase != retainedPhaseZero && criticalPhase != retainedPhaseOne) {
      throw new IllegalArgumentException("critical phase must select one retained EOS root");
    }
    this.template = template.clone();
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.incipientPhase = incipientPhase;
    this.criticalPhase = criticalPhase;
  }

  /** Sets Newton, curvature, third-phase stationarity, and distinctness gates. */
  public HydrocarbonWaterRetainedPairCriticalEndpointSolver setNumericalControls(int maximumIterations,
      double criticalityTolerance, double finiteDifferenceStep, double curvatureFiniteDifferenceStep,
      double stationarityTolerance, double minimumDistinctPhaseDistance) {
    if (maximumIterations < 1 || !positive(criticalityTolerance) || !positive(finiteDifferenceStep)
        || !positive(curvatureFiniteDifferenceStep) || !positive(stationarityTolerance)
        || !positive(minimumDistinctPhaseDistance)) {
      throw new IllegalArgumentException("invalid retained-pair critical endpoint controls");
    }
    this.maximumIterations = maximumIterations;
    this.criticalityTolerance = criticalityTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    this.curvatureFiniteDifferenceStep = curvatureFiniteDifferenceStep;
    this.stationarityTolerance = stationarityTolerance;
    this.minimumDistinctPhaseDistance = minimumDistinctPhaseDistance;
    return this;
  }

  /** Solves from a corrected boundary root whose retained phases are approaching coalescence. */
  public Result solve(TwoToThreePhaseBoundaryPointSolver.Result initial) {
    validateInitial(initial);
    double[] variables = new double[] {Math.log(initial.getTemperatureK()), Math.log(initial.getPressureBara())};
    double[] lowerBounds = new double[] {Math.log(Math.max(50.0, 0.8 * initial.getTemperatureK())),
        Math.log(Math.max(1.0e-6, 0.5 * initial.getPressureBara()))};
    double[] upperBounds = new double[] {Math.log(Math.min(2500.0, 1.2 * initial.getTemperatureK())),
        Math.log(Math.min(1.0e6, 1.5 * initial.getPressureBara()))};
    Evaluation evaluation;
    try {
      evaluation = evaluate(variables);
    } catch (RuntimeException error) {
      return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase, criticalPhase, error.getMessage());
    }
    double initialResidual = evaluation.maximumResidual;
    int iterations = 0;
    double conditionNumber = Double.NaN;
    String failureMessage = null;
    while (iterations < maximumIterations && evaluation.maximumResidual > criticalityTolerance) {
      iterations++;
      Matrix jacobian;
      try {
        jacobian = numericalJacobian(variables, evaluation, lowerBounds, upperBounds);
        conditionNumber = jacobian.cond();
      } catch (RuntimeException error) {
        failureMessage = "retained-pair critical Jacobian evaluation failed: " + error.getMessage();
        break;
      }
      Matrix correction;
      try {
        correction = jacobian.solve(new Matrix(evaluation.residual, 2));
      } catch (RuntimeException error) {
        failureMessage = "retained-pair critical Jacobian solve failed: " + error.getMessage();
        break;
      }
      if (!Double.isFinite(correction.normInf())) {
        failureMessage = "retained-pair critical Newton correction is non-finite";
        break;
      }
      limitCorrection(correction);
      boolean accepted = false;
      double damping = 1.0;
      for (int lineSearch = 0; lineSearch < 18; lineSearch++) {
        double[] trialVariables = variables.clone();
        for (int index = 0; index < 2; index++) {
          trialVariables[index] = Math.max(lowerBounds[index],
              Math.min(upperBounds[index], variables[index] - damping * correction.get(index, 0)));
        }
        try {
          Evaluation trial = evaluate(trialVariables);
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
        failureMessage = "line search could not reduce the retained-pair criticality residual";
        break;
      }
    }

    boolean criticalityConverged = evaluation.maximumResidual <= criticalityTolerance;
    if (failureMessage == null && !criticalityConverged) {
      failureMessage = "maximum retained-pair critical endpoint iteration count reached";
    }
    double temperatureK = Math.exp(variables[0]);
    double pressureBara = Math.exp(variables[1]);
    IncipientPhaseStationaryPointSolver.Result incipient = new IncipientPhaseStationaryPointSolver(
        homogeneousSystem(temperatureK, pressureBara), incipientPhase)
        .setNumericalControls(100, stationarityTolerance, 2.0e-5).solve(initial.getIncipientComposition());
    double phaseDistance = incipient.isConverged()
        ? compositionDistance(overallComposition(), incipient.getComposition())
        : Double.POSITIVE_INFINITY;
    boolean gradientConverged = evaluation.criticality.getMaximumGradient() <= 10.0 * criticalityTolerance;
    boolean incipientAtZeroTpd = incipient.isConverged() && !incipient.isTrivial()
        && incipient.getPhysicalPhase() == incipientPhase
        && Math.abs(incipient.getTangentPlaneDistance()) <= stationarityTolerance
        && incipient.getStationarityResidual() <= stationarityTolerance
        && phaseDistance >= minimumDistinctPhaseDistance;
    boolean physicalEndpoint = criticalityConverged && gradientConverged && incipientAtZeroTpd;
    if (failureMessage == null && criticalityConverged && !gradientConverged) {
      failureMessage = "criticality equations converged but the homogeneous TPD gradient did not";
    }
    if (failureMessage == null && criticalityConverged && !incipientAtZeroTpd) {
      failureMessage = "retained-pair critical point does not coincide with a distinct zero-TPD third phase";
    }
    return new Result(retainedPhaseZero, retainedPhaseOne, incipientPhase, criticalPhase, temperatureK, pressureBara,
        evaluation.criticality, incipient, phaseDistance, initialResidual, iterations, conditionNumber,
        criticalityConverged, gradientConverged, incipientAtZeroTpd, physicalEndpoint, failureMessage);
  }

  private Matrix numericalJacobian(double[] variables, Evaluation base, double[] lowerBounds, double[] upperBounds) {
    Matrix jacobian = new Matrix(2, 2);
    for (int column = 0; column < 2; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(variables[column]));
      double[] plus = variables.clone();
      double[] minus = variables.clone();
      plus[column] = Math.min(upperBounds[column], plus[column] + step);
      minus[column] = Math.max(lowerBounds[column], minus[column] - step);
      double denominator = plus[column] - minus[column];
      Evaluation plusEvaluation = evaluate(plus);
      Evaluation minusEvaluation;
      try {
        minusEvaluation = evaluate(minus);
      } catch (RuntimeException error) {
        minusEvaluation = base;
        denominator = plus[column] - variables[column];
      }
      if (!(denominator > 0.0)) {
        throw new IllegalStateException("finite-difference variable is pinned at a local bound");
      }
      for (int row = 0; row < 2; row++) {
        jacobian.set(row, column, (plusEvaluation.residual[row] - minusEvaluation.residual[row]) / denominator);
      }
    }
    return jacobian;
  }

  private Evaluation evaluate(double[] variables) {
    double temperatureK = Math.exp(variables[0]);
    double pressureBara = Math.exp(variables[1]);
    double[] overall = overallComposition();
    TwoToThreePhaseArcLengthCorrector.State homogeneous = TwoToThreePhaseArcLengthCorrector.State.create(
        retainedPhaseZero, retainedPhaseOne, incipientPhase, temperatureK, pressureBara, 0.5, overall, overall,
        overall);
    IncipientPhaseCurvatureAnalyzer.Result criticality = new IncipientPhaseCurvatureAnalyzer(template, criticalPhase,
        criticalPhase).setFiniteDifferenceStep(curvatureFiniteDifferenceStep).analyze(homogeneous);
    double[] residual = new double[] {criticality.getMinimumEigenvalue(), criticality.getThirdDirectionalDerivative()};
    double maximumResidual = Math.max(Math.abs(residual[0]), Math.abs(residual[1]));
    if (!Double.isFinite(maximumResidual)) {
      throw new IllegalStateException("retained-pair criticality residual is non-finite");
    }
    return new Evaluation(criticality, residual, maximumResidual);
  }

  private SystemInterface homogeneousSystem(double temperatureK, double pressureBara) {
    SystemInterface working = template.clone();
    working.setMultiPhaseCheck(false);
    working.setNumberOfPhases(1);
    PhaseInterface phaseTemplate = template.getPhase(0);
    working.setPhase(phaseTemplate.clone(), 0);
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    working.setBeta(0, 1.0);
    working.setPhaseType(0, toPhaseType(criticalPhase));
    double[] overall = overallComposition();
    for (int componentIndex = 0; componentIndex < overall.length; componentIndex++) {
      working.getPhase(0).getComponent(componentIndex).setx(overall[componentIndex]);
    }
    working.getPhase(0).normalize();
    working.setPhaseType(0, toPhaseType(criticalPhase));
    working.init(1, 0);
    return working;
  }

  private double[] overallComposition() {
    double[] composition = new double[template.getPhase(0).getNumberOfComponents()];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(),
          MINIMUM_COMPOSITION);
      total += composition[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }

  private void validateInitial(TwoToThreePhaseBoundaryPointSolver.Result initial) {
    if (initial == null || !initial.isConverged() || initial.getRetainedPhaseZero() != retainedPhaseZero
        || initial.getRetainedPhaseOne() != retainedPhaseOne || initial.getIncipientPhase() != incipientPhase) {
      throw new IllegalArgumentException("initial boundary root does not match the retained-pair endpoint topology");
    }
  }

  private static void limitCorrection(Matrix correction) {
    double scale = 1.0;
    double[] maximum = new double[] {0.08, 0.15};
    for (int index = 0; index < 2; index++) {
      double magnitude = Math.abs(correction.get(index, 0));
      if (magnitude > maximum[index]) {
        scale = Math.min(scale, maximum[index] / magnitude);
      }
    }
    correction.timesEquals(scale);
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
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

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  private static final class Evaluation {
    private final IncipientPhaseCurvatureAnalyzer.Result criticality;
    private final double[] residual;
    private final double maximumResidual;

    private Evaluation(IncipientPhaseCurvatureAnalyzer.Result criticality, double[] residual, double maximumResidual) {
      this.criticality = criticality;
      this.residual = residual;
      this.maximumResidual = maximumResidual;
    }
  }

  /** Immutable retained-pair criticality and independent third-phase evidence. */
  public static final class Result {
    private final CandidatePhase retainedPhaseZero;
    private final CandidatePhase retainedPhaseOne;
    private final CandidatePhase incipientPhase;
    private final CandidatePhase criticalPhase;
    private final double temperatureK;
    private final double pressureBara;
    private final IncipientPhaseCurvatureAnalyzer.Result criticality;
    private final IncipientPhaseStationaryPointSolver.Result incipient;
    private final double incipientPhaseDistance;
    private final double initialCriticalityResidual;
    private final int iterations;
    private final double jacobianConditionNumber;
    private final boolean criticalityConverged;
    private final boolean gradientConverged;
    private final boolean incipientAtZeroTpd;
    private final boolean physicalEndpoint;
    private final String failureMessage;

    private Result(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase,
        CandidatePhase criticalPhase, double temperatureK, double pressureBara,
        IncipientPhaseCurvatureAnalyzer.Result criticality, IncipientPhaseStationaryPointSolver.Result incipient,
        double incipientPhaseDistance, double initialCriticalityResidual, int iterations,
        double jacobianConditionNumber, boolean criticalityConverged, boolean gradientConverged,
        boolean incipientAtZeroTpd, boolean physicalEndpoint, String failureMessage) {
      this.retainedPhaseZero = retainedPhaseZero;
      this.retainedPhaseOne = retainedPhaseOne;
      this.incipientPhase = incipientPhase;
      this.criticalPhase = criticalPhase;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.criticality = criticality;
      this.incipient = incipient;
      this.incipientPhaseDistance = incipientPhaseDistance;
      this.initialCriticalityResidual = initialCriticalityResidual;
      this.iterations = iterations;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.criticalityConverged = criticalityConverged;
      this.gradientConverged = gradientConverged;
      this.incipientAtZeroTpd = incipientAtZeroTpd;
      this.physicalEndpoint = physicalEndpoint;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase incipientPhase, CandidatePhase criticalPhase, String failureMessage) {
      return new Result(retainedPhaseZero, retainedPhaseOne, incipientPhase, criticalPhase, Double.NaN, Double.NaN,
          null, null, Double.POSITIVE_INFINITY, Double.NaN, 0, Double.NaN, false, false, false, false, failureMessage);
    }

    public boolean isPhysicalEndpoint() {
      return physicalEndpoint && failureMessage == null;
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

    public CandidatePhase getCriticalPhase() {
      return criticalPhase;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double getMinimumEigenvalue() {
      return criticality == null ? Double.NaN : criticality.getMinimumEigenvalue();
    }

    public double getThirdDirectionalDerivative() {
      return criticality == null ? Double.NaN : criticality.getThirdDirectionalDerivative();
    }

    public double getHomogeneousTpdGradient() {
      return criticality == null ? Double.NaN : criticality.getMaximumGradient();
    }

    public double getIncipientTangentPlaneDistance() {
      return incipient == null ? Double.NaN : incipient.getTangentPlaneDistance();
    }

    public double getIncipientStationarityResidual() {
      return incipient == null ? Double.NaN : incipient.getStationarityResidual();
    }

    public double getIncipientPhaseDistance() {
      return incipientPhaseDistance;
    }

    public double getInitialCriticalityResidual() {
      return initialCriticalityResidual;
    }

    public int getIterations() {
      return iterations;
    }

    public double getJacobianConditionNumber() {
      return jacobianConditionNumber;
    }

    public boolean isCriticalityConverged() {
      return criticalityConverged;
    }

    public boolean isGradientConverged() {
      return gradientConverged;
    }

    public boolean isIncipientAtZeroTpd() {
      return incipientAtZeroTpd;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
