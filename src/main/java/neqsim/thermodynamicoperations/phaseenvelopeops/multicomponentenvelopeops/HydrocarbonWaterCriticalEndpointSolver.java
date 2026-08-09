package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import Jama.Matrix;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Solves a two-phase critical endpoint from retained-phase equilibrium plus Michelsen criticality conditions.
 *
 * <p>
 * At every trial T/P the retained two phases are fully re-flashed. The phase selected by {@code criticalPhase} then
 * satisfies two independent mixture-criticality equations: the smallest homogeneous TPD-Hessian eigenvalue is zero and
 * the third directional derivative along its null eigenvector is zero. A converged mathematical candidate is only
 * promoted to a physical critical endpoint when an independently solved incipient-phase stationary point coalesces with
 * the critical retained phase.
 * </p>
 */
public final class HydrocarbonWaterCriticalEndpointSolver {
  private final SystemInterface template;
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase incipientPhase;
  private final CandidatePhase criticalPhase;
  private int maximumIterations = 20;
  private double criticalityTolerance = 1.0e-7;
  private double flashResidualTolerance = 1.0e-9;
  private double finiteDifferenceStep = 1.0e-3;
  private double curvatureFiniteDifferenceStep = 2.0e-4;
  private double stationarityTolerance = 1.0e-7;
  private double phaseCoalescenceTolerance = 1.0e-4;

  /** Creates a critical endpoint solver for one explicit two-to-three-phase boundary topology. */
  public HydrocarbonWaterCriticalEndpointSolver(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase, CandidatePhase criticalPhase) {
    if (template == null || retainedPhaseZero == null || retainedPhaseOne == null || incipientPhase == null
        || criticalPhase == null) {
      throw new IllegalArgumentException("template and all phase identities are required");
    }
    if (retainedPhaseZero == retainedPhaseOne || incipientPhase == retainedPhaseOne
        || incipientPhase == retainedPhaseZero) {
      throw new IllegalArgumentException("the retained and incipient phase families must be distinct");
    }
    if (criticalPhase != retainedPhaseZero && criticalPhase != retainedPhaseOne) {
      throw new IllegalArgumentException("critical phase must be one of the retained phases");
    }
    this.template = template.clone();
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.incipientPhase = incipientPhase;
    this.criticalPhase = criticalPhase;
  }

  /** Sets outer Newton, flash, numerical derivative, stationarity, and coalescence tolerances. */
  public HydrocarbonWaterCriticalEndpointSolver setNumericalControls(int maximumIterations, double criticalityTolerance,
      double flashResidualTolerance, double finiteDifferenceStep, double curvatureFiniteDifferenceStep,
      double stationarityTolerance, double phaseCoalescenceTolerance) {
    if (maximumIterations < 1 || !positive(criticalityTolerance) || !positive(flashResidualTolerance)
        || !positive(finiteDifferenceStep) || !positive(curvatureFiniteDifferenceStep)
        || !positive(stationarityTolerance) || !positive(phaseCoalescenceTolerance)) {
      throw new IllegalArgumentException("invalid critical endpoint numerical controls");
    }
    this.maximumIterations = maximumIterations;
    this.criticalityTolerance = criticalityTolerance;
    this.flashResidualTolerance = flashResidualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    this.curvatureFiniteDifferenceStep = curvatureFiniteDifferenceStep;
    this.stationarityTolerance = stationarityTolerance;
    this.phaseCoalescenceTolerance = phaseCoalescenceTolerance;
    return this;
  }

  /** Solves from a nearby corrected two-to-three-phase boundary state. */
  public Result solve(TwoToThreePhaseArcLengthCorrector.State initial) {
    validateInitial(initial);
    SpecifiedTwoPhaseFlashSolver flashSolver = new SpecifiedTwoPhaseFlashSolver(template, retainedPhaseZero,
        retainedPhaseOne).setNumericalControls(80, flashResidualTolerance, 2.0e-5);
    double[] variables = new double[] { Math.log(initial.getTemperatureK()), Math.log(initial.getPressureBara()) };
    FlashSeed seed = new FlashSeed(initial.getBeta(), initial.getPhaseZeroComposition(),
        initial.getPhaseOneComposition());
    Evaluation evaluation;
    try {
      evaluation = evaluate(flashSolver, variables, seed);
    } catch (RuntimeException error) {
      return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase, criticalPhase, error.getMessage());
    }
    double initialResidual = evaluation.maximumCriticalityResidual;
    int iterations = 0;
    double conditionNumber = Double.NaN;
    String failureMessage = null;
    while (iterations < maximumIterations && evaluation.maximumCriticalityResidual > criticalityTolerance) {
      iterations++;
      Matrix jacobian;
      try {
        jacobian = numericalJacobian(flashSolver, variables, evaluation, seed);
        conditionNumber = jacobian.cond();
      } catch (RuntimeException error) {
        failureMessage = "critical endpoint Jacobian evaluation failed: " + error.getMessage();
        break;
      }
      Matrix correction;
      try {
        correction = jacobian.solve(new Matrix(evaluation.residual, 2));
      } catch (RuntimeException error) {
        failureMessage = "critical endpoint Jacobian solve failed: " + error.getMessage();
        break;
      }
      if (!Double.isFinite(correction.normInf())) {
        failureMessage = "critical endpoint Newton correction is non-finite";
        break;
      }
      limitCorrection(correction);
      boolean accepted = false;
      double damping = 1.0;
      for (int lineSearch = 0; lineSearch < 18; lineSearch++) {
        double[] trialVariables = new double[] { variables[0] - damping * correction.get(0, 0),
            variables[1] - damping * correction.get(1, 0) };
        clamp(trialVariables);
        try {
          Evaluation trial = evaluate(flashSolver, trialVariables, new FlashSeed(evaluation.flash.getBeta(),
              evaluation.flash.getPhaseZeroComposition(), evaluation.flash.getPhaseOneComposition()));
          if (trial.maximumCriticalityResidual < evaluation.maximumCriticalityResidual) {
            variables = trialVariables;
            evaluation = trial;
            seed = new FlashSeed(trial.flash.getBeta(), trial.flash.getPhaseZeroComposition(),
                trial.flash.getPhaseOneComposition());
            accepted = true;
            break;
          }
        } catch (RuntimeException error) {
          // Retry with a shorter line-search step.
        }
        damping *= 0.5;
      }
      if (!accepted) {
        failureMessage = "line search could not reduce the criticality residual";
        break;
      }
    }
    boolean criticalityConverged = evaluation.maximumCriticalityResidual <= criticalityTolerance;
    if (!criticalityConverged && failureMessage == null) {
      failureMessage = "maximum critical endpoint iteration count reached";
    }

    IncipientPhaseStationaryPointSolver.Result coalescence = new IncipientPhaseStationaryPointSolver(
        flashSolver.toThermodynamicSystem(evaluation.flash), incipientPhase)
        .setNumericalControls(80, stationarityTolerance, 2.0e-5).solve(criticalComposition(evaluation.flash));
    double phaseDistance = coalescence.isConverged()
        ? compositionDistance(criticalComposition(evaluation.flash), coalescence.getComposition())
        : Double.POSITIVE_INFINITY;
    boolean coalesced = coalescence.isConverged() && phaseDistance <= phaseCoalescenceTolerance
        && Math.abs(coalescence.getTangentPlaneDistance()) <= stationarityTolerance;
    boolean gradientConverged = evaluation.criticality.getMaximumGradient() <= 10.0 * criticalityTolerance;
    boolean physicalEndpoint = criticalityConverged && gradientConverged && coalesced;
    if (failureMessage == null && criticalityConverged && !gradientConverged) {
      failureMessage = "criticality equations converged but the homogeneous TPD gradient did not";
    }
    if (failureMessage == null && criticalityConverged && !coalesced) {
      failureMessage = "criticality equations converged but the incipient phase did not coalesce";
    }
    return new Result(retainedPhaseZero, retainedPhaseOne, incipientPhase, criticalPhase, evaluation.flash,
        evaluation.criticality, coalescence, phaseDistance, initialResidual, iterations, conditionNumber,
        criticalityConverged, gradientConverged, coalesced, physicalEndpoint, failureMessage);
  }

  private Matrix numericalJacobian(SpecifiedTwoPhaseFlashSolver flashSolver, double[] variables, Evaluation base,
      FlashSeed seed) {
    Matrix jacobian = new Matrix(2, 2);
    FlashSeed warm = new FlashSeed(base.flash.getBeta(), base.flash.getPhaseZeroComposition(),
        base.flash.getPhaseOneComposition());
    for (int column = 0; column < 2; column++) {
      double step = finiteDifferenceStep * Math.max(1.0, Math.abs(variables[column]));
      double[] plus = variables.clone();
      double[] minus = variables.clone();
      plus[column] += step;
      minus[column] -= step;
      clamp(plus);
      clamp(minus);
      Evaluation plusEvaluation = evaluate(flashSolver, plus, warm);
      Evaluation minusEvaluation;
      double denominator = plus[column] - minus[column];
      try {
        minusEvaluation = evaluate(flashSolver, minus, warm);
      } catch (RuntimeException error) {
        minusEvaluation = base;
        denominator = plus[column] - variables[column];
      }
      for (int row = 0; row < 2; row++) {
        jacobian.set(row, column, (plusEvaluation.residual[row] - minusEvaluation.residual[row]) / denominator);
      }
    }
    return jacobian;
  }

  private Evaluation evaluate(SpecifiedTwoPhaseFlashSolver flashSolver, double[] variables, FlashSeed seed) {
    double temperatureK = Math.exp(variables[0]);
    double pressureBara = Math.exp(variables[1]);
    SpecifiedTwoPhaseFlashSolver.Result flash = flashSolver.solve(temperatureK, pressureBara, seed.beta,
        seed.phaseZeroComposition, seed.phaseOneComposition);
    if (!flash.isConverged()) {
      throw new IllegalStateException("specified retained two-phase flash failed: " + flash.getFailureMessage());
    }
    double[] criticalComposition = criticalComposition(flash);
    TwoToThreePhaseArcLengthCorrector.State criticalState = TwoToThreePhaseArcLengthCorrector.State.create(
        retainedPhaseZero, retainedPhaseOne, incipientPhase, temperatureK, pressureBara, flash.getBeta(),
        flash.getPhaseZeroComposition(), flash.getPhaseOneComposition(), criticalComposition);
    IncipientPhaseCurvatureAnalyzer.Result criticality = new IncipientPhaseCurvatureAnalyzer(template, criticalPhase,
        criticalPhase).setFiniteDifferenceStep(curvatureFiniteDifferenceStep).analyze(criticalState);
    double[] residual = new double[] { criticality.getMinimumEigenvalue(),
        criticality.getThirdDirectionalDerivative() };
    double maximumResidual = Math.max(Math.abs(residual[0]), Math.abs(residual[1]));
    if (!Double.isFinite(maximumResidual)) {
      throw new IllegalStateException("criticality residual is non-finite");
    }
    return new Evaluation(flash, criticality, residual, maximumResidual);
  }

  private double[] criticalComposition(SpecifiedTwoPhaseFlashSolver.Result flash) {
    return criticalPhase == retainedPhaseZero ? flash.getPhaseZeroComposition() : flash.getPhaseOneComposition();
  }

  private void validateInitial(TwoToThreePhaseArcLengthCorrector.State initial) {
    if (initial == null || initial.getRetainedPhaseZero() != retainedPhaseZero
        || initial.getRetainedPhaseOne() != retainedPhaseOne || initial.getIncipientPhase() != incipientPhase) {
      throw new IllegalArgumentException("initial state does not match the critical endpoint topology");
    }
  }

  private static void limitCorrection(Matrix correction) {
    double scale = 1.0;
    double[] maximum = new double[] { 0.08, 0.15 };
    for (int index = 0; index < 2; index++) {
      if (Math.abs(correction.get(index, 0)) > maximum[index]) {
        scale = Math.min(scale, maximum[index] / Math.abs(correction.get(index, 0)));
      }
    }
    correction.timesEquals(scale);
  }

  private static void clamp(double[] variables) {
    variables[0] = Math.max(Math.log(50.0), Math.min(Math.log(2500.0), variables[0]));
    variables[1] = Math.max(Math.log(1.0e-6), Math.min(Math.log(1.0e6), variables[1]));
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  private static final class FlashSeed {
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;

    private FlashSeed(double beta, double[] phaseZeroComposition, double[] phaseOneComposition) {
      this.beta = beta;
      this.phaseZeroComposition = phaseZeroComposition.clone();
      this.phaseOneComposition = phaseOneComposition.clone();
    }
  }

  private static final class Evaluation {
    private final SpecifiedTwoPhaseFlashSolver.Result flash;
    private final IncipientPhaseCurvatureAnalyzer.Result criticality;
    private final double[] residual;
    private final double maximumCriticalityResidual;

    private Evaluation(SpecifiedTwoPhaseFlashSolver.Result flash, IncipientPhaseCurvatureAnalyzer.Result criticality,
        double[] residual, double maximumCriticalityResidual) {
      this.flash = flash;
      this.criticality = criticality;
      this.residual = residual;
      this.maximumCriticalityResidual = maximumCriticalityResidual;
    }
  }

  /** Immutable critical-equation, stationary-point, and phase-coalescence evidence. */
  public static final class Result {
    private final CandidatePhase retainedPhaseZero;
    private final CandidatePhase retainedPhaseOne;
    private final CandidatePhase incipientPhase;
    private final CandidatePhase criticalPhase;
    private final SpecifiedTwoPhaseFlashSolver.Result flash;
    private final IncipientPhaseCurvatureAnalyzer.Result criticality;
    private final IncipientPhaseStationaryPointSolver.Result coalescence;
    private final double phaseDistance;
    private final double initialCriticalityResidual;
    private final int iterations;
    private final double jacobianConditionNumber;
    private final boolean criticalityConverged;
    private final boolean gradientConverged;
    private final boolean coalesced;
    private final boolean physicalEndpoint;
    private final String failureMessage;

    private Result(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase,
        CandidatePhase criticalPhase, SpecifiedTwoPhaseFlashSolver.Result flash,
        IncipientPhaseCurvatureAnalyzer.Result criticality, IncipientPhaseStationaryPointSolver.Result coalescence,
        double phaseDistance, double initialCriticalityResidual, int iterations, double jacobianConditionNumber,
        boolean criticalityConverged, boolean gradientConverged, boolean coalesced, boolean physicalEndpoint,
        String failureMessage) {
      this.retainedPhaseZero = retainedPhaseZero;
      this.retainedPhaseOne = retainedPhaseOne;
      this.incipientPhase = incipientPhase;
      this.criticalPhase = criticalPhase;
      this.flash = flash;
      this.criticality = criticality;
      this.coalescence = coalescence;
      this.phaseDistance = phaseDistance;
      this.initialCriticalityResidual = initialCriticalityResidual;
      this.iterations = iterations;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.criticalityConverged = criticalityConverged;
      this.gradientConverged = gradientConverged;
      this.coalesced = coalesced;
      this.physicalEndpoint = physicalEndpoint;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase incipientPhase, CandidatePhase criticalPhase, String failureMessage) {
      return new Result(retainedPhaseZero, retainedPhaseOne, incipientPhase, criticalPhase, null, null, null,
          Double.POSITIVE_INFINITY, Double.NaN, 0, Double.NaN, false, false, false, false, failureMessage);
    }

    public boolean isPhysicalEndpoint() {
      return physicalEndpoint && failureMessage == null;
    }

    public boolean isCriticalityConverged() {
      return criticalityConverged;
    }

    public boolean isGradientConverged() {
      return gradientConverged;
    }

    public boolean isCoalesced() {
      return coalesced;
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
      return flash == null ? Double.NaN : flash.getTemperatureK();
    }

    public double getPressureBara() {
      return flash == null ? Double.NaN : flash.getPressureBara();
    }

    public double getFlashResidual() {
      return flash == null ? Double.NaN : flash.getMaximumResidual();
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

    public double getPhaseDistance() {
      return phaseDistance;
    }

    public double getIncipientTangentPlaneDistance() {
      return coalescence == null ? Double.NaN : coalescence.getTangentPlaneDistance();
    }

    public double getIncipientStationarityResidual() {
      return coalescence == null ? Double.NaN : coalescence.getStationarityResidual();
    }

    /** Converts a fully accepted endpoint to quality-gate evidence with the independently coalesced composition. */
    public TwoToThreePhaseBoundaryQualityGate.EvidencePoint toEvidencePoint() {
      if (!isPhysicalEndpoint()) {
        throw new IllegalStateException("only a physical critical endpoint can become boundary evidence");
      }
      return new TwoToThreePhaseBoundaryQualityGate.EvidencePoint(flash.getTemperatureK(), flash.getPressureBara(),
          flash.getBeta(), flash.getPhaseZeroComposition(), flash.getPhaseOneComposition(),
          coalescence.getComposition(), flash.getMaximumResidual(), coalescence.getTangentPlaneDistance(),
          coalescence.getStationarityResidual());
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

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
