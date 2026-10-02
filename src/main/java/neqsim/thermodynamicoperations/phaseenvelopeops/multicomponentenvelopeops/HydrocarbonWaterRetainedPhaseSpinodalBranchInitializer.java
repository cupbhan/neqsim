package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Builds the first distinct three-retained-phase state born at a retained-phase spinodal. */
public final class HydrocarbonWaterRetainedPhaseSpinodalBranchInitializer {
  private final SystemInterface template;
  private final CandidatePhase bifurcatingPhase;
  private double pressureOffsetBara = 5.0;
  private int boundaryMaximumIterations = 100;
  private double boundaryResidualTolerance = 1.0e-8;
  private double boundaryFiniteDifferenceStep = 2.0e-5;
  private int flashMaximumIterations = 200;
  private double flashResidualTolerance = 1.0e-9;
  private double flashFiniteDifferenceStep = 2.0e-5;
  private double minimumDuplicatePhaseDistance = 1.0e-4;

  /** Creates an initializer for the retained phase that loses local composition stability. */
  public HydrocarbonWaterRetainedPhaseSpinodalBranchInitializer(SystemInterface template,
      CandidatePhase bifurcatingPhase) {
    if (template == null || bifurcatingPhase == null) {
      throw new IllegalArgumentException("template and bifurcating phase are required");
    }
    this.template = template.clone();
    this.bifurcatingPhase = bifurcatingPhase;
  }

  /** Sets the offset into the unstable side and the boundary/flash numerical controls. */
  public HydrocarbonWaterRetainedPhaseSpinodalBranchInitializer setNumericalControls(double pressureOffsetBara,
      int boundaryMaximumIterations, double boundaryResidualTolerance, double boundaryFiniteDifferenceStep,
      int flashMaximumIterations, double flashResidualTolerance, double flashFiniteDifferenceStep,
      double minimumDuplicatePhaseDistance) {
    if (!positive(pressureOffsetBara) || boundaryMaximumIterations < 1 || !positive(boundaryResidualTolerance)
        || !positive(boundaryFiniteDifferenceStep) || flashMaximumIterations < 1 || !positive(flashResidualTolerance)
        || !positive(flashFiniteDifferenceStep) || !positive(minimumDuplicatePhaseDistance)) {
      throw new IllegalArgumentException("invalid retained-phase branch initialization controls");
    }
    this.pressureOffsetBara = pressureOffsetBara;
    this.boundaryMaximumIterations = boundaryMaximumIterations;
    this.boundaryResidualTolerance = boundaryResidualTolerance;
    this.boundaryFiniteDifferenceStep = boundaryFiniteDifferenceStep;
    this.flashMaximumIterations = flashMaximumIterations;
    this.flashResidualTolerance = flashResidualTolerance;
    this.flashFiniteDifferenceStep = flashFiniteDifferenceStep;
    this.minimumDuplicatePhaseDistance = minimumDuplicatePhaseDistance;
    return this;
  }

  /**
   * Moves a short distance into the unstable side, perturbs both signs of the null mode, and corrects a distinct
   * branch.
   */
  public Result initialize(HydrocarbonWaterRetainedPhaseSpinodalSolver.Result spinodal) {
    if (spinodal == null || !spinodal.isConverged()) {
      throw new IllegalArgumentException("a converged retained-phase spinodal is required");
    }
    HydrocarbonWaterRetainedPhaseSpinodalSolver.CurvatureState negativeSide = negativeSide(spinodal);
    HydrocarbonWaterRetainedPhaseSpinodalSolver.CurvatureState positiveSide = negativeSide == spinodal
        .getLowerPressure() ? spinodal.getUpperPressure() : spinodal.getLowerPressure();
    if (negativeSide == null || !Double.isFinite(negativeSide.getMinimumEigenvalue())
        || negativeSide.getMinimumEigenvalue() >= 0.0) {
      return Result.failure(Collections.<Attempt>emptyList(), "spinodal does not expose a negative-eigenvalue side");
    }
    TwoToThreePhaseArcLengthCorrector.State spinodalState = spinodal.getSpinodal().getState();
    if (bifurcatingPhase != spinodalState.getRetainedPhaseZero()
        && bifurcatingPhase != spinodalState.getRetainedPhaseOne()) {
      throw new IllegalArgumentException("configured bifurcating phase is not retained at the supplied spinodal");
    }
    double pressureDirection = Math.copySign(1.0,
        negativeSide.getState().getPressureBara() - positiveSide.getState().getPressureBara());
    double targetPressureBara = spinodalState.getPressureBara() + pressureDirection * pressureOffsetBara;
    if (!(targetPressureBara > 1.0e-6)) {
      return Result.failure(Collections.<Attempt>emptyList(), "unstable-side pressure target is non-physical");
    }
    TwoToThreePhaseArcLengthCorrector boundaryCorrector = new TwoToThreePhaseArcLengthCorrector(template,
        spinodalState.getRetainedPhaseZero(), spinodalState.getRetainedPhaseOne(), spinodalState.getIncipientPhase())
        .setNumericalControls(boundaryMaximumIterations, boundaryResidualTolerance, boundaryFiniteDifferenceStep);
    TwoToThreePhaseArcLengthCorrector.Result correctedBoundary = boundaryCorrector
        .correctAtPressure(negativeSide.getState(), targetPressureBara);
    if (!correctedBoundary.isConverged()) {
      return Result.failure(Collections.<Attempt>emptyList(),
          "unstable-side boundary correction failed: " + correctedBoundary.getFailureMessage());
    }
    TwoToThreePhaseArcLengthCorrector.State boundaryState = correctedBoundary.getState();
    double[] bifurcatingComposition = retainedComposition(boundaryState, bifurcatingPhase);
    TwoToThreePhaseArcLengthCorrector.State homogeneous = TwoToThreePhaseArcLengthCorrector.State.create(
        boundaryState.getRetainedPhaseZero(), boundaryState.getRetainedPhaseOne(), bifurcatingPhase,
        boundaryState.getTemperatureK(), boundaryState.getPressureBara(), boundaryState.getBeta(),
        boundaryState.getPhaseZeroComposition(), boundaryState.getPhaseOneComposition(), bifurcatingComposition);
    IncipientPhaseStationarityJacobianAnalyzer.Result mode = new IncipientPhaseStationarityJacobianAnalyzer(template,
        bifurcatingPhase, bifurcatingPhase).setFiniteDifferenceStep(2.0e-4).analyze(homogeneous);
    if (!Double.isFinite(mode.getBifurcationEigenvalue()) || mode.getBifurcationEigenvalue() >= 0.0) {
      return Result.failure(Collections.<Attempt>emptyList(),
          "target boundary did not remain on the negative retained-phase eigenvalue side");
    }

    CandidatePhase otherPhase = bifurcatingPhase == boundaryState.getRetainedPhaseZero()
        ? boundaryState.getRetainedPhaseOne()
        : boundaryState.getRetainedPhaseZero();
    double[] otherComposition = retainedComposition(boundaryState, otherPhase);
    double bifurcatingFraction = retainedFraction(boundaryState, bifurcatingPhase);
    double otherFraction = 1.0 - bifurcatingFraction;
    List<Attempt> attempts = new ArrayList<Attempt>();
    SpecifiedThreePhaseFlashSolver.Result accepted = null;
    Candidate aqueousCandidate = null;
    for (double split : new double[] {0.05, 0.01, 0.10}) {
      for (double amplitude : new double[] {0.20, 0.10, 0.40}) {
        for (double sign : new double[] {-1.0, 1.0}) {
          double firstAmplitude = -sign * amplitude * split / (1.0 - split);
          double secondAmplitude = sign * amplitude;
          double[] firstComposition = perturbAlongMode(bifurcatingComposition, mode, firstAmplitude);
          double[] secondComposition = perturbAlongMode(bifurcatingComposition, mode, secondAmplitude);
          double[] fractions = new double[] {otherFraction, bifurcatingFraction * (1.0 - split),
              bifurcatingFraction * split};
          SpecifiedThreePhaseGibbsSeedPreconditioner.Result preconditioned = new SpecifiedThreePhaseGibbsSeedPreconditioner(
              template, otherPhase, bifurcatingPhase, bifurcatingPhase).setNumericalControls(4, 1.0e-10)
              .precondition(boundaryState.getTemperatureK(), boundaryState.getPressureBara(), fractions,
                  otherComposition, firstComposition, secondComposition);
          SpecifiedThreePhaseFlashSolver.Result flash = null;
          if (preconditioned.isConverged()) {
            flash = new SpecifiedThreePhaseFlashSolver(template, otherPhase, bifurcatingPhase, bifurcatingPhase)
                .setNumericalControls(flashMaximumIterations, flashResidualTolerance, flashFiniteDifferenceStep)
                .setCoincidentPhaseDeflation(1, 2, 1.0, 1.0).solve(boundaryState.getTemperatureK(),
                    boundaryState.getPressureBara(), preconditioned.getPhaseFractions(),
                    preconditioned.getPhaseComposition(0), preconditioned.getPhaseComposition(1),
                    preconditioned.getPhaseComposition(2));
          }
          boolean usable = usable(flash);
          attempts.add(new Attempt(split, firstAmplitude, secondAmplitude, preconditioned, flash, usable));
          if (!usable) {
            continue;
          }
          accepted = flash;
          try {
            SystemInterface threePhaseSystem = new SpecifiedThreePhaseFlashSolver(template, otherPhase,
                bifurcatingPhase, bifurcatingPhase).toThermodynamicSystem(flash);
            aqueousCandidate = new IncipientPhaseStabilityAnalyzer(threePhaseSystem).setMaximumIterations(500)
                .setDampingFactor(0.2)
                .analyzeCandidate(boundaryState.getIncipientPhase(), boundaryState.getIncipientComposition());
          } catch (RuntimeException error) {
            aqueousCandidate = null;
          }
          break;
        }
        if (accepted != null) {
          break;
        }
      }
      if (accepted != null) {
        break;
      }
    }
    if (accepted == null) {
      return new Result(boundaryState, mode, otherPhase, attempts, null, null, false,
          "no distinct three-retained-phase branch seed converged");
    }
    return new Result(boundaryState, mode, otherPhase, attempts, accepted, aqueousCandidate, true, null);
  }

  private boolean usable(SpecifiedThreePhaseFlashSolver.Result flash) {
    if (flash == null || !flash.isConverged() || flash.getPhaseCompositionDistance(1, 2) < minimumDuplicatePhaseDistance
        || flash.getThermodynamicMaximumResidual() > 10.0 * flashResidualTolerance
        || flash.getMaterialBalanceResidual() > 10.0 * flashResidualTolerance
        || !matchesConfiguredFamily(flash.getPhaseComposition(1))
        || !matchesConfiguredFamily(flash.getPhaseComposition(2))) {
      return false;
    }
    for (double fraction : flash.getPhaseFractions()) {
      if (!Double.isFinite(fraction) || fraction <= 1.0e-12) {
        return false;
      }
    }
    return true;
  }

  private boolean matchesConfiguredFamily(double[] composition) {
    if (!template.getPhase(0).hasComponent("water")) {
      return true;
    }
    int waterIndex = template.getPhase(0).getComponent("water").getComponentNumber();
    double waterMoleFraction = composition[waterIndex];
    if (bifurcatingPhase == CandidatePhase.OIL) {
      return waterMoleFraction < 0.5;
    }
    if (bifurcatingPhase == CandidatePhase.AQUEOUS) {
      return waterMoleFraction >= 0.5;
    }
    return true;
  }

  private static HydrocarbonWaterRetainedPhaseSpinodalSolver.CurvatureState negativeSide(
      HydrocarbonWaterRetainedPhaseSpinodalSolver.Result spinodal) {
    if (spinodal.getLowerPressure().getMinimumEigenvalue() < 0.0) {
      return spinodal.getLowerPressure();
    }
    if (spinodal.getUpperPressure().getMinimumEigenvalue() < 0.0) {
      return spinodal.getUpperPressure();
    }
    return null;
  }

  private static double retainedFraction(TwoToThreePhaseArcLengthCorrector.State state, CandidatePhase phase) {
    if (phase == state.getRetainedPhaseZero()) {
      return state.getBeta();
    }
    if (phase == state.getRetainedPhaseOne()) {
      return 1.0 - state.getBeta();
    }
    throw new IllegalArgumentException("phase is not retained in the supplied boundary state");
  }

  private static double[] retainedComposition(TwoToThreePhaseArcLengthCorrector.State state, CandidatePhase phase) {
    if (phase == state.getRetainedPhaseZero()) {
      return state.getPhaseZeroComposition();
    }
    if (phase == state.getRetainedPhaseOne()) {
      return state.getPhaseOneComposition();
    }
    throw new IllegalArgumentException("phase is not retained in the supplied boundary state");
  }

  private static double[] perturbAlongMode(double[] base, IncipientPhaseStationarityJacobianAnalyzer.Result mode,
      double amplitude) {
    int referenceIndex = mode.getReferenceComponentIndex();
    double[] eigenvector = mode.getBifurcationEigenvector();
    double[] logarithms = new double[base.length];
    double reference = Math.max(base[referenceIndex], 1.0e-100);
    int coordinateIndex = 0;
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < base.length; componentIndex++) {
      if (componentIndex != referenceIndex) {
        logarithms[componentIndex] = Math.log(Math.max(base[componentIndex], 1.0e-100) / reference)
            + amplitude * eigenvector[coordinateIndex++];
      }
      maximum = Math.max(maximum, logarithms[componentIndex]);
    }
    double[] composition = new double[base.length];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = Math.exp(Math.max(-700.0, logarithms[componentIndex] - maximum));
      total += composition[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** One signed null-mode seed attempt and its two numerical stages. */
  public static final class Attempt {
    private final double splitFraction;
    private final double firstModeAmplitude;
    private final double secondModeAmplitude;
    private final SpecifiedThreePhaseGibbsSeedPreconditioner.Result preconditioned;
    private final SpecifiedThreePhaseFlashSolver.Result flash;
    private final boolean usable;

    private Attempt(double splitFraction, double firstModeAmplitude, double secondModeAmplitude,
        SpecifiedThreePhaseGibbsSeedPreconditioner.Result preconditioned, SpecifiedThreePhaseFlashSolver.Result flash,
        boolean usable) {
      this.splitFraction = splitFraction;
      this.firstModeAmplitude = firstModeAmplitude;
      this.secondModeAmplitude = secondModeAmplitude;
      this.preconditioned = preconditioned;
      this.flash = flash;
      this.usable = usable;
    }

    public double getSplitFraction() {
      return splitFraction;
    }

    public double getFirstModeAmplitude() {
      return firstModeAmplitude;
    }

    public double getSecondModeAmplitude() {
      return secondModeAmplitude;
    }

    public SpecifiedThreePhaseGibbsSeedPreconditioner.Result getPreconditioned() {
      return preconditioned;
    }

    public SpecifiedThreePhaseFlashSolver.Result getFlash() {
      return flash;
    }

    public boolean isUsable() {
      return usable;
    }
  }

  /** Immutable first post-spinodal branch state and diagnostics. */
  public static final class Result {
    private final TwoToThreePhaseArcLengthCorrector.State oldBoundaryState;
    private final IncipientPhaseStationarityJacobianAnalyzer.Result bifurcationMode;
    private final CandidatePhase otherPhase;
    private final List<Attempt> attempts;
    private final SpecifiedThreePhaseFlashSolver.Result branchState;
    private final Candidate formerIncipientPhaseCandidate;
    private final boolean initialized;
    private final String failureMessage;

    private Result(TwoToThreePhaseArcLengthCorrector.State oldBoundaryState,
        IncipientPhaseStationarityJacobianAnalyzer.Result bifurcationMode, CandidatePhase otherPhase,
        List<Attempt> attempts, SpecifiedThreePhaseFlashSolver.Result branchState,
        Candidate formerIncipientPhaseCandidate, boolean initialized, String failureMessage) {
      this.oldBoundaryState = oldBoundaryState;
      this.bifurcationMode = bifurcationMode;
      this.otherPhase = otherPhase;
      this.attempts = Collections.unmodifiableList(new ArrayList<Attempt>(attempts));
      this.branchState = branchState;
      this.formerIncipientPhaseCandidate = formerIncipientPhaseCandidate;
      this.initialized = initialized;
      this.failureMessage = failureMessage;
    }

    private static Result failure(List<Attempt> attempts, String failureMessage) {
      return new Result(null, null, null, attempts, null, null, false, failureMessage);
    }

    public boolean isInitialized() {
      return initialized && failureMessage == null && branchState != null && branchState.isConverged();
    }

    public TwoToThreePhaseArcLengthCorrector.State getOldBoundaryState() {
      return oldBoundaryState;
    }

    public IncipientPhaseStationarityJacobianAnalyzer.Result getBifurcationMode() {
      return bifurcationMode;
    }

    public CandidatePhase getOtherPhase() {
      return otherPhase;
    }

    public List<Attempt> getAttempts() {
      return attempts;
    }

    public SpecifiedThreePhaseFlashSolver.Result getBranchState() {
      return branchState;
    }

    public Candidate getFormerIncipientPhaseCandidate() {
      return formerIncipientPhaseCandidate;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
