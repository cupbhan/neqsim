package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.SeedEvaluation;

/** Rejects zero-TPD roots whose retained two-phase state has a deeper non-trivial instability. */
public final class HydrocarbonWaterBoundaryGlobalStabilityGate {
  private final SystemInterface template;
  private double globalTangentPlaneTolerance = 1.0e-7;
  private double targetTangentPlaneTolerance = 1.0e-7;
  private double targetCompositionTolerance = 1.0e-4;
  private double retainedPhaseCurvatureStep = 2.0e-4;
  private double retainedPhaseEigenvalueTolerance = 1.0e-6;
  private double retainedPhaseResidualTolerance = 1.0e-8;

  /** Creates a non-destructive global stability gate for one configured EOS and mixing rule. */
  public HydrocarbonWaterBoundaryGlobalStabilityGate(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets global-negative-TPD, target-zero-TPD, and target-composition tolerances. */
  public HydrocarbonWaterBoundaryGlobalStabilityGate setTolerances(double globalTangentPlaneTolerance,
      double targetTangentPlaneTolerance, double targetCompositionTolerance) {
    if (!positive(globalTangentPlaneTolerance) || !positive(targetTangentPlaneTolerance)
        || !positive(targetCompositionTolerance)) {
      throw new IllegalArgumentException("invalid global boundary stability tolerances");
    }
    this.globalTangentPlaneTolerance = globalTangentPlaneTolerance;
    this.targetTangentPlaneTolerance = targetTangentPlaneTolerance;
    this.targetCompositionTolerance = targetCompositionTolerance;
    return this;
  }

  /** Sets stationary-Jacobian step, minimum-eigenvalue, and homogeneous-residual tolerances for retained phases. */
  public HydrocarbonWaterBoundaryGlobalStabilityGate setRetainedPhaseCurvatureControls(double curvatureStep,
      double eigenvalueTolerance, double residualTolerance) {
    if (!positive(curvatureStep) || !positive(eigenvalueTolerance) || !positive(residualTolerance)) {
      throw new IllegalArgumentException("invalid retained-phase curvature controls");
    }
    this.retainedPhaseCurvatureStep = curvatureStep;
    this.retainedPhaseEigenvalueTolerance = eigenvalueTolerance;
    this.retainedPhaseResidualTolerance = residualTolerance;
    return this;
  }

  /**
   * Re-flashes the retained phases and searches all physical trial families plus the corrected incipient seed.
   *
   * @param root converged specified two-phase plus zero-third-phase-TPD root
   * @return immutable global stability evidence
   */
  public Result evaluate(TwoToThreePhaseBoundaryPointSolver.Result root) {
    return evaluate(root, Collections.<StationarySeed>emptyList());
  }

  /**
   * Re-flashes and evaluates generic trials, the target root, and explicitly continued stationary-branch seeds.
   *
   * @param root converged specified two-phase plus zero-third-phase-TPD root
   * @param stationarySeeds compositions retained from distinct stationary branches at adjacent continuation points
   * @return immutable global stability evidence
   */
  public Result evaluate(TwoToThreePhaseBoundaryPointSolver.Result root, List<StationarySeed> stationarySeeds) {
    if (root == null || !root.isConverged()) {
      throw new IllegalArgumentException("a converged strict boundary root is required");
    }
    if (stationarySeeds == null) {
      throw new IllegalArgumentException("stationarySeeds must not be null");
    }
    // The independent re-flash must be materially tighter than the zero-TPD gate. Otherwise a
    // small retained-phase fugacity mismatch is amplified when the incipient stationary point is
    // reconstructed and a valid boundary root can appear to have a slightly negative TPD.
    SpecifiedTwoPhaseFlashSolver flashSolver = new SpecifiedTwoPhaseFlashSolver(template, root.getRetainedPhaseZero(),
        root.getRetainedPhaseOne()).setNumericalControls(160, 1.0e-11, 2.0e-5);
    SpecifiedTwoPhaseFlashSolver.Result retained = flashSolver.solve(root.getTemperatureK(), root.getPressureBara(),
        root.getBeta(), root.getPhaseZeroComposition(), root.getPhaseOneComposition());
    if (!retained.isConverged()) {
      return Result.failure(root, retained.getFailureMessage());
    }

    SystemInterface reference = flashSolver.toThermodynamicSystem(retained);
    IncipientPhaseStabilityAnalyzer analyzer = new IncipientPhaseStabilityAnalyzer(reference).setMaximumIterations(500)
        .setTolerances(1.0e-9, -globalTangentPlaneTolerance).setDampingFactor(0.2);
    IncipientPhaseStabilityAnalyzer.Result generic = analyzer.analyze();
    List<Candidate> trials = new ArrayList<Candidate>(generic.getTrials());
    SeedEvaluation targetSeedEvaluation = analyzer.evaluateCandidateSeed(root.getIncipientPhase(),
        root.getIncipientComposition());
    SeedEvaluation targetCorrectorSlotSeedEvaluation = analyzer.evaluateCandidateSeed(root.getIncipientPhase(),
        root.getIncipientComposition(), 2);
    TwoToThreePhaseArcLengthCorrector.State continuationState = TwoToThreePhaseArcLengthCorrector.State.create(
        root.getRetainedPhaseZero(), root.getRetainedPhaseOne(), root.getIncipientPhase(), root.getTemperatureK(),
        root.getPressureBara(), root.getBeta(), root.getPhaseZeroComposition(), root.getPhaseOneComposition(),
        root.getIncipientComposition());
    TwoToThreePhaseArcLengthCorrector.EquationResiduals replayedEquationResiduals = new TwoToThreePhaseArcLengthCorrector(
        template, root.getRetainedPhaseZero(), root.getRetainedPhaseOne(), root.getIncipientPhase())
        .replayEquationResiduals(continuationState);
    Candidate target = analyzer.analyzeCandidate(root.getIncipientPhase(), root.getIncipientComposition());
    trials.add(target);
    List<Candidate> seededCandidates = new ArrayList<Candidate>();
    for (StationarySeed seed : stationarySeeds) {
      if (seed == null || seed.getComposition().length != root.getPhaseZeroComposition().length) {
        throw new IllegalArgumentException("every stationary seed must match the boundary component count");
      }
      Candidate seeded = analyzer.analyzeCandidate(seed.getPhase(), seed.getComposition());
      seededCandidates.add(seeded);
      trials.add(seeded);
    }
    List<RetainedPhaseLocalStability> retainedPhaseStability = new ArrayList<RetainedPhaseLocalStability>();
    retainedPhaseStability.add(retainedPhaseStability(root, root.getRetainedPhaseZero()));
    retainedPhaseStability.add(retainedPhaseStability(root, root.getRetainedPhaseOne()));

    double minimumNonTrivialTpd = Double.POSITIVE_INFINITY;
    Candidate mostUnstable = null;
    for (Candidate candidate : trials) {
      if (!candidate.isConverged() || candidate.isTrivial() || !Double.isFinite(candidate.getTangentPlaneDistance())) {
        continue;
      }
      if (candidate.getTangentPlaneDistance() < minimumNonTrivialTpd) {
        minimumNonTrivialTpd = candidate.getTangentPlaneDistance();
        mostUnstable = candidate;
      }
    }
    double targetDistance = target.isConverged()
        ? compositionDistance(target.getComposition(), root.getIncipientComposition())
        : Double.POSITIVE_INFINITY;
    boolean targetMatches = target.isConverged() && !target.isTrivial() && target.getPhase() == root.getIncipientPhase()
        && targetDistance <= targetCompositionTolerance
        && Math.abs(target.getTangentPlaneDistance()) <= targetTangentPlaneTolerance;
    boolean retainedPhasesLocallyStable = true;
    RetainedPhaseLocalStability unstableRetainedPhase = null;
    for (RetainedPhaseLocalStability retainedPhase : retainedPhaseStability) {
      if (!retainedPhase.isAccepted()) {
        retainedPhasesLocallyStable = false;
        if (unstableRetainedPhase == null
            || retainedPhase.getMinimumEigenvalue() < unstableRetainedPhase.getMinimumEigenvalue()) {
          unstableRetainedPhase = retainedPhase;
        }
      }
    }
    boolean globallyStable = targetMatches && minimumNonTrivialTpd >= -globalTangentPlaneTolerance
        && retainedPhasesLocallyStable;
    String failureMessage = null;
    if (!targetMatches) {
      failureMessage = "corrected incipient root was not recovered by the independent stability solve";
    } else if (!retainedPhasesLocallyStable) {
      failureMessage = unstableRetainedPhase.getFailureMessage() == null
          ? "retained " + unstableRetainedPhase.getPhase() + " phase has a negative homogeneous TPD eigenvalue"
          : unstableRetainedPhase.getFailureMessage();
    } else if (!globallyStable) {
      failureMessage = "retained two-phase state has a deeper negative-TPD stationary point";
    }
    return new Result(root, retained, generic, trials, seededCandidates, retainedPhaseStability, targetSeedEvaluation,
        targetCorrectorSlotSeedEvaluation, replayedEquationResiduals, target, mostUnstable, minimumNonTrivialTpd,
        targetDistance, targetMatches, globallyStable, failureMessage);
  }

  private RetainedPhaseLocalStability retainedPhaseStability(TwoToThreePhaseBoundaryPointSolver.Result root,
      CandidatePhase phase) {
    try {
      double[] homogeneousComposition = phase == root.getRetainedPhaseZero() ? root.getPhaseZeroComposition()
          : root.getPhaseOneComposition();
      TwoToThreePhaseArcLengthCorrector.State homogeneous = TwoToThreePhaseArcLengthCorrector.State.create(
          root.getRetainedPhaseZero(), root.getRetainedPhaseOne(), phase, root.getTemperatureK(),
          root.getPressureBara(), root.getBeta(), root.getPhaseZeroComposition(), root.getPhaseOneComposition(),
          homogeneousComposition);
      IncipientPhaseStationarityJacobianAnalyzer.Result curvature = new IncipientPhaseStationarityJacobianAnalyzer(
          template, phase, phase).setFiniteDifferenceStep(retainedPhaseCurvatureStep).analyze(homogeneous);
      boolean accepted = Double.isFinite(curvature.getBifurcationEigenvalue())
          && curvature.getBifurcationEigenvalue() >= -retainedPhaseEigenvalueTolerance
          && Double.isFinite(curvature.getMaximumResidual())
          && curvature.getMaximumResidual() <= retainedPhaseResidualTolerance;
      String failureMessage = null;
      if (!curvature.isFiniteEvaluation() || !Double.isFinite(curvature.getMaximumResidual())) {
        failureMessage = "retained " + phase + " phase curvature is non-finite";
      } else if (!curvature.hasRealBifurcationMode()) {
        failureMessage = "retained " + phase
            + " phase has no real static-bifurcation eigenmode in the finite stationarity spectrum";
      } else if (curvature.getMaximumResidual() > retainedPhaseResidualTolerance) {
        failureMessage = "retained " + phase + " phase stationary self-residual exceeds tolerance";
      }
      return new RetainedPhaseLocalStability(phase, curvature, accepted, failureMessage);
    } catch (RuntimeException error) {
      return new RetainedPhaseLocalStability(phase, null, false,
          "retained " + phase + " phase curvature evaluation failed: " + error.getMessage());
    }
  }

  private static double compositionDistance(double[] first, double[] second) {
    if (first.length != second.length) {
      return Double.POSITIVE_INFINITY;
    }
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Explicit composition seed preserving the identity of one stationary branch between continuation points. */
  public static final class StationarySeed {
    private final CandidatePhase phase;
    private final double[] composition;

    public StationarySeed(CandidatePhase phase, double[] composition) {
      if (phase == null || composition == null || composition.length == 0) {
        throw new IllegalArgumentException("stationary phase and composition are required");
      }
      this.phase = phase;
      this.composition = normalized(composition);
    }

    public static StationarySeed from(Candidate candidate) {
      if (candidate == null || !candidate.isConverged()) {
        throw new IllegalArgumentException("a converged stationary candidate is required");
      }
      return new StationarySeed(candidate.getPhase(), candidate.getComposition());
    }

    public CandidatePhase getPhase() {
      return phase;
    }

    public double[] getComposition() {
      return composition.clone();
    }

    private static double[] normalized(double[] values) {
      double[] result = values.clone();
      double total = 0.0;
      for (double value : result) {
        if (!Double.isFinite(value) || value < 0.0) {
          throw new IllegalArgumentException("stationary composition contains an invalid value");
        }
        total += value;
      }
      if (!(total > 0.0) || !Double.isFinite(total)) {
        throw new IllegalArgumentException("stationary composition cannot be normalized");
      }
      for (int index = 0; index < result.length; index++) {
        result[index] /= total;
      }
      return result;
    }
  }

  /** Independent homogeneous TPD-Hessian evidence for one retained phase. */
  public static final class RetainedPhaseLocalStability {
    private final CandidatePhase phase;
    private final IncipientPhaseStationarityJacobianAnalyzer.Result curvature;
    private final boolean accepted;
    private final String failureMessage;

    private RetainedPhaseLocalStability(CandidatePhase phase,
        IncipientPhaseStationarityJacobianAnalyzer.Result curvature, boolean accepted, String failureMessage) {
      this.phase = phase;
      this.curvature = curvature;
      this.accepted = accepted;
      this.failureMessage = failureMessage;
    }

    public CandidatePhase getPhase() {
      return phase;
    }

    public IncipientPhaseStationarityJacobianAnalyzer.Result getCurvature() {
      return curvature;
    }

    public double getMinimumEigenvalue() {
      return curvature == null ? Double.NEGATIVE_INFINITY : curvature.getBifurcationEigenvalue();
    }

    public boolean isAccepted() {
      return accepted && failureMessage == null;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }

  /** Immutable global stability result retaining every generic and target-seeded trial. */
  public static final class Result {
    private final TwoToThreePhaseBoundaryPointSolver.Result boundaryRoot;
    private final SpecifiedTwoPhaseFlashSolver.Result retainedFlash;
    private final IncipientPhaseStabilityAnalyzer.Result genericAnalysis;
    private final List<Candidate> trials;
    private final List<Candidate> seededCandidates;
    private final List<RetainedPhaseLocalStability> retainedPhaseStability;
    private final SeedEvaluation targetSeedEvaluation;
    private final SeedEvaluation targetCorrectorSlotSeedEvaluation;
    private final TwoToThreePhaseArcLengthCorrector.EquationResiduals replayedEquationResiduals;
    private final Candidate targetCandidate;
    private final Candidate mostUnstableCandidate;
    private final double minimumNonTrivialTangentPlaneDistance;
    private final double targetCompositionDistance;
    private final boolean targetMatches;
    private final boolean globallyStable;
    private final String failureMessage;

    private Result(TwoToThreePhaseBoundaryPointSolver.Result boundaryRoot,
        SpecifiedTwoPhaseFlashSolver.Result retainedFlash, IncipientPhaseStabilityAnalyzer.Result genericAnalysis,
        List<Candidate> trials, List<Candidate> seededCandidates,
        List<RetainedPhaseLocalStability> retainedPhaseStability, SeedEvaluation targetSeedEvaluation,
        SeedEvaluation targetCorrectorSlotSeedEvaluation,
        TwoToThreePhaseArcLengthCorrector.EquationResiduals replayedEquationResiduals, Candidate targetCandidate,
        Candidate mostUnstableCandidate, double minimumNonTrivialTangentPlaneDistance, double targetCompositionDistance,
        boolean targetMatches, boolean globallyStable, String failureMessage) {
      this.boundaryRoot = boundaryRoot;
      this.retainedFlash = retainedFlash;
      this.genericAnalysis = genericAnalysis;
      this.trials = Collections.unmodifiableList(new ArrayList<Candidate>(trials));
      this.seededCandidates = Collections.unmodifiableList(new ArrayList<Candidate>(seededCandidates));
      this.retainedPhaseStability = Collections
          .unmodifiableList(new ArrayList<RetainedPhaseLocalStability>(retainedPhaseStability));
      this.targetSeedEvaluation = targetSeedEvaluation;
      this.targetCorrectorSlotSeedEvaluation = targetCorrectorSlotSeedEvaluation;
      this.replayedEquationResiduals = replayedEquationResiduals;
      this.targetCandidate = targetCandidate;
      this.mostUnstableCandidate = mostUnstableCandidate;
      this.minimumNonTrivialTangentPlaneDistance = minimumNonTrivialTangentPlaneDistance;
      this.targetCompositionDistance = targetCompositionDistance;
      this.targetMatches = targetMatches;
      this.globallyStable = globallyStable;
      this.failureMessage = failureMessage;
    }

    private static Result failure(TwoToThreePhaseBoundaryPointSolver.Result root, String failureMessage) {
      return new Result(root, null, null, Collections.<Candidate>emptyList(), Collections.<Candidate>emptyList(),
          Collections.<RetainedPhaseLocalStability>emptyList(), null, null, null, null, null, Double.NEGATIVE_INFINITY,
          Double.POSITIVE_INFINITY, false, false, failureMessage);
    }

    public boolean isAccepted() {
      return globallyStable && failureMessage == null;
    }

    public TwoToThreePhaseBoundaryPointSolver.Result getBoundaryRoot() {
      return boundaryRoot;
    }

    public SpecifiedTwoPhaseFlashSolver.Result getRetainedFlash() {
      return retainedFlash;
    }

    public IncipientPhaseStabilityAnalyzer.Result getGenericAnalysis() {
      return genericAnalysis;
    }

    public List<Candidate> getTrials() {
      return trials;
    }

    /** @return candidates evaluated from stationary identities propagated by the continuation driver */
    public List<Candidate> getSeededCandidates() {
      return seededCandidates;
    }

    /** @return homogeneous TPD-Hessian evidence for both retained phases */
    public List<RetainedPhaseLocalStability> getRetainedPhaseStability() {
      return retainedPhaseStability;
    }

    public Candidate getTargetCandidate() {
      return targetCandidate;
    }

    /** @return direct non-iterative evidence for the corrected incipient composition */
    public SeedEvaluation getTargetSeedEvaluation() {
      return targetSeedEvaluation;
    }

    /** @return the same direct seed evaluated in the third phase slot used by the strict corrector */
    public SeedEvaluation getTargetCorrectorSlotSeedEvaluation() {
      return targetCorrectorSlotSeedEvaluation;
    }

    /** @return cold-start replay residual of the strict boundary equations */
    public double getReplayedBoundaryResidual() {
      return replayedEquationResiduals == null ? Double.NaN
          : replayedEquationResiduals.getThermodynamicMaximumResidual();
    }

    /** @return cold-start replay residuals split by physical equation group */
    public TwoToThreePhaseArcLengthCorrector.EquationResiduals getReplayedEquationResiduals() {
      return replayedEquationResiduals;
    }

    public Candidate getMostUnstableCandidate() {
      return mostUnstableCandidate;
    }

    public double getMinimumNonTrivialTangentPlaneDistance() {
      return minimumNonTrivialTangentPlaneDistance;
    }

    public double getTargetCompositionDistance() {
      return targetCompositionDistance;
    }

    public boolean isTargetMatched() {
      return targetMatches;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
