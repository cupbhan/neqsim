package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Attempts a retained-phase branch switch after a globally rejected two-to-three-phase boundary correction.
 *
 * <p>
 * A lower-TPD stationary point belonging to an already retained phase family can mean either that the specified
 * two-phase flash converged to the wrong local branch, or that an additional same-family liquid phase is required. This
 * solver first tests the less expansive explanation. It replaces the matching retained-phase seed, scans independent
 * phase-fraction and composition-blend starts, and then re-solves the requested third-phase zero-TPD boundary on every
 * distinct retained branch recovered. If no globally stable corrected boundary is found, the full evidence remains
 * available to a subsequent multi-liquid-phase solver; no guessed branch is promoted.
 * </p>
 */
public final class HydrocarbonWaterRetainedPhaseBranchSwitcher {
  private final SystemInterface template;
  private int maximumFlashIterations = 120;
  private double flashResidualTolerance = 1.0e-9;
  private double finiteDifferenceStep = 2.0e-5;
  private int boundaryScanIntervals = 32;
  private int maximumBoundaryBisections = 80;
  private double temperatureToleranceK = 1.0e-5;
  private double tangentPlaneTolerance = 1.0e-8;
  private double distinctBranchCompositionDistance = 1.0e-4;

  /** Creates a non-destructive switcher for one configured EOS and mixing rule. */
  public HydrocarbonWaterRetainedPhaseBranchSwitcher(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets the specified two-phase flash controls used by every branch-switch start. */
  public HydrocarbonWaterRetainedPhaseBranchSwitcher setFlashControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    if (maximumIterations < 1 || !positive(residualTolerance) || !positive(finiteDifferenceStep)) {
      throw new IllegalArgumentException("invalid retained-phase branch-switch flash controls");
    }
    this.maximumFlashIterations = maximumIterations;
    this.flashResidualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /** Sets fixed-pressure water-onset correction controls and the minimum distinct retained-branch distance. */
  public HydrocarbonWaterRetainedPhaseBranchSwitcher setBoundaryControls(int scanIntervals, int maximumBisections,
      double temperatureToleranceK, double tangentPlaneTolerance, double distinctBranchCompositionDistance) {
    if (scanIntervals < 2 || maximumBisections < 1 || !positive(temperatureToleranceK)
        || !positive(tangentPlaneTolerance) || !positive(distinctBranchCompositionDistance)) {
      throw new IllegalArgumentException("invalid retained-phase branch-switch boundary controls");
    }
    this.boundaryScanIntervals = scanIntervals;
    this.maximumBoundaryBisections = maximumBisections;
    this.temperatureToleranceK = temperatureToleranceK;
    this.tangentPlaneTolerance = tangentPlaneTolerance;
    this.distinctBranchCompositionDistance = distinctBranchCompositionDistance;
    return this;
  }

  /**
   * Attempts to replace the destabilized retained phase and correct the target boundary at the same pressure.
   *
   * @param rejected globally rejected, independently target-matched boundary correction
   * @param minimumTemperatureK lower correction temperature
   * @param maximumTemperatureK upper correction temperature
   * @return immutable multistart and corrected-root evidence
   */
  public Result switchAndCorrect(HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejected,
      double minimumTemperatureK, double maximumTemperatureK) {
    return switchAndCorrect(rejected, rejected == null ? null : rejected.getMostUnstableCandidate(),
        minimumTemperatureK, maximumTemperatureK);
  }

  /**
   * Attempts the branch switch with an explicitly continued same-family stationary point.
   *
   * <p>
   * The global stability multistart can miss a narrow secondary stationary branch close to a retained-phase spinodal. A
   * composition-continuous branch tracker may therefore supply stronger evidence than the generic seed inventory. The
   * explicit candidate is still required to be converged, non-trivial, lower-TPD, and in a retained phase family; it
   * never bypasses the subsequent flash, boundary, or global-topology gates.
   * </p>
   *
   * @param rejected globally rejected, independently target-matched boundary correction
   * @param explicitCandidate composition-continuous lower-TPD retained-family stationary point
   * @param minimumTemperatureK lower correction temperature
   * @param maximumTemperatureK upper correction temperature
   * @return immutable multistart and corrected-root evidence
   */
  public Result switchAndCorrect(HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejected,
      Candidate explicitCandidate, double minimumTemperatureK, double maximumTemperatureK) {
    if (rejected == null || rejected.getBoundaryRoot() == null || !Double.isFinite(minimumTemperatureK)
        || !Double.isFinite(maximumTemperatureK) || minimumTemperatureK < 50.0
        || maximumTemperatureK <= minimumTemperatureK) {
      throw new IllegalArgumentException("rejected boundary evidence and a finite temperature interval are required");
    }
    TwoToThreePhaseBoundaryPointSolver.Result root = rejected.getBoundaryRoot();
    Candidate candidate = explicitCandidate;
    if (rejected.isAccepted() || !rejected.isTargetMatched() || candidate == null || !candidate.isConverged()
        || candidate.isTrivial() || !(candidate.getTangentPlaneDistance() < 0.0)) {
      return Result.failure(rejected, null, Collections.<Attempt>emptyList(), Collections.<Correction>emptyList(),
          "global-stability evidence does not contain a lower-TPD non-trivial stationary branch");
    }
    CandidatePhase replacementPhase = candidate.getPhase();
    boolean replacePhaseZero = replacementPhase == root.getRetainedPhaseZero();
    boolean replacePhaseOne = replacementPhase == root.getRetainedPhaseOne();
    if (!replacePhaseZero && !replacePhaseOne) {
      return Result.failure(rejected, replacementPhase, Collections.<Attempt>emptyList(),
          Collections.<Correction>emptyList(), "the destabilizing phase is not one of the retained phase families");
    }

    double[] originalZero = root.getPhaseZeroComposition();
    double[] originalOne = root.getPhaseOneComposition();
    double[] alternative = candidate.getComposition();
    double[] blendFractions = new double[] {1.0, 0.75, 0.5, 0.25};
    double[] betaSeeds = distinctBetaSeeds(root.getBeta());
    SpecifiedTwoPhaseFlashSolver flashSolver = new SpecifiedTwoPhaseFlashSolver(template, root.getRetainedPhaseZero(),
        root.getRetainedPhaseOne())
        .setNumericalControls(maximumFlashIterations, flashResidualTolerance, finiteDifferenceStep);
    List<Attempt> attempts = new ArrayList<Attempt>();
    List<SpecifiedTwoPhaseFlashSolver.Result> distinctBranches = new ArrayList<SpecifiedTwoPhaseFlashSolver.Result>();
    for (double blendFraction : blendFractions) {
      double[] zeroSeed = replacePhaseZero ? blend(originalZero, alternative, blendFraction) : originalZero;
      double[] oneSeed = replacePhaseOne ? blend(originalOne, alternative, blendFraction) : originalOne;
      for (double betaSeed : betaSeeds) {
        SpecifiedTwoPhaseFlashSolver.Result flash = flashSolver.solve(root.getTemperatureK(), root.getPressureBara(),
            betaSeed, zeroSeed, oneSeed);
        double replacementDistance = flash.isConverged()
            ? compositionDistance(replacePhaseZero ? flash.getPhaseZeroComposition() : flash.getPhaseOneComposition(),
                replacePhaseZero ? originalZero : originalOne)
            : Double.NaN;
        double retainedPairDistance = flash.isConverged()
            ? compositionDistance(flash.getPhaseZeroComposition(), flash.getPhaseOneComposition())
            : Double.NaN;
        boolean distinct = flash.isConverged() && replacementDistance > distinctBranchCompositionDistance
            && retainedPairDistance > distinctBranchCompositionDistance;
        attempts.add(new Attempt(blendFraction, betaSeed, replacementDistance, distinct, flash));
        if (distinct && !containsEquivalentBranch(distinctBranches, flash)) {
          distinctBranches.add(flash);
        }
      }
    }

    List<Correction> corrections = new ArrayList<Correction>();
    List<HydrocarbonWaterBoundaryEndpointClassifier.Result> accepted = new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>();
    for (SpecifiedTwoPhaseFlashSolver.Result branch : distinctBranches) {
      TwoToThreePhaseBoundaryPointSolver.RootSet roots = new TwoToThreePhaseBoundaryPointSolver(template,
          root.getRetainedPhaseZero(), root.getRetainedPhaseOne(), root.getIncipientPhase())
          .setNumericalControls(boundaryScanIntervals, maximumBoundaryBisections, temperatureToleranceK,
              tangentPlaneTolerance)
          .solveAll(root.getPressureBara(), minimumTemperatureK, maximumTemperatureK, branch.getBeta(),
              branch.getPhaseZeroComposition(), branch.getPhaseOneComposition(),
              Collections.singletonList(root.getIncipientComposition()));
      List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications = new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>();
      HydrocarbonWaterBoundaryEndpointClassifier classifier = new HydrocarbonWaterBoundaryEndpointClassifier(template);
      for (TwoToThreePhaseBoundaryPointSolver.Result correctedRoot : roots.getRoots()) {
        HydrocarbonWaterBoundaryEndpointClassifier.Result classification = classifier.classify(correctedRoot);
        classifications.add(classification);
        if (classification.isOrdinaryBoundaryPoint()) {
          accepted.add(classification);
        }
      }
      corrections.add(new Correction(branch, roots, classifications));
    }

    String failureMessage = null;
    if (distinctBranches.isEmpty()) {
      failureMessage = "multistart specified flash returned no retained branch distinct from the rejected branch";
    } else if (corrections.stream().noneMatch(correction -> correction.rootSet.hasConvergedRoots())) {
      failureMessage = "distinct retained branches produced no corrected target-phase zero-TPD root";
    } else if (accepted.isEmpty()) {
      failureMessage = "corrected roots on the switched retained branches did not pass global topology gates";
    }
    return new Result(rejected, replacementPhase, attempts, corrections, accepted, distinctBranches.size(),
        failureMessage);
  }

  private static double[] distinctBetaSeeds(double sourceBeta) {
    double[] candidates = new double[] {sourceBeta, 0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98};
    List<Double> distinct = new ArrayList<Double>();
    for (double candidate : candidates) {
      boolean duplicate = false;
      for (double current : distinct) {
        duplicate |= Math.abs(candidate - current) <= 1.0e-10;
      }
      if (!duplicate) {
        distinct.add(candidate);
      }
    }
    double[] values = new double[distinct.size()];
    for (int index = 0; index < values.length; index++) {
      values[index] = distinct.get(index);
    }
    return values;
  }

  private static double[] blend(double[] original, double[] alternative, double alternativeFraction) {
    double[] result = new double[original.length];
    double sum = 0.0;
    for (int index = 0; index < result.length; index++) {
      result[index] = Math.max(0.0,
          (1.0 - alternativeFraction) * original[index] + alternativeFraction * alternative[index]);
      sum += result[index];
    }
    for (int index = 0; index < result.length; index++) {
      result[index] /= sum;
    }
    return result;
  }

  private boolean containsEquivalentBranch(List<SpecifiedTwoPhaseFlashSolver.Result> branches,
      SpecifiedTwoPhaseFlashSolver.Result candidate) {
    for (SpecifiedTwoPhaseFlashSolver.Result branch : branches) {
      double distance = Math.max(
          compositionDistance(branch.getPhaseZeroComposition(), candidate.getPhaseZeroComposition()),
          compositionDistance(branch.getPhaseOneComposition(), candidate.getPhaseOneComposition()));
      if (distance <= distinctBranchCompositionDistance && Math.abs(branch.getBeta() - candidate.getBeta()) <= 1.0e-5) {
        return true;
      }
    }
    return false;
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

  /** One specified-flash start and its auditable convergence result. */
  public static final class Attempt {
    private final double alternativeCompositionFraction;
    private final double betaSeed;
    private final double replacementCompositionDistance;
    private final boolean distinctBranch;
    private final SpecifiedTwoPhaseFlashSolver.Result flash;

    private Attempt(double alternativeCompositionFraction, double betaSeed, double replacementCompositionDistance,
        boolean distinctBranch, SpecifiedTwoPhaseFlashSolver.Result flash) {
      this.alternativeCompositionFraction = alternativeCompositionFraction;
      this.betaSeed = betaSeed;
      this.replacementCompositionDistance = replacementCompositionDistance;
      this.distinctBranch = distinctBranch;
      this.flash = flash;
    }

    public double getAlternativeCompositionFraction() {
      return alternativeCompositionFraction;
    }

    public double getBetaSeed() {
      return betaSeed;
    }

    public double getReplacementCompositionDistance() {
      return replacementCompositionDistance;
    }

    public boolean isDistinctBranch() {
      return distinctBranch;
    }

    public SpecifiedTwoPhaseFlashSolver.Result getFlash() {
      return flash;
    }
  }

  /** Corrected target-phase boundary evidence for one distinct retained branch. */
  public static final class Correction {
    private final SpecifiedTwoPhaseFlashSolver.Result retainedBranch;
    private final TwoToThreePhaseBoundaryPointSolver.RootSet rootSet;
    private final List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications;

    private Correction(SpecifiedTwoPhaseFlashSolver.Result retainedBranch,
        TwoToThreePhaseBoundaryPointSolver.RootSet rootSet,
        List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications) {
      this.retainedBranch = retainedBranch;
      this.rootSet = rootSet;
      this.classifications = Collections
          .unmodifiableList(new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>(classifications));
    }

    public SpecifiedTwoPhaseFlashSolver.Result getRetainedBranch() {
      return retainedBranch;
    }

    public TwoToThreePhaseBoundaryPointSolver.RootSet getRootSet() {
      return rootSet;
    }

    public List<HydrocarbonWaterBoundaryEndpointClassifier.Result> getClassifications() {
      return classifications;
    }
  }

  /** Immutable branch-switch decision retaining every attempted start and corrected root. */
  public static final class Result {
    private final HydrocarbonWaterBoundaryGlobalStabilityGate.Result sourceEvidence;
    private final CandidatePhase replacementPhase;
    private final List<Attempt> attempts;
    private final List<Correction> corrections;
    private final List<HydrocarbonWaterBoundaryEndpointClassifier.Result> acceptedBoundaryRoots;
    private final int distinctRetainedBranchCount;
    private final String failureMessage;

    private Result(HydrocarbonWaterBoundaryGlobalStabilityGate.Result sourceEvidence, CandidatePhase replacementPhase,
        List<Attempt> attempts, List<Correction> corrections,
        List<HydrocarbonWaterBoundaryEndpointClassifier.Result> acceptedBoundaryRoots, int distinctRetainedBranchCount,
        String failureMessage) {
      this.sourceEvidence = sourceEvidence;
      this.replacementPhase = replacementPhase;
      this.attempts = Collections.unmodifiableList(new ArrayList<Attempt>(attempts));
      this.corrections = Collections.unmodifiableList(new ArrayList<Correction>(corrections));
      this.acceptedBoundaryRoots = Collections
          .unmodifiableList(new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>(acceptedBoundaryRoots));
      this.distinctRetainedBranchCount = distinctRetainedBranchCount;
      this.failureMessage = failureMessage;
    }

    private static Result failure(HydrocarbonWaterBoundaryGlobalStabilityGate.Result sourceEvidence,
        CandidatePhase replacementPhase, List<Attempt> attempts, List<Correction> corrections, String failureMessage) {
      return new Result(sourceEvidence, replacementPhase, attempts, corrections,
          Collections.<HydrocarbonWaterBoundaryEndpointClassifier.Result>emptyList(), 0, failureMessage);
    }

    public boolean hasGloballyStableBoundaryRoot() {
      return !acceptedBoundaryRoots.isEmpty();
    }

    public HydrocarbonWaterBoundaryGlobalStabilityGate.Result getSourceEvidence() {
      return sourceEvidence;
    }

    public CandidatePhase getReplacementPhase() {
      return replacementPhase;
    }

    public List<Attempt> getAttempts() {
      return attempts;
    }

    public List<Correction> getCorrections() {
      return corrections;
    }

    public List<HydrocarbonWaterBoundaryEndpointClassifier.Result> getAcceptedBoundaryRoots() {
      return acceptedBoundaryRoots;
    }

    public int getDistinctRetainedBranchCount() {
      return distinctRetainedBranchCount;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
