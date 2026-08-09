package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Tracks one explicitly seeded secondary TPD stationary branch along an existing PT boundary.
 *
 * <p>
 * Generic Wilson gas/oil/water stability seeds are not sufficient when two stationary points of the same physical
 * family coexist. This tracker starts from a known terminal stationary composition and walks backwards over ordered
 * boundary roots. The converged composition at one root is the only seed used at the preceding root. A composition
 * jump, a phase-family change, a trivial retained-phase solution, or a failed retained flash terminates the tracked
 * segment instead of silently switching stationary-point identity.
 * </p>
 */
public final class HydrocarbonWaterSecondaryStationaryBranchTracker {
  private final SystemInterface template;
  private int flashMaximumIterations = 100;
  private int stabilityMaximumIterations = 500;
  private double convergenceTolerance = 1.0e-9;
  private double finiteDifferenceStep = 2.0e-5;
  private double dampingFactor = 0.2;
  private double maximumCompositionJump = 0.15;
  private double zeroTangentPlaneTolerance = 1.0e-7;
  private double targetCompositionTolerance = 1.0e-4;

  /** Creates a non-destructive branch tracker for one configured EOS and mixing rule. */
  public HydrocarbonWaterSecondaryStationaryBranchTracker(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets retained-flash, stationary-iteration, composition-continuity, and zero-TPD controls. */
  public HydrocarbonWaterSecondaryStationaryBranchTracker setNumericalControls(int flashMaximumIterations,
      int stabilityMaximumIterations, double convergenceTolerance, double finiteDifferenceStep, double dampingFactor,
      double maximumCompositionJump, double zeroTangentPlaneTolerance) {
    if (flashMaximumIterations < 5 || stabilityMaximumIterations < 5 || !positive(convergenceTolerance)
        || !positive(finiteDifferenceStep) || !positive(dampingFactor) || dampingFactor > 1.0
        || !positive(maximumCompositionJump) || !positive(zeroTangentPlaneTolerance)) {
      throw new IllegalArgumentException("invalid secondary-stationary-branch numerical controls");
    }
    this.flashMaximumIterations = flashMaximumIterations;
    this.stabilityMaximumIterations = stabilityMaximumIterations;
    this.convergenceTolerance = convergenceTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    this.dampingFactor = dampingFactor;
    this.maximumCompositionJump = maximumCompositionJump;
    this.zeroTangentPlaneTolerance = zeroTangentPlaneTolerance;
    return this;
  }

  /** Sets the identity gate used to stop when the secondary branch merges with the boundary target branch. */
  public HydrocarbonWaterSecondaryStationaryBranchTracker setTargetCompositionTolerance(
      double targetCompositionTolerance) {
    if (!positive(targetCompositionTolerance)) {
      throw new IllegalArgumentException("target-composition tolerance must be positive");
    }
    this.targetCompositionTolerance = targetCompositionTolerance;
    return this;
  }

  /**
   * Tracks a known stationary phase backwards over boundary roots ordered from early to late continuation state.
   *
   * @param orderedBoundaryRoots converged roots in boundary-continuation order
   * @param stationaryPhase physical family of the secondary stationary branch
   * @param terminalCompositionSeed known composition at the last supplied root
   * @return immutable samples in the same order as the input roots
   */
  public Result trackBackward(List<TwoToThreePhaseBoundaryPointSolver.Result> orderedBoundaryRoots,
      CandidatePhase stationaryPhase, double[] terminalCompositionSeed) {
    if (orderedBoundaryRoots == null || orderedBoundaryRoots.isEmpty() || stationaryPhase == null
        || terminalCompositionSeed == null) {
      throw new IllegalArgumentException("ordered roots, stationary phase, and terminal composition are required");
    }
    validateRoots(orderedBoundaryRoots, terminalCompositionSeed.length);

    List<Sample> reverseSamples = new ArrayList<Sample>();
    double[] seed = normalized(terminalCompositionSeed);
    String terminationMessage = null;
    for (int rootIndex = orderedBoundaryRoots.size() - 1; rootIndex >= 0; rootIndex--) {
      TwoToThreePhaseBoundaryPointSolver.Result root = orderedBoundaryRoots.get(rootIndex);
      Sample sample = evaluate(root, stationaryPhase, seed);
      reverseSamples.add(sample);
      if (!sample.isTracked()) {
        terminationMessage = sample.getFailureMessage();
        break;
      }
      seed = sample.getCandidate().getComposition();
    }
    Collections.reverse(reverseSamples);

    List<ZeroTpdBracket> brackets = new ArrayList<ZeroTpdBracket>();
    List<BranchLimitBracket> branchLimits = new ArrayList<BranchLimitBracket>();
    List<TargetBranchMergeBracket> targetBranchMerges = new ArrayList<TargetBranchMergeBracket>();
    for (int index = 1; index < reverseSamples.size(); index++) {
      Sample lower = reverseSamples.get(index - 1);
      Sample upper = reverseSamples.get(index);
      if (lower.isTracked() && upper.isTracked()
          && bracketsZeroTpd(lower.getTangentPlaneDistance(), upper.getTangentPlaneDistance())) {
        brackets.add(new ZeroTpdBracket(lower, upper));
      }
      if (lower.isCollapsedStationaryBranch() && upper.isTracked()) {
        branchLimits.add(new BranchLimitBracket(lower, upper));
      }
      if (lower.isCollapsedTargetBranch() && upper.isTracked()) {
        targetBranchMerges.add(new TargetBranchMergeBracket(lower, upper));
      }
    }
    return new Result(stationaryPhase, reverseSamples, brackets, branchLimits, targetBranchMerges, terminationMessage);
  }

  private Sample evaluate(TwoToThreePhaseBoundaryPointSolver.Result root, CandidatePhase stationaryPhase,
      double[] seed) {
    SpecifiedTwoPhaseFlashSolver flashSolver = new SpecifiedTwoPhaseFlashSolver(template, root.getRetainedPhaseZero(),
        root.getRetainedPhaseOne())
        .setNumericalControls(flashMaximumIterations, convergenceTolerance, finiteDifferenceStep);
    SpecifiedTwoPhaseFlashSolver.Result retained = flashSolver.solve(root.getTemperatureK(), root.getPressureBara(),
        root.getBeta(), root.getPhaseZeroComposition(), root.getPhaseOneComposition());
    if (!retained.isConverged()) {
      return Sample.failure(root, retained, null, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, false,
          "retained two-phase flash failed while tracking secondary stationary branch: "
              + retained.getFailureMessage());
    }

    SystemInterface reference = flashSolver.toThermodynamicSystem(retained);
    Candidate candidate = new IncipientPhaseStabilityAnalyzer(reference)
        .setMaximumIterations(stabilityMaximumIterations)
        .setTolerances(convergenceTolerance, -zeroTangentPlaneTolerance).setDampingFactor(dampingFactor)
        .analyzeCandidate(stationaryPhase, seed);
    double compositionJump = candidate.isConverged() ? compositionDistance(seed, candidate.getComposition())
        : Double.POSITIVE_INFINITY;
    double targetCompositionDistance = candidate.isConverged()
        ? compositionDistance(candidate.getComposition(), root.getIncipientComposition())
        : Double.POSITIVE_INFINITY;
    String failureMessage = null;
    boolean collapsedTargetBranch = false;
    if (!candidate.isConverged()) {
      failureMessage = "secondary stationary solve failed: " + candidate.getFailureMessage();
    } else if (candidate.getPhase() != stationaryPhase) {
      failureMessage = "secondary stationary solve changed physical family from " + stationaryPhase + " to "
          + candidate.getPhase();
    } else if (candidate.isTrivial()) {
      failureMessage = "secondary stationary branch collapsed onto an existing retained phase";
    } else if (!Double.isFinite(candidate.getTangentPlaneDistance())) {
      failureMessage = "secondary stationary branch returned a non-finite tangent-plane distance";
    } else if (targetCompositionDistance <= targetCompositionTolerance) {
      collapsedTargetBranch = true;
      failureMessage = "secondary stationary branch merged with the boundary target incipient branch";
    } else if (compositionJump > maximumCompositionJump) {
      failureMessage = "secondary stationary composition continuity gate rejected jump=" + compositionJump;
    }
    return failureMessage == null
        ? Sample.tracked(root, retained, candidate, compositionJump, targetCompositionDistance)
        : Sample.failure(root, retained, candidate, compositionJump, targetCompositionDistance, collapsedTargetBranch,
            failureMessage);
  }

  private void validateRoots(List<TwoToThreePhaseBoundaryPointSolver.Result> roots, int componentCount) {
    CandidatePhase retainedZero = null;
    CandidatePhase retainedOne = null;
    CandidatePhase incipient = null;
    for (TwoToThreePhaseBoundaryPointSolver.Result root : roots) {
      if (root == null || !root.isConverged() || root.getPhaseZeroComposition().length != componentCount
          || root.getPhaseOneComposition().length != componentCount
          || root.getIncipientComposition().length != componentCount) {
        throw new IllegalArgumentException(
            "every tracked boundary root must be converged with compatible compositions");
      }
      if (retainedZero == null) {
        retainedZero = root.getRetainedPhaseZero();
        retainedOne = root.getRetainedPhaseOne();
        incipient = root.getIncipientPhase();
      } else if (root.getRetainedPhaseZero() != retainedZero || root.getRetainedPhaseOne() != retainedOne
          || root.getIncipientPhase() != incipient) {
        throw new IllegalArgumentException("all roots must belong to the same physical boundary family");
      }
    }
  }

  private boolean bracketsZeroTpd(double first, double second) {
    if (!Double.isFinite(first) || !Double.isFinite(second)) {
      return false;
    }
    if (Math.abs(first) <= zeroTangentPlaneTolerance || Math.abs(second) <= zeroTangentPlaneTolerance) {
      return true;
    }
    return Math.copySign(1.0, first) != Math.copySign(1.0, second);
  }

  private static double[] normalized(double[] composition) {
    double[] result = composition.clone();
    double sum = 0.0;
    for (double value : result) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException("terminal stationary composition contains an invalid value");
      }
      sum += value;
    }
    if (!(sum > 0.0) || !Double.isFinite(sum)) {
      throw new IllegalArgumentException("terminal stationary composition cannot be normalized");
    }
    for (int index = 0; index < result.length; index++) {
      result[index] /= sum;
    }
    return result;
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

  /** One explicit stationary-branch evaluation at one boundary root. */
  public static final class Sample {
    private final TwoToThreePhaseBoundaryPointSolver.Result boundaryRoot;
    private final SpecifiedTwoPhaseFlashSolver.Result retainedFlash;
    private final Candidate candidate;
    private final double compositionJump;
    private final double targetCompositionDistance;
    private final boolean collapsedTargetBranch;
    private final boolean tracked;
    private final String failureMessage;

    private Sample(TwoToThreePhaseBoundaryPointSolver.Result boundaryRoot,
        SpecifiedTwoPhaseFlashSolver.Result retainedFlash, Candidate candidate, double compositionJump,
        double targetCompositionDistance, boolean collapsedTargetBranch, boolean tracked, String failureMessage) {
      this.boundaryRoot = boundaryRoot;
      this.retainedFlash = retainedFlash;
      this.candidate = candidate;
      this.compositionJump = compositionJump;
      this.targetCompositionDistance = targetCompositionDistance;
      this.collapsedTargetBranch = collapsedTargetBranch;
      this.tracked = tracked;
      this.failureMessage = failureMessage;
    }

    private static Sample tracked(TwoToThreePhaseBoundaryPointSolver.Result root,
        SpecifiedTwoPhaseFlashSolver.Result retained, Candidate candidate, double compositionJump,
        double targetCompositionDistance) {
      return new Sample(root, retained, candidate, compositionJump, targetCompositionDistance, false, true, null);
    }

    private static Sample failure(TwoToThreePhaseBoundaryPointSolver.Result root,
        SpecifiedTwoPhaseFlashSolver.Result retained, Candidate candidate, double compositionJump,
        double targetCompositionDistance, boolean collapsedTargetBranch, String failureMessage) {
      return new Sample(root, retained, candidate, compositionJump, targetCompositionDistance, collapsedTargetBranch,
          false, failureMessage);
    }

    public TwoToThreePhaseBoundaryPointSolver.Result getBoundaryRoot() {
      return boundaryRoot;
    }

    public SpecifiedTwoPhaseFlashSolver.Result getRetainedFlash() {
      return retainedFlash;
    }

    public Candidate getCandidate() {
      return candidate;
    }

    public double getCompositionJump() {
      return compositionJump;
    }

    public double getTargetCompositionDistance() {
      return targetCompositionDistance;
    }

    public double getTangentPlaneDistance() {
      return candidate == null ? Double.NaN : candidate.getTangentPlaneDistance();
    }

    public boolean isTracked() {
      return tracked;
    }

    /** @return true when the explicitly continued branch has merged with an existing retained phase */
    public boolean isCollapsedStationaryBranch() {
      return !tracked && candidate != null && candidate.isConverged() && candidate.isTrivial()
          && Double.isFinite(candidate.getTangentPlaneDistance());
    }

    /** @return true when the secondary identity merged with the boundary's current incipient-phase identity */
    public boolean isCollapsedTargetBranch() {
      return collapsedTargetBranch;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }

  /** Adjacent collapsed and non-trivial samples enclosing a stationary-branch bifurcation. */
  public static final class BranchLimitBracket {
    private final Sample collapsed;
    private final Sample nonTrivial;

    private BranchLimitBracket(Sample collapsed, Sample nonTrivial) {
      this.collapsed = collapsed;
      this.nonTrivial = nonTrivial;
    }

    public Sample getCollapsed() {
      return collapsed;
    }

    public Sample getNonTrivial() {
      return nonTrivial;
    }
  }

  /** Adjacent merged and distinct samples enclosing a secondary-to-target stationary-branch junction. */
  public static final class TargetBranchMergeBracket {
    private final Sample merged;
    private final Sample distinct;

    private TargetBranchMergeBracket(Sample merged, Sample distinct) {
      this.merged = merged;
      this.distinct = distinct;
    }

    public Sample getMerged() {
      return merged;
    }

    public Sample getDistinct() {
      return distinct;
    }
  }

  /** Adjacent samples enclosing or touching zero TPD on the same stationary branch. */
  public static final class ZeroTpdBracket {
    private final Sample first;
    private final Sample second;

    private ZeroTpdBracket(Sample first, Sample second) {
      this.first = first;
      this.second = second;
    }

    public Sample getFirst() {
      return first;
    }

    public Sample getSecond() {
      return second;
    }
  }

  /** Immutable backward-tracking evidence in the original boundary order. */
  public static final class Result {
    private final CandidatePhase stationaryPhase;
    private final List<Sample> samples;
    private final List<ZeroTpdBracket> zeroTpdBrackets;
    private final List<BranchLimitBracket> branchLimitBrackets;
    private final List<TargetBranchMergeBracket> targetBranchMergeBrackets;
    private final String terminationMessage;

    private Result(CandidatePhase stationaryPhase, List<Sample> samples, List<ZeroTpdBracket> zeroTpdBrackets,
        List<BranchLimitBracket> branchLimitBrackets, List<TargetBranchMergeBracket> targetBranchMergeBrackets,
        String terminationMessage) {
      this.stationaryPhase = stationaryPhase;
      this.samples = Collections.unmodifiableList(new ArrayList<Sample>(samples));
      this.zeroTpdBrackets = Collections.unmodifiableList(new ArrayList<ZeroTpdBracket>(zeroTpdBrackets));
      this.branchLimitBrackets = Collections.unmodifiableList(new ArrayList<BranchLimitBracket>(branchLimitBrackets));
      this.targetBranchMergeBrackets = Collections
          .unmodifiableList(new ArrayList<TargetBranchMergeBracket>(targetBranchMergeBrackets));
      this.terminationMessage = terminationMessage;
    }

    public CandidatePhase getStationaryPhase() {
      return stationaryPhase;
    }

    public List<Sample> getSamples() {
      return samples;
    }

    public List<ZeroTpdBracket> getZeroTpdBrackets() {
      return zeroTpdBrackets;
    }

    public List<BranchLimitBracket> getBranchLimitBrackets() {
      return branchLimitBrackets;
    }

    public List<TargetBranchMergeBracket> getTargetBranchMergeBrackets() {
      return targetBranchMergeBrackets;
    }

    public boolean reachesFirstBoundaryRoot() {
      return samples.size() > 0 && samples.get(0).isTracked() && terminationMessage == null;
    }

    public String getTerminationMessage() {
      return terminationMessage;
    }
  }
}
