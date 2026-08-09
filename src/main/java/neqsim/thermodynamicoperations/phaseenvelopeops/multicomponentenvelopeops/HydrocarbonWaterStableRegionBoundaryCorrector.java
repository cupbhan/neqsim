package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterStableRegionTransitionScanner.StableState;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Corrects a stable adjacent-region bracket to globally stable zero-TPD boundary roots. */
public final class HydrocarbonWaterStableRegionBoundaryCorrector {
  private final SystemInterface template;
  private int scanIntervals = 8;
  private int maximumBisections = 80;
  private double temperatureToleranceK = 1.0e-7;
  private double tangentPlaneTolerance = 1.0e-8;
  private double maximumExpansionHalfWidthK = 25.0;
  private int maximumIntervalAttempts = Integer.MAX_VALUE;
  private boolean useThreePhaseSideRetainedSeedFallback;

  /** Creates a non-destructive stable-bracket boundary corrector. */
  public HydrocarbonWaterStableRegionBoundaryCorrector(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets strict fixed-pressure correction controls. */
  public HydrocarbonWaterStableRegionBoundaryCorrector setNumericalControls(int scanIntervals, int maximumBisections,
      double temperatureToleranceK, double tangentPlaneTolerance) {
    if (scanIntervals < 2 || maximumBisections < 1 || !positive(temperatureToleranceK)
        || !positive(tangentPlaneTolerance)) {
      throw new IllegalArgumentException("invalid stable-bracket boundary controls");
    }
    this.scanIntervals = scanIntervals;
    this.maximumBisections = maximumBisections;
    this.temperatureToleranceK = temperatureToleranceK;
    this.tangentPlaneTolerance = tangentPlaneTolerance;
    return this;
  }

  /** Sets the largest auditable temperature half-width used when the raw topology bracket is numerically displaced. */
  public HydrocarbonWaterStableRegionBoundaryCorrector setMaximumExpansionHalfWidthK(
      double maximumExpansionHalfWidthK) {
    if (!positive(maximumExpansionHalfWidthK)) {
      throw new IllegalArgumentException("maximum correction expansion must be positive");
    }
    this.maximumExpansionHalfWidthK = maximumExpansionHalfWidthK;
    return this;
  }

  /**
   * Limits the number of successively wider intervals considered for one already-proven stable-region adjacency.
   */
  public HydrocarbonWaterStableRegionBoundaryCorrector setMaximumIntervalAttempts(int maximumIntervalAttempts) {
    if (maximumIntervalAttempts < 1) {
      throw new IllegalArgumentException("maximum interval-attempt count must be positive");
    }
    this.maximumIntervalAttempts = maximumIntervalAttempts;
    return this;
  }

  /**
   * Enables the expensive retained-GOW reverse seed as a diagnostic fallback after the ordinary two-phase seed fails.
   *
   * <p>
   * This is disabled by default because equivalent retained branches double the stability-search cost without adding
   * information. Model-homotopy seeding is the preferred recovery path for a proven TP-flash branch switch.
   * </p>
   */
  public HydrocarbonWaterStableRegionBoundaryCorrector setUseThreePhaseSideRetainedSeedFallback(boolean enabled) {
    this.useThreePhaseSideRetainedSeedFallback = enabled;
    return this;
  }

  /**
   * Uses retained-phase compositions from the stable 2P side and the third-phase composition from the stable GOW side.
   */
  public Result correct(TransitionBracket bracket) {
    if (bracket == null) {
      throw new IllegalArgumentException("a stable-region transition bracket is required");
    }
    BoundaryFamily family = bracket.getFamily();
    StableState retained = bracket.getTwoPhaseState();
    StableState threePhase = bracket.getThreePhaseState();
    CandidatePhase phaseZero = family.getRetainedPhaseZero();
    CandidatePhase phaseOne = family.getRetainedPhaseOne();
    CandidatePhase incipient = family.getIncipientPhase();
    if (!retained.hasPhase(phaseZero) || !retained.hasPhase(phaseOne) || retained.hasPhase(incipient)
        || !threePhase.hasPhase(phaseZero) || !threePhase.hasPhase(phaseOne) || !threePhase.hasPhase(incipient)) {
      throw new IllegalArgumentException("stable bracket phase inventory does not match its topology family");
    }
    double retainedFraction = retained.getPhaseFraction(phaseZero) + retained.getPhaseFraction(phaseOne);
    double beta = retained.getPhaseFraction(phaseZero) / retainedFraction;
    if (!(beta > 0.0 && beta < 1.0) || !Double.isFinite(beta)) {
      return Result.failure(bracket, "stable two-phase side has no finite retained-phase split");
    }
    double threePhaseRetainedFraction = threePhase.getPhaseFraction(phaseZero) + threePhase.getPhaseFraction(phaseOne);
    double threePhaseBeta = threePhase.getPhaseFraction(phaseZero) / threePhaseRetainedFraction;
    if (!(threePhaseBeta > 0.0 && threePhaseBeta < 1.0) || !Double.isFinite(threePhaseBeta)) {
      return Result.failure(bracket, "stable three-phase side has no finite retained-phase split");
    }
    List<double[]> incipientSeeds = Collections.singletonList(threePhase.getPhaseComposition(incipient));
    HydrocarbonWaterBoundaryEndpointClassifier classifier = new HydrocarbonWaterBoundaryEndpointClassifier(template)
        .setTolerances(1.0e-6, 1.0e-5, Math.max(tangentPlaneTolerance, 1.0e-8));
    SeedCorrection twoPhaseSeed = correctFromSeed(bracket, phaseZero, phaseOne, incipient, beta,
        retained.getPhaseComposition(phaseZero), retained.getPhaseComposition(phaseOne), incipientSeeds, classifier,
        "TWO_PHASE_SIDE");
    if (twoPhaseSeed.acceptedOrdinaryRootCount > 0) {
      return twoPhaseSeed.toResult(bracket);
    }
    if (!useThreePhaseSideRetainedSeedFallback) {
      return twoPhaseSeed.toResult(bracket);
    }

    // A stable multiphase TP solver can switch between disconnected local minima at the apparent
    // phase-count transition. In that case the retained phases on the two-phase side are not on the
    // branch that bounds the observed GOW state. Use the retained GOW compositions as an independent
    // reverse-continuation seed, but subject every resulting root to the same endpoint and global
    // stability classifiers. Numerical convergence from this fallback never bypasses a quality gate.
    SeedCorrection threePhaseSeed = correctFromSeed(bracket, phaseZero, phaseOne, incipient, threePhaseBeta,
        threePhase.getPhaseComposition(phaseZero), threePhase.getPhaseComposition(phaseOne), incipientSeeds, classifier,
        "THREE_PHASE_SIDE_RETAINED");
    if (threePhaseSeed.acceptedOrdinaryRootCount > 0) {
      return SeedCorrection.combine(twoPhaseSeed, threePhaseSeed).toResult(bracket);
    }
    return SeedCorrection.combine(twoPhaseSeed, threePhaseSeed).toResult(bracket);
  }

  private SeedCorrection correctFromSeed(TransitionBracket bracket, CandidatePhase phaseZero, CandidatePhase phaseOne,
      CandidatePhase incipient, double beta, double[] phaseZeroComposition, double[] phaseOneComposition,
      List<double[]> incipientSeeds, HydrocarbonWaterBoundaryEndpointClassifier classifier, String seedSource) {
    List<IntervalAttempt> attempts = new ArrayList<IntervalAttempt>();
    List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classificationCache = new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>();
    int classificationEvaluationCount = 0;
    int classificationReuseCount = 0;
    double midpoint = bracket.getMidpointTemperatureK();
    double halfWidth = 0.5 * bracket.getTemperatureWidthK();
    int expansionAttemptCount = 0;
    while (true) {
      double minimumTemperatureK = expansionAttemptCount == 0 ? bracket.getLowerTemperatureState().getTemperatureK()
          : Math.max(50.0, midpoint - halfWidth);
      double maximumTemperatureK = expansionAttemptCount == 0 ? bracket.getUpperTemperatureState().getTemperatureK()
          : Math.min(2500.0, midpoint + halfWidth);
      expansionAttemptCount++;
      TwoToThreePhaseBoundaryPointSolver.RootSet roots = new TwoToThreePhaseBoundaryPointSolver(template, phaseZero,
          phaseOne, incipient)
          .setNumericalControls(scanIntervals, maximumBisections, temperatureToleranceK, tangentPlaneTolerance)
          .solveAll(bracket.getPressureBara(), minimumTemperatureK, maximumTemperatureK, beta, phaseZeroComposition,
              phaseOneComposition, incipientSeeds);
      List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications = new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>();
      int accepted = 0;
      for (TwoToThreePhaseBoundaryPointSolver.Result root : roots.getRoots()) {
        HydrocarbonWaterBoundaryEndpointClassifier.Result classification = cachedClassification(root,
            classificationCache);
        if (classification == null) {
          classification = classifier.classify(root);
          classificationCache.add(classification);
          classificationEvaluationCount++;
        } else {
          classificationReuseCount++;
        }
        classifications.add(classification);
        accepted += classification.isOrdinaryBoundaryPoint() ? 1 : 0;
      }
      attempts.add(
          new IntervalAttempt(seedSource, minimumTemperatureK, maximumTemperatureK, roots, classifications, accepted));
      if (accepted > 0) {
        return new SeedCorrection(roots, classifications, accepted, attempts, null, classificationEvaluationCount,
            classificationReuseCount);
      }
      if (!classifications.isEmpty()) {
        return new SeedCorrection(roots, classifications, 0, attempts,
            seedSource + ": local root was rejected by the independent endpoint/global-stability gate",
            classificationEvaluationCount, classificationReuseCount);
      }
      if (expansionAttemptCount >= maximumIntervalAttempts) {
        String failure = roots.getFailureMessage() == null
            ? "stable topology bracket reached the controlled interval-attempt limit without a strict root"
            : roots.getFailureMessage();
        return new SeedCorrection(roots, classifications, 0, attempts, seedSource + ": " + failure,
            classificationEvaluationCount, classificationReuseCount);
      }
      if (halfWidth >= maximumExpansionHalfWidthK || minimumTemperatureK <= 50.0 && maximumTemperatureK >= 2500.0) {
        String failure = roots.getFailureMessage() == null
            ? "stable topology bracket produced no globally stable ordinary root"
            : roots.getFailureMessage();
        return new SeedCorrection(roots, classifications, 0, attempts, seedSource + ": " + failure,
            classificationEvaluationCount, classificationReuseCount);
      }
      halfWidth = Math.min(maximumExpansionHalfWidthK, Math.max(0.025, 4.0 * halfWidth));
    }
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  private static HydrocarbonWaterBoundaryEndpointClassifier.Result cachedClassification(
      TwoToThreePhaseBoundaryPointSolver.Result root, List<HydrocarbonWaterBoundaryEndpointClassifier.Result> cache) {
    for (HydrocarbonWaterBoundaryEndpointClassifier.Result cached : cache) {
      TwoToThreePhaseBoundaryPointSolver.Result previous = cached.getBoundaryRoot();
      if (previous.getRetainedPhaseZero() == root.getRetainedPhaseZero()
          && previous.getRetainedPhaseOne() == root.getRetainedPhaseOne()
          && previous.getIncipientPhase() == root.getIncipientPhase()
          && Math.abs(previous.getTemperatureK() - root.getTemperatureK()) <= 1.0e-5
          && Math.abs(previous.getPressureBara() - root.getPressureBara()) <= 1.0e-8
          && Math.abs(previous.getBeta() - root.getBeta()) <= 1.0e-7
          && compositionDistance(previous.getPhaseZeroComposition(), root.getPhaseZeroComposition()) <= 1.0e-6
          && compositionDistance(previous.getPhaseOneComposition(), root.getPhaseOneComposition()) <= 1.0e-6
          && compositionDistance(previous.getIncipientComposition(), root.getIncipientComposition()) <= 1.0e-6) {
        return cached;
      }
    }
    return null;
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

  private static final class SeedCorrection {
    private final TwoToThreePhaseBoundaryPointSolver.RootSet rootSet;
    private final List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications;
    private final int acceptedOrdinaryRootCount;
    private final List<IntervalAttempt> intervalAttempts;
    private final String failureMessage;
    private final int classificationEvaluationCount;
    private final int classificationReuseCount;

    private SeedCorrection(TwoToThreePhaseBoundaryPointSolver.RootSet rootSet,
        List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications, int acceptedOrdinaryRootCount,
        List<IntervalAttempt> intervalAttempts, String failureMessage, int classificationEvaluationCount,
        int classificationReuseCount) {
      this.rootSet = rootSet;
      this.classifications = new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>(classifications);
      this.acceptedOrdinaryRootCount = acceptedOrdinaryRootCount;
      this.intervalAttempts = new ArrayList<IntervalAttempt>(intervalAttempts);
      this.failureMessage = failureMessage;
      this.classificationEvaluationCount = classificationEvaluationCount;
      this.classificationReuseCount = classificationReuseCount;
    }

    private static SeedCorrection combine(SeedCorrection first, SeedCorrection second) {
      List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications = new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>();
      classifications.addAll(first.classifications);
      classifications.addAll(second.classifications);
      List<IntervalAttempt> attempts = new ArrayList<IntervalAttempt>();
      attempts.addAll(first.intervalAttempts);
      attempts.addAll(second.intervalAttempts);
      SeedCorrection accepted = second.acceptedOrdinaryRootCount > 0 ? second
          : first.acceptedOrdinaryRootCount > 0 ? first : second;
      String failure = accepted.acceptedOrdinaryRootCount > 0 ? null
          : joinFailures(first.failureMessage, second.failureMessage);
      return new SeedCorrection(accepted.rootSet, classifications,
          first.acceptedOrdinaryRootCount + second.acceptedOrdinaryRootCount, attempts, failure,
          first.classificationEvaluationCount + second.classificationEvaluationCount,
          first.classificationReuseCount + second.classificationReuseCount);
    }

    private Result toResult(TransitionBracket bracket) {
      return new Result(bracket, rootSet, classifications, acceptedOrdinaryRootCount, intervalAttempts, failureMessage,
          classificationEvaluationCount, classificationReuseCount);
    }

    private static String joinFailures(String first, String second) {
      if (first == null) {
        return second;
      }
      if (second == null) {
        return first;
      }
      return first + "; " + second;
    }
  }

  /** Immutable correction evidence for one stable topology bracket. */
  public static final class Result {
    private final TransitionBracket bracket;
    private final TwoToThreePhaseBoundaryPointSolver.RootSet rootSet;
    private final List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications;
    private final List<HydrocarbonWaterBoundaryEndpointClassifier.Result> allClassifications;
    private final int acceptedOrdinaryRootCount;
    private final List<IntervalAttempt> intervalAttempts;
    private final String failureMessage;
    private final int classificationEvaluationCount;
    private final int classificationReuseCount;

    private Result(TransitionBracket bracket, TwoToThreePhaseBoundaryPointSolver.RootSet rootSet,
        List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications, int acceptedOrdinaryRootCount,
        List<IntervalAttempt> intervalAttempts, String failureMessage, int classificationEvaluationCount,
        int classificationReuseCount) {
      this.bracket = bracket;
      this.rootSet = rootSet;
      this.classifications = Collections
          .unmodifiableList(new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>(classifications));
      this.allClassifications = Collections
          .unmodifiableList(distinctClassifications(intervalAttempts, classifications));
      this.acceptedOrdinaryRootCount = acceptedOrdinaryRootCount;
      this.intervalAttempts = Collections.unmodifiableList(new ArrayList<IntervalAttempt>(intervalAttempts));
      this.failureMessage = failureMessage;
      this.classificationEvaluationCount = classificationEvaluationCount;
      this.classificationReuseCount = classificationReuseCount;
    }

    private static Result failure(TransitionBracket bracket, String failureMessage) {
      return new Result(bracket, null, Collections.<HydrocarbonWaterBoundaryEndpointClassifier.Result>emptyList(), 0,
          Collections.<IntervalAttempt>emptyList(), failureMessage, 0, 0);
    }

    public boolean isConverged() {
      return acceptedOrdinaryRootCount > 0 && failureMessage == null;
    }

    public TransitionBracket getBracket() {
      return bracket;
    }

    public TwoToThreePhaseBoundaryPointSolver.RootSet getRootSet() {
      return rootSet;
    }

    public List<HydrocarbonWaterBoundaryEndpointClassifier.Result> getClassifications() {
      return classifications;
    }

    /** @return every distinct classification produced by all controlled interval expansions */
    public List<HydrocarbonWaterBoundaryEndpointClassifier.Result> getAllClassifications() {
      return allClassifications;
    }

    public int getAcceptedOrdinaryRootCount() {
      return acceptedOrdinaryRootCount;
    }

    public List<IntervalAttempt> getIntervalAttempts() {
      return intervalAttempts;
    }

    public String getFailureMessage() {
      return failureMessage;
    }

    public int getClassificationEvaluationCount() {
      return classificationEvaluationCount;
    }

    public int getClassificationReuseCount() {
      return classificationReuseCount;
    }

    private static List<HydrocarbonWaterBoundaryEndpointClassifier.Result> distinctClassifications(
        List<IntervalAttempt> attempts, List<HydrocarbonWaterBoundaryEndpointClassifier.Result> fallback) {
      List<HydrocarbonWaterBoundaryEndpointClassifier.Result> distinct = new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>();
      if (attempts.isEmpty()) {
        distinct.addAll(fallback);
        return distinct;
      }
      for (IntervalAttempt attempt : attempts) {
        for (HydrocarbonWaterBoundaryEndpointClassifier.Result candidate : attempt.getClassifications()) {
          boolean duplicate = false;
          for (HydrocarbonWaterBoundaryEndpointClassifier.Result existing : distinct) {
            TwoToThreePhaseBoundaryPointSolver.Result first = existing.getBoundaryRoot();
            TwoToThreePhaseBoundaryPointSolver.Result second = candidate.getBoundaryRoot();
            if (existing.getClassification() == candidate.getClassification()
                && Math.abs(first.getTemperatureK() - second.getTemperatureK()) <= 1.0e-7
                && Math.abs(first.getPressureBara() - second.getPressureBara()) <= 1.0e-7) {
              duplicate = true;
              break;
            }
          }
          if (!duplicate) {
            distinct.add(candidate);
          }
        }
      }
      return distinct;
    }
  }

  /** One preserved strict-root search interval used during controlled bracket expansion. */
  public static final class IntervalAttempt {
    private final String seedSource;
    private final double minimumTemperatureK;
    private final double maximumTemperatureK;
    private final TwoToThreePhaseBoundaryPointSolver.RootSet rootSet;
    private final List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications;
    private final int acceptedOrdinaryRootCount;

    private IntervalAttempt(String seedSource, double minimumTemperatureK, double maximumTemperatureK,
        TwoToThreePhaseBoundaryPointSolver.RootSet rootSet,
        List<HydrocarbonWaterBoundaryEndpointClassifier.Result> classifications, int acceptedOrdinaryRootCount) {
      this.seedSource = seedSource;
      this.minimumTemperatureK = minimumTemperatureK;
      this.maximumTemperatureK = maximumTemperatureK;
      this.rootSet = rootSet;
      this.classifications = Collections
          .unmodifiableList(new ArrayList<HydrocarbonWaterBoundaryEndpointClassifier.Result>(classifications));
      this.acceptedOrdinaryRootCount = acceptedOrdinaryRootCount;
    }

    public String getSeedSource() {
      return seedSource;
    }

    public double getMinimumTemperatureK() {
      return minimumTemperatureK;
    }

    public double getMaximumTemperatureK() {
      return maximumTemperatureK;
    }

    public TwoToThreePhaseBoundaryPointSolver.RootSet getRootSet() {
      return rootSet;
    }

    public List<HydrocarbonWaterBoundaryEndpointClassifier.Result> getClassifications() {
      return classifications;
    }

    public int getAcceptedOrdinaryRootCount() {
      return acceptedOrdinaryRootCount;
    }
  }
}
