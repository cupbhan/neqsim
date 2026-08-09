package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBranchSeedScanner.Seed;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.BoundaryDefinition;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Discovers and clusters rigorously corrected pressure-anchor points for all hydrocarbon-water 2P-to-3P boundaries. */
public final class HydrocarbonWaterBoundaryAnchorDiscoverer {
  private static final double MINIMUM_ORDINARY_PHASE_DISTANCE = 1.0e-5;
  private final SystemInterface template;
  private int scanIntervals = 80;
  private int maximumBisections = 80;
  private double temperatureToleranceK = 1.0e-5;
  private double tangentPlaneTolerance = 1.0e-8;
  private int stableScanMaximumBisections = 80;
  private double stableScanTemperatureToleranceK = 1.0e-4;
  private int maximumCorrectionIntervalAttempts = Integer.MAX_VALUE;
  private double maximumCompositionJump = 0.35;
  private double maximumLogPressureJump = 1.5;
  private double maximumRelativeTemperatureJump = 0.5;

  /** Creates a non-destructive global anchor discoverer. */
  public HydrocarbonWaterBoundaryAnchorDiscoverer(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets fixed-pressure correction controls applied identically to every topology family. */
  public HydrocarbonWaterBoundaryAnchorDiscoverer setCorrectionControls(int scanIntervals, int maximumBisections,
      double temperatureToleranceK, double tangentPlaneTolerance) {
    if (scanIntervals < 2 || maximumBisections < 1 || !Double.isFinite(temperatureToleranceK)
        || temperatureToleranceK <= 0.0 || !Double.isFinite(tangentPlaneTolerance) || tangentPlaneTolerance <= 0.0) {
      throw new IllegalArgumentException("invalid boundary-anchor correction controls");
    }
    this.scanIntervals = scanIntervals;
    this.maximumBisections = maximumBisections;
    this.temperatureToleranceK = temperatureToleranceK;
    this.tangentPlaneTolerance = tangentPlaneTolerance;
    return this;
  }

  /**
   * Sets the deliberately looser topology-bracket refinement controls used before strict zero-TPD correction.
   *
   * <p>
   * A stable-region bracket only proves which two phase regions are adjacent. It is not an engineering boundary point;
   * the subsequent strict corrector recomputes the zero-TPD root to its own tolerance. Keeping these controls separate
   * avoids spending most of a characterized-fluid audit refining a TP phase-count switch to a precision that is then
   * discarded.
   * </p>
   */
  public HydrocarbonWaterBoundaryAnchorDiscoverer setStableScanControls(int maximumBisections,
      double temperatureToleranceK) {
    if (maximumBisections < 1 || !Double.isFinite(temperatureToleranceK) || temperatureToleranceK <= 0.0) {
      throw new IllegalArgumentException("invalid stable-region scan controls");
    }
    this.stableScanMaximumBisections = maximumBisections;
    this.stableScanTemperatureToleranceK = temperatureToleranceK;
    return this;
  }

  /** Sets the controlled expansion-attempt limit used after each stable-region transition bracket is proven. */
  public HydrocarbonWaterBoundaryAnchorDiscoverer setMaximumCorrectionIntervalAttempts(
      int maximumCorrectionIntervalAttempts) {
    if (maximumCorrectionIntervalAttempts < 1) {
      throw new IllegalArgumentException("maximum correction interval-attempt count must be positive");
    }
    this.maximumCorrectionIntervalAttempts = maximumCorrectionIntervalAttempts;
    return this;
  }

  /** Sets the one-to-one composition/PT continuity thresholds used only to cluster corrected anchors. */
  public HydrocarbonWaterBoundaryAnchorDiscoverer setClusteringControls(double maximumCompositionJump,
      double maximumLogPressureJump, double maximumRelativeTemperatureJump) {
    if (!Double.isFinite(maximumCompositionJump) || maximumCompositionJump <= 0.0
        || !Double.isFinite(maximumLogPressureJump) || maximumLogPressureJump <= 0.0
        || !Double.isFinite(maximumRelativeTemperatureJump) || maximumRelativeTemperatureJump <= 0.0) {
      throw new IllegalArgumentException("invalid boundary-anchor clustering controls");
    }
    this.maximumCompositionJump = maximumCompositionJump;
    this.maximumLogPressureJump = maximumLogPressureJump;
    this.maximumRelativeTemperatureJump = maximumRelativeTemperatureJump;
    return this;
  }

  /** Discovers all three physical two-phase-to-GOW topology families. */
  public Result discover(double[] seedTemperaturesK, double[] anchorPressuresBara, double minimumTemperatureK,
      double maximumTemperatureK) {
    return discover(seedTemperaturesK, anchorPressuresBara, minimumTemperatureK, maximumTemperatureK,
        EnumSet.allOf(BoundaryFamily.class));
  }

  /**
   * Discovers physical anchors only from explicitly adjacent stable 2P/GOW topology brackets.
   *
   * <p>
   * This is the preferred path for high-water heavy-fluid systems. The stable TP flash first proves which two regions
   * meet. The narrow bracket is then corrected by the same specified-two-phase plus zero-third-phase-TPD equations and
   * must pass the independent global stability gate before it can become an ordinary continuation anchor.
   * </p>
   */
  public StableDiscoveryResult discoverFromStableRegionTransitions(double[] temperaturesK,
      double[] anchorPressuresBara) {
    validateInputs(temperaturesK, anchorPressuresBara, minimum(temperaturesK), maximum(temperaturesK),
        EnumSet.allOf(BoundaryFamily.class));
    long stableScanStarted = System.currentTimeMillis();
    HydrocarbonWaterStableRegionTransitionScanner.Result stableScan = new HydrocarbonWaterStableRegionTransitionScanner(
        template).setNumericalControls(stableScanMaximumBisections, stableScanTemperatureToleranceK, 1.0e-9)
        .scan(temperaturesK, anchorPressuresBara);
    long stableScanTimeMs = System.currentTimeMillis() - stableScanStarted;
    HydrocarbonWaterStableRegionBoundaryCorrector corrector = new HydrocarbonWaterStableRegionBoundaryCorrector(
        template)
        .setNumericalControls(Math.max(8, scanIntervals / 4), maximumBisections,
            Math.min(temperatureToleranceK, 1.0e-7), tangentPlaneTolerance)
        .setMaximumIntervalAttempts(maximumCorrectionIntervalAttempts);
    List<AnchorPoint> correctedAnchors = new ArrayList<AnchorPoint>();
    List<EndpointCandidate> endpointCandidates = new ArrayList<EndpointCandidate>();
    List<CorrectionAttempt> attempts = new ArrayList<CorrectionAttempt>();
    List<HydrocarbonWaterStableRegionBoundaryCorrector.Result> corrections = new ArrayList<HydrocarbonWaterStableRegionBoundaryCorrector.Result>();
    List<Long> correctionTimesMs = new ArrayList<Long>();
    long correctionStarted = System.currentTimeMillis();
    for (HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket : stableScan.getBrackets()) {
      long oneCorrectionStarted = System.currentTimeMillis();
      HydrocarbonWaterStableRegionBoundaryCorrector.Result correction = corrector.correct(bracket);
      correctionTimesMs.add(System.currentTimeMillis() - oneCorrectionStarted);
      corrections.add(correction);
      int ordinary = 0;
      int withheld = 0;
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : correction.getAllClassifications()) {
        if (classification.isOrdinaryBoundaryPoint()) {
          correctedAnchors.add(new AnchorPoint(bracket.getFamily(), classification));
          ordinary++;
        } else {
          endpointCandidates.add(new EndpointCandidate(bracket.getFamily(), classification));
          withheld++;
        }
      }
      attempts.add(CorrectionAttempt.fromStableBracket(bracket, correction, ordinary, withheld));
    }
    List<Branch> branches = cluster(correctedAnchors);
    int preservedAnchors = branches.stream().mapToInt(branch -> branch.points.size()).sum();
    if (preservedAnchors != correctedAnchors.size()) {
      throw new IllegalStateException("stable-transition boundary-anchor preservation failed");
    }
    Result discovery = new Result(branches, endpointCandidates, attempts, stableScan.getFlashEvaluations(),
        correctedAnchors.size());
    return new StableDiscoveryResult(discovery, stableScan, corrections, stableScanTimeMs,
        System.currentTimeMillis() - correctionStarted, correctionTimesMs);
  }

  /**
   * Discovers selected topology families without ever promoting an uncorrected TP-grid point to a boundary anchor.
   *
   * <p>
   * At each pressure the TP grid supplies same-state phase-composition seeds. Each family is then recomputed by a
   * specified retained two-phase flash plus a zero-TPD third-phase condition. All strict roots are preserved and
   * clustered one-to-one across pressure anchors. An unmatched root starts a new branch, so isolated branches remain in
   * the result instead of being silently discarded.
   * </p>
   */
  public Result discover(double[] seedTemperaturesK, double[] anchorPressuresBara, double minimumTemperatureK,
      double maximumTemperatureK, Set<BoundaryFamily> requestedFamilies) {
    validateInputs(seedTemperaturesK, anchorPressuresBara, minimumTemperatureK, maximumTemperatureK, requestedFamilies);
    List<AnchorPoint> correctedAnchors = new ArrayList<AnchorPoint>();
    List<EndpointCandidate> endpointCandidates = new ArrayList<EndpointCandidate>();
    List<CorrectionAttempt> attempts = new ArrayList<CorrectionAttempt>();
    HydrocarbonWaterBoundaryEndpointClassifier endpointClassifier = new HydrocarbonWaterBoundaryEndpointClassifier(
        template).setTolerances(1.0e-6, MINIMUM_ORDINARY_PHASE_DISTANCE, tangentPlaneTolerance);
    int gridEvaluations = 0;
    for (BoundaryFamily family : BoundaryFamily.values()) {
      if (!requestedFamilies.contains(family)) {
        continue;
      }
      HydrocarbonWaterBoundarySeedCorrector corrector = new HydrocarbonWaterBoundarySeedCorrector(template,
          family.retainedPhaseZero, family.retainedPhaseOne, family.incipientPhase)
          .setNumericalControls(scanIntervals, maximumBisections, temperatureToleranceK, tangentPlaneTolerance);
      for (double pressureBara : sortedUnique(anchorPressuresBara)) {
        HydrocarbonWaterBranchSeedScanner.Result scan = new HydrocarbonWaterBranchSeedScanner(template,
            family.incipientPhase).scan(seedTemperaturesK, new double[] { pressureBara });
        gridEvaluations += scan.getEvaluatedPointCount();
        List<Seed> threePhaseSeeds = new ArrayList<Seed>();
        for (Seed seed : scan.getSeeds()) {
          if (seed.hasPhase(family.retainedPhaseZero) && seed.hasPhase(family.retainedPhaseOne)
              && seed.hasPhase(family.incipientPhase)) {
            threePhaseSeeds.add(seed);
          }
        }
        if (threePhaseSeeds.isEmpty()) {
          attempts.add(CorrectionAttempt.noSeed(family, pressureBara, scan.getFailures().size()));
          continue;
        }
        TwoToThreePhaseBoundaryPointSolver.RootSet roots = corrector.correctAtSeedPressure(threePhaseSeeds,
            minimumTemperatureK, maximumTemperatureK);
        int ordinaryAnchorCount = 0;
        int endpointCandidateCount = 0;
        for (TwoToThreePhaseBoundaryPointSolver.Result root : roots.getRoots()) {
          HydrocarbonWaterBoundaryEndpointClassifier.Result classification = endpointClassifier.classify(root);
          if (classification.isOrdinaryBoundaryPoint()) {
            correctedAnchors.add(new AnchorPoint(family, classification));
            ordinaryAnchorCount++;
          } else {
            endpointCandidates.add(new EndpointCandidate(family, classification));
            endpointCandidateCount++;
          }
        }
        attempts.add(CorrectionAttempt.from(family, pressureBara, threePhaseSeeds.size(), scan.getFailures().size(),
            roots, ordinaryAnchorCount, endpointCandidateCount));
      }
    }
    List<Branch> branches = cluster(correctedAnchors);
    int preservedAnchors = 0;
    for (Branch branch : branches) {
      preservedAnchors += branch.points.size();
    }
    if (preservedAnchors != correctedAnchors.size()) {
      throw new IllegalStateException("corrected boundary-anchor preservation failed");
    }
    return new Result(branches, endpointCandidates, attempts, gridEvaluations, correctedAnchors.size());
  }

  /** Clusters independently corrected roots while preserving every input anchor exactly once. */
  public List<Branch> clusterCorrectedAnchors(List<AnchorPoint> anchors) {
    if (anchors == null) {
      throw new IllegalArgumentException("corrected anchor list is required");
    }
    for (AnchorPoint anchor : anchors) {
      if (anchor == null) {
        throw new IllegalArgumentException("corrected anchor list must not contain null entries");
      }
    }
    List<Branch> branches = cluster(new ArrayList<AnchorPoint>(anchors));
    int preserved = 0;
    for (Branch branch : branches) {
      preserved += branch.points.size();
    }
    if (preserved != anchors.size()) {
      throw new IllegalStateException("corrected boundary-anchor preservation failed");
    }
    return branches;
  }

  /**
   * Re-clusters a discovery result after an independent local-topology audit supplies additional globally stable
   * ordinary anchors.
   *
   * <p>
   * Only strict {@link AnchorPoint} instances are accepted. Existing anchors, endpoint candidates, correction attempts,
   * and discovery evaluation counts are preserved; PT/composition duplicates are removed before the complete inventory
   * is clustered again.
   */
  public Result augmentWithIndependentAnchors(Result discovery, List<AnchorPoint> additionalAnchors) {
    if (discovery == null || additionalAnchors == null) {
      throw new IllegalArgumentException("discovery result and independent anchors are required");
    }
    List<AnchorPoint> merged = new ArrayList<AnchorPoint>();
    for (Branch branch : discovery.getBranches()) {
      merged.addAll(branch.getPoints());
    }
    for (AnchorPoint additional : additionalAnchors) {
      if (additional == null) {
        throw new IllegalArgumentException("independent anchor inventory must not contain null entries");
      }
      if (!containsEquivalentAnchor(merged, additional)) {
        merged.add(additional);
      }
    }
    List<Branch> branches = clusterCorrectedAnchors(merged);
    return new Result(branches, discovery.endpointCandidates, discovery.attempts, discovery.gridEvaluationCount,
        merged.size());
  }

  private static boolean containsEquivalentAnchor(List<AnchorPoint> anchors, AnchorPoint candidate) {
    for (AnchorPoint existing : anchors) {
      double relativeTemperature = Math.abs(existing.temperatureK - candidate.temperatureK)
          / Math.max(1.0, candidate.temperatureK);
      double logPressure = Math.abs(Math.log(existing.pressureBara / candidate.pressureBara));
      if (existing.family == candidate.family && relativeTemperature <= 1.0e-7 && logPressure <= 1.0e-7
          && compositionJump(existing, candidate) <= 1.0e-6) {
        return true;
      }
    }
    return false;
  }

  private List<Branch> cluster(List<AnchorPoint> anchors) {
    List<MutableBranch> mutable = new ArrayList<MutableBranch>();
    for (BoundaryFamily family : BoundaryFamily.values()) {
      List<AnchorPoint> familyAnchors = new ArrayList<AnchorPoint>();
      for (AnchorPoint anchor : anchors) {
        if (anchor.family == family) {
          familyAnchors.add(anchor);
        }
      }
      familyAnchors.sort(
          Comparator.comparingDouble(AnchorPoint::getPressureBara).thenComparingDouble(AnchorPoint::getTemperatureK));
      int start = 0;
      while (start < familyAnchors.size()) {
        double pressureBara = familyAnchors.get(start).pressureBara;
        int end = start + 1;
        while (end < familyAnchors.size()
            && Math.abs(familyAnchors.get(end).pressureBara - pressureBara) <= 1.0e-10 * Math.max(1.0, pressureBara)) {
          end++;
        }
        assignPressureLevel(family, familyAnchors.subList(start, end), mutable);
        start = end;
      }
    }
    List<Branch> frozen = new ArrayList<Branch>();
    int familyBranchIndex = 0;
    BoundaryFamily previousFamily = null;
    for (MutableBranch branch : mutable) {
      if (branch.family != previousFamily) {
        previousFamily = branch.family;
        familyBranchIndex = 0;
      }
      frozen.add(new Branch(branch.family.getCode() + "#" + familyBranchIndex++, branch.family, branch.points));
    }
    return frozen;
  }

  private void assignPressureLevel(BoundaryFamily family, List<AnchorPoint> current, List<MutableBranch> branches) {
    List<Link> links = new ArrayList<Link>();
    for (int branchIndex = 0; branchIndex < branches.size(); branchIndex++) {
      MutableBranch branch = branches.get(branchIndex);
      if (branch.family != family || branch.points.isEmpty()) {
        continue;
      }
      AnchorPoint previous = branch.points.get(branch.points.size() - 1);
      if (!(previous.pressureBara < current.get(0).pressureBara)) {
        continue;
      }
      for (int pointIndex = 0; pointIndex < current.size(); pointIndex++) {
        AnchorPoint candidate = current.get(pointIndex);
        double compositionJump = compositionJump(previous, candidate);
        double logPressureJump = Math.abs(Math.log(candidate.pressureBara / previous.pressureBara));
        double relativeTemperatureJump = Math.abs(candidate.temperatureK - previous.temperatureK)
            / Math.max(1.0, Math.min(candidate.temperatureK, previous.temperatureK));
        if (compositionJump <= maximumCompositionJump && logPressureJump <= maximumLogPressureJump
            && relativeTemperatureJump <= maximumRelativeTemperatureJump) {
          links.add(new Link(branchIndex, pointIndex,
              compositionJump + 0.1 * logPressureJump + 0.1 * relativeTemperatureJump));
        }
      }
    }
    links.sort(Comparator.comparingDouble(link -> link.score));
    Set<Integer> assignedBranches = new HashSet<Integer>();
    Set<Integer> assignedPoints = new HashSet<Integer>();
    for (Link link : links) {
      if (assignedBranches.add(link.branchIndex) && assignedPoints.add(link.pointIndex)) {
        branches.get(link.branchIndex).points.add(current.get(link.pointIndex));
      }
    }
    for (int pointIndex = 0; pointIndex < current.size(); pointIndex++) {
      if (!assignedPoints.contains(pointIndex)) {
        MutableBranch branch = new MutableBranch(family);
        branch.points.add(current.get(pointIndex));
        branches.add(branch);
      }
    }
  }

  private static double compositionJump(AnchorPoint first, AnchorPoint second) {
    return Math.max(compositionDistance(first.phaseZeroComposition, second.phaseZeroComposition),
        Math.max(compositionDistance(first.phaseOneComposition, second.phaseOneComposition),
            compositionDistance(first.incipientComposition, second.incipientComposition)));
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }

  private static double[] sortedUnique(double[] values) {
    double[] sorted = values.clone();
    java.util.Arrays.sort(sorted);
    List<Double> unique = new ArrayList<Double>();
    for (double value : sorted) {
      if (unique.isEmpty() || Math.abs(value - unique.get(unique.size() - 1)) > 1.0e-12 * Math.max(1.0, value)) {
        unique.add(value);
      }
    }
    double[] result = new double[unique.size()];
    for (int index = 0; index < result.length; index++) {
      result[index] = unique.get(index);
    }
    return result;
  }

  private static double minimum(double[] values) {
    if (values == null || values.length == 0) {
      return Double.NaN;
    }
    double minimum = Double.POSITIVE_INFINITY;
    for (double value : values) {
      minimum = Math.min(minimum, value);
    }
    return minimum;
  }

  private static double maximum(double[] values) {
    if (values == null || values.length == 0) {
      return Double.NaN;
    }
    double maximum = Double.NEGATIVE_INFINITY;
    for (double value : values) {
      maximum = Math.max(maximum, value);
    }
    return maximum;
  }

  private static void validateInputs(double[] temperaturesK, double[] pressuresBara, double minimumTemperatureK,
      double maximumTemperatureK, Set<BoundaryFamily> requestedFamilies) {
    if (temperaturesK == null || temperaturesK.length == 0 || pressuresBara == null || pressuresBara.length == 0
        || requestedFamilies == null || requestedFamilies.isEmpty()) {
      throw new IllegalArgumentException("seed grids and at least one boundary family are required");
    }
    if (!Double.isFinite(minimumTemperatureK) || !Double.isFinite(maximumTemperatureK) || minimumTemperatureK < 50.0
        || maximumTemperatureK <= minimumTemperatureK) {
      throw new IllegalArgumentException("invalid boundary correction temperature interval");
    }
    for (BoundaryFamily family : requestedFamilies) {
      if (family == null) {
        throw new IllegalArgumentException("requested boundary family must not be null");
      }
    }
  }

  /** The three possible two-phase-to-GOW boundary families. */
  public enum BoundaryFamily {
    GO_TO_GOW(CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS,
        new BoundaryDefinition(Region.GAS_OIL, Region.GAS_OIL_AQUEOUS)),
    GW_TO_GOW(CandidatePhase.GAS, CandidatePhase.AQUEOUS, CandidatePhase.OIL,
        new BoundaryDefinition(Region.GAS_AQUEOUS, Region.GAS_OIL_AQUEOUS)),
    OW_TO_GOW(CandidatePhase.OIL, CandidatePhase.AQUEOUS, CandidatePhase.GAS,
        new BoundaryDefinition(Region.OIL_AQUEOUS, Region.GAS_OIL_AQUEOUS));

    private final CandidatePhase retainedPhaseZero;
    private final CandidatePhase retainedPhaseOne;
    private final CandidatePhase incipientPhase;
    private final BoundaryDefinition definition;

    BoundaryFamily(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase,
        BoundaryDefinition definition) {
      this.retainedPhaseZero = retainedPhaseZero;
      this.retainedPhaseOne = retainedPhaseOne;
      this.incipientPhase = incipientPhase;
      this.definition = definition;
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

    public BoundaryDefinition getDefinition() {
      return definition;
    }

    public String getCode() {
      return definition.getCode();
    }
  }

  /** One strict fixed-pressure boundary root; no raw TP-grid coordinates are retained here. */
  public static final class AnchorPoint {
    private final BoundaryFamily family;
    private final TwoToThreePhaseBoundaryPointSolver.Result boundaryRoot;
    private final double temperatureK;
    private final double pressureBara;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double[] incipientComposition;
    private final double retainedFlashResidual;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;
    private final HydrocarbonWaterBoundaryGlobalStabilityGate.Result globalStabilityResult;

    private AnchorPoint(BoundaryFamily family, HydrocarbonWaterBoundaryEndpointClassifier.Result classification) {
      if (classification == null || !classification.isOrdinaryBoundaryPoint()
          || classification.getGlobalStabilityResult() == null
          || !classification.getGlobalStabilityResult().isAccepted()) {
        throw new IllegalArgumentException("a globally stable ordinary boundary classification is required");
      }
      TwoToThreePhaseBoundaryPointSolver.Result root = classification.getBoundaryRoot();
      if (family == null || root.getRetainedPhaseZero() != family.retainedPhaseZero
          || root.getRetainedPhaseOne() != family.retainedPhaseOne
          || root.getIncipientPhase() != family.incipientPhase) {
        throw new IllegalArgumentException("strict boundary root does not match the requested topology family");
      }
      double minimumPhaseDistance = Math.min(
          compositionDistance(root.getPhaseZeroComposition(), root.getPhaseOneComposition()),
          Math.min(compositionDistance(root.getPhaseZeroComposition(), root.getIncipientComposition()),
              compositionDistance(root.getPhaseOneComposition(), root.getIncipientComposition())));
      if (minimumPhaseDistance <= MINIMUM_ORDINARY_PHASE_DISTANCE) {
        throw new IllegalArgumentException(
            "a duplicate-phase endpoint candidate cannot become an ordinary boundary anchor");
      }
      this.family = family;
      this.boundaryRoot = root;
      this.temperatureK = root.getTemperatureK();
      this.pressureBara = root.getPressureBara();
      this.beta = root.getBeta();
      this.phaseZeroComposition = root.getPhaseZeroComposition();
      this.phaseOneComposition = root.getPhaseOneComposition();
      this.incipientComposition = root.getIncipientComposition();
      this.retainedFlashResidual = root.getRetainedFlashResidual();
      this.tangentPlaneDistance = root.getTangentPlaneDistance();
      this.stationarityResidual = root.getStationarityResidual();
      this.globalStabilityResult = classification.getGlobalStabilityResult();
    }

    /** Independently classifies a fixed-pressure root before converting it to a globally clusterable anchor. */
    public static AnchorPoint from(BoundaryFamily family, TwoToThreePhaseBoundaryPointSolver.Result root,
        SystemInterface template) {
      return new AnchorPoint(family, new HydrocarbonWaterBoundaryEndpointClassifier(template).classify(root));
    }

    /** Creates an anchor without repeating a classification already produced by an independent strict corrector. */
    public static AnchorPoint fromClassification(BoundaryFamily family,
        HydrocarbonWaterBoundaryEndpointClassifier.Result classification) {
      return new AnchorPoint(family, classification);
    }

    public BoundaryFamily getFamily() {
      return family;
    }

    /** @return immutable strict root from which this anchor was accepted */
    public TwoToThreePhaseBoundaryPointSolver.Result getBoundaryRoot() {
      return boundaryRoot;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double getBeta() {
      return beta;
    }

    public double[] getPhaseZeroComposition() {
      return phaseZeroComposition.clone();
    }

    public double[] getPhaseOneComposition() {
      return phaseOneComposition.clone();
    }

    public double[] getIncipientComposition() {
      return incipientComposition.clone();
    }

    public double getRetainedFlashResidual() {
      return retainedFlashResidual;
    }

    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    public double getStationarityResidual() {
      return stationarityResidual;
    }

    public HydrocarbonWaterBoundaryGlobalStabilityGate.Result getGlobalStabilityResult() {
      return globalStabilityResult;
    }

    /** Converts this strict root to immutable quality-gate evidence. */
    public TwoToThreePhaseBoundaryQualityGate.EvidencePoint toEvidencePoint() {
      return new TwoToThreePhaseBoundaryQualityGate.EvidencePoint(temperatureK, pressureBara, beta,
          phaseZeroComposition, phaseOneComposition, incipientComposition, retainedFlashResidual, tangentPlaneDistance,
          stationarityResidual);
    }

    /** Converts this anchor to a pseudo-arclength continuation seed. */
    public TwoToThreePhaseArcLengthCorrector.State toContinuationState() {
      return TwoToThreePhaseArcLengthCorrector.State.create(family.retainedPhaseZero, family.retainedPhaseOne,
          family.incipientPhase, temperatureK, pressureBara, beta, phaseZeroComposition, phaseOneComposition,
          incipientComposition);
    }
  }

  /** Corrected zero-TPD root withheld from ordinary continuation pending or carrying endpoint classification. */
  public static final class EndpointCandidate {
    private final BoundaryFamily family;
    private final HydrocarbonWaterBoundaryEndpointClassifier.Result classification;

    private EndpointCandidate(BoundaryFamily family, HydrocarbonWaterBoundaryEndpointClassifier.Result classification) {
      if (family == null || classification == null || classification.isOrdinaryBoundaryPoint()) {
        throw new IllegalArgumentException("a non-ordinary endpoint classification is required");
      }
      this.family = family;
      this.classification = classification;
    }

    public BoundaryFamily getFamily() {
      return family;
    }

    public HydrocarbonWaterBoundaryEndpointClassifier.Result getClassification() {
      return classification;
    }
  }

  /** Ordered corrected anchors tentatively belonging to one composition-continuous branch. */
  public static final class Branch {
    private final String identifier;
    private final BoundaryFamily family;
    private final List<AnchorPoint> points;

    private Branch(String identifier, BoundaryFamily family, List<AnchorPoint> points) {
      this.identifier = identifier;
      this.family = family;
      this.points = Collections.unmodifiableList(new ArrayList<AnchorPoint>(points));
    }

    public String getIdentifier() {
      return identifier;
    }

    public BoundaryFamily getFamily() {
      return family;
    }

    public List<AnchorPoint> getPoints() {
      return points;
    }

    public boolean isIsolated() {
      return points.size() == 1;
    }

    public boolean canSeedContinuation() {
      return points.size() >= 2;
    }
  }

  /** Auditable result of one family/pressure seed and correction attempt. */
  public static final class CorrectionAttempt {
    private final BoundaryFamily family;
    private final double pressureBara;
    private final int threePhaseSeedCount;
    private final int gridFailureCount;
    private final int correctedRootCount;
    private final int ordinaryAnchorCount;
    private final int endpointCandidateCount;
    private final int flashEvaluations;
    private final int stabilityEvaluations;
    private final String failureMessage;

    private CorrectionAttempt(BoundaryFamily family, double pressureBara, int threePhaseSeedCount, int gridFailureCount,
        int correctedRootCount, int ordinaryAnchorCount, int endpointCandidateCount, int flashEvaluations,
        int stabilityEvaluations, String failureMessage) {
      this.family = family;
      this.pressureBara = pressureBara;
      this.threePhaseSeedCount = threePhaseSeedCount;
      this.gridFailureCount = gridFailureCount;
      this.correctedRootCount = correctedRootCount;
      this.ordinaryAnchorCount = ordinaryAnchorCount;
      this.endpointCandidateCount = endpointCandidateCount;
      this.flashEvaluations = flashEvaluations;
      this.stabilityEvaluations = stabilityEvaluations;
      this.failureMessage = failureMessage;
    }

    private static CorrectionAttempt noSeed(BoundaryFamily family, double pressureBara, int gridFailureCount) {
      return new CorrectionAttempt(family, pressureBara, 0, gridFailureCount, 0, 0, 0, 0, 0,
          "NO_THREE_PHASE_GRID_SEED");
    }

    private static CorrectionAttempt from(BoundaryFamily family, double pressureBara, int seedCount,
        int gridFailureCount, TwoToThreePhaseBoundaryPointSolver.RootSet roots, int ordinaryAnchorCount,
        int endpointCandidateCount) {
      return new CorrectionAttempt(family, pressureBara, seedCount, gridFailureCount, roots.getRoots().size(),
          ordinaryAnchorCount, endpointCandidateCount, roots.getFlashEvaluations(), roots.getStabilityEvaluations(),
          roots.getFailureMessage());
    }

    private static CorrectionAttempt fromStableBracket(
        HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket,
        HydrocarbonWaterStableRegionBoundaryCorrector.Result correction, int ordinaryAnchorCount,
        int endpointCandidateCount) {
      TwoToThreePhaseBoundaryPointSolver.RootSet roots = correction.getRootSet();
      return new CorrectionAttempt(bracket.getFamily(), bracket.getPressureBara(), 1, 0,
          roots == null ? 0 : roots.getRoots().size(), ordinaryAnchorCount, endpointCandidateCount,
          roots == null ? 0 : roots.getFlashEvaluations(), roots == null ? 0 : roots.getStabilityEvaluations(),
          correction.getFailureMessage());
    }

    public BoundaryFamily getFamily() {
      return family;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public int getThreePhaseSeedCount() {
      return threePhaseSeedCount;
    }

    public int getGridFailureCount() {
      return gridFailureCount;
    }

    public int getCorrectedRootCount() {
      return correctedRootCount;
    }

    public int getOrdinaryAnchorCount() {
      return ordinaryAnchorCount;
    }

    public int getEndpointCandidateCount() {
      return endpointCandidateCount;
    }

    public int getFlashEvaluations() {
      return flashEvaluations;
    }

    public int getStabilityEvaluations() {
      return stabilityEvaluations;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }

  /** Immutable global discovery result preserving successful, failed, and isolated topology evidence. */
  public static final class Result {
    private final List<Branch> branches;
    private final List<EndpointCandidate> endpointCandidates;
    private final List<CorrectionAttempt> attempts;
    private final int gridEvaluationCount;
    private final int correctedAnchorCount;

    private Result(List<Branch> branches, List<EndpointCandidate> endpointCandidates, List<CorrectionAttempt> attempts,
        int gridEvaluationCount, int correctedAnchorCount) {
      this.branches = Collections.unmodifiableList(new ArrayList<Branch>(branches));
      this.endpointCandidates = Collections.unmodifiableList(new ArrayList<EndpointCandidate>(endpointCandidates));
      this.attempts = Collections.unmodifiableList(new ArrayList<CorrectionAttempt>(attempts));
      this.gridEvaluationCount = gridEvaluationCount;
      this.correctedAnchorCount = correctedAnchorCount;
    }

    public List<Branch> getBranches() {
      return branches;
    }

    public List<CorrectionAttempt> getAttempts() {
      return attempts;
    }

    /** @return corrected special-point candidates excluded from ordinary branch clustering */
    public List<EndpointCandidate> getEndpointCandidates() {
      return endpointCandidates;
    }

    public int getEndpointCandidateCount() {
      return endpointCandidates.size();
    }

    public int getGridEvaluationCount() {
      return gridEvaluationCount;
    }

    public int getCorrectedAnchorCount() {
      return correctedAnchorCount;
    }

    public int getIsolatedBranchCount() {
      int isolated = 0;
      for (Branch branch : branches) {
        if (branch.isIsolated()) {
          isolated++;
        }
      }
      return isolated;
    }
  }

  /** Stable-region scan, strict corrections, and ordinary discovery kept together for full auditability. */
  public static final class StableDiscoveryResult {
    private final Result discovery;
    private final HydrocarbonWaterStableRegionTransitionScanner.Result stableRegionScan;
    private final List<HydrocarbonWaterStableRegionBoundaryCorrector.Result> corrections;
    private final long stableScanTimeMs;
    private final long correctionTimeMs;
    private final List<Long> correctionTimesMs;

    private StableDiscoveryResult(Result discovery,
        HydrocarbonWaterStableRegionTransitionScanner.Result stableRegionScan,
        List<HydrocarbonWaterStableRegionBoundaryCorrector.Result> corrections, long stableScanTimeMs,
        long correctionTimeMs, List<Long> correctionTimesMs) {
      this.discovery = discovery;
      this.stableRegionScan = stableRegionScan;
      this.corrections = Collections
          .unmodifiableList(new ArrayList<HydrocarbonWaterStableRegionBoundaryCorrector.Result>(corrections));
      this.stableScanTimeMs = stableScanTimeMs;
      this.correctionTimeMs = correctionTimeMs;
      this.correctionTimesMs = Collections.unmodifiableList(new ArrayList<Long>(correctionTimesMs));
    }

    public Result getDiscovery() {
      return discovery;
    }

    public HydrocarbonWaterStableRegionTransitionScanner.Result getStableRegionScan() {
      return stableRegionScan;
    }

    public List<HydrocarbonWaterStableRegionBoundaryCorrector.Result> getCorrections() {
      return corrections;
    }

    public long getStableScanTimeMs() {
      return stableScanTimeMs;
    }

    public long getCorrectionTimeMs() {
      return correctionTimeMs;
    }

    public List<Long> getCorrectionTimesMs() {
      return correctionTimesMs;
    }
  }

  private static final class MutableBranch {
    private final BoundaryFamily family;
    private final List<AnchorPoint> points = new ArrayList<AnchorPoint>();

    private MutableBranch(BoundaryFamily family) {
      this.family = family;
    }
  }

  private static final class Link {
    private final int branchIndex;
    private final int pointIndex;
    private final double score;

    private Link(int branchIndex, int pointIndex, double score) {
      this.branchIndex = branchIndex;
      this.pointIndex = pointIndex;
      this.score = score;
    }
  }
}
