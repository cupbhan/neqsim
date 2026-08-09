package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryEndpointClassifier.Classification;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket;

/**
 * Classifies an isolated strict boundary root from independently corrected local topology evidence.
 *
 * <p>
 * A finite-fraction zero-TPD root is not promoted to a curve merely because a Newton corrector converged. The
 * classifier scans a symmetric local PT neighborhood, requires stable 2P/GOW adjacency, strictly corrects every
 * matching bracket, reapplies the endpoint and global-stability gates, and finally checks phase-composition continuity
 * with the supplied root. It reports support for a branch or turning point but never performs continuation itself.
 */
public final class HydrocarbonWaterIsolatedBoundaryRootClassifier {
  private final SystemInterface template;
  private int pressureOffsetCountPerSide = 3;
  private double maximumLogPressureOffset = 0.12;
  private int temperatureProbeCount = 25;
  private double relativeTemperatureHalfWidth = 0.12;
  private double minimumTemperatureHalfWidthK = 15.0;
  private double maximumCompositionJump = 0.35;
  private double minimumSuccessfulGridFraction = 0.90;

  /** Creates a non-destructive local topology classifier. */
  public HydrocarbonWaterIsolatedBoundaryRootClassifier(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets the symmetric log-pressure and temperature neighborhood sampled around the root. */
  public HydrocarbonWaterIsolatedBoundaryRootClassifier setNeighborhoodControls(int pressureOffsetCountPerSide,
      double maximumLogPressureOffset, int temperatureProbeCount, double relativeTemperatureHalfWidth,
      double minimumTemperatureHalfWidthK) {
    if (pressureOffsetCountPerSide < 1 || !positive(maximumLogPressureOffset) || temperatureProbeCount < 5
        || temperatureProbeCount % 2 == 0 || !positive(relativeTemperatureHalfWidth)
        || !positive(minimumTemperatureHalfWidthK)) {
      throw new IllegalArgumentException("invalid isolated-root neighborhood controls");
    }
    this.pressureOffsetCountPerSide = pressureOffsetCountPerSide;
    this.maximumLogPressureOffset = maximumLogPressureOffset;
    this.temperatureProbeCount = temperatureProbeCount;
    this.relativeTemperatureHalfWidth = relativeTemperatureHalfWidth;
    this.minimumTemperatureHalfWidthK = minimumTemperatureHalfWidthK;
    return this;
  }

  /** Sets phase-composition continuity and minimum stable-grid-success gates. */
  public HydrocarbonWaterIsolatedBoundaryRootClassifier setQualityControls(double maximumCompositionJump,
      double minimumSuccessfulGridFraction) {
    if (!positive(maximumCompositionJump) || maximumCompositionJump > 2.0
        || !Double.isFinite(minimumSuccessfulGridFraction) || minimumSuccessfulGridFraction <= 0.0
        || minimumSuccessfulGridFraction > 1.0) {
      throw new IllegalArgumentException("invalid isolated-root quality controls");
    }
    this.maximumCompositionJump = maximumCompositionJump;
    this.minimumSuccessfulGridFraction = minimumSuccessfulGridFraction;
    return this;
  }

  /**
   * Audits one globally stable ordinary anchor without assuming that it belongs to a branch.
   *
   * @param anchor isolated strict anchor produced by the global discovery stage
   * @return immutable local topology classification and all corrected neighboring-root evidence
   */
  public Result classify(AnchorPoint anchor) {
    if (anchor == null) {
      throw new IllegalArgumentException("an isolated strict anchor is required");
    }
    TwoToThreePhaseBoundaryPointSolver.Result root = anchor.getBoundaryRoot();
    HydrocarbonWaterBoundaryEndpointClassifier.Result endpoint = new HydrocarbonWaterBoundaryEndpointClassifier(
        template).classify(root);
    if (!endpoint.isOrdinaryBoundaryPoint()) {
      RootTopology topology = endpoint.getClassification() == Classification.METASTABLE_BOUNDARY_ROOT
          ? RootTopology.METASTABLE_BOUNDARY_ROOT
          : RootTopology.NON_ORDINARY_SPECIAL_POINT;
      return Result.nonOrdinary(anchor.getFamily(), endpoint, topology);
    }

    double[] pressures = pressureProbes(root.getPressureBara());
    double[] temperatures = temperatureProbes(root.getTemperatureK());
    HydrocarbonWaterStableRegionTransitionScanner.Result scan = new HydrocarbonWaterStableRegionTransitionScanner(
        template).setNumericalControls(80, 0.05, 1.0e-9).scan(temperatures, pressures);
    List<BracketCorrection> corrections = new ArrayList<BracketCorrection>();
    List<TwoToThreePhaseBoundaryPointSolver.Result> matchingRoots = new ArrayList<TwoToThreePhaseBoundaryPointSolver.Result>();
    int matchingBracketCount = 0;
    int rejectedCompositionCount = 0;
    for (TransitionBracket bracket : scan.getBrackets()) {
      if (bracket.getFamily() != anchor.getFamily()) {
        continue;
      }
      matchingBracketCount++;
      HydrocarbonWaterStableRegionBoundaryCorrector.Result correction = new HydrocarbonWaterStableRegionBoundaryCorrector(
          template).setNumericalControls(24, 80, 1.0e-7, 1.0e-8).setMaximumExpansionHalfWidthK(5.0)
          .setMaximumIntervalAttempts(3).correct(bracket);
      int acceptedContinuousRoots = 0;
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result candidate : correction.getAllClassifications()) {
        if (!candidate.isOrdinaryBoundaryPoint()) {
          continue;
        }
        TwoToThreePhaseBoundaryPointSolver.Result candidateRoot = candidate.getBoundaryRoot();
        if (maximumPhaseCompositionJump(root, candidateRoot) <= maximumCompositionJump) {
          addDistinct(matchingRoots, candidateRoot);
          acceptedContinuousRoots++;
        } else {
          rejectedCompositionCount++;
        }
      }
      corrections.add(new BracketCorrection(bracket, correction, acceptedContinuousRoots));
    }

    double successfulGridFraction = (double) scan.getStates().size() / (pressures.length * temperatures.length);
    boolean lowerPressureSupport = hasPressureSupport(matchingRoots, root.getPressureBara(), -1);
    boolean higherPressureSupport = hasPressureSupport(matchingRoots, root.getPressureBara(), 1);
    boolean turningPointPattern = hasTurningPointPattern(matchingRoots, root);
    int nonCentralRootCount = countNonCentralRoots(matchingRoots, root);
    RootTopology topology;
    String diagnostic;
    if (successfulGridFraction < minimumSuccessfulGridFraction) {
      topology = RootTopology.INCONCLUSIVE;
      diagnostic = "local stable-flash grid did not meet the configured success fraction";
    } else if (lowerPressureSupport && higherPressureSupport) {
      topology = RootTopology.REGULAR_BRANCH_SUPPORTED;
      diagnostic = "independently corrected composition-continuous roots support both pressure directions";
    } else if (turningPointPattern) {
      topology = RootTopology.PRESSURE_TURNING_POINT_CANDIDATE;
      diagnostic = "two composition-continuous roots occur on one pressure side and straddle the anchor temperature";
    } else if (nonCentralRootCount > 0) {
      topology = RootTopology.ONE_SIDED_BRANCH_SUPPORTED;
      diagnostic = "independently corrected composition-continuous roots support only one pressure direction";
    } else {
      topology = RootTopology.ISOLATED_ORDINARY_ROOT;
      diagnostic = "zero-TPD root is ordinary and globally stable, but no independent neighboring strict root was found";
    }
    return new Result(anchor.getFamily(), endpoint, topology, scan, corrections, matchingRoots, matchingBracketCount,
        rejectedCompositionCount, successfulGridFraction, lowerPressureSupport, higherPressureSupport, diagnostic);
  }

  private double[] pressureProbes(double pressureBara) {
    double[] probes = new double[2 * pressureOffsetCountPerSide + 1];
    for (int index = -pressureOffsetCountPerSide; index <= pressureOffsetCountPerSide; index++) {
      double fraction = (double) index / pressureOffsetCountPerSide;
      probes[index + pressureOffsetCountPerSide] = pressureBara * Math.exp(fraction * maximumLogPressureOffset);
    }
    return probes;
  }

  private double[] temperatureProbes(double temperatureK) {
    double halfWidth = Math.max(minimumTemperatureHalfWidthK, relativeTemperatureHalfWidth * temperatureK);
    double minimum = Math.max(50.0, temperatureK - halfWidth);
    double maximum = Math.min(2500.0, temperatureK + halfWidth);
    double[] probes = new double[temperatureProbeCount];
    for (int index = 0; index < temperatureProbeCount; index++) {
      probes[index] = minimum + (maximum - minimum) * index / (temperatureProbeCount - 1);
    }
    return probes;
  }

  private static boolean hasPressureSupport(List<TwoToThreePhaseBoundaryPointSolver.Result> roots,
      double anchorPressureBara, int direction) {
    double tolerance = 1.0e-7 * Math.max(1.0, anchorPressureBara);
    for (TwoToThreePhaseBoundaryPointSolver.Result root : roots) {
      if (direction < 0 && root.getPressureBara() < anchorPressureBara - tolerance
          || direction > 0 && root.getPressureBara() > anchorPressureBara + tolerance) {
        return true;
      }
    }
    return false;
  }

  private static int countNonCentralRoots(List<TwoToThreePhaseBoundaryPointSolver.Result> roots,
      TwoToThreePhaseBoundaryPointSolver.Result anchor) {
    int count = 0;
    double pressureTolerance = 1.0e-7 * Math.max(1.0, anchor.getPressureBara());
    for (TwoToThreePhaseBoundaryPointSolver.Result root : roots) {
      if (Math.abs(root.getPressureBara() - anchor.getPressureBara()) > pressureTolerance) {
        count++;
      }
    }
    return count;
  }

  private static boolean hasTurningPointPattern(List<TwoToThreePhaseBoundaryPointSolver.Result> roots,
      TwoToThreePhaseBoundaryPointSolver.Result anchor) {
    double pressureTolerance = 1.0e-7 * Math.max(1.0, anchor.getPressureBara());
    double temperatureTolerance = 1.0e-5 * Math.max(1.0, anchor.getTemperatureK());
    for (int firstIndex = 0; firstIndex < roots.size(); firstIndex++) {
      TwoToThreePhaseBoundaryPointSolver.Result first = roots.get(firstIndex);
      if (Math.abs(first.getPressureBara() - anchor.getPressureBara()) <= pressureTolerance) {
        continue;
      }
      for (int secondIndex = firstIndex + 1; secondIndex < roots.size(); secondIndex++) {
        TwoToThreePhaseBoundaryPointSolver.Result second = roots.get(secondIndex);
        if (Math.abs(first.getPressureBara() - second.getPressureBara()) > pressureTolerance
            || Math.signum(first.getPressureBara() - anchor.getPressureBara()) != Math
                .signum(second.getPressureBara() - anchor.getPressureBara())) {
          continue;
        }
        boolean firstBelow = first.getTemperatureK() < anchor.getTemperatureK() - temperatureTolerance;
        boolean firstAbove = first.getTemperatureK() > anchor.getTemperatureK() + temperatureTolerance;
        boolean secondBelow = second.getTemperatureK() < anchor.getTemperatureK() - temperatureTolerance;
        boolean secondAbove = second.getTemperatureK() > anchor.getTemperatureK() + temperatureTolerance;
        if (firstBelow && secondAbove || firstAbove && secondBelow) {
          return true;
        }
      }
    }
    return false;
  }

  private static double maximumPhaseCompositionJump(TwoToThreePhaseBoundaryPointSolver.Result first,
      TwoToThreePhaseBoundaryPointSolver.Result second) {
    if (first.getRetainedPhaseZero() != second.getRetainedPhaseZero()
        || first.getRetainedPhaseOne() != second.getRetainedPhaseOne()
        || first.getIncipientPhase() != second.getIncipientPhase()) {
      return Double.POSITIVE_INFINITY;
    }
    return Math.max(compositionDistance(first.getPhaseZeroComposition(), second.getPhaseZeroComposition()),
        Math.max(compositionDistance(first.getPhaseOneComposition(), second.getPhaseOneComposition()),
            compositionDistance(first.getIncipientComposition(), second.getIncipientComposition())));
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

  private static void addDistinct(List<TwoToThreePhaseBoundaryPointSolver.Result> roots,
      TwoToThreePhaseBoundaryPointSolver.Result candidate) {
    for (TwoToThreePhaseBoundaryPointSolver.Result existing : roots) {
      double relativeTemperature = Math.abs(existing.getTemperatureK() - candidate.getTemperatureK())
          / Math.max(1.0, candidate.getTemperatureK());
      double logPressure = Math.abs(Math.log(existing.getPressureBara() / candidate.getPressureBara()));
      if (relativeTemperature <= 1.0e-7 && logPressure <= 1.0e-7
          && maximumPhaseCompositionJump(existing, candidate) <= 1.0e-6) {
        return;
      }
    }
    roots.add(candidate);
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Mutually exclusive evidence-level outcomes; none automatically promotes a boundary curve. */
  public enum RootTopology {
    REGULAR_BRANCH_SUPPORTED, ONE_SIDED_BRANCH_SUPPORTED, PRESSURE_TURNING_POINT_CANDIDATE, ISOLATED_ORDINARY_ROOT,
    NON_ORDINARY_SPECIAL_POINT, METASTABLE_BOUNDARY_ROOT, INCONCLUSIVE
  }

  /** One stable-region bracket and its independent strict-correction audit. */
  public static final class BracketCorrection {
    private final TransitionBracket bracket;
    private final HydrocarbonWaterStableRegionBoundaryCorrector.Result correction;
    private final int acceptedCompositionContinuousRootCount;

    private BracketCorrection(TransitionBracket bracket,
        HydrocarbonWaterStableRegionBoundaryCorrector.Result correction, int acceptedCompositionContinuousRootCount) {
      this.bracket = bracket;
      this.correction = correction;
      this.acceptedCompositionContinuousRootCount = acceptedCompositionContinuousRootCount;
    }

    public TransitionBracket getBracket() {
      return bracket;
    }

    public HydrocarbonWaterStableRegionBoundaryCorrector.Result getCorrection() {
      return correction;
    }

    public int getAcceptedCompositionContinuousRootCount() {
      return acceptedCompositionContinuousRootCount;
    }
  }

  /** Immutable local-topology evidence for one isolated strict root. */
  public static final class Result {
    private final BoundaryFamily family;
    private final HydrocarbonWaterBoundaryEndpointClassifier.Result endpointClassification;
    private final RootTopology topology;
    private final HydrocarbonWaterStableRegionTransitionScanner.Result stableRegionScan;
    private final List<BracketCorrection> bracketCorrections;
    private final List<TwoToThreePhaseBoundaryPointSolver.Result> matchingRoots;
    private final int matchingBracketCount;
    private final int rejectedCompositionCount;
    private final double successfulGridFraction;
    private final boolean lowerPressureSupport;
    private final boolean higherPressureSupport;
    private final String diagnostic;

    private Result(BoundaryFamily family, HydrocarbonWaterBoundaryEndpointClassifier.Result endpointClassification,
        RootTopology topology, HydrocarbonWaterStableRegionTransitionScanner.Result stableRegionScan,
        List<BracketCorrection> bracketCorrections, List<TwoToThreePhaseBoundaryPointSolver.Result> matchingRoots,
        int matchingBracketCount, int rejectedCompositionCount, double successfulGridFraction,
        boolean lowerPressureSupport, boolean higherPressureSupport, String diagnostic) {
      this.family = family;
      this.endpointClassification = endpointClassification;
      this.topology = topology;
      this.stableRegionScan = stableRegionScan;
      this.bracketCorrections = Collections.unmodifiableList(new ArrayList<BracketCorrection>(bracketCorrections));
      this.matchingRoots = Collections
          .unmodifiableList(new ArrayList<TwoToThreePhaseBoundaryPointSolver.Result>(matchingRoots));
      this.matchingBracketCount = matchingBracketCount;
      this.rejectedCompositionCount = rejectedCompositionCount;
      this.successfulGridFraction = successfulGridFraction;
      this.lowerPressureSupport = lowerPressureSupport;
      this.higherPressureSupport = higherPressureSupport;
      this.diagnostic = diagnostic;
    }

    private static Result nonOrdinary(BoundaryFamily family,
        HydrocarbonWaterBoundaryEndpointClassifier.Result endpointClassification, RootTopology topology) {
      return new Result(family, endpointClassification, topology, null, Collections.<BracketCorrection>emptyList(),
          Collections.<TwoToThreePhaseBoundaryPointSolver.Result>emptyList(), 0, 0, 0.0, false, false,
          endpointClassification.getFailureMessage());
    }

    public BoundaryFamily getFamily() {
      return family;
    }

    public HydrocarbonWaterBoundaryEndpointClassifier.Result getEndpointClassification() {
      return endpointClassification;
    }

    public RootTopology getTopology() {
      return topology;
    }

    public HydrocarbonWaterStableRegionTransitionScanner.Result getStableRegionScan() {
      return stableRegionScan;
    }

    public List<BracketCorrection> getBracketCorrections() {
      return bracketCorrections;
    }

    public List<TwoToThreePhaseBoundaryPointSolver.Result> getMatchingRoots() {
      return matchingRoots;
    }

    public int getMatchingBracketCount() {
      return matchingBracketCount;
    }

    public int getRejectedCompositionCount() {
      return rejectedCompositionCount;
    }

    public double getSuccessfulGridFraction() {
      return successfulGridFraction;
    }

    public boolean hasLowerPressureSupport() {
      return lowerPressureSupport;
    }

    public boolean hasHigherPressureSupport() {
      return higherPressureSupport;
    }

    public String getDiagnostic() {
      return diagnostic;
    }

    /** @return true when independent local evidence proves at least one non-central branch point */
    public boolean hasIndependentBranchSupport() {
      return topology == RootTopology.REGULAR_BRANCH_SUPPORTED || topology == RootTopology.ONE_SIDED_BRANCH_SUPPORTED
          || topology == RootTopology.PRESSURE_TURNING_POINT_CANDIDATE;
    }
  }
}
