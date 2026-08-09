package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket;

/**
 * Continues a regular hydrocarbon-water boundary with log pressure as the local parameter.
 *
 * <p>
 * Every accepted point is independently rediscovered as a stable 2P/GOW adjacency, corrected with the specified
 * retained-two-phase plus zero-third-phase-TPD equations, and passed through the endpoint/global-stability classifier.
 * Pressure is only a regular local parameter; repeated step collapse is reported as a pressure-parameter limit so the
 * caller can switch to pseudo-arclength.
 */
public final class HydrocarbonWaterRegularPressureBoundaryTracer {
  private final SystemInterface template;
  private final BoundaryFamily family;
  private double initialLogPressureStep = 0.04;
  private double minimumLogPressureStep = 0.0025;
  private double maximumLogPressureStep = 0.12;
  private int maximumRetriesPerPoint = 6;
  private int temperatureProbeCount = 11;
  private double minimumTemperatureHalfWidthK = 8.0;
  private double relativeTemperatureHalfWidth = 0.04;
  private double maximumCompositionJump = 0.35;
  private double maximumPredictionErrorK = 30.0;
  private boolean domainBoundsEnabled;
  private double minimumTemperatureK;
  private double maximumTemperatureK;
  private double minimumPressureBara;
  private double maximumPressureBara;

  /** Creates a non-destructive regular-pressure tracer for one explicit topology family. */
  public HydrocarbonWaterRegularPressureBoundaryTracer(SystemInterface template, BoundaryFamily family) {
    if (template == null || family == null) {
      throw new IllegalArgumentException("thermodynamic template and boundary family are required");
    }
    this.template = template.clone();
    this.family = family;
  }

  /** Sets adaptive logarithmic pressure-step and retry controls. */
  public HydrocarbonWaterRegularPressureBoundaryTracer setStepControls(double initialLogPressureStep,
      double minimumLogPressureStep, double maximumLogPressureStep, int maximumRetriesPerPoint) {
    if (!positive(initialLogPressureStep) || !positive(minimumLogPressureStep) || !positive(maximumLogPressureStep)
        || minimumLogPressureStep > initialLogPressureStep || initialLogPressureStep > maximumLogPressureStep
        || maximumRetriesPerPoint < 1) {
      throw new IllegalArgumentException("invalid regular-pressure continuation controls");
    }
    this.initialLogPressureStep = initialLogPressureStep;
    this.minimumLogPressureStep = minimumLogPressureStep;
    this.maximumLogPressureStep = maximumLogPressureStep;
    this.maximumRetriesPerPoint = maximumRetriesPerPoint;
    return this;
  }

  /** Sets the local temperature scan and branch-continuity gates. */
  public HydrocarbonWaterRegularPressureBoundaryTracer setCorrectionControls(int temperatureProbeCount,
      double minimumTemperatureHalfWidthK, double relativeTemperatureHalfWidth, double maximumCompositionJump,
      double maximumPredictionErrorK) {
    if (temperatureProbeCount < 5 || temperatureProbeCount % 2 == 0 || !positive(minimumTemperatureHalfWidthK)
        || !positive(relativeTemperatureHalfWidth) || !positive(maximumCompositionJump)
        || !positive(maximumPredictionErrorK)) {
      throw new IllegalArgumentException("invalid regular-pressure correction controls");
    }
    this.temperatureProbeCount = temperatureProbeCount;
    this.minimumTemperatureHalfWidthK = minimumTemperatureHalfWidthK;
    this.relativeTemperatureHalfWidth = relativeTemperatureHalfWidth;
    this.maximumCompositionJump = maximumCompositionJump;
    this.maximumPredictionErrorK = maximumPredictionErrorK;
    return this;
  }

  /** Sets the declared PT domain; the first accepted crossing point is retained as evidence. */
  public HydrocarbonWaterRegularPressureBoundaryTracer setDomainBounds(double minimumTemperatureK,
      double maximumTemperatureK, double minimumPressureBara, double maximumPressureBara) {
    if (!Double.isFinite(minimumTemperatureK) || !Double.isFinite(maximumTemperatureK)
        || maximumTemperatureK <= minimumTemperatureK || !positive(minimumPressureBara)
        || !positive(maximumPressureBara) || maximumPressureBara <= minimumPressureBara) {
      throw new IllegalArgumentException("invalid regular-pressure PT domain bounds");
    }
    this.minimumTemperatureK = minimumTemperatureK;
    this.maximumTemperatureK = maximumTemperatureK;
    this.minimumPressureBara = minimumPressureBara;
    this.maximumPressureBara = maximumPressureBara;
    this.domainBoundsEnabled = true;
    return this;
  }

  /** Continues from two ordered strict anchors in their pressure direction. */
  public Result trace(AnchorPoint previous, AnchorPoint current, int additionalPointCount) {
    validateSeeds(previous, current, additionalPointCount);
    List<AnchorPoint> points = new ArrayList<AnchorPoint>();
    points.add(previous);
    points.add(current);
    List<Attempt> attempts = new ArrayList<Attempt>();
    double direction = Math.signum(Math.log(current.getPressureBara() / previous.getPressureBara()));
    double secantStep = Math.abs(Math.log(current.getPressureBara() / previous.getPressureBara()));
    double logPressureStep = Math.max(minimumLogPressureStep,
        Math.min(maximumLogPressureStep, Math.min(initialLogPressureStep, secantStep)));
    TerminationReason terminationReason = null;
    String failureMessage = null;

    while (points.size() < additionalPointCount + 2) {
      AnchorPoint first = points.get(points.size() - 2);
      AnchorPoint second = points.get(points.size() - 1);
      boolean accepted = false;
      for (int retry = 0; retry < maximumRetriesPerPoint; retry++) {
        double targetLogPressure = Math.log(second.getPressureBara()) + direction * logPressureStep;
        double targetPressure = Math.exp(targetLogPressure);
        if (domainBoundsEnabled) {
          targetPressure = Math.max(minimumPressureBara, Math.min(maximumPressureBara, targetPressure));
          targetLogPressure = Math.log(targetPressure);
        }
        if (Math.abs(targetLogPressure - Math.log(second.getPressureBara())) < 1.0e-12) {
          terminationReason = TerminationReason.DOMAIN_EXIT;
          break;
        }
        double predictedTemperature = predictTemperature(first, second, targetLogPressure);
        Attempt attempt = correct(second, targetPressure, predictedTemperature, logPressureStep);
        attempts.add(attempt);
        if (attempt.acceptedAnchor != null) {
          points.add(attempt.acceptedAnchor);
          accepted = true;
          double predictionError = Math.abs(attempt.acceptedAnchor.getTemperatureK() - predictedTemperature);
          double compositionJump = compositionJump(second, attempt.acceptedAnchor);
          if (predictionError < 0.25 * maximumPredictionErrorK && compositionJump < 0.25 * maximumCompositionJump) {
            logPressureStep = Math.min(maximumLogPressureStep, 1.25 * logPressureStep);
          }
          if (domainBoundsEnabled && outsideDomain(attempt.acceptedAnchor)) {
            terminationReason = TerminationReason.DOMAIN_EXIT;
          }
          break;
        }
        failureMessage = attempt.failureMessage;
        logPressureStep *= 0.5;
        if (logPressureStep < minimumLogPressureStep) {
          terminationReason = TerminationReason.PRESSURE_PARAMETER_LIMIT;
          break;
        }
      }
      if (terminationReason != null || !accepted) {
        if (terminationReason == null) {
          terminationReason = TerminationReason.CORRECTION_RETRY_LIMIT;
        }
        break;
      }
    }
    boolean completed = points.size() == additionalPointCount + 2;
    if (completed) {
      terminationReason = TerminationReason.REQUESTED_POINT_COUNT;
      failureMessage = null;
    }
    return new Result(points, attempts, completed, logPressureStep, terminationReason, failureMessage);
  }

  private Attempt correct(AnchorPoint current, double pressureBara, double predictedTemperatureK,
      double logPressureStep) {
    double halfWidth = Math.max(minimumTemperatureHalfWidthK,
        relativeTemperatureHalfWidth * Math.max(50.0, predictedTemperatureK));
    double lower = Math.max(50.0, predictedTemperatureK - halfWidth);
    double upper = Math.min(2500.0, predictedTemperatureK + halfWidth);
    double[] temperatures = linearGrid(lower, upper, temperatureProbeCount);
    HydrocarbonWaterStableRegionTransitionScanner.Result scan = new HydrocarbonWaterStableRegionTransitionScanner(
        template).setNumericalControls(60, 0.05, 1.0e-9).scan(temperatures, new double[] { pressureBara });
    AnchorPoint best = null;
    double bestScore = Double.POSITIVE_INFINITY;
    int matchingBracketCount = 0;
    int strictRootCount = 0;
    String failure = null;
    for (TransitionBracket bracket : scan.getBrackets()) {
      if (bracket.getFamily() != family) {
        continue;
      }
      matchingBracketCount++;
      HydrocarbonWaterStableRegionBoundaryCorrector.Result correction = new HydrocarbonWaterStableRegionBoundaryCorrector(
          template).setNumericalControls(24, 80, 1.0e-7, 1.0e-8).setMaximumExpansionHalfWidthK(5.0)
          .setMaximumIntervalAttempts(3).correct(bracket);
      failure = correction.getFailureMessage();
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : correction.getAllClassifications()) {
        if (!classification.isOrdinaryBoundaryPoint()) {
          continue;
        }
        strictRootCount++;
        AnchorPoint candidate = AnchorPoint.fromClassification(family, classification);
        double compositionJump = compositionJump(current, candidate);
        double predictionError = Math.abs(candidate.getTemperatureK() - predictedTemperatureK);
        if (compositionJump > maximumCompositionJump || predictionError > maximumPredictionErrorK) {
          continue;
        }
        double score = compositionJump + predictionError / Math.max(1.0, halfWidth);
        if (score < bestScore) {
          best = candidate;
          bestScore = score;
        }
      }
    }
    if (best == null && failure == null) {
      failure = matchingBracketCount == 0 ? "no matching stable-region adjacency at the trial pressure"
          : "matching stable-region brackets produced no composition-continuous ordinary root";
    }
    return new Attempt(pressureBara, predictedTemperatureK, logPressureStep, scan.getFlashEvaluations(),
        scan.getFailures().size(), matchingBracketCount, strictRootCount, best, best == null ? failure : null);
  }

  private void validateSeeds(AnchorPoint previous, AnchorPoint current, int additionalPointCount) {
    if (previous == null || current == null || previous.getFamily() != family || current.getFamily() != family
        || additionalPointCount < 1 || previous.getPressureBara() == current.getPressureBara()) {
      throw new IllegalArgumentException("two ordered family-matching anchors and a positive point count are required");
    }
  }

  private boolean outsideDomain(AnchorPoint point) {
    return point.getTemperatureK() < minimumTemperatureK || point.getTemperatureK() > maximumTemperatureK
        || point.getPressureBara() < minimumPressureBara || point.getPressureBara() > maximumPressureBara;
  }

  private static double predictTemperature(AnchorPoint previous, AnchorPoint current, double targetLogPressure) {
    double previousLogPressure = Math.log(previous.getPressureBara());
    double currentLogPressure = Math.log(current.getPressureBara());
    double denominator = currentLogPressure - previousLogPressure;
    return current.getTemperatureK() + (current.getTemperatureK() - previous.getTemperatureK()) / denominator
        * (targetLogPressure - currentLogPressure);
  }

  private static double[] linearGrid(double minimum, double maximum, int count) {
    double[] grid = new double[count];
    for (int index = 0; index < count; index++) {
      grid[index] = minimum + (maximum - minimum) * index / (count - 1);
    }
    return grid;
  }

  private static double compositionJump(AnchorPoint first, AnchorPoint second) {
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

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Explicit reason why regular-pressure continuation stopped. */
  public enum TerminationReason {
    REQUESTED_POINT_COUNT, DOMAIN_EXIT, PRESSURE_PARAMETER_LIMIT, CORRECTION_RETRY_LIMIT
  }

  /** One auditable fixed-pressure topology scan and strict correction attempt. */
  public static final class Attempt {
    private final double pressureBara;
    private final double predictedTemperatureK;
    private final double logPressureStep;
    private final int stableFlashEvaluations;
    private final int stableScanFailures;
    private final int matchingBracketCount;
    private final int strictOrdinaryRootCount;
    private final AnchorPoint acceptedAnchor;
    private final String failureMessage;

    private Attempt(double pressureBara, double predictedTemperatureK, double logPressureStep,
        int stableFlashEvaluations, int stableScanFailures, int matchingBracketCount, int strictOrdinaryRootCount,
        AnchorPoint acceptedAnchor, String failureMessage) {
      this.pressureBara = pressureBara;
      this.predictedTemperatureK = predictedTemperatureK;
      this.logPressureStep = logPressureStep;
      this.stableFlashEvaluations = stableFlashEvaluations;
      this.stableScanFailures = stableScanFailures;
      this.matchingBracketCount = matchingBracketCount;
      this.strictOrdinaryRootCount = strictOrdinaryRootCount;
      this.acceptedAnchor = acceptedAnchor;
      this.failureMessage = failureMessage;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double getPredictedTemperatureK() {
      return predictedTemperatureK;
    }

    public double getLogPressureStep() {
      return logPressureStep;
    }

    public int getStableFlashEvaluations() {
      return stableFlashEvaluations;
    }

    public int getStableScanFailures() {
      return stableScanFailures;
    }

    public int getMatchingBracketCount() {
      return matchingBracketCount;
    }

    public int getStrictOrdinaryRootCount() {
      return strictOrdinaryRootCount;
    }

    public AnchorPoint getAcceptedAnchor() {
      return acceptedAnchor;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }

  /** Immutable regular-pressure continuation evidence. */
  public static final class Result {
    private final List<AnchorPoint> points;
    private final List<Attempt> attempts;
    private final boolean completedRequestedPoints;
    private final double finalLogPressureStep;
    private final TerminationReason terminationReason;
    private final String failureMessage;

    private Result(List<AnchorPoint> points, List<Attempt> attempts, boolean completedRequestedPoints,
        double finalLogPressureStep, TerminationReason terminationReason, String failureMessage) {
      this.points = Collections.unmodifiableList(new ArrayList<AnchorPoint>(points));
      this.attempts = Collections.unmodifiableList(new ArrayList<Attempt>(attempts));
      this.completedRequestedPoints = completedRequestedPoints;
      this.finalLogPressureStep = finalLogPressureStep;
      this.terminationReason = terminationReason;
      this.failureMessage = failureMessage;
    }

    public List<AnchorPoint> getPoints() {
      return points;
    }

    public List<Attempt> getAttempts() {
      return attempts;
    }

    public boolean hasCompletedRequestedPoints() {
      return completedRequestedPoints;
    }

    public double getFinalLogPressureStep() {
      return finalLogPressureStep;
    }

    public TerminationReason getTerminationReason() {
      return terminationReason;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
