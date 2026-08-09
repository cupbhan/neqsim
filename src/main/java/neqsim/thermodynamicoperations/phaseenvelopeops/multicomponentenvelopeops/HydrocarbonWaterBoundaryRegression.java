package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.BoundaryDefinition;

/** Same-topology, same-pressure regression for non-reactive hydrocarbon-water PT boundaries. */
public final class HydrocarbonWaterBoundaryRegression {
  private static final double PRESSURE_INTERSECTION_TOLERANCE = 1.0e-10;

  /** Development-reference provenance. Production calculations must not require either source. */
  public enum ReferenceSource {
    PAPER, PVTSIM, NEQSIM
  }

  /** Compares one candidate curve with one independently labelled reference curve. */
  public Result compare(Curve reference, Curve candidate, AcceptanceCriteria criteria) {
    if (reference == null || candidate == null) {
      throw new IllegalArgumentException("reference and candidate curves are required");
    }
    List<String> violations = contractViolations(reference, candidate);
    if (!violations.isEmpty()) {
      return Result.ineligible(violations, criteria != null);
    }

    DirectionalMetrics referenceToCandidate = directional(reference.getPoints(), candidate.getPoints());
    DirectionalMetrics candidateToReference = directional(candidate.getPoints(), reference.getPoints());
    CriticalPointMetrics criticalPointMetrics = criticalPointMetrics(reference.getCriticalPoint(),
        candidate.getCriticalPoint());
    Metrics metrics = new Metrics(referenceToCandidate, candidateToReference, criticalPointMetrics);

    if (criteria == null) {
      return Result.pending(metrics);
    }
    addThresholdViolations(metrics, criteria, violations);
    return Result.evaluated(metrics, violations);
  }

  private static List<String> contractViolations(Curve reference, Curve candidate) {
    List<String> violations = new ArrayList<String>();
    if (!reference.isTopologyResolved()) {
      violations.add("REFERENCE_TOPOLOGY_UNRESOLVED");
    }
    if (!candidate.isTopologyResolved()) {
      violations.add("CANDIDATE_TOPOLOGY_UNRESOLVED");
    }
    if (!reference.getBoundaryDefinition().getCode().equals(candidate.getBoundaryDefinition().getCode())) {
      violations.add("BOUNDARY_TOPOLOGY_MISMATCH");
    }
    ModelFingerprint referenceModel = reference.getModelFingerprint();
    ModelFingerprint candidateModel = candidate.getModelFingerprint();
    if (!sameToken(referenceModel.getEquationOfState(), candidateModel.getEquationOfState())) {
      violations.add("EQUATION_OF_STATE_MISMATCH");
    }
    if (!sameToken(referenceModel.getVolumeTranslation(), candidateModel.getVolumeTranslation())) {
      violations.add("VOLUME_TRANSLATION_MISMATCH");
    }
    if (!sameToken(referenceModel.getMixingRule(), candidateModel.getMixingRule())) {
      violations.add("MIXING_RULE_MISMATCH");
    }
    if (!referenceModel.isParameterFingerprintResolved() || !candidateModel.isParameterFingerprintResolved()) {
      violations.add("PARAMETER_FINGERPRINT_UNRESOLVED");
    } else if (!referenceModel.getParameterFingerprint().equals(candidateModel.getParameterFingerprint())) {
      violations.add("PARAMETER_FINGERPRINT_MISMATCH");
    }
    if (!referenceModel.isCompositionFingerprintResolved() || !candidateModel.isCompositionFingerprintResolved()) {
      violations.add("COMPOSITION_FINGERPRINT_UNRESOLVED");
    } else if (!referenceModel.getCompositionFingerprint().equals(candidateModel.getCompositionFingerprint())) {
      violations.add("COMPOSITION_FINGERPRINT_MISMATCH");
    }
    if (uniquePoints(reference.getPoints()).size() < 2) {
      violations.add("REFERENCE_CURVE_INSUFFICIENT");
    }
    if (uniquePoints(candidate.getPoints()).size() < 2) {
      violations.add("CANDIDATE_CURVE_INSUFFICIENT");
    }
    return violations;
  }

  private static void addThresholdViolations(Metrics metrics, AcceptanceCriteria criteria, List<String> violations) {
    if (metrics.getReferenceCoverageFraction() < criteria.getMinimumReferenceCoverageFraction()) {
      violations.add("REFERENCE_PRESSURE_COVERAGE_BELOW_THRESHOLD");
    }
    if (metrics.getCandidateCoverageFraction() < criteria.getMinimumCandidateCoverageFraction()) {
      violations.add("CANDIDATE_PRESSURE_COVERAGE_BELOW_THRESHOLD");
    }
    if (metrics.getMaximumRelativeTemperatureDifference() > criteria.getMaximumRelativeTemperatureDifference()) {
      violations.add("MAXIMUM_SAME_PRESSURE_TEMPERATURE_DIFFERENCE_ABOVE_THRESHOLD");
    }
    if (metrics.getRmsRelativeTemperatureDifference() > criteria.getMaximumRmsRelativeTemperatureDifference()) {
      violations.add("RMS_SAME_PRESSURE_TEMPERATURE_DIFFERENCE_ABOVE_THRESHOLD");
    }
    if (metrics.getP95RelativeTemperatureDifference() > criteria.getMaximumP95RelativeTemperatureDifference()) {
      violations.add("P95_SAME_PRESSURE_TEMPERATURE_DIFFERENCE_ABOVE_THRESHOLD");
    }
    if (metrics.getSymmetricMaximumRelativeTemperatureDifference() > criteria
        .getMaximumSymmetricRelativeTemperatureDifference()) {
      violations.add("SYMMETRIC_SAME_PRESSURE_TEMPERATURE_DIFFERENCE_ABOVE_THRESHOLD");
    }
    CriticalPointMetrics critical = metrics.getCriticalPointMetrics();
    if (criteria.isCriticalPointRequired() && !critical.isAvailable()) {
      violations.add("CRITICAL_POINT_COMPARISON_MISSING");
    } else if (critical.isAvailable()) {
      if (critical.getRelativeTemperatureDifference() > criteria.getMaximumCriticalTemperatureDifference()) {
        violations.add("CRITICAL_TEMPERATURE_DIFFERENCE_ABOVE_THRESHOLD");
      }
      if (critical.getRelativePressureDifference() > criteria.getMaximumCriticalPressureDifference()) {
        violations.add("CRITICAL_PRESSURE_DIFFERENCE_ABOVE_THRESHOLD");
      }
    }
  }

  private static DirectionalMetrics directional(List<Point> source, List<Point> target) {
    List<Point> sourcePoints = uniquePoints(source);
    List<Point> targetPoints = uniquePoints(target);
    List<Double> absoluteErrors = new ArrayList<Double>();
    List<Double> relativeErrors = new ArrayList<Double>();
    for (Point point : sourcePoints) {
      double matchedTemperature = closestTemperatureAtPressure(targetPoints, point.getPressureBara(),
          point.getTemperatureK());
      if (!Double.isFinite(matchedTemperature)) {
        continue;
      }
      double absolute = Math.abs(matchedTemperature - point.getTemperatureK());
      absoluteErrors.add(absolute);
      relativeErrors.add(absolute / point.getTemperatureK());
    }
    return new DirectionalMetrics(sourcePoints.size(), absoluteErrors, relativeErrors);
  }

  private static double closestTemperatureAtPressure(List<Point> curve, double pressureBara,
      double targetTemperatureK) {
    double bestTemperatureK = Double.NaN;
    double bestError = Double.POSITIVE_INFINITY;
    for (int index = 1; index < curve.size(); index++) {
      Point first = curve.get(index - 1);
      Point second = curve.get(index);
      double firstPressure = first.getPressureBara();
      double secondPressure = second.getPressureBara();
      double tolerance = PRESSURE_INTERSECTION_TOLERANCE
          * Math.max(1.0, Math.max(Math.abs(firstPressure), Math.abs(secondPressure)));
      double candidateTemperatureK;
      if (Math.abs(secondPressure - firstPressure) <= tolerance) {
        if (Math.abs(pressureBara - firstPressure) > tolerance) {
          continue;
        }
        double minimumTemperatureK = Math.min(first.getTemperatureK(), second.getTemperatureK());
        double maximumTemperatureK = Math.max(first.getTemperatureK(), second.getTemperatureK());
        candidateTemperatureK = Math.max(minimumTemperatureK, Math.min(maximumTemperatureK, targetTemperatureK));
      } else {
        double fraction = (pressureBara - firstPressure) / (secondPressure - firstPressure);
        if (fraction < -PRESSURE_INTERSECTION_TOLERANCE || fraction > 1.0 + PRESSURE_INTERSECTION_TOLERANCE) {
          continue;
        }
        fraction = Math.max(0.0, Math.min(1.0, fraction));
        candidateTemperatureK = first.getTemperatureK()
            + fraction * (second.getTemperatureK() - first.getTemperatureK());
      }
      double error = Math.abs(candidateTemperatureK - targetTemperatureK);
      if (error < bestError) {
        bestError = error;
        bestTemperatureK = candidateTemperatureK;
      }
    }
    return bestTemperatureK;
  }

  private static CriticalPointMetrics criticalPointMetrics(Point reference, Point candidate) {
    if (reference == null || candidate == null) {
      return CriticalPointMetrics.missing();
    }
    double temperatureDifferenceK = Math.abs(candidate.getTemperatureK() - reference.getTemperatureK());
    double pressureDifferenceBara = Math.abs(candidate.getPressureBara() - reference.getPressureBara());
    return CriticalPointMetrics.available(temperatureDifferenceK, pressureDifferenceBara,
        temperatureDifferenceK / reference.getTemperatureK(), pressureDifferenceBara / reference.getPressureBara());
  }

  private static List<Point> uniquePoints(List<Point> points) {
    List<Point> unique = new ArrayList<Point>();
    if (points == null) {
      return unique;
    }
    for (Point point : points) {
      if (point == null) {
        continue;
      }
      if (unique.isEmpty() || !samePoint(unique.get(unique.size() - 1), point)) {
        unique.add(point);
      }
    }
    if (unique.size() > 2 && samePoint(unique.get(0), unique.get(unique.size() - 1))) {
      unique.remove(unique.size() - 1);
    }
    return unique;
  }

  private static boolean samePoint(Point first, Point second) {
    return Math.abs(first.getTemperatureK() - second.getTemperatureK()) <= 1.0e-9
        && Math.abs(first.getPressureBara() - second.getPressureBara()) <= 1.0e-9;
  }

  private static boolean sameToken(String first, String second) {
    return normalizedToken(first).equals(normalizedToken(second));
  }

  private static String normalizedToken(String value) {
    return value.trim().toLowerCase(Locale.ROOT).replaceAll("[^a-z0-9]", "");
  }

  private static double maximum(List<Double> values) {
    double maximum = 0.0;
    for (double value : values) {
      maximum = Math.max(maximum, value);
    }
    return values.isEmpty() ? Double.POSITIVE_INFINITY : maximum;
  }

  private static double rms(List<Double> values) {
    if (values.isEmpty()) {
      return Double.POSITIVE_INFINITY;
    }
    double sum = 0.0;
    for (double value : values) {
      sum += value * value;
    }
    return Math.sqrt(sum / values.size());
  }

  private static double percentile95(List<Double> values) {
    if (values.isEmpty()) {
      return Double.POSITIVE_INFINITY;
    }
    List<Double> sorted = new ArrayList<Double>(values);
    Collections.sort(sorted, Comparator.naturalOrder());
    int index = Math.max(0, (int) Math.ceil(0.95 * sorted.size()) - 1);
    return sorted.get(index);
  }

  private static void requireText(String value, String label) {
    if (value == null || value.trim().isEmpty()) {
      throw new IllegalArgumentException(label + " is required");
    }
  }

  /** One PT point in absolute units. */
  public static final class Point {
    private final double temperatureK;
    private final double pressureBara;

    public Point(double temperatureK, double pressureBara) {
      if (!Double.isFinite(temperatureK) || temperatureK <= 0.0 || !Double.isFinite(pressureBara)
          || pressureBara <= 0.0) {
        throw new IllegalArgumentException("benchmark PT coordinates must be positive and finite");
      }
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }
  }

  /** Exact model-and-composition identity required before numeric errors may be calculated. */
  public static final class ModelFingerprint {
    private final String equationOfState;
    private final String volumeTranslation;
    private final String mixingRule;
    private final String compositionFingerprint;
    private final String parameterFingerprint;

    public ModelFingerprint(String equationOfState, String volumeTranslation, String mixingRule,
        String compositionFingerprint, String parameterFingerprint) {
      requireText(equationOfState, "equation of state");
      requireText(volumeTranslation, "volume translation");
      requireText(mixingRule, "mixing rule");
      requireText(compositionFingerprint, "composition fingerprint");
      requireText(parameterFingerprint, "parameter fingerprint");
      this.equationOfState = equationOfState;
      this.volumeTranslation = volumeTranslation;
      this.mixingRule = mixingRule;
      this.compositionFingerprint = compositionFingerprint;
      this.parameterFingerprint = parameterFingerprint;
    }

    public String getEquationOfState() {
      return equationOfState;
    }

    public String getVolumeTranslation() {
      return volumeTranslation;
    }

    public String getMixingRule() {
      return mixingRule;
    }

    public String getCompositionFingerprint() {
      return compositionFingerprint;
    }

    public boolean isCompositionFingerprintResolved() {
      return isResolvedFingerprint(compositionFingerprint);
    }

    public String getParameterFingerprint() {
      return parameterFingerprint;
    }

    public boolean isParameterFingerprintResolved() {
      return isResolvedFingerprint(parameterFingerprint);
    }

    private static boolean isResolvedFingerprint(String fingerprint) {
      String normalized = normalizedToken(fingerprint);
      return !normalized.equals("unknown") && !normalized.equals("unresolved") && !normalized.equals("notproven");
    }
  }

  /** One auditable source curve with an independently resolved topology. */
  public static final class Curve {
    private final String identifier;
    private final ReferenceSource source;
    private final ModelFingerprint modelFingerprint;
    private final BoundaryDefinition boundaryDefinition;
    private final boolean topologyResolved;
    private final String topologyEvidence;
    private final List<Point> points;
    private final Point criticalPoint;

    public Curve(String identifier, ReferenceSource source, ModelFingerprint modelFingerprint,
        BoundaryDefinition boundaryDefinition, boolean topologyResolved, String topologyEvidence, List<Point> points,
        Point criticalPoint) {
      requireText(identifier, "curve identifier");
      requireText(topologyEvidence, "topology evidence");
      if (source == null || modelFingerprint == null || boundaryDefinition == null || points == null) {
        throw new IllegalArgumentException("curve source, model, topology, and points are required");
      }
      this.identifier = identifier;
      this.source = source;
      this.modelFingerprint = modelFingerprint;
      this.boundaryDefinition = boundaryDefinition;
      this.topologyResolved = topologyResolved;
      this.topologyEvidence = topologyEvidence;
      this.points = Collections.unmodifiableList(new ArrayList<Point>(points));
      this.criticalPoint = criticalPoint;
    }

    public String getIdentifier() {
      return identifier;
    }

    public ReferenceSource getSource() {
      return source;
    }

    public ModelFingerprint getModelFingerprint() {
      return modelFingerprint;
    }

    public BoundaryDefinition getBoundaryDefinition() {
      return boundaryDefinition;
    }

    public boolean isTopologyResolved() {
      return topologyResolved;
    }

    public String getTopologyEvidence() {
      return topologyEvidence;
    }

    public List<Point> getPoints() {
      return points;
    }

    public Point getCriticalPoint() {
      return criticalPoint;
    }
  }

  /** Frozen pass/fail thresholds. Passing is impossible while this object is absent. */
  public static final class AcceptanceCriteria {
    private final double minimumReferenceCoverageFraction;
    private final double minimumCandidateCoverageFraction;
    private final double maximumRelativeTemperatureDifference;
    private final double maximumRmsRelativeTemperatureDifference;
    private final double maximumP95RelativeTemperatureDifference;
    private final double maximumSymmetricRelativeTemperatureDifference;
    private final boolean criticalPointRequired;
    private final double maximumCriticalTemperatureDifference;
    private final double maximumCriticalPressureDifference;

    public AcceptanceCriteria(double minimumReferenceCoverageFraction, double minimumCandidateCoverageFraction,
        double maximumRelativeTemperatureDifference, double maximumRmsRelativeTemperatureDifference,
        double maximumP95RelativeTemperatureDifference, double maximumSymmetricRelativeTemperatureDifference,
        boolean criticalPointRequired, double maximumCriticalTemperatureDifference,
        double maximumCriticalPressureDifference) {
      if (!fraction(minimumReferenceCoverageFraction) || !fraction(minimumCandidateCoverageFraction)
          || !positive(maximumRelativeTemperatureDifference) || !positive(maximumRmsRelativeTemperatureDifference)
          || !positive(maximumP95RelativeTemperatureDifference)
          || !positive(maximumSymmetricRelativeTemperatureDifference) || !positive(maximumCriticalTemperatureDifference)
          || !positive(maximumCriticalPressureDifference)) {
        throw new IllegalArgumentException("invalid hydrocarbon-water regression acceptance criteria");
      }
      this.minimumReferenceCoverageFraction = minimumReferenceCoverageFraction;
      this.minimumCandidateCoverageFraction = minimumCandidateCoverageFraction;
      this.maximumRelativeTemperatureDifference = maximumRelativeTemperatureDifference;
      this.maximumRmsRelativeTemperatureDifference = maximumRmsRelativeTemperatureDifference;
      this.maximumP95RelativeTemperatureDifference = maximumP95RelativeTemperatureDifference;
      this.maximumSymmetricRelativeTemperatureDifference = maximumSymmetricRelativeTemperatureDifference;
      this.criticalPointRequired = criticalPointRequired;
      this.maximumCriticalTemperatureDifference = maximumCriticalTemperatureDifference;
      this.maximumCriticalPressureDifference = maximumCriticalPressureDifference;
    }

    private static boolean fraction(double value) {
      return Double.isFinite(value) && value >= 0.0 && value <= 1.0;
    }

    private static boolean positive(double value) {
      return Double.isFinite(value) && value > 0.0;
    }

    public double getMinimumReferenceCoverageFraction() {
      return minimumReferenceCoverageFraction;
    }

    public double getMinimumCandidateCoverageFraction() {
      return minimumCandidateCoverageFraction;
    }

    public double getMaximumRelativeTemperatureDifference() {
      return maximumRelativeTemperatureDifference;
    }

    public double getMaximumRmsRelativeTemperatureDifference() {
      return maximumRmsRelativeTemperatureDifference;
    }

    public double getMaximumP95RelativeTemperatureDifference() {
      return maximumP95RelativeTemperatureDifference;
    }

    public double getMaximumSymmetricRelativeTemperatureDifference() {
      return maximumSymmetricRelativeTemperatureDifference;
    }

    public boolean isCriticalPointRequired() {
      return criticalPointRequired;
    }

    public double getMaximumCriticalTemperatureDifference() {
      return maximumCriticalTemperatureDifference;
    }

    public double getMaximumCriticalPressureDifference() {
      return maximumCriticalPressureDifference;
    }
  }

  /** Same-pressure statistics in one directed source-to-target projection. */
  public static final class DirectionalMetrics {
    private final int sourcePointCount;
    private final int matchedPointCount;
    private final double coverageFraction;
    private final double maximumAbsoluteTemperatureDifferenceK;
    private final double rmsAbsoluteTemperatureDifferenceK;
    private final double p95AbsoluteTemperatureDifferenceK;
    private final double maximumRelativeTemperatureDifference;
    private final double rmsRelativeTemperatureDifference;
    private final double p95RelativeTemperatureDifference;

    private DirectionalMetrics(int sourcePointCount, List<Double> absoluteErrors, List<Double> relativeErrors) {
      this.sourcePointCount = sourcePointCount;
      this.matchedPointCount = absoluteErrors.size();
      this.coverageFraction = sourcePointCount == 0 ? 0.0 : (double) matchedPointCount / sourcePointCount;
      this.maximumAbsoluteTemperatureDifferenceK = maximum(absoluteErrors);
      this.rmsAbsoluteTemperatureDifferenceK = rms(absoluteErrors);
      this.p95AbsoluteTemperatureDifferenceK = percentile95(absoluteErrors);
      this.maximumRelativeTemperatureDifference = maximum(relativeErrors);
      this.rmsRelativeTemperatureDifference = rms(relativeErrors);
      this.p95RelativeTemperatureDifference = percentile95(relativeErrors);
    }

    public int getSourcePointCount() {
      return sourcePointCount;
    }

    public int getMatchedPointCount() {
      return matchedPointCount;
    }

    public double getCoverageFraction() {
      return coverageFraction;
    }

    public double getMaximumAbsoluteTemperatureDifferenceK() {
      return maximumAbsoluteTemperatureDifferenceK;
    }

    public double getRmsAbsoluteTemperatureDifferenceK() {
      return rmsAbsoluteTemperatureDifferenceK;
    }

    public double getP95AbsoluteTemperatureDifferenceK() {
      return p95AbsoluteTemperatureDifferenceK;
    }

    public double getMaximumRelativeTemperatureDifference() {
      return maximumRelativeTemperatureDifference;
    }

    public double getRmsRelativeTemperatureDifference() {
      return rmsRelativeTemperatureDifference;
    }

    public double getP95RelativeTemperatureDifference() {
      return p95RelativeTemperatureDifference;
    }
  }

  /** Optional critical-point difference, kept separate from curve-shape statistics. */
  public static final class CriticalPointMetrics {
    private final boolean available;
    private final double absoluteTemperatureDifferenceK;
    private final double absolutePressureDifferenceBara;
    private final double relativeTemperatureDifference;
    private final double relativePressureDifference;

    private CriticalPointMetrics(boolean available, double absoluteTemperatureDifferenceK,
        double absolutePressureDifferenceBara, double relativeTemperatureDifference,
        double relativePressureDifference) {
      this.available = available;
      this.absoluteTemperatureDifferenceK = absoluteTemperatureDifferenceK;
      this.absolutePressureDifferenceBara = absolutePressureDifferenceBara;
      this.relativeTemperatureDifference = relativeTemperatureDifference;
      this.relativePressureDifference = relativePressureDifference;
    }

    private static CriticalPointMetrics missing() {
      return new CriticalPointMetrics(false, Double.NaN, Double.NaN, Double.NaN, Double.NaN);
    }

    private static CriticalPointMetrics available(double absoluteTemperatureDifferenceK,
        double absolutePressureDifferenceBara, double relativeTemperatureDifference,
        double relativePressureDifference) {
      return new CriticalPointMetrics(true, absoluteTemperatureDifferenceK, absolutePressureDifferenceBara,
          relativeTemperatureDifference, relativePressureDifference);
    }

    public boolean isAvailable() {
      return available;
    }

    public double getAbsoluteTemperatureDifferenceK() {
      return absoluteTemperatureDifferenceK;
    }

    public double getAbsolutePressureDifferenceBara() {
      return absolutePressureDifferenceBara;
    }

    public double getRelativeTemperatureDifference() {
      return relativeTemperatureDifference;
    }

    public double getRelativePressureDifference() {
      return relativePressureDifference;
    }
  }

  /** Bidirectional curve metrics; the primary error is PVTsim/paper reference to candidate. */
  public static final class Metrics {
    private final DirectionalMetrics referenceToCandidate;
    private final DirectionalMetrics candidateToReference;
    private final CriticalPointMetrics criticalPointMetrics;

    private Metrics(DirectionalMetrics referenceToCandidate, DirectionalMetrics candidateToReference,
        CriticalPointMetrics criticalPointMetrics) {
      this.referenceToCandidate = referenceToCandidate;
      this.candidateToReference = candidateToReference;
      this.criticalPointMetrics = criticalPointMetrics;
    }

    public DirectionalMetrics getReferenceToCandidate() {
      return referenceToCandidate;
    }

    public DirectionalMetrics getCandidateToReference() {
      return candidateToReference;
    }

    public double getReferenceCoverageFraction() {
      return referenceToCandidate.getCoverageFraction();
    }

    public double getCandidateCoverageFraction() {
      return candidateToReference.getCoverageFraction();
    }

    public double getMaximumRelativeTemperatureDifference() {
      return referenceToCandidate.getMaximumRelativeTemperatureDifference();
    }

    public double getRmsRelativeTemperatureDifference() {
      return referenceToCandidate.getRmsRelativeTemperatureDifference();
    }

    public double getP95RelativeTemperatureDifference() {
      return referenceToCandidate.getP95RelativeTemperatureDifference();
    }

    public double getSymmetricMaximumRelativeTemperatureDifference() {
      return Math.max(referenceToCandidate.getMaximumRelativeTemperatureDifference(),
          candidateToReference.getMaximumRelativeTemperatureDifference());
    }

    public CriticalPointMetrics getCriticalPointMetrics() {
      return criticalPointMetrics;
    }
  }

  /** Contract, threshold, and final benchmark verdict. */
  public static final class Result {
    private final boolean comparisonEligible;
    private final boolean thresholdsFrozen;
    private final boolean benchmarkPending;
    private final boolean accepted;
    private final Metrics metrics;
    private final List<String> violations;

    private Result(boolean comparisonEligible, boolean thresholdsFrozen, boolean benchmarkPending, boolean accepted,
        Metrics metrics, List<String> violations) {
      this.comparisonEligible = comparisonEligible;
      this.thresholdsFrozen = thresholdsFrozen;
      this.benchmarkPending = benchmarkPending;
      this.accepted = accepted;
      this.metrics = metrics;
      this.violations = Collections.unmodifiableList(new ArrayList<String>(violations));
    }

    private static Result ineligible(List<String> violations, boolean thresholdsFrozen) {
      return new Result(false, thresholdsFrozen, true, false, null, violations);
    }

    private static Result pending(Metrics metrics) {
      return new Result(true, false, true, false, metrics, Collections.<String>emptyList());
    }

    private static Result evaluated(Metrics metrics, List<String> violations) {
      return new Result(true, true, false, violations.isEmpty(), metrics, violations);
    }

    public boolean isComparisonEligible() {
      return comparisonEligible;
    }

    public boolean isThresholdsFrozen() {
      return thresholdsFrozen;
    }

    public boolean isBenchmarkPending() {
      return benchmarkPending;
    }

    public boolean isAccepted() {
      return accepted;
    }

    public Metrics getMetrics() {
      return metrics;
    }

    public List<String> getViolations() {
      return violations;
    }
  }
}
