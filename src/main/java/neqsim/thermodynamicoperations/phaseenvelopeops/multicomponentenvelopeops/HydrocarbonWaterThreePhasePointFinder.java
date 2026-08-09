package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Discovers and strictly corrects points where two phases become incipient from one mother phase.
 *
 * <p>
 * A three-phase point is never inferred from a change in the number of phases returned by a TP flash. The finder first
 * traces the two independent single-phase zero-TPD boundaries over pressure. It then brackets intersections of their
 * temperatures and passes the two independently calculated incipient compositions to {@link ThreePhasePointSolver}.
 * Only a converged, distinct-phase coupled solution is returned as a physical point.
 * </p>
 */
public final class HydrocarbonWaterThreePhasePointFinder {
  private final SystemInterface template;
  private final CandidatePhase motherPhase;
  private final CandidatePhase firstIncipientPhase;
  private final CandidatePhase secondIncipientPhase;
  private int pressureScanPoints = 25;
  private int temperatureScanIntervals = 40;
  private int maximumPressureRefinements = 40;
  private int maximumCoupledIterations = 60;
  private double pressureLogTolerance = 1.0e-6;
  private double temperatureCoincidenceToleranceK = 1.0e-3;
  private double maximumCandidateTemperatureGapK = 5.0;
  private double tangentPlaneTolerance = 1.0e-8;
  private double coupledResidualTolerance = 1.0e-7;

  /** Creates a finder for two distinct incipient phases from one mother phase. */
  public HydrocarbonWaterThreePhasePointFinder(SystemInterface template, CandidatePhase motherPhase,
      CandidatePhase firstIncipientPhase, CandidatePhase secondIncipientPhase) {
    if (template == null || motherPhase == null || firstIncipientPhase == null || secondIncipientPhase == null) {
      throw new IllegalArgumentException("template and all phase families must be specified");
    }
    if (motherPhase == firstIncipientPhase || motherPhase == secondIncipientPhase
        || firstIncipientPhase == secondIncipientPhase) {
      throw new IllegalArgumentException("mother and incipient phase families must be distinct");
    }
    this.template = template.clone();
    this.motherPhase = motherPhase;
    this.firstIncipientPhase = firstIncipientPhase;
    this.secondIncipientPhase = secondIncipientPhase;
  }

  /**
   * Sets scan, one-dimensional intersection, and coupled-corrector controls.
   *
   * @return this finder
   */
  public HydrocarbonWaterThreePhasePointFinder setNumericalControls(int pressureScanPoints,
      int temperatureScanIntervals, int maximumPressureRefinements, int maximumCoupledIterations,
      double pressureLogTolerance, double temperatureCoincidenceToleranceK, double maximumCandidateTemperatureGapK,
      double tangentPlaneTolerance, double coupledResidualTolerance) {
    if (pressureScanPoints < 3 || temperatureScanIntervals < 2 || maximumPressureRefinements < 1
        || maximumCoupledIterations < 1 || !isPositiveFinite(pressureLogTolerance)
        || !isPositiveFinite(temperatureCoincidenceToleranceK) || !isPositiveFinite(maximumCandidateTemperatureGapK)
        || !isPositiveFinite(tangentPlaneTolerance) || !isPositiveFinite(coupledResidualTolerance)) {
      throw new IllegalArgumentException("invalid three-phase finder numerical controls");
    }
    this.pressureScanPoints = pressureScanPoints;
    this.temperatureScanIntervals = temperatureScanIntervals;
    this.maximumPressureRefinements = maximumPressureRefinements;
    this.maximumCoupledIterations = maximumCoupledIterations;
    this.pressureLogTolerance = pressureLogTolerance;
    this.temperatureCoincidenceToleranceK = temperatureCoincidenceToleranceK;
    this.maximumCandidateTemperatureGapK = maximumCandidateTemperatureGapK;
    this.tangentPlaneTolerance = tangentPlaneTolerance;
    this.coupledResidualTolerance = coupledResidualTolerance;
    return this;
  }

  /**
   * Scans a log-pressure domain and returns every strictly corrected simultaneous-incidence point found.
   *
   * @param minimumPressureBara lower pressure bound in bara
   * @param maximumPressureBara upper pressure bound in bara
   * @param minimumTemperatureK lower boundary-search temperature in kelvin
   * @param maximumTemperatureK upper boundary-search temperature in kelvin
   * @return immutable scan result with full diagnostics, including a valid no-point outcome
   */
  public Result find(double minimumPressureBara, double maximumPressureBara, double minimumTemperatureK,
      double maximumTemperatureK) {
    validateDomain(minimumPressureBara, maximumPressureBara, minimumTemperatureK, maximumTemperatureK);
    List<Sample> samples = new ArrayList<Sample>();
    double minimumLogPressure = Math.log(minimumPressureBara);
    double maximumLogPressure = Math.log(maximumPressureBara);
    for (int pointIndex = 0; pointIndex < pressureScanPoints; pointIndex++) {
      double fraction = (double) pointIndex / (pressureScanPoints - 1);
      double pressureBara = Math.exp(minimumLogPressure + fraction * (maximumLogPressure - minimumLogPressure));
      samples.add(evaluateAtPressure(pressureBara, minimumTemperatureK, maximumTemperatureK));
    }

    List<Bracket> brackets = discoverBrackets(samples);
    List<Sample> candidateSeeds = discoverNearCoincidentSamples(samples, brackets);
    List<CoupledTrial> coupledTrials = new ArrayList<CoupledTrial>();
    List<ThreePhasePointSolver.Result> physicalPoints = new ArrayList<ThreePhasePointSolver.Result>();
    for (Bracket bracket : brackets) {
      RefinedIntersection intersection = refine(bracket, minimumTemperatureK, maximumTemperatureK);
      samples.addAll(intersection.additionalSamples);
      CoupledTrial trial = correct(intersection.bestSample, "SIGN_CHANGE");
      coupledTrials.add(trial);
      if (trial.point.isConverged() && insideDomain(trial.point, minimumPressureBara, maximumPressureBara,
          minimumTemperatureK, maximumTemperatureK)) {
        addDistinctPhysicalPoint(physicalPoints, trial.point);
      }
    }
    for (Sample candidate : candidateSeeds) {
      CoupledTrial trial = correct(candidate, "NEAR_COINCIDENCE");
      coupledTrials.add(trial);
      if (trial.point.isConverged() && insideDomain(trial.point, minimumPressureBara, maximumPressureBara,
          minimumTemperatureK, maximumTemperatureK)) {
        addDistinctPhysicalPoint(physicalPoints, trial.point);
      }
    }
    Collections.sort(samples, Comparator.comparingDouble(Sample::getPressureBara));
    Collections.sort(physicalPoints, Comparator.comparingDouble(ThreePhasePointSolver.Result::getPressureBara));
    String diagnostic = physicalPoints.isEmpty() ? failureDiagnostic(samples, brackets, coupledTrials) : null;
    return new Result(motherPhase, firstIncipientPhase, secondIncipientPhase, samples, coupledTrials, physicalPoints,
        diagnostic);
  }

  private Sample evaluateAtPressure(double pressureBara, double minimumTemperatureK, double maximumTemperatureK) {
    IncipientPhaseBoundaryPointSolver.Result first = boundarySolver(firstIncipientPhase).solve(pressureBara,
        minimumTemperatureK, maximumTemperatureK);
    IncipientPhaseBoundaryPointSolver.Result second = boundarySolver(secondIncipientPhase).solve(pressureBara,
        minimumTemperatureK, maximumTemperatureK);
    return new Sample(pressureBara, first, second);
  }

  private IncipientPhaseBoundaryPointSolver boundarySolver(CandidatePhase incipientPhase) {
    return new IncipientPhaseBoundaryPointSolver(template, motherPhase, incipientPhase).setNumericalControls(
        temperatureScanIntervals, 80, Math.max(1.0e-5, 0.1 * temperatureCoincidenceToleranceK), tangentPlaneTolerance);
  }

  private List<Bracket> discoverBrackets(List<Sample> samples) {
    List<Bracket> brackets = new ArrayList<Bracket>();
    Sample previous = null;
    for (Sample current : samples) {
      if (!current.isUsable()) {
        previous = null;
        continue;
      }
      if (previous != null && hasSignChange(previous.temperatureGapK, current.temperatureGapK)) {
        brackets.add(new Bracket(previous, current));
      }
      previous = current;
    }
    return brackets;
  }

  private List<Sample> discoverNearCoincidentSamples(List<Sample> samples, List<Bracket> brackets) {
    List<Sample> candidates = new ArrayList<Sample>();
    for (int index = 0; index < samples.size(); index++) {
      Sample sample = samples.get(index);
      if (!sample.isUsable() || Math.abs(sample.temperatureGapK) > maximumCandidateTemperatureGapK
          || belongsToBracket(sample, brackets)) {
        continue;
      }
      boolean localMinimum = index == 0 || !samples.get(index - 1).isUsable()
          || Math.abs(sample.temperatureGapK) <= Math.abs(samples.get(index - 1).temperatureGapK);
      localMinimum &= index + 1 == samples.size() || !samples.get(index + 1).isUsable()
          || Math.abs(sample.temperatureGapK) <= Math.abs(samples.get(index + 1).temperatureGapK);
      if (localMinimum || Math.abs(sample.temperatureGapK) <= temperatureCoincidenceToleranceK) {
        candidates.add(sample);
      }
    }
    return candidates;
  }

  private RefinedIntersection refine(Bracket bracket, double minimumTemperatureK, double maximumTemperatureK) {
    Sample lower = bracket.lower;
    Sample upper = bracket.upper;
    Sample best = better(lower, upper);
    List<Sample> additional = new ArrayList<Sample>();
    int refinements = 0;
    while (refinements < maximumPressureRefinements
        && Math.log(upper.pressureBara / lower.pressureBara) > pressureLogTolerance
        && Math.abs(best.temperatureGapK) > temperatureCoincidenceToleranceK) {
      refinements++;
      double middlePressure = Math.sqrt(lower.pressureBara * upper.pressureBara);
      Sample middle = evaluateAtPressure(middlePressure, minimumTemperatureK, maximumTemperatureK);
      additional.add(middle);
      if (!middle.isUsable()) {
        break;
      }
      best = better(best, middle);
      if (hasSignChange(lower.temperatureGapK, middle.temperatureGapK)) {
        upper = middle;
      } else {
        lower = middle;
      }
    }
    return new RefinedIntersection(best, additional, refinements);
  }

  private CoupledTrial correct(Sample seed, String source) {
    double temperatureK = 0.5 * (seed.firstBoundary.getTemperatureK() + seed.secondBoundary.getTemperatureK());
    ThreePhasePointSolver.Result point = new ThreePhasePointSolver(template, motherPhase, firstIncipientPhase,
        secondIncipientPhase).setNumericalControls(maximumCoupledIterations, coupledResidualTolerance, 2.0e-5)
        .solve(temperatureK, seed.pressureBara, seed.firstBoundary.getIncipientComposition(),
            seed.secondBoundary.getIncipientComposition());
    return new CoupledTrial(source, seed, point);
  }

  private static void addDistinctPhysicalPoint(List<ThreePhasePointSolver.Result> points,
      ThreePhasePointSolver.Result candidate) {
    for (ThreePhasePointSolver.Result existing : points) {
      double relativeTemperature = Math.abs(existing.getTemperatureK() - candidate.getTemperatureK())
          / Math.max(1.0, candidate.getTemperatureK());
      double logPressure = Math.abs(Math.log(existing.getPressureBara() / candidate.getPressureBara()));
      if (relativeTemperature <= 1.0e-5 && logPressure <= 1.0e-5) {
        return;
      }
    }
    points.add(candidate);
  }

  private static boolean insideDomain(ThreePhasePointSolver.Result point, double minimumPressureBara,
      double maximumPressureBara, double minimumTemperatureK, double maximumTemperatureK) {
    return point.getPressureBara() >= minimumPressureBara && point.getPressureBara() <= maximumPressureBara
        && point.getTemperatureK() >= minimumTemperatureK && point.getTemperatureK() <= maximumTemperatureK;
  }

  private static String failureDiagnostic(List<Sample> samples, List<Bracket> brackets,
      List<CoupledTrial> coupledTrials) {
    int usable = 0;
    double minimumGap = Double.POSITIVE_INFINITY;
    for (Sample sample : samples) {
      if (sample.isUsable()) {
        usable++;
        minimumGap = Math.min(minimumGap, Math.abs(sample.temperatureGapK));
      }
    }
    if (usable == 0) {
      return "neither pair of single-phase incipient boundaries was jointly traceable over the scan domain";
    }
    if (brackets.isEmpty() && coupledTrials.isEmpty()) {
      return "no intersection of the two independently traced incipient boundaries; minimum temperature gap="
          + minimumGap + " K";
    }
    return "intersection candidates were found but every coupled three-phase correction failed the residual, "
        + "distinct-phase, or requested-domain gate; candidates=" + coupledTrials.size() + ", minimum gap=" + minimumGap
        + " K";
  }

  private static boolean belongsToBracket(Sample sample, List<Bracket> brackets) {
    for (Bracket bracket : brackets) {
      if (sample == bracket.lower || sample == bracket.upper) {
        return true;
      }
    }
    return false;
  }

  private static Sample better(Sample first, Sample second) {
    return Math.abs(first.temperatureGapK) <= Math.abs(second.temperatureGapK) ? first : second;
  }

  private static boolean hasSignChange(double first, double second) {
    return Double.isFinite(first) && Double.isFinite(second)
        && (first <= 0.0 && second >= 0.0 || first >= 0.0 && second <= 0.0);
  }

  private static boolean isPositiveFinite(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  private static void validateDomain(double minimumPressureBara, double maximumPressureBara, double minimumTemperatureK,
      double maximumTemperatureK) {
    if (!isPositiveFinite(minimumPressureBara) || !isPositiveFinite(maximumPressureBara)
        || maximumPressureBara <= minimumPressureBara || !Double.isFinite(minimumTemperatureK)
        || minimumTemperatureK < 50.0 || !Double.isFinite(maximumTemperatureK)
        || maximumTemperatureK <= minimumTemperatureK) {
      throw new IllegalArgumentException("invalid three-phase point search domain");
    }
  }

  private static final class Bracket {
    private final Sample lower;
    private final Sample upper;

    private Bracket(Sample first, Sample second) {
      this.lower = first.pressureBara <= second.pressureBara ? first : second;
      this.upper = first.pressureBara <= second.pressureBara ? second : first;
    }
  }

  private static final class RefinedIntersection {
    private final Sample bestSample;
    private final List<Sample> additionalSamples;
    @SuppressWarnings("unused")
    private final int refinements;

    private RefinedIntersection(Sample bestSample, List<Sample> additionalSamples, int refinements) {
      this.bestSample = bestSample;
      this.additionalSamples = additionalSamples;
      this.refinements = refinements;
    }
  }

  /** One pressure sample containing both independently corrected incipient boundaries. */
  public static final class Sample {
    private final double pressureBara;
    private final IncipientPhaseBoundaryPointSolver.Result firstBoundary;
    private final IncipientPhaseBoundaryPointSolver.Result secondBoundary;
    private final double temperatureGapK;

    private Sample(double pressureBara, IncipientPhaseBoundaryPointSolver.Result firstBoundary,
        IncipientPhaseBoundaryPointSolver.Result secondBoundary) {
      this.pressureBara = pressureBara;
      this.firstBoundary = firstBoundary;
      this.secondBoundary = secondBoundary;
      this.temperatureGapK = firstBoundary.isConverged() && secondBoundary.isConverged()
          ? firstBoundary.getTemperatureK() - secondBoundary.getTemperatureK()
          : Double.NaN;
    }

    /** @return pressure in bara */
    public double getPressureBara() {
      return pressureBara;
    }

    /** @return first incipient-boundary solve */
    public IncipientPhaseBoundaryPointSolver.Result getFirstBoundary() {
      return firstBoundary;
    }

    /** @return second incipient-boundary solve */
    public IncipientPhaseBoundaryPointSolver.Result getSecondBoundary() {
      return secondBoundary;
    }

    /** @return first-boundary temperature minus second-boundary temperature in kelvin */
    public double getTemperatureGapK() {
      return temperatureGapK;
    }

    /** @return true when both independent boundary points passed their zero-TPD gates */
    public boolean isUsable() {
      return firstBoundary.isConverged() && secondBoundary.isConverged() && Double.isFinite(temperatureGapK);
    }
  }

  /** One strict coupled correction attempt and the independent seed that produced it. */
  public static final class CoupledTrial {
    private final String source;
    private final Sample seed;
    private final ThreePhasePointSolver.Result point;

    private CoupledTrial(String source, Sample seed, ThreePhasePointSolver.Result point) {
      this.source = source;
      this.seed = seed;
      this.point = point;
    }

    public String getSource() {
      return source;
    }

    public Sample getSeed() {
      return seed;
    }

    public ThreePhasePointSolver.Result getPoint() {
      return point;
    }
  }

  /** Immutable full-domain search result. An empty physical-point list is a valid negative topology result. */
  public static final class Result {
    private final CandidatePhase motherPhase;
    private final CandidatePhase firstIncipientPhase;
    private final CandidatePhase secondIncipientPhase;
    private final List<Sample> samples;
    private final List<CoupledTrial> coupledTrials;
    private final List<ThreePhasePointSolver.Result> physicalPoints;
    private final String diagnostic;

    private Result(CandidatePhase motherPhase, CandidatePhase firstIncipientPhase, CandidatePhase secondIncipientPhase,
        List<Sample> samples, List<CoupledTrial> coupledTrials, List<ThreePhasePointSolver.Result> physicalPoints,
        String diagnostic) {
      this.motherPhase = motherPhase;
      this.firstIncipientPhase = firstIncipientPhase;
      this.secondIncipientPhase = secondIncipientPhase;
      this.samples = Collections.unmodifiableList(new ArrayList<Sample>(samples));
      this.coupledTrials = Collections.unmodifiableList(new ArrayList<CoupledTrial>(coupledTrials));
      this.physicalPoints = Collections.unmodifiableList(new ArrayList<ThreePhasePointSolver.Result>(physicalPoints));
      this.diagnostic = diagnostic;
    }

    /** @return true when at least one physical simultaneous-incidence point was strictly corrected */
    public boolean hasPhysicalPoint() {
      return !physicalPoints.isEmpty();
    }

    public CandidatePhase getMotherPhase() {
      return motherPhase;
    }

    public CandidatePhase getFirstIncipientPhase() {
      return firstIncipientPhase;
    }

    public CandidatePhase getSecondIncipientPhase() {
      return secondIncipientPhase;
    }

    public List<Sample> getSamples() {
      return samples;
    }

    public List<CoupledTrial> getCoupledTrials() {
      return coupledTrials;
    }

    public List<ThreePhasePointSolver.Result> getPhysicalPoints() {
      return physicalPoints;
    }

    public String getDiagnostic() {
      return diagnostic;
    }

    /** @return smallest absolute independent-boundary temperature difference, or positive infinity */
    public double getMinimumTemperatureGapK() {
      double minimum = Double.POSITIVE_INFINITY;
      for (Sample sample : samples) {
        if (sample.isUsable()) {
          minimum = Math.min(minimum, Math.abs(sample.temperatureGapK));
        }
      }
      return minimum;
    }
  }
}
