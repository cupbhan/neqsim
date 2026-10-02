package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Globally searches simultaneous-incidence three-phase points from every possible mother-phase family. */
public final class HydrocarbonWaterThreePhaseTopologyFinder {
  private final SystemInterface template;
  private int pressureScanPoints = 25;
  private int temperatureScanIntervals = 40;
  private int maximumPressureRefinements = 40;
  private int maximumCoupledIterations = 60;
  private double pressureLogTolerance = 1.0e-6;
  private double temperatureCoincidenceToleranceK = 1.0e-3;
  private double maximumCandidateTemperatureGapK = 5.0;
  private double tangentPlaneTolerance = 1.0e-8;
  private double coupledResidualTolerance = 1.0e-7;

  /** Creates a non-destructive global three-phase topology finder. */
  public HydrocarbonWaterThreePhaseTopologyFinder(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets identical numerical controls for all three mother-phase searches. */
  public HydrocarbonWaterThreePhaseTopologyFinder setNumericalControls(int pressureScanPoints,
      int temperatureScanIntervals, int maximumPressureRefinements, int maximumCoupledIterations,
      double pressureLogTolerance, double temperatureCoincidenceToleranceK, double maximumCandidateTemperatureGapK,
      double tangentPlaneTolerance, double coupledResidualTolerance) {
    if (pressureScanPoints < 3 || temperatureScanIntervals < 2 || maximumPressureRefinements < 1
        || maximumCoupledIterations < 1 || !positive(pressureLogTolerance)
        || !positive(temperatureCoincidenceToleranceK) || !positive(maximumCandidateTemperatureGapK)
        || !positive(tangentPlaneTolerance) || !positive(coupledResidualTolerance)) {
      throw new IllegalArgumentException("invalid global three-phase topology controls");
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
   * Searches gas-, oil-, and aqueous-mother simultaneous-incidence topologies over one common PT domain.
   *
   * <p>
   * Every mother search and its no-point diagnostic are retained. Physical coupled points are deduplicated only after
   * all searches finish, so a failure in one topology can never be hidden by a successful point in another.
   * </p>
   */
  public Result find(double minimumPressureBara, double maximumPressureBara, double minimumTemperatureK,
      double maximumTemperatureK) {
    List<MotherSearch> searches = new ArrayList<MotherSearch>();
    List<PhysicalPoint> physicalPoints = new ArrayList<PhysicalPoint>();
    for (CandidatePhase mother : CandidatePhase.values()) {
      CandidatePhase[] incipient = otherPhases(mother);
      HydrocarbonWaterThreePhasePointFinder.Result search = new HydrocarbonWaterThreePhasePointFinder(template, mother,
          incipient[0], incipient[1])
          .setNumericalControls(pressureScanPoints, temperatureScanIntervals, maximumPressureRefinements,
              maximumCoupledIterations, pressureLogTolerance, temperatureCoincidenceToleranceK,
              maximumCandidateTemperatureGapK, tangentPlaneTolerance, coupledResidualTolerance)
          .find(minimumPressureBara, maximumPressureBara, minimumTemperatureK, maximumTemperatureK);
      MotherSearch motherSearch = new MotherSearch(mother, incipient[0], incipient[1], search);
      searches.add(motherSearch);
      for (ThreePhasePointSolver.Result point : search.getPhysicalPoints()) {
        addDistinct(physicalPoints, new PhysicalPoint(motherSearch, point));
      }
    }
    return new Result(searches, physicalPoints);
  }

  private static CandidatePhase[] otherPhases(CandidatePhase mother) {
    switch (mother) {
    case GAS:
      return new CandidatePhase[] {CandidatePhase.OIL, CandidatePhase.AQUEOUS};
    case OIL:
      return new CandidatePhase[] {CandidatePhase.GAS, CandidatePhase.AQUEOUS};
    case AQUEOUS:
      return new CandidatePhase[] {CandidatePhase.GAS, CandidatePhase.OIL};
    default:
      throw new IllegalArgumentException("unsupported mother phase " + mother);
    }
  }

  private static void addDistinct(List<PhysicalPoint> points, PhysicalPoint candidate) {
    for (PhysicalPoint existing : points) {
      double relativeTemperature = Math.abs(existing.point.getTemperatureK() - candidate.point.getTemperatureK())
          / Math.max(1.0, candidate.point.getTemperatureK());
      double logPressure = Math.abs(Math.log(existing.point.getPressureBara() / candidate.point.getPressureBara()));
      if (relativeTemperature <= 1.0e-5 && logPressure <= 1.0e-5) {
        return;
      }
    }
    points.add(candidate);
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /** Complete diagnostics for one explicit mother-phase topology. */
  public static final class MotherSearch {
    private final CandidatePhase motherPhase;
    private final CandidatePhase firstIncipientPhase;
    private final CandidatePhase secondIncipientPhase;
    private final HydrocarbonWaterThreePhasePointFinder.Result result;

    private MotherSearch(CandidatePhase motherPhase, CandidatePhase firstIncipientPhase,
        CandidatePhase secondIncipientPhase, HydrocarbonWaterThreePhasePointFinder.Result result) {
      this.motherPhase = motherPhase;
      this.firstIncipientPhase = firstIncipientPhase;
      this.secondIncipientPhase = secondIncipientPhase;
      this.result = result;
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

    public HydrocarbonWaterThreePhasePointFinder.Result getResult() {
      return result;
    }
  }

  /** One globally deduplicated physical point retaining the mother search that found it. */
  public static final class PhysicalPoint {
    private final MotherSearch source;
    private final ThreePhasePointSolver.Result point;

    private PhysicalPoint(MotherSearch source, ThreePhasePointSolver.Result point) {
      this.source = source;
      this.point = point;
    }

    public MotherSearch getSource() {
      return source;
    }

    public ThreePhasePointSolver.Result getPoint() {
      return point;
    }
  }

  /** Immutable global search retaining all mother results and every distinct physical point. */
  public static final class Result {
    private final List<MotherSearch> searches;
    private final List<PhysicalPoint> physicalPoints;

    private Result(List<MotherSearch> searches, List<PhysicalPoint> physicalPoints) {
      this.searches = Collections.unmodifiableList(new ArrayList<MotherSearch>(searches));
      this.physicalPoints = Collections.unmodifiableList(new ArrayList<PhysicalPoint>(physicalPoints));
    }

    public List<MotherSearch> getSearches() {
      return searches;
    }

    public List<PhysicalPoint> getPhysicalPoints() {
      return physicalPoints;
    }

    public boolean hasPhysicalPoint() {
      return !physicalPoints.isEmpty();
    }
  }
}
