package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Finds one single-to-two-phase incipient boundary point at fixed pressure from a zero-TPD condition. */
public final class IncipientPhaseBoundaryPointSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;

  private final SystemInterface template;
  private final CandidatePhase motherPhase;
  private final CandidatePhase incipientPhase;
  private int scanIntervals = 32;
  private int maximumBisections = 80;
  private double temperatureToleranceK = 1.0e-5;
  private double tangentPlaneTolerance = 1.0e-8;

  /**
   * Creates a fixed-pressure incipient-boundary point solver.
   *
   * @param template configured EOS/mixing-rule system
   * @param motherPhase fixed-composition mother phase
   * @param incipientPhase distinct phase whose TPD is driven to zero
   */
  public IncipientPhaseBoundaryPointSolver(SystemInterface template, CandidatePhase motherPhase,
      CandidatePhase incipientPhase) {
    if (template == null || motherPhase == null || incipientPhase == null) {
      throw new IllegalArgumentException("template, mother phase, and incipient phase must be specified");
    }
    if (motherPhase == incipientPhase) {
      throw new IllegalArgumentException("mother and incipient phase must be distinct");
    }
    this.template = template.clone();
    this.motherPhase = motherPhase;
    this.incipientPhase = incipientPhase;
  }

  /**
   * Sets scanning and bisection controls.
   *
   * @param scanIntervals temperature intervals used to discover a non-trivial sign change
   * @param maximumBisections maximum root bisections
   * @param temperatureToleranceK final temperature bracket width
   * @param tangentPlaneTolerance final absolute TPD tolerance
   * @return this solver
   */
  public IncipientPhaseBoundaryPointSolver setNumericalControls(int scanIntervals, int maximumBisections,
      double temperatureToleranceK, double tangentPlaneTolerance) {
    if (scanIntervals < 2 || maximumBisections < 1 || !Double.isFinite(temperatureToleranceK)
        || temperatureToleranceK <= 0.0 || !Double.isFinite(tangentPlaneTolerance) || tangentPlaneTolerance <= 0.0) {
      throw new IllegalArgumentException("invalid incipient-boundary numerical controls");
    }
    this.scanIntervals = scanIntervals;
    this.maximumBisections = maximumBisections;
    this.temperatureToleranceK = temperatureToleranceK;
    this.tangentPlaneTolerance = tangentPlaneTolerance;
    return this;
  }

  /**
   * Locates the first non-trivial TPD sign change in the supplied temperature interval and refines it.
   *
   * @param pressureBara fixed pressure in bara
   * @param minimumTemperatureK scan lower temperature in kelvin
   * @param maximumTemperatureK scan upper temperature in kelvin
   * @return immutable boundary point or explicit failure diagnostic
   */
  public Result solve(double pressureBara, double minimumTemperatureK, double maximumTemperatureK) {
    if (!Double.isFinite(pressureBara) || pressureBara <= 0.0 || !Double.isFinite(minimumTemperatureK)
        || !Double.isFinite(maximumTemperatureK) || minimumTemperatureK < 50.0
        || maximumTemperatureK <= minimumTemperatureK) {
      throw new IllegalArgumentException("invalid fixed-pressure boundary search domain");
    }

    Trial lower = null;
    Trial upper = null;
    Trial previous = null;
    int evaluations = 0;
    for (int scanIndex = 0; scanIndex <= scanIntervals; scanIndex++) {
      double fraction = (double) scanIndex / scanIntervals;
      double temperatureK = minimumTemperatureK + fraction * (maximumTemperatureK - minimumTemperatureK);
      Trial current = previous == null ? evaluateFromIndependentSeeds(temperatureK, pressureBara)
          : evaluate(temperatureK, pressureBara, previous.composition);
      if (!current.isUsable() && previous != null) {
        current = evaluateFromIndependentSeeds(temperatureK, pressureBara);
      }
      evaluations++;
      if (!current.isUsable()) {
        previous = null;
        continue;
      }
      if (previous != null && previous.isUsable()
          && hasSignChange(previous.tangentPlaneDistance, current.tangentPlaneDistance)) {
        lower = previous.temperatureK <= current.temperatureK ? previous : current;
        upper = previous.temperatureK <= current.temperatureK ? current : previous;
        break;
      }
      previous = current;
    }
    if (lower == null || upper == null) {
      return Result.failure(motherPhase, incipientPhase, pressureBara, evaluations,
          "no non-trivial TPD sign change found in the temperature interval");
    }

    Trial best = Math.abs(lower.tangentPlaneDistance) <= Math.abs(upper.tangentPlaneDistance) ? lower : upper;
    int bisections = 0;
    while (bisections < maximumBisections && upper.temperatureK - lower.temperatureK > temperatureToleranceK
        && Math.abs(best.tangentPlaneDistance) > tangentPlaneTolerance) {
      bisections++;
      double midpointTemperature = 0.5 * (lower.temperatureK + upper.temperatureK);
      Trial seedSource = midpointTemperature - lower.temperatureK <= upper.temperatureK - midpointTemperature ? lower
          : upper;
      Trial midpoint = evaluate(midpointTemperature, pressureBara, seedSource.composition);
      evaluations++;
      if (!midpoint.isUsable()) {
        midpoint = evaluateFromIndependentSeeds(midpointTemperature, pressureBara);
        evaluations++;
      }
      if (!midpoint.isUsable()) {
        return Result.failure(motherPhase, incipientPhase, pressureBara, evaluations,
            "non-trivial incipient stationary point was lost during bisection");
      }
      if (Math.abs(midpoint.tangentPlaneDistance) < Math.abs(best.tangentPlaneDistance)) {
        best = midpoint;
      }
      if (hasSignChange(lower.tangentPlaneDistance, midpoint.tangentPlaneDistance)) {
        upper = midpoint;
      } else {
        lower = midpoint;
      }
    }
    double secantDenominator = upper.tangentPlaneDistance - lower.tangentPlaneDistance;
    if (Double.isFinite(secantDenominator) && Math.abs(secantDenominator) > 1.0e-20) {
      double secantTemperature = (lower.temperatureK * upper.tangentPlaneDistance
          - upper.temperatureK * lower.tangentPlaneDistance) / secantDenominator;
      if (secantTemperature > lower.temperatureK && secantTemperature < upper.temperatureK) {
        Trial secant = evaluate(secantTemperature, pressureBara,
            Math.abs(secantTemperature - lower.temperatureK) <= Math.abs(upper.temperatureK - secantTemperature)
                ? lower.composition
                : upper.composition);
        evaluations++;
        if (secant.isUsable() && Math.abs(secant.tangentPlaneDistance) < Math.abs(best.tangentPlaneDistance)) {
          best = secant;
        }
      }
    }
    boolean converged = upper.temperatureK - lower.temperatureK <= temperatureToleranceK
        || Math.abs(best.tangentPlaneDistance) <= tangentPlaneTolerance;
    return new Result(motherPhase, incipientPhase, best.temperatureK, pressureBara, best.composition,
        best.tangentPlaneDistance, best.stationarityResidual, upper.temperatureK - lower.temperatureK, evaluations,
        bisections, converged, converged ? null : "maximum bisection count reached");
  }

  private Trial evaluate(double temperatureK, double pressureBara, double[] compositionSeed) {
    Trial successiveSubstitution = evaluateBySuccessiveSubstitution(temperatureK, pressureBara, compositionSeed);
    if (successiveSubstitution.isUsable()) {
      return successiveSubstitution;
    }
    Trial newton = evaluateByNewton(temperatureK, pressureBara, compositionSeed);
    if (newton.isUsable()) {
      return newton;
    }
    if (newton.converged && !successiveSubstitution.converged) {
      return newton;
    }
    return successiveSubstitution;
  }

  private Trial evaluateFromIndependentSeeds(double temperatureK, double pressureBara) {
    java.util.List<double[]> seeds = independentSeeds(temperatureK, pressureBara);
    Trial bestUsable = null;
    Trial bestFailure = null;
    for (double[] seed : seeds) {
      Trial trial = evaluate(temperatureK, pressureBara, seed);
      if (trial.isUsable() && (bestUsable == null || trial.tangentPlaneDistance < bestUsable.tangentPlaneDistance)) {
        bestUsable = trial;
      } else if (bestFailure == null || trial.stationarityResidual < bestFailure.stationarityResidual) {
        bestFailure = trial;
      }
    }
    return bestUsable == null ? bestFailure : bestUsable;
  }

  private java.util.List<double[]> independentSeeds(double temperatureK, double pressureBara) {
    java.util.List<double[]> seeds = new java.util.ArrayList<double[]>();
    addDistinctSeed(seeds, gridPhaseSeed(temperatureK, pressureBara));
    addDistinctSeed(seeds, physicalSeed(temperatureK, pressureBara));
    addDistinctSeed(seeds, heavyEnrichedSeed());
    addDistinctSeed(seeds, heaviestHydrocarbonSeed());
    addDistinctSeed(seeds, overallPerturbedTowardHeaviest(0.05));
    addDistinctSeed(seeds, overallPerturbedTowardHeaviest(0.25));
    addDistinctSeed(seeds, overallComposition());
    return seeds;
  }

  private double[] gridPhaseSeed(double temperatureK, double pressureBara) {
    try {
      SystemInterface grid = template.clone();
      grid.setTemperature(temperatureK);
      grid.setPressure(pressureBara);
      grid.setMultiPhaseCheck(true);
      grid.setMaxNumberOfPhases(3);
      new ThermodynamicOperations(grid).TPflash();
      grid.init(1);
      for (int phaseIndex = 0; phaseIndex < grid.getNumberOfPhases(); phaseIndex++) {
        PhaseType phaseType = grid.getPhase(phaseIndex).getType();
        boolean matches = incipientPhase == CandidatePhase.GAS && phaseType == PhaseType.GAS
            || incipientPhase == CandidatePhase.OIL && (phaseType == PhaseType.OIL || phaseType == PhaseType.LIQUID)
            || incipientPhase == CandidatePhase.AQUEOUS && phaseType == PhaseType.AQUEOUS;
        if (!matches) {
          continue;
        }
        double[] composition = new double[grid.getPhase(phaseIndex).getNumberOfComponents()];
        for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
          composition[componentIndex] = Math.max(grid.getPhase(phaseIndex).getComponent(componentIndex).getx(),
              MINIMUM_COMPOSITION);
        }
        normalize(composition);
        return composition;
      }
    } catch (RuntimeException error) {
      // A coarse TP flash is seed evidence only; strict zero-TPD correction below remains authoritative.
    }
    return null;
  }

  private double[] overallComposition() {
    double[] composition = new double[template.getPhase(0).getNumberOfComponents()];
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(),
          MINIMUM_COMPOSITION);
    }
    normalize(composition);
    return composition;
  }

  private double[] heavyEnrichedSeed() {
    double[] seed = overallComposition();
    for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
      ComponentInterface component = template.getPhase(0).getComponent(componentIndex);
      if (component.getComponentName().equalsIgnoreCase("water")) {
        seed[componentIndex] *= 1.0e-12;
      } else if (component.isHydrocarbon()) {
        double exponent = Math.max(-4.0, Math.min(4.0, (component.getTC() - 300.0) / 75.0));
        seed[componentIndex] *= Math.exp(exponent);
      } else {
        seed[componentIndex] *= 0.1;
      }
      seed[componentIndex] = Math.max(seed[componentIndex], MINIMUM_COMPOSITION);
    }
    normalize(seed);
    return seed;
  }

  private double[] heaviestHydrocarbonSeed() {
    double[] seed = overallComposition();
    int heaviest = heaviestHydrocarbonIndex();
    for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
      seed[componentIndex] = Math.max(seed[componentIndex] * 1.0e-6, MINIMUM_COMPOSITION);
    }
    if (heaviest >= 0) {
      seed[heaviest] = 1.0;
    }
    normalize(seed);
    return seed;
  }

  private double[] overallPerturbedTowardHeaviest(double fraction) {
    double[] seed = overallComposition();
    int heaviest = heaviestHydrocarbonIndex();
    for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
      seed[componentIndex] = Math.max((1.0 - fraction) * seed[componentIndex], MINIMUM_COMPOSITION);
    }
    if (heaviest >= 0) {
      seed[heaviest] += fraction;
    }
    normalize(seed);
    return seed;
  }

  private int heaviestHydrocarbonIndex() {
    int index = -1;
    double maximumCriticalTemperature = Double.NEGATIVE_INFINITY;
    for (int componentIndex = 0; componentIndex < template.getPhase(0).getNumberOfComponents(); componentIndex++) {
      ComponentInterface component = template.getPhase(0).getComponent(componentIndex);
      if (component.isHydrocarbon() && component.getTC() > maximumCriticalTemperature) {
        maximumCriticalTemperature = component.getTC();
        index = componentIndex;
      }
    }
    return index;
  }

  private static void addDistinctSeed(java.util.List<double[]> seeds, double[] candidate) {
    if (candidate == null) {
      return;
    }
    for (double[] seed : seeds) {
      double distance = 0.0;
      for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
        distance += Math.abs(seed[componentIndex] - candidate[componentIndex]);
      }
      if (distance <= 1.0e-7) {
        return;
      }
    }
    seeds.add(candidate);
  }

  private Trial evaluateBySuccessiveSubstitution(double temperatureK, double pressureBara, double[] compositionSeed) {
    try {
      SystemInterface mother = createMother(temperatureK, pressureBara);
      Candidate candidate = new IncipientPhaseStabilityAnalyzer(mother).setMaximumIterations(300).setDampingFactor(0.35)
          .analyzeCandidate(incipientPhase, compositionSeed);
      return new Trial(temperatureK, candidate.getComposition(), candidate.getTangentPlaneDistance(),
          candidate.getStationarityResidual(), candidate.isConverged(), candidate.isTrivial(), candidate.getPhase());
    } catch (RuntimeException error) {
      return new Trial(temperatureK, new double[0], Double.NaN, Double.POSITIVE_INFINITY, false, false, null);
    }
  }

  private Trial evaluateByNewton(double temperatureK, double pressureBara, double[] compositionSeed) {
    try {
      IncipientPhaseStationaryPointSolver.Result candidate = new IncipientPhaseStationaryPointSolver(
          createMother(temperatureK, pressureBara), incipientPhase).setNumericalControls(80, 1.0e-9, 2.0e-5)
          .solve(compositionSeed);
      return new Trial(temperatureK, candidate.getComposition(), candidate.getTangentPlaneDistance(),
          candidate.getStationarityResidual(), candidate.isConverged(), candidate.isTrivial(),
          candidate.getPhysicalPhase());
    } catch (RuntimeException error) {
      return new Trial(temperatureK, new double[0], Double.NaN, Double.POSITIVE_INFINITY, false, false, null);
    }
  }

  private SystemInterface createMother(double temperatureK, double pressureBara) {
    SystemInterface mother = template.clone();
    mother.setNumberOfPhases(1);
    mother.setTemperature(temperatureK);
    mother.setPressure(pressureBara);
    mother.setPhaseType(0, toPhaseType(motherPhase));
    for (int componentIndex = 0; componentIndex < mother.getPhase(0).getNumberOfComponents(); componentIndex++) {
      mother.getPhase(0).getComponent(componentIndex)
          .setx(Math.max(mother.getPhase(0).getComponent(componentIndex).getz(), MINIMUM_COMPOSITION));
    }
    mother.getPhase(0).normalize();
    mother.init(1);
    return mother;
  }

  private double[] physicalSeed(double temperatureK, double pressureBara) {
    double[] composition = new double[template.getPhase(0).getNumberOfComponents()];
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      ComponentInterface component = template.getPhase(0).getComponent(componentIndex);
      double z = Math.max(component.getz(), MINIMUM_COMPOSITION);
      double wilsonK = component.getPC() / pressureBara
          * Math.exp(5.373 * (1.0 + component.getAcentricFactor()) * (1.0 - component.getTC() / temperatureK));
      wilsonK = Math.max(1.0e-20, Math.min(1.0e20, wilsonK));
      if (incipientPhase == CandidatePhase.GAS) {
        composition[componentIndex] = z * wilsonK;
      } else if (incipientPhase == CandidatePhase.OIL) {
        composition[componentIndex] = z / wilsonK;
      } else if (component.getComponentName().equalsIgnoreCase("water")) {
        composition[componentIndex] = Math.max(z, 0.99);
      } else if (component.isHydrocarbon()) {
        composition[componentIndex] = z * 1.0e-8;
      } else {
        composition[componentIndex] = z;
      }
      composition[componentIndex] = Math.max(composition[componentIndex], MINIMUM_COMPOSITION);
    }
    normalize(composition);
    return composition;
  }

  private static boolean hasSignChange(double first, double second) {
    return Double.isFinite(first) && Double.isFinite(second)
        && (first <= 0.0 && second >= 0.0 || first >= 0.0 && second <= 0.0);
  }

  private static PhaseType toPhaseType(CandidatePhase phase) {
    switch (phase) {
    case GAS:
      return PhaseType.GAS;
    case OIL:
      return PhaseType.OIL;
    case AQUEOUS:
      return PhaseType.AQUEOUS;
    default:
      throw new IllegalArgumentException("unsupported phase family " + phase);
    }
  }

  private static void normalize(double[] composition) {
    double total = 0.0;
    for (double value : composition) {
      total += value;
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalArgumentException("composition cannot be normalized");
    }
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] /= total;
    }
  }

  private final class Trial {
    private final double temperatureK;
    private final double[] composition;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;
    private final boolean converged;
    private final boolean trivial;
    private final CandidatePhase physicalPhase;

    private Trial(double temperatureK, double[] composition, double tangentPlaneDistance, double stationarityResidual,
        boolean converged, boolean trivial, CandidatePhase physicalPhase) {
      this.temperatureK = temperatureK;
      this.composition = composition;
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.stationarityResidual = stationarityResidual;
      this.converged = converged;
      this.trivial = trivial;
      this.physicalPhase = physicalPhase;
    }

    private boolean isUsable() {
      return converged && !trivial && physicalPhase == incipientPhase && Double.isFinite(tangentPlaneDistance);
    }
  }

  /** Immutable fixed-pressure zero-TPD boundary point. */
  public static final class Result {
    private final CandidatePhase motherPhase;
    private final CandidatePhase incipientPhase;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] incipientComposition;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;
    private final double temperatureBracketWidthK;
    private final int evaluations;
    private final int bisections;
    private final boolean converged;
    private final String failureMessage;

    private Result(CandidatePhase motherPhase, CandidatePhase incipientPhase, double temperatureK, double pressureBara,
        double[] incipientComposition, double tangentPlaneDistance, double stationarityResidual,
        double temperatureBracketWidthK, int evaluations, int bisections, boolean converged, String failureMessage) {
      this.motherPhase = motherPhase;
      this.incipientPhase = incipientPhase;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.incipientComposition = incipientComposition.clone();
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.stationarityResidual = stationarityResidual;
      this.temperatureBracketWidthK = temperatureBracketWidthK;
      this.evaluations = evaluations;
      this.bisections = bisections;
      this.converged = converged;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase motherPhase, CandidatePhase incipientPhase, double pressureBara,
        int evaluations, String failureMessage) {
      return new Result(motherPhase, incipientPhase, Double.NaN, pressureBara, new double[0], Double.NaN,
          Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, evaluations, 0, false, failureMessage);
    }

    /** @return true when a non-trivial zero-TPD root reached the requested tolerance */
    public boolean isConverged() {
      return converged && failureMessage == null;
    }

    public CandidatePhase getMotherPhase() {
      return motherPhase;
    }

    public CandidatePhase getIncipientPhase() {
      return incipientPhase;
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double[] getIncipientComposition() {
      return incipientComposition.clone();
    }

    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    public double getStationarityResidual() {
      return stationarityResidual;
    }

    public double getTemperatureBracketWidthK() {
      return temperatureBracketWidthK;
    }

    public int getEvaluations() {
      return evaluations;
    }

    public int getBisections() {
      return bisections;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
