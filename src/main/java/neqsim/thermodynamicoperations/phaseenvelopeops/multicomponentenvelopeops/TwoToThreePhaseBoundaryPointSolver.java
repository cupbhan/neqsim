package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Finds a specified two-phase to three-phase boundary point at fixed pressure.
 *
 * <p>
 * The retained two phases are solved by {@link SpecifiedTwoPhaseFlashSolver}; the third phase is followed as a
 * non-trivial Michelsen stationary point. A root is accepted only when the specified two-phase flash converges and the
 * requested third phase has zero tangent-plane distance. This prevents an ordinary multiphase flash from silently
 * replacing the aqueous phase by an oil phase while the boundary is being traced.
 * </p>
 */
public final class TwoToThreePhaseBoundaryPointSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;

  private final SystemInterface template;
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase incipientPhase;
  private int scanIntervals = 32;
  private int maximumBisections = 80;
  private double temperatureToleranceK = 1.0e-5;
  private double tangentPlaneTolerance = 1.0e-8;

  /**
   * Creates a fixed-pressure two-to-three-phase boundary solver.
   *
   * @param template configured EOS and mixing-rule system
   * @param retainedPhaseZero first retained phase; its mole fraction is beta
   * @param retainedPhaseOne second retained phase
   * @param incipientPhase third phase whose TPD is driven to zero
   */
  public TwoToThreePhaseBoundaryPointSolver(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase) {
    if (template == null || retainedPhaseZero == null || retainedPhaseOne == null || incipientPhase == null) {
      throw new IllegalArgumentException("template and all three phase families must be specified");
    }
    if (retainedPhaseZero == retainedPhaseOne || retainedPhaseZero == incipientPhase
        || retainedPhaseOne == incipientPhase) {
      throw new IllegalArgumentException("the three physical phase families must be distinct");
    }
    this.template = template.clone();
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.incipientPhase = incipientPhase;
  }

  /** Sets scan, bisection, and root tolerances. */
  public TwoToThreePhaseBoundaryPointSolver setNumericalControls(int scanIntervals, int maximumBisections,
      double temperatureToleranceK, double tangentPlaneTolerance) {
    if (scanIntervals < 2 || maximumBisections < 1 || !Double.isFinite(temperatureToleranceK)
        || temperatureToleranceK <= 0.0 || !Double.isFinite(tangentPlaneTolerance) || tangentPlaneTolerance <= 0.0) {
      throw new IllegalArgumentException("invalid two-to-three-phase boundary controls");
    }
    this.scanIntervals = scanIntervals;
    this.maximumBisections = maximumBisections;
    this.temperatureToleranceK = temperatureToleranceK;
    this.tangentPlaneTolerance = tangentPlaneTolerance;
    return this;
  }

  /**
   * Locates one fixed-pressure two-to-three-phase onset.
   *
   * @param pressureBara fixed pressure in bara
   * @param minimumTemperatureK lower scan temperature
   * @param maximumTemperatureK upper scan temperature
   * @param betaSeed retained-phase-zero fraction at the upper-temperature start
   * @param phaseZeroCompositionSeed retained-phase-zero composition seed
   * @param phaseOneCompositionSeed retained-phase-one composition seed
   * @return immutable boundary point or an explicit failure result
   */
  public Result solve(double pressureBara, double minimumTemperatureK, double maximumTemperatureK, double betaSeed,
      double[] phaseZeroCompositionSeed, double[] phaseOneCompositionSeed) {
    RootSet rootSet = solveRoots(pressureBara, minimumTemperatureK, maximumTemperatureK, betaSeed,
        phaseZeroCompositionSeed, phaseOneCompositionSeed, false, Collections.<double[]>emptyList());
    if (!rootSet.roots.isEmpty()) {
      return rootSet.roots.get(0);
    }
    return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase, pressureBara, rootSet.flashEvaluations,
        rootSet.stabilityEvaluations, rootSet.failureMessage);
  }

  /**
   * Locates every resolved fixed-pressure onset on the continuously traceable third-phase stationary branch.
   *
   * <p>
   * Unlike {@link #solve(double, double, double, double, double[], double[])}, this method does not stop after the
   * first TPD sign change. All disjoint zero brackets found over the requested temperature interval are refined and
   * returned in increasing-temperature order. Gaps where the requested stationary point is not usable break continuity,
   * so a sign change is never inferred across a missing segment.
   * </p>
   *
   * @param pressureBara fixed pressure in bara
   * @param minimumTemperatureK lower scan temperature
   * @param maximumTemperatureK upper scan temperature
   * @param betaSeed retained-phase-zero fraction at the upper-temperature start
   * @param phaseZeroCompositionSeed retained-phase-zero composition seed
   * @param phaseOneCompositionSeed retained-phase-one composition seed
   * @return all converged roots plus scan diagnostics
   */
  public RootSet solveAll(double pressureBara, double minimumTemperatureK, double maximumTemperatureK, double betaSeed,
      double[] phaseZeroCompositionSeed, double[] phaseOneCompositionSeed) {
    return solveRoots(pressureBara, minimumTemperatureK, maximumTemperatureK, betaSeed, phaseZeroCompositionSeed,
        phaseOneCompositionSeed, true, Collections.<double[]>emptyList());
  }

  /**
   * Locates all fixed-pressure roots and adds externally discovered incipient-composition seeds.
   *
   * <p>
   * The additional seeds are intended for a coarse global TP stability/flash scanner. They let an isolated TPD
   * stationary branch enter the rigorous zero-TPD continuation without treating grid points as final envelope data.
   * Every returned root is still independently corrected and must pass the same residual gates.
   * </p>
   */
  public RootSet solveAll(double pressureBara, double minimumTemperatureK, double maximumTemperatureK, double betaSeed,
      double[] phaseZeroCompositionSeed, double[] phaseOneCompositionSeed, List<double[]> additionalIncipientSeeds) {
    if (additionalIncipientSeeds == null) {
      throw new IllegalArgumentException("additional incipient seeds must not be null");
    }
    return solveRoots(pressureBara, minimumTemperatureK, maximumTemperatureK, betaSeed, phaseZeroCompositionSeed,
        phaseOneCompositionSeed, true, additionalIncipientSeeds);
  }

  private RootSet solveRoots(double pressureBara, double minimumTemperatureK, double maximumTemperatureK,
      double betaSeed, double[] phaseZeroCompositionSeed, double[] phaseOneCompositionSeed, boolean multistart,
      List<double[]> additionalIncipientSeeds) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    validateInputs(pressureBara, minimumTemperatureK, maximumTemperatureK, betaSeed, phaseZeroCompositionSeed,
        phaseOneCompositionSeed, componentCount);

    SpecifiedTwoPhaseFlashSolver flashSolver = new SpecifiedTwoPhaseFlashSolver(template, retainedPhaseZero,
        retainedPhaseOne).setNumericalControls(80, 1.0e-9, 2.0e-5);
    EvaluationCounter counter = new EvaluationCounter();
    List<RetainedState> descending = new ArrayList<RetainedState>();
    double currentBeta = betaSeed;
    double[] currentZero = phaseZeroCompositionSeed.clone();
    double[] currentOne = phaseOneCompositionSeed.clone();
    for (int scanIndex = scanIntervals; scanIndex >= 0; scanIndex--) {
      double fraction = (double) scanIndex / scanIntervals;
      double temperatureK = minimumTemperatureK + fraction * (maximumTemperatureK - minimumTemperatureK);
      SpecifiedTwoPhaseFlashSolver.Result flash = flashSolver.solve(temperatureK, pressureBara, currentBeta,
          currentZero, currentOne);
      counter.flashEvaluations++;
      if (!flash.isConverged()) {
        continue;
      }
      RetainedState retained = new RetainedState(flash);
      descending.add(retained);
      currentBeta = flash.getBeta();
      currentZero = flash.getPhaseZeroComposition();
      currentOne = flash.getPhaseOneComposition();
    }
    if (descending.size() < 2) {
      return RootSet.failure(pressureBara, counter, 0, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
          "specified retained two-phase branch could not be traced");
    }

    Collections.reverse(descending);
    List<StabilityBracket> brackets = new ArrayList<StabilityBracket>();
    List<StabilityState> exactRoots = new ArrayList<StabilityState>();
    List<StabilityTrack> stabilityTracks = stabilityTracks(multistart, additionalIncipientSeeds, componentCount);
    int usableStabilityStates = 0;
    double minimumTangentPlaneDistance = Double.POSITIVE_INFINITY;
    double maximumTangentPlaneDistance = Double.NEGATIVE_INFINITY;
    for (RetainedState retained : descending) {
      for (StabilityTrack track : stabilityTracks) {
        double[] restartSeed = seedFor(track, retained, pressureBara);
        double[] requestedSeed = track.explicitSeed != null || track.previous == null ? restartSeed
            : track.previous.incipientComposition;
        StabilityState current = evaluateStability(flashSolver, retained, requestedSeed, multistart);
        counter.stabilityEvaluations++;
        if (!current.isUsable(incipientPhase) && track.previous != null) {
          StabilityState restarted = evaluateStability(flashSolver, retained, restartSeed, multistart);
          counter.stabilityEvaluations++;
          if (restarted.isUsable(incipientPhase)) {
            current = restarted;
          }
        }
        if (!current.isUsable(incipientPhase)) {
          track.previous = null;
          continue;
        }
        usableStabilityStates++;
        minimumTangentPlaneDistance = Math.min(minimumTangentPlaneDistance, current.tangentPlaneDistance);
        maximumTangentPlaneDistance = Math.max(maximumTangentPlaneDistance, current.tangentPlaneDistance);
        if (Math.abs(current.tangentPlaneDistance) <= tangentPlaneTolerance) {
          exactRoots.add(current);
        }
        if (track.previous != null
            && hasStrictSignChange(track.previous.tangentPlaneDistance, current.tangentPlaneDistance)) {
          StabilityState lower = track.previous.temperatureK <= current.temperatureK ? track.previous : current;
          StabilityState upper = track.previous.temperatureK <= current.temperatureK ? current : track.previous;
          brackets.add(new StabilityBracket(lower, upper));
        }
        track.previous = current;
      }
    }

    List<Result> roots = new ArrayList<Result>();
    for (StabilityState exactRoot : exactRoots) {
      addDistinctRoot(roots, resultFromState(exactRoot, 0.0, counter, 0));
    }
    for (StabilityBracket bracket : brackets) {
      Result root = refineBracket(flashSolver, bracket.lower, bracket.upper, counter, multistart);
      if (root.isConverged()) {
        addDistinctRoot(roots, root);
      }
    }
    Collections.sort(roots, (first, second) -> Double.compare(first.temperatureK, second.temperatureK));
    String failureMessage = null;
    if (roots.isEmpty()) {
      failureMessage = String.format(Locale.ROOT,
          "no non-trivial third-phase TPD sign change found on the retained two-phase branch; "
              + "usable=%d, TPD range=[%.8g, %.8g]",
          usableStabilityStates, minimumTangentPlaneDistance, maximumTangentPlaneDistance);
    }
    return new RootSet(pressureBara, roots, counter.flashEvaluations, counter.stabilityEvaluations,
        usableStabilityStates, minimumTangentPlaneDistance, maximumTangentPlaneDistance, failureMessage);
  }

  private Result refineBracket(SpecifiedTwoPhaseFlashSolver flashSolver, StabilityState initialLower,
      StabilityState initialUpper, EvaluationCounter counter, boolean robustStationarySolver) {
    double pressureBara = initialLower.retained.flashResult.getPressureBara();
    StabilityState lower = initialLower;
    StabilityState upper = initialUpper;
    StabilityState best = Math.abs(lower.tangentPlaneDistance) <= Math.abs(upper.tangentPlaneDistance) ? lower : upper;
    int bisections = 0;
    while (bisections < maximumBisections && upper.temperatureK - lower.temperatureK > temperatureToleranceK
        && Math.abs(best.tangentPlaneDistance) > tangentPlaneTolerance) {
      bisections++;
      double midpointTemperatureK = 0.5 * (lower.temperatureK + upper.temperatureK);
      StabilityState seedSource = midpointTemperatureK - lower.temperatureK <= upper.temperatureK - midpointTemperatureK
          ? lower
          : upper;
      RetainedState retained = solveRetained(flashSolver, midpointTemperatureK, pressureBara, seedSource);
      counter.flashEvaluations++;
      if (retained == null) {
        return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase,
            initialLower.retained.flashResult.getPressureBara(), counter.flashEvaluations, counter.stabilityEvaluations,
            "specified retained two-phase flash was lost during bisection");
      }
      StabilityState midpoint = evaluateStability(flashSolver, retained, seedSource.incipientComposition,
          robustStationarySolver);
      counter.stabilityEvaluations++;
      if (!midpoint.isUsable(incipientPhase)) {
        midpoint = evaluateStability(flashSolver, retained, physicalSeed(midpointTemperatureK, pressureBara),
            robustStationarySolver);
        counter.stabilityEvaluations++;
      }
      if (!midpoint.isUsable(incipientPhase)) {
        return Result.failure(retainedPhaseZero, retainedPhaseOne, incipientPhase,
            initialLower.retained.flashResult.getPressureBara(), counter.flashEvaluations, counter.stabilityEvaluations,
            "non-trivial third-phase stationary point was lost during bisection");
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
      double secantTemperatureK = (lower.temperatureK * upper.tangentPlaneDistance
          - upper.temperatureK * lower.tangentPlaneDistance) / secantDenominator;
      if (secantTemperatureK > lower.temperatureK && secantTemperatureK < upper.temperatureK) {
        StabilityState seedSource = secantTemperatureK - lower.temperatureK <= upper.temperatureK - secantTemperatureK
            ? lower
            : upper;
        RetainedState retained = solveRetained(flashSolver, secantTemperatureK, pressureBara, seedSource);
        counter.flashEvaluations++;
        if (retained != null) {
          StabilityState secant = evaluateStability(flashSolver, retained, seedSource.incipientComposition,
              robustStationarySolver);
          counter.stabilityEvaluations++;
          if (secant.isUsable(incipientPhase)
              && Math.abs(secant.tangentPlaneDistance) < Math.abs(best.tangentPlaneDistance)) {
            best = secant;
          }
        }
      }
    }

    boolean converged = Math.abs(best.tangentPlaneDistance) <= tangentPlaneTolerance;
    Result result = resultFromState(best, upper.temperatureK - lower.temperatureK, counter, bisections);
    if (converged) {
      return result;
    }
    return new Result(retainedPhaseZero, retainedPhaseOne, incipientPhase, result.temperatureK, result.pressureBara,
        result.beta, result.phaseZeroComposition, result.phaseOneComposition, result.incipientComposition,
        result.retainedFlashResidual, result.tangentPlaneDistance, result.stationarityResidual,
        result.temperatureBracketWidthK, result.flashEvaluations, result.stabilityEvaluations, result.bisections, false,
        "maximum two-to-three-phase bisection count reached");
  }

  private Result resultFromState(StabilityState state, double bracketWidthK, EvaluationCounter counter,
      int bisections) {
    return new Result(retainedPhaseZero, retainedPhaseOne, incipientPhase, state.temperatureK,
        state.retained.flashResult.getPressureBara(), state.retained.beta, state.retained.phaseZeroComposition,
        state.retained.phaseOneComposition, state.incipientComposition, state.retained.maximumResidual,
        state.tangentPlaneDistance, state.stationarityResidual, bracketWidthK, counter.flashEvaluations,
        counter.stabilityEvaluations, bisections, true, null);
  }

  private void addDistinctRoot(List<Result> roots, Result candidate) {
    double temperatureTolerance = Math.max(10.0 * temperatureToleranceK, 1.0e-4);
    for (int rootIndex = 0; rootIndex < roots.size(); rootIndex++) {
      Result current = roots.get(rootIndex);
      if (Math.abs(current.temperatureK - candidate.temperatureK) <= temperatureTolerance) {
        if (Math.abs(candidate.tangentPlaneDistance) < Math.abs(current.tangentPlaneDistance)) {
          roots.set(rootIndex, candidate);
        }
        return;
      }
    }
    roots.add(candidate);
  }

  private static List<StabilityTrack> stabilityTracks(boolean multistart, List<double[]> additionalIncipientSeeds,
      int componentCount) {
    List<StabilityTrack> tracks = new ArrayList<StabilityTrack>();
    tracks.add(new StabilityTrack(SeedKind.PHYSICAL));
    if (multistart) {
      tracks.add(new StabilityTrack(SeedKind.RETAINED_HEAVY_05));
      tracks.add(new StabilityTrack(SeedKind.OVERALL));
      tracks.add(new StabilityTrack(SeedKind.HEAVY_ENRICHED));
      tracks.add(new StabilityTrack(SeedKind.HEAVIEST_COMPONENT));
    }
    for (double[] additionalSeed : additionalIncipientSeeds) {
      if (additionalSeed == null || additionalSeed.length != componentCount) {
        throw new IllegalArgumentException("each additional incipient seed must contain one value per component");
      }
      double[] normalizedSeed = additionalSeed.clone();
      for (double value : normalizedSeed) {
        if (!Double.isFinite(value) || value < 0.0) {
          throw new IllegalArgumentException("additional incipient seed contains an invalid value");
        }
      }
      normalize(normalizedSeed);
      tracks.add(new StabilityTrack(normalizedSeed));
    }
    return tracks;
  }

  private double[] seedFor(StabilityTrack track, RetainedState retained, double pressureBara) {
    if (track.explicitSeed != null) {
      return track.explicitSeed.clone();
    }
    switch (track.seedKind) {
    case PHYSICAL:
      return physicalSeed(retained.temperatureK, pressureBara);
    case RETAINED_HEAVY_05:
      return phaseZeroPerturbedTowardHeaviest(retained.phaseZeroComposition, 0.05);
    case OVERALL:
      return overallComposition();
    case HEAVY_ENRICHED:
      return heavyEnrichedSeed(retained.phaseZeroComposition);
    case HEAVIEST_COMPONENT:
      return heaviestComponentSeed();
    default:
      throw new IllegalArgumentException("unsupported stability seed family " + track.seedKind);
    }
  }

  private double[] overallComposition() {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    double[] composition = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(),
          MINIMUM_COMPOSITION);
    }
    normalize(composition);
    return composition;
  }

  private double[] heavyEnrichedSeed(double[] baseComposition) {
    double[] seed = baseComposition.clone();
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

  private double[] heaviestComponentSeed() {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    int heaviestHydrocarbonIndex = -1;
    double maximumCriticalTemperature = Double.NEGATIVE_INFINITY;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      ComponentInterface component = template.getPhase(0).getComponent(componentIndex);
      if (component.isHydrocarbon() && component.getTC() > maximumCriticalTemperature) {
        maximumCriticalTemperature = component.getTC();
        heaviestHydrocarbonIndex = componentIndex;
      }
    }
    double[] seed = overallComposition();
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      seed[componentIndex] = Math.max(seed[componentIndex] * 1.0e-6, MINIMUM_COMPOSITION);
    }
    if (heaviestHydrocarbonIndex >= 0) {
      seed[heaviestHydrocarbonIndex] = 1.0;
    }
    normalize(seed);
    return seed;
  }

  private double[] phaseZeroPerturbedTowardHeaviest(double[] phaseZeroComposition, double fraction) {
    double[] seed = phaseZeroComposition.clone();
    int heaviestHydrocarbonIndex = heaviestHydrocarbonIndex();
    for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
      seed[componentIndex] *= 1.0 - fraction;
      seed[componentIndex] = Math.max(seed[componentIndex], MINIMUM_COMPOSITION);
    }
    if (heaviestHydrocarbonIndex >= 0) {
      seed[heaviestHydrocarbonIndex] += fraction;
    }
    normalize(seed);
    return seed;
  }

  private int heaviestHydrocarbonIndex() {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    int heaviestHydrocarbonIndex = -1;
    double maximumCriticalTemperature = Double.NEGATIVE_INFINITY;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      ComponentInterface component = template.getPhase(0).getComponent(componentIndex);
      if (component.isHydrocarbon() && component.getTC() > maximumCriticalTemperature) {
        maximumCriticalTemperature = component.getTC();
        heaviestHydrocarbonIndex = componentIndex;
      }
    }
    return heaviestHydrocarbonIndex;
  }

  private RetainedState solveRetained(SpecifiedTwoPhaseFlashSolver flashSolver, double temperatureK,
      double pressureBara, StabilityState seed) {
    SpecifiedTwoPhaseFlashSolver.Result flash = flashSolver.solve(temperatureK, pressureBara, seed.retained.beta,
        seed.retained.phaseZeroComposition, seed.retained.phaseOneComposition);
    return flash.isConverged() ? new RetainedState(flash) : null;
  }

  private StabilityState evaluateStability(SpecifiedTwoPhaseFlashSolver flashSolver, RetainedState retained,
      double[] incipientCompositionSeed, boolean robustStationarySolver) {
    StabilityState successiveSubstitution = evaluateStabilityBySuccessiveSubstitution(flashSolver, retained,
        incipientCompositionSeed);
    if (!robustStationarySolver) {
      return successiveSubstitution;
    }
    StabilityState newton = evaluateStabilityByNewton(flashSolver, retained, incipientCompositionSeed);
    boolean successiveUsable = successiveSubstitution.isUsable(incipientPhase);
    boolean newtonUsable = newton.isUsable(incipientPhase);
    if (successiveUsable && newtonUsable) {
      double successiveDistance = compositionDistance(successiveSubstitution.incipientComposition,
          incipientCompositionSeed);
      double newtonDistance = compositionDistance(newton.incipientComposition, incipientCompositionSeed);
      return newtonDistance < successiveDistance ? newton : successiveSubstitution;
    }
    if (newtonUsable) {
      return newton;
    }
    return successiveSubstitution;
  }

  private StabilityState evaluateStabilityBySuccessiveSubstitution(SpecifiedTwoPhaseFlashSolver flashSolver,
      RetainedState retained, double[] incipientCompositionSeed) {
    try {
      IncipientPhaseStabilityAnalyzer analyzer = new IncipientPhaseStabilityAnalyzer(
          flashSolver.toThermodynamicSystem(retained.flashResult)).setMaximumIterations(500).setDampingFactor(0.2);
      Candidate candidate = analyzer.analyzeCandidate(incipientPhase, incipientCompositionSeed);
      return new StabilityState(retained, candidate.getComposition(), candidate.getTangentPlaneDistance(),
          candidate.getStationarityResidual(), candidate.isConverged(), candidate.isTrivial(), candidate.getPhase());
    } catch (RuntimeException error) {
      return new StabilityState(retained, new double[0], Double.NaN, Double.POSITIVE_INFINITY, false, false, null);
    }
  }

  private StabilityState evaluateStabilityByNewton(SpecifiedTwoPhaseFlashSolver flashSolver, RetainedState retained,
      double[] incipientCompositionSeed) {
    try {
      IncipientPhaseStationaryPointSolver.Result candidate = new IncipientPhaseStationaryPointSolver(
          flashSolver.toThermodynamicSystem(retained.flashResult), incipientPhase)
          .setNumericalControls(80, 1.0e-9, 2.0e-5).solve(incipientCompositionSeed);
      return new StabilityState(retained, candidate.getComposition(), candidate.getTangentPlaneDistance(),
          candidate.getStationarityResidual(), candidate.isConverged(), candidate.isTrivial(),
          candidate.getPhysicalPhase());
    } catch (RuntimeException error) {
      return new StabilityState(retained, new double[0], Double.NaN, Double.POSITIVE_INFINITY, false, false, null);
    }
  }

  private static double compositionDistance(double[] first, double[] second) {
    if (first.length != second.length) {
      return Double.POSITIVE_INFINITY;
    }
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

  private double[] physicalSeed(double temperatureK, double pressureBara) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    double[] seed = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      ComponentInterface component = template.getPhase(0).getComponent(componentIndex);
      double z = Math.max(component.getz(), MINIMUM_COMPOSITION);
      double wilsonK = component.getPC() / pressureBara
          * Math.exp(5.373 * (1.0 + component.getAcentricFactor()) * (1.0 - component.getTC() / temperatureK));
      wilsonK = Math.max(1.0e-20, Math.min(1.0e20, wilsonK));
      if (incipientPhase == CandidatePhase.GAS) {
        seed[componentIndex] = z * wilsonK;
      } else if (incipientPhase == CandidatePhase.OIL) {
        seed[componentIndex] = component.getComponentName().equalsIgnoreCase("water") ? z * 1.0e-12 : z / wilsonK;
      } else if (component.getComponentName().equalsIgnoreCase("water")) {
        seed[componentIndex] = Math.max(z, 0.99);
      } else if (component.isHydrocarbon()) {
        seed[componentIndex] = z * 1.0e-8;
      } else {
        seed[componentIndex] = z;
      }
      seed[componentIndex] = Math.max(seed[componentIndex], MINIMUM_COMPOSITION);
    }
    normalize(seed);
    return seed;
  }

  private static void validateInputs(double pressureBara, double minimumTemperatureK, double maximumTemperatureK,
      double betaSeed, double[] phaseZeroCompositionSeed, double[] phaseOneCompositionSeed, int componentCount) {
    if (!Double.isFinite(pressureBara) || pressureBara <= 0.0 || !Double.isFinite(minimumTemperatureK)
        || minimumTemperatureK < 50.0 || !Double.isFinite(maximumTemperatureK)
        || maximumTemperatureK <= minimumTemperatureK || !Double.isFinite(betaSeed) || betaSeed <= 0.0
        || betaSeed >= 1.0 || phaseZeroCompositionSeed == null || phaseOneCompositionSeed == null
        || phaseZeroCompositionSeed.length != componentCount || phaseOneCompositionSeed.length != componentCount) {
      throw new IllegalArgumentException("invalid fixed-pressure two-to-three-phase boundary input");
    }
  }

  private static boolean hasSignChange(double first, double second) {
    return Double.isFinite(first) && Double.isFinite(second)
        && (first <= 0.0 && second >= 0.0 || first >= 0.0 && second <= 0.0);
  }

  private static boolean hasStrictSignChange(double first, double second) {
    return Double.isFinite(first) && Double.isFinite(second)
        && (first < 0.0 && second > 0.0 || first > 0.0 && second < 0.0);
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

  private static final class RetainedState {
    private final SpecifiedTwoPhaseFlashSolver.Result flashResult;
    private final double temperatureK;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double maximumResidual;

    private RetainedState(SpecifiedTwoPhaseFlashSolver.Result flashResult) {
      this.flashResult = flashResult;
      this.temperatureK = flashResult.getTemperatureK();
      this.beta = flashResult.getBeta();
      this.phaseZeroComposition = flashResult.getPhaseZeroComposition();
      this.phaseOneComposition = flashResult.getPhaseOneComposition();
      this.maximumResidual = flashResult.getMaximumResidual();
    }
  }

  private static final class StabilityBracket {
    private final StabilityState lower;
    private final StabilityState upper;

    private StabilityBracket(StabilityState lower, StabilityState upper) {
      this.lower = lower;
      this.upper = upper;
    }
  }

  private static final class EvaluationCounter {
    private int flashEvaluations;
    private int stabilityEvaluations;
  }

  private enum SeedKind {
    PHYSICAL, RETAINED_HEAVY_05, OVERALL, HEAVY_ENRICHED, HEAVIEST_COMPONENT
  }

  private static final class StabilityTrack {
    private final SeedKind seedKind;
    private final double[] explicitSeed;
    private StabilityState previous;

    private StabilityTrack(SeedKind seedKind) {
      this.seedKind = seedKind;
      this.explicitSeed = null;
    }

    private StabilityTrack(double[] explicitSeed) {
      this.seedKind = null;
      this.explicitSeed = explicitSeed.clone();
    }
  }

  private static final class StabilityState {
    private final RetainedState retained;
    private final double temperatureK;
    private final double[] incipientComposition;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;
    private final boolean converged;
    private final boolean trivial;
    private final CandidatePhase physicalPhase;

    private StabilityState(RetainedState retained, double[] incipientComposition, double tangentPlaneDistance,
        double stationarityResidual, boolean converged, boolean trivial, CandidatePhase physicalPhase) {
      this.retained = retained;
      this.temperatureK = retained.temperatureK;
      this.incipientComposition = incipientComposition.clone();
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.stationarityResidual = stationarityResidual;
      this.converged = converged;
      this.trivial = trivial;
      this.physicalPhase = physicalPhase;
    }

    private boolean isUsable(CandidatePhase expectedPhase) {
      return converged && !trivial && physicalPhase == expectedPhase && Double.isFinite(tangentPlaneDistance);
    }
  }

  /** Immutable set of every resolved root and the full fixed-pressure scan diagnostics. */
  public static final class RootSet {
    private final double pressureBara;
    private final List<Result> roots;
    private final int flashEvaluations;
    private final int stabilityEvaluations;
    private final int usableStabilityStates;
    private final double minimumTangentPlaneDistance;
    private final double maximumTangentPlaneDistance;
    private final String failureMessage;

    private RootSet(double pressureBara, List<Result> roots, int flashEvaluations, int stabilityEvaluations,
        int usableStabilityStates, double minimumTangentPlaneDistance, double maximumTangentPlaneDistance,
        String failureMessage) {
      this.pressureBara = pressureBara;
      this.roots = Collections.unmodifiableList(new ArrayList<Result>(roots));
      this.flashEvaluations = flashEvaluations;
      this.stabilityEvaluations = stabilityEvaluations;
      this.usableStabilityStates = usableStabilityStates;
      this.minimumTangentPlaneDistance = minimumTangentPlaneDistance;
      this.maximumTangentPlaneDistance = maximumTangentPlaneDistance;
      this.failureMessage = failureMessage;
    }

    private static RootSet failure(double pressureBara, EvaluationCounter counter, int usableStabilityStates,
        double minimumTangentPlaneDistance, double maximumTangentPlaneDistance, String failureMessage) {
      return new RootSet(pressureBara, Collections.<Result>emptyList(), counter.flashEvaluations,
          counter.stabilityEvaluations, usableStabilityStates, minimumTangentPlaneDistance, maximumTangentPlaneDistance,
          failureMessage);
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public List<Result> getRoots() {
      return roots;
    }

    public boolean hasConvergedRoots() {
      return !roots.isEmpty();
    }

    public int getFlashEvaluations() {
      return flashEvaluations;
    }

    public int getStabilityEvaluations() {
      return stabilityEvaluations;
    }

    public int getUsableStabilityStates() {
      return usableStabilityStates;
    }

    public double getMinimumTangentPlaneDistance() {
      return minimumTangentPlaneDistance;
    }

    public double getMaximumTangentPlaneDistance() {
      return maximumTangentPlaneDistance;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }

  /** Immutable fixed-pressure specified two-to-three-phase boundary point. */
  public static final class Result {
    private final CandidatePhase retainedPhaseZero;
    private final CandidatePhase retainedPhaseOne;
    private final CandidatePhase incipientPhase;
    private final double temperatureK;
    private final double pressureBara;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double[] incipientComposition;
    private final double retainedFlashResidual;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;
    private final double temperatureBracketWidthK;
    private final int flashEvaluations;
    private final int stabilityEvaluations;
    private final int bisections;
    private final boolean converged;
    private final String failureMessage;

    private Result(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase,
        double temperatureK, double pressureBara, double beta, double[] phaseZeroComposition,
        double[] phaseOneComposition, double[] incipientComposition, double retainedFlashResidual,
        double tangentPlaneDistance, double stationarityResidual, double temperatureBracketWidthK, int flashEvaluations,
        int stabilityEvaluations, int bisections, boolean converged, String failureMessage) {
      this.retainedPhaseZero = retainedPhaseZero;
      this.retainedPhaseOne = retainedPhaseOne;
      this.incipientPhase = incipientPhase;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.beta = beta;
      this.phaseZeroComposition = phaseZeroComposition.clone();
      this.phaseOneComposition = phaseOneComposition.clone();
      this.incipientComposition = incipientComposition.clone();
      this.retainedFlashResidual = retainedFlashResidual;
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.stationarityResidual = stationarityResidual;
      this.temperatureBracketWidthK = temperatureBracketWidthK;
      this.flashEvaluations = flashEvaluations;
      this.stabilityEvaluations = stabilityEvaluations;
      this.bisections = bisections;
      this.converged = converged;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase incipientPhase, double pressureBara, int flashEvaluations, int stabilityEvaluations,
        String failureMessage) {
      return new Result(retainedPhaseZero, retainedPhaseOne, incipientPhase, Double.NaN, pressureBara, Double.NaN,
          new double[0], new double[0], new double[0], Double.NaN, Double.NaN, Double.POSITIVE_INFINITY,
          Double.POSITIVE_INFINITY, flashEvaluations, stabilityEvaluations, 0, false, failureMessage);
    }

    /**
     * Converts an already converged pseudo-arclength state into the common boundary-state evidence contract.
     *
     * <p>
     * This does not claim a fixed-pressure solve. It only exposes the identical PT, beta, and phase-composition state
     * to independent classifiers such as the global stability gate.
     * </p>
     */
    public static Result fromContinuationState(TwoToThreePhaseArcLengthCorrector.State state,
        double thermodynamicMaximumResidual) {
      if (state == null || !Double.isFinite(thermodynamicMaximumResidual) || thermodynamicMaximumResidual < 0.0) {
        throw new IllegalArgumentException("a finite converged continuation state and residual are required");
      }
      return new Result(state.getRetainedPhaseZero(), state.getRetainedPhaseOne(), state.getIncipientPhase(),
          state.getTemperatureK(), state.getPressureBara(), state.getBeta(), state.getPhaseZeroComposition(),
          state.getPhaseOneComposition(), state.getIncipientComposition(), thermodynamicMaximumResidual,
          thermodynamicMaximumResidual, thermodynamicMaximumResidual, 0.0, 0, 0, 0, true, null);
    }

    /** @return true when both retained equilibrium and zero third-phase TPD satisfy the requested tolerances */
    public boolean isConverged() {
      return converged && failureMessage == null;
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

    public double getTemperatureBracketWidthK() {
      return temperatureBracketWidthK;
    }

    public int getFlashEvaluations() {
      return flashEvaluations;
    }

    public int getStabilityEvaluations() {
      return stabilityEvaluations;
    }

    public int getBisections() {
      return bisections;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
