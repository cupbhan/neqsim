package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Corrects a post-spinodal three-retained-phase/fourth-incipient-phase boundary point.
 *
 * <p>
 * The retained topology is typically {@code GAS/OIL/OIL}, while an aqueous fourth phase remains incipient. Temperature
 * and pressure are solved together with the three retained phase compositions and fractions. A signed composition
 * separation along the retained-phase spinodal eigenvector selects a non-coincident branch and closes the otherwise
 * one-dimensional boundary equations.
 * </p>
 */
public final class HydrocarbonWaterThreeToFourPhaseBoundaryPointSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;

  private final SystemInterface template;
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase bifurcatingPhase;
  private final CandidatePhase incipientPhase;
  private int maximumIterations = 240;
  private double residualTolerance = 1.0e-8;
  private double finiteDifferenceStep = 1.0e-5;
  private double initialTrustRadius = 0.25;
  private double maximumTrustRadius = 40.0;

  /** Creates a specified post-spinodal boundary corrector. */
  public HydrocarbonWaterThreeToFourPhaseBoundaryPointSolver(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase bifurcatingPhase, CandidatePhase incipientPhase) {
    if (template == null || retainedPhaseZero == null || retainedPhaseOne == null || bifurcatingPhase == null
        || incipientPhase == null) {
      throw new IllegalArgumentException("template and four phase identities are required");
    }
    if (retainedPhaseZero == retainedPhaseOne
        || bifurcatingPhase != retainedPhaseZero && bifurcatingPhase != retainedPhaseOne
        || incipientPhase == retainedPhaseZero || incipientPhase == retainedPhaseOne) {
      throw new IllegalArgumentException("invalid retained, bifurcating, or incipient phase topology");
    }
    this.template = template.clone();
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.bifurcatingPhase = bifurcatingPhase;
    this.incipientPhase = incipientPhase;
  }

  /** Sets trust-region and finite-difference controls. */
  public HydrocarbonWaterThreeToFourPhaseBoundaryPointSolver setNumericalControls(int maximumIterations,
      double residualTolerance, double finiteDifferenceStep, double initialTrustRadius, double maximumTrustRadius) {
    if (maximumIterations < 1 || !positive(residualTolerance) || !positive(finiteDifferenceStep)
        || !positive(initialTrustRadius) || !positive(maximumTrustRadius) || initialTrustRadius > maximumTrustRadius) {
      throw new IllegalArgumentException("invalid three-to-four-phase boundary controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    this.initialTrustRadius = initialTrustRadius;
    this.maximumTrustRadius = maximumTrustRadius;
    return this;
  }

  /**
   * Corrects one finite separation from a retained-phase spinodal.
   *
   * @param spinodalState converged two-retained-phase/third-incipient boundary state
   * @param mode retained-phase spinodal eigenmode
   * @param targetSeparation signed log-composition separation along {@code mode}
   * @param bifurcatingFractionSeed positive total fraction seed for the new retained phase
   * @return corrected boundary point and numerical diagnostics
   */
  public Result solve(TwoToThreePhaseArcLengthCorrector.State spinodalState,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode, double targetSeparation, double bifurcatingFractionSeed) {
    validateInputs(spinodalState, mode, targetSeparation, bifurcatingFractionSeed);
    final int componentCount = componentCount();
    final int coordinateCount = componentCount - 1;
    InitialState initial = initialState(spinodalState, mode, targetSeparation, bifurcatingFractionSeed);
    double[] variables = initial.variables;
    int fractionStart = 3 * coordinateCount;
    int incipientStart = fractionStart + 2;
    int temperatureIndex = incipientStart + coordinateCount;
    int pressureIndex = temperatureIndex + 1;
    double[] lowerBounds = new double[variables.length];
    double[] upperBounds = new double[variables.length];
    double[] scales = new double[variables.length];
    for (int index = 0; index < variables.length; index++) {
      if (index < fractionStart) {
        lowerBounds[index] = -80.0;
        upperBounds[index] = 80.0;
        scales[index] = 1.0;
      } else if (index < incipientStart) {
        lowerBounds[index] = -32.0;
        upperBounds[index] = 32.0;
        scales[index] = 1.0;
      } else if (index < temperatureIndex) {
        lowerBounds[index] = -80.0;
        upperBounds[index] = 80.0;
        scales[index] = 1.0;
      } else if (index == temperatureIndex) {
        lowerBounds[index] = Math.log(50.0);
        upperBounds[index] = Math.log(2500.0);
        scales[index] = 0.03;
      } else {
        lowerBounds[index] = Math.log(1.0e-4);
        upperBounds[index] = Math.log(1.0e6);
        scales[index] = 0.05;
      }
    }
    final IncipientPhaseStationarityJacobianAnalyzer.Result fixedMode = mode;
    final double fixedSeparation = targetSeparation;
    final int[] fixedReferenceIndices = initial.referenceIndices;
    final int fixedMaterialBalanceOmittedIndex = largestComponentIndex(overallComposition());
    ScaledTrustRegionLeastSquaresSolver trustRegion = new ScaledTrustRegionLeastSquaresSolver(trial -> evaluate(trial,
        fixedReferenceIndices, fixedMaterialBalanceOmittedIndex, fixedMode, fixedSeparation).residual)
        .setNumericalControls(maximumIterations, residualTolerance, finiteDifferenceStep, initialTrustRadius,
            maximumTrustRadius)
        .setVariableScales(scales).setBounds(lowerBounds, upperBounds);
    ScaledTrustRegionLeastSquaresSolver.Result numerical = trustRegion.solve(variables);
    Evaluation evaluation;
    try {
      evaluation = evaluate(numerical.getVariables(), fixedReferenceIndices, fixedMaterialBalanceOmittedIndex, mode,
          targetSeparation);
    } catch (RuntimeException error) {
      return Result.failure(retainedPhaseZero, retainedPhaseOne, bifurcatingPhase, incipientPhase, componentCount,
          "final post-spinodal boundary evaluation failed: " + error.getMessage());
    }
    boolean distinct = compositionDistance(evaluation.duplicatedRetainedComposition,
        evaluation.bifurcatingComposition) > 1.0e-5;
    boolean finiteFractions = evaluation.phaseFractions[0] > 1.0e-12 && evaluation.phaseFractions[1] > 1.0e-12
        && evaluation.phaseFractions[2] > 1.0e-12;
    boolean physicalIdentity = phaseIdentity(evaluation.phaseCompositions, evaluation.incipientComposition);
    String failureMessage = numerical.getFailureMessage();
    if (failureMessage == null && !distinct) {
      failureMessage = "post-spinodal correction returned coincident retained compositions";
    }
    if (failureMessage == null && !finiteFractions) {
      failureMessage = "post-spinodal correction returned a vanishing retained phase";
    }
    if (failureMessage == null && !physicalIdentity) {
      failureMessage = "post-spinodal correction changed a requested physical phase family";
    }
    return new Result(retainedPhaseZero, retainedPhaseOne, bifurcatingPhase, incipientPhase, evaluation.temperatureK,
        evaluation.pressureBara, evaluation.phaseFractions, evaluation.phaseCompositions,
        evaluation.incipientComposition, evaluation.residual, numerical.getInitialMaximumResidual(),
        evaluation.maximumResidual, evaluation.equilibriumMaximumResidual, evaluation.incipientMaximumResidual,
        evaluation.separationResidual, materialBalanceResidual(evaluation.phaseFractions, evaluation.phaseCompositions),
        numerical.getIterations(), numerical.getAcceptedSteps(), numerical.getRejectedSteps(),
        numerical.getFinalTrustRadius(), numerical.getJacobianConditionNumber(), numerical.isConverged(), distinct,
        finiteFractions, physicalIdentity, failureMessage);
  }

  private InitialState initialState(TwoToThreePhaseArcLengthCorrector.State state,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode, double targetSeparation, double bifurcatingFractionSeed) {
    int componentCount = componentCount();
    int coordinateCount = componentCount - 1;
    int fractionStart = 3 * coordinateCount;
    int incipientStart = fractionStart + 2;
    double[] phaseZero = normalized(state.getPhaseZeroComposition());
    double[] phaseOne = normalized(state.getPhaseOneComposition());
    double parentFraction = bifurcatingPhase == retainedPhaseZero ? state.getBeta() : 1.0 - state.getBeta();
    double splitShare = bifurcatingFractionSeed / parentFraction;
    double[] parent = bifurcatingPhase == retainedPhaseZero ? phaseZero : phaseOne;
    double[] retainedSplit = perturbAlongMode(parent, mode, -splitShare * targetSeparation);
    double[] bifurcating = perturbAlongMode(parent, mode, (1.0 - splitShare) * targetSeparation);
    if (bifurcatingPhase == retainedPhaseZero) {
      phaseZero = retainedSplit;
    } else {
      phaseOne = retainedSplit;
    }
    double[] incipient = normalized(state.getIncipientComposition());
    int[] referenceIndices = new int[] { largestComponentIndex(phaseZero), largestComponentIndex(phaseOne),
        largestComponentIndex(bifurcating), largestComponentIndex(incipient) };
    double[] variables = new double[4 * componentCount];
    writeLogRatioCoordinates(variables, 0, phaseZero, referenceIndices[0]);
    writeLogRatioCoordinates(variables, coordinateCount, phaseOne, referenceIndices[1]);
    writeLogRatioCoordinates(variables, 2 * coordinateCount, bifurcating, referenceIndices[2]);
    writeLogRatioCoordinates(variables, incipientStart, incipient, referenceIndices[3]);
    double fractionZero = state.getBeta();
    double fractionOne = 1.0 - state.getBeta();
    if (bifurcatingPhase == retainedPhaseZero) {
      fractionZero -= bifurcatingFractionSeed;
    } else {
      fractionOne -= bifurcatingFractionSeed;
    }
    variables[fractionStart] = Math.log(fractionOne / fractionZero);
    variables[fractionStart + 1] = Math.log(bifurcatingFractionSeed / fractionZero);
    variables[incipientStart + coordinateCount] = Math.log(state.getTemperatureK());
    variables[incipientStart + coordinateCount + 1] = Math.log(state.getPressureBara());
    return new InitialState(variables, referenceIndices);
  }

  private Evaluation evaluate(double[] variables, int[] referenceIndices, int materialBalanceOmittedIndex,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode, double targetSeparation) {
    int componentCount = componentCount();
    int coordinateCount = componentCount - 1;
    int fractionStart = 3 * coordinateCount;
    int incipientStart = fractionStart + 2;
    int temperatureIndex = incipientStart + coordinateCount;
    int pressureIndex = temperatureIndex + 1;
    double temperatureK = Math.exp(variables[temperatureIndex]);
    double pressureBara = Math.exp(variables[pressureIndex]);
    double[] fractions = softmaxThree(variables[fractionStart], variables[fractionStart + 1]);
    double[] overall = overallComposition();
    double[][] compositions = new double[][] {
        compositionFromLogRatios(variables, 0, referenceIndices[0], componentCount),
        compositionFromLogRatios(variables, coordinateCount, referenceIndices[1], componentCount),
        compositionFromLogRatios(variables, 2 * coordinateCount, referenceIndices[2], componentCount) };
    double[] incipientComposition = compositionFromLogRatios(variables, incipientStart, referenceIndices[3],
        componentCount);
    SystemInterface working = createWorkingSystem(temperatureK, pressureBara, fractions, compositions,
        incipientComposition);
    double[] residual = new double[variables.length];
    double equilibriumMaximumResidual = 0.0;
    int phaseOneFugacityStart = 0;
    int phaseTwoFugacityStart = componentCount;
    int materialBalanceStart = 2 * componentCount;
    int incipientFugacityStart = materialBalanceStart + componentCount - 1;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double referenceLogFugacity = chemicalPotential(working, 0, componentIndex);
      residual[phaseOneFugacityStart + componentIndex] = chemicalPotential(working, 1, componentIndex)
          - referenceLogFugacity;
      residual[phaseTwoFugacityStart + componentIndex] = chemicalPotential(working, 2, componentIndex)
          - referenceLogFugacity;
      residual[incipientFugacityStart + componentIndex] = chemicalPotential(working, 3, componentIndex)
          - referenceLogFugacity;
      equilibriumMaximumResidual = Math.max(equilibriumMaximumResidual,
          Math.abs(residual[phaseOneFugacityStart + componentIndex]));
      equilibriumMaximumResidual = Math.max(equilibriumMaximumResidual,
          Math.abs(residual[phaseTwoFugacityStart + componentIndex]));
    }
    int materialResidualIndex = materialBalanceStart;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      if (componentIndex == materialBalanceOmittedIndex) {
        continue;
      }
      double reconstructed = fractions[0] * compositions[0][componentIndex]
          + fractions[1] * compositions[1][componentIndex] + fractions[2] * compositions[2][componentIndex];
      residual[materialResidualIndex] = reconstructed - overall[componentIndex];
      equilibriumMaximumResidual = Math.max(equilibriumMaximumResidual, Math.abs(residual[materialResidualIndex++]));
    }
    double incipientMaximumResidual = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      incipientMaximumResidual = Math.max(incipientMaximumResidual,
          Math.abs(residual[incipientFugacityStart + componentIndex]));
    }
    double[] duplicated = bifurcatingPhase == retainedPhaseZero ? compositions[0] : compositions[1];
    double separationResidual = modeSeparation(duplicated, compositions[2], mode) - targetSeparation;
    residual[variables.length - 1] = separationResidual;
    double maximumResidual = Math.max(equilibriumMaximumResidual,
        Math.max(incipientMaximumResidual, Math.abs(separationResidual)));
    if (!Double.isFinite(maximumResidual)) {
      throw new IllegalStateException("post-spinodal boundary residual is non-finite");
    }
    return new Evaluation(residual, maximumResidual, equilibriumMaximumResidual, incipientMaximumResidual,
        separationResidual, temperatureK, pressureBara, fractions, compositions, incipientComposition, duplicated,
        compositions[2]);
  }

  private SystemInterface createWorkingSystem(double temperatureK, double pressureBara, double[] fractions,
      double[][] compositions, double[] incipientComposition) {
    SystemInterface working = template.clone();
    working.setMultiPhaseCheck(false);
    working.setMaxNumberOfPhases(4);
    working.setNumberOfPhases(4);
    PhaseInterface phaseTemplate = template.getPhase(0);
    for (int phaseIndex = 0; phaseIndex < 4; phaseIndex++) {
      working.getPhases()[phaseIndex] = phaseTemplate.clone();
    }
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    CandidatePhase[] slots = new CandidatePhase[] { retainedPhaseZero, retainedPhaseOne, bifurcatingPhase,
        incipientPhase };
    for (int phaseIndex = 0; phaseIndex < 4; phaseIndex++) {
      working.setPhaseType(phaseIndex, toPhaseType(slots[phaseIndex]));
      working.setBeta(phaseIndex, phaseIndex < 3 ? fractions[phaseIndex] * (1.0 - 1.0e-14) : 1.0e-14);
      double[] composition = phaseIndex < 3 ? compositions[phaseIndex] : incipientComposition;
      for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
        working.getPhase(phaseIndex).getComponent(componentIndex).setx(composition[componentIndex]);
      }
      working.getPhase(phaseIndex).normalize();
      working.setPhaseType(phaseIndex, toPhaseType(slots[phaseIndex]));
      working.init(1, phaseIndex);
    }
    return working;
  }

  private boolean phaseIdentity(double[][] retainedCompositions, double[] incipientComposition) {
    int waterIndex = waterComponentIndex();
    if (waterIndex < 0) {
      return true;
    }
    CandidatePhase[] slots = new CandidatePhase[] { retainedPhaseZero, retainedPhaseOne, bifurcatingPhase };
    for (int phaseIndex = 0; phaseIndex < slots.length; phaseIndex++) {
      double water = retainedCompositions[phaseIndex][waterIndex];
      if (slots[phaseIndex] == CandidatePhase.OIL && water >= 0.5
          || slots[phaseIndex] == CandidatePhase.AQUEOUS && water < 0.5) {
        return false;
      }
    }
    double incipientWater = incipientComposition[waterIndex];
    if (incipientPhase == CandidatePhase.OIL) {
      return incipientWater < 0.5;
    }
    if (incipientPhase == CandidatePhase.AQUEOUS) {
      return incipientWater >= 0.5;
    }
    return true;
  }

  private void validateInputs(TwoToThreePhaseArcLengthCorrector.State state,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode, double targetSeparation, double bifurcatingFractionSeed) {
    int componentCount = componentCount();
    if (state == null) {
      throw new IllegalArgumentException("invalid three-to-four-phase boundary input");
    }
    double parentFraction = bifurcatingPhase == retainedPhaseZero ? state.getBeta() : 1.0 - state.getBeta();
    if (state.getRetainedPhaseZero() != retainedPhaseZero || state.getRetainedPhaseOne() != retainedPhaseOne
        || state.getIncipientPhase() != incipientPhase || mode == null
        || mode.getBifurcationEigenvector().length != componentCount - 1 || !Double.isFinite(targetSeparation)
        || Math.abs(targetSeparation) < 1.0e-8 || !positive(bifurcatingFractionSeed)
        || bifurcatingFractionSeed >= parentFraction) {
      throw new IllegalArgumentException("invalid three-to-four-phase boundary input");
    }
  }

  private int componentCount() {
    return template.getPhase(0).getNumberOfComponents();
  }

  private int waterComponentIndex() {
    for (int componentIndex = 0; componentIndex < componentCount(); componentIndex++) {
      if (template.getPhase(0).getComponent(componentIndex).getComponentName().equalsIgnoreCase("water")) {
        return componentIndex;
      }
    }
    return -1;
  }

  private double[] overallComposition() {
    double[] result = new double[componentCount()];
    for (int componentIndex = 0; componentIndex < result.length; componentIndex++) {
      result[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(), MINIMUM_COMPOSITION);
    }
    return normalized(result);
  }

  private double materialBalanceResidual(double[] fractions, double[][] compositions) {
    double[] overall = overallComposition();
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < overall.length; componentIndex++) {
      double reconstructed = 0.0;
      for (int phaseIndex = 0; phaseIndex < 3; phaseIndex++) {
        reconstructed += fractions[phaseIndex] * compositions[phaseIndex][componentIndex];
      }
      maximum = Math.max(maximum, Math.abs(reconstructed - overall[componentIndex]));
    }
    return maximum;
  }

  private static double[] perturbAlongMode(double[] base, IncipientPhaseStationarityJacobianAnalyzer.Result mode,
      double amplitude) {
    int referenceIndex = mode.getReferenceComponentIndex();
    double[] eigenvector = mode.getBifurcationEigenvector();
    double[] logarithms = new double[base.length];
    double reference = Math.max(base[referenceIndex], MINIMUM_COMPOSITION);
    int coordinateIndex = 0;
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < base.length; componentIndex++) {
      if (componentIndex != referenceIndex) {
        logarithms[componentIndex] = Math.log(Math.max(base[componentIndex], MINIMUM_COMPOSITION) / reference)
            + amplitude * eigenvector[coordinateIndex++];
      }
      maximum = Math.max(maximum, logarithms[componentIndex]);
    }
    double[] result = new double[base.length];
    for (int componentIndex = 0; componentIndex < result.length; componentIndex++) {
      result[componentIndex] = Math.exp(Math.max(-80.0, logarithms[componentIndex] - maximum));
    }
    return normalized(result);
  }

  private static double modeSeparation(double[] retained, double[] bifurcating,
      IncipientPhaseStationarityJacobianAnalyzer.Result mode) {
    int referenceIndex = mode.getReferenceComponentIndex();
    double retainedReference = Math.max(retained[referenceIndex], MINIMUM_COMPOSITION);
    double bifurcatingReference = Math.max(bifurcating[referenceIndex], MINIMUM_COMPOSITION);
    double[] eigenvector = mode.getBifurcationEigenvector();
    double separation = 0.0;
    int coordinateIndex = 0;
    for (int componentIndex = 0; componentIndex < retained.length; componentIndex++) {
      if (componentIndex == referenceIndex) {
        continue;
      }
      double difference = Math.log(Math.max(bifurcating[componentIndex], MINIMUM_COMPOSITION) / bifurcatingReference)
          - Math.log(Math.max(retained[componentIndex], MINIMUM_COMPOSITION) / retainedReference);
      separation += eigenvector[coordinateIndex++] * difference;
    }
    return separation;
  }

  private static double[] softmaxThree(double firstLogit, double secondLogit) {
    double maximum = Math.max(0.0, Math.max(firstLogit, secondLogit));
    double zero = Math.exp(-maximum);
    double one = Math.exp(firstLogit - maximum);
    double two = Math.exp(secondLogit - maximum);
    double total = zero + one + two;
    return new double[] { zero / total, one / total, two / total };
  }

  private static void writeLogRatioCoordinates(double[] variables, int start, double[] composition,
      int referenceIndex) {
    double reference = Math.max(composition[referenceIndex], MINIMUM_COMPOSITION);
    int coordinateIndex = 0;
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      if (componentIndex != referenceIndex) {
        variables[start + coordinateIndex++] = Math
            .log(Math.max(composition[componentIndex], MINIMUM_COMPOSITION) / reference);
      }
    }
  }

  private static double[] compositionFromLogRatios(double[] variables, int start, int referenceIndex,
      int componentCount) {
    double[] logarithms = new double[componentCount];
    int coordinateIndex = 0;
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      if (componentIndex != referenceIndex) {
        logarithms[componentIndex] = variables[start + coordinateIndex++];
      }
      maximum = Math.max(maximum, logarithms[componentIndex]);
    }
    double[] composition = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = Math.exp(Math.max(-80.0, logarithms[componentIndex] - maximum));
    }
    return normalized(composition);
  }

  private static int largestComponentIndex(double[] composition) {
    int result = 0;
    for (int componentIndex = 1; componentIndex < composition.length; componentIndex++) {
      if (composition[componentIndex] > composition[result]) {
        result = componentIndex;
      }
    }
    return result;
  }

  private static double chemicalPotential(SystemInterface system, int phaseIndex, int componentIndex) {
    return Math.log(Math.max(system.getPhase(phaseIndex).getComponent(componentIndex).getx(), MINIMUM_COMPOSITION))
        + system.getPhase(phaseIndex).getComponent(componentIndex).getLogFugacityCoefficient();
  }

  private static double logSumExp(double[] values, int start, int end) {
    double maximum = Double.NEGATIVE_INFINITY;
    for (int index = start; index < end; index++) {
      maximum = Math.max(maximum, values[index]);
    }
    double total = 0.0;
    for (int index = start; index < end; index++) {
      total += Math.exp(values[index] - maximum);
    }
    return maximum + Math.log(total);
  }

  private static double[] normalizedExponentials(double[] values, int start, int end, double logSum) {
    double[] result = new double[end - start];
    for (int index = start; index < end; index++) {
      result[index - start] = Math.exp(values[index] - logSum);
    }
    return result;
  }

  private static double[] normalized(double[] values) {
    double[] result = values.clone();
    double total = 0.0;
    for (double value : result) {
      total += value;
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalStateException("composition cannot be normalized");
    }
    for (int index = 0; index < result.length; index++) {
      result[index] /= total;
    }
    return result;
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
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

  private static double[][] copy(double[][] matrix) {
    double[][] result = new double[matrix.length][];
    for (int row = 0; row < matrix.length; row++) {
      result[row] = matrix[row].clone();
    }
    return result;
  }

  private static final class InitialState {
    private final double[] variables;
    private final int[] referenceIndices;

    private InitialState(double[] variables, int[] referenceIndices) {
      this.variables = variables;
      this.referenceIndices = referenceIndices;
    }
  }

  private static final class Evaluation {
    private final double[] residual;
    private final double maximumResidual;
    private final double equilibriumMaximumResidual;
    private final double incipientMaximumResidual;
    private final double separationResidual;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] phaseFractions;
    private final double[][] phaseCompositions;
    private final double[] incipientComposition;
    private final double[] duplicatedRetainedComposition;
    private final double[] bifurcatingComposition;

    private Evaluation(double[] residual, double maximumResidual, double equilibriumMaximumResidual,
        double incipientMaximumResidual, double separationResidual, double temperatureK, double pressureBara,
        double[] phaseFractions, double[][] phaseCompositions, double[] incipientComposition,
        double[] duplicatedRetainedComposition, double[] bifurcatingComposition) {
      this.residual = residual;
      this.maximumResidual = maximumResidual;
      this.equilibriumMaximumResidual = equilibriumMaximumResidual;
      this.incipientMaximumResidual = incipientMaximumResidual;
      this.separationResidual = separationResidual;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.phaseFractions = phaseFractions;
      this.phaseCompositions = phaseCompositions;
      this.incipientComposition = incipientComposition;
      this.duplicatedRetainedComposition = duplicatedRetainedComposition;
      this.bifurcatingComposition = bifurcatingComposition;
    }
  }

  /** Immutable corrected three-to-four-phase boundary point. */
  public static final class Result {
    private final CandidatePhase retainedPhaseZero;
    private final CandidatePhase retainedPhaseOne;
    private final CandidatePhase bifurcatingPhase;
    private final CandidatePhase incipientPhase;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] phaseFractions;
    private final double[][] phaseCompositions;
    private final double[] incipientComposition;
    private final double[] residual;
    private final double initialMaximumResidual;
    private final double maximumResidual;
    private final double equilibriumMaximumResidual;
    private final double incipientMaximumResidual;
    private final double separationResidual;
    private final double materialBalanceResidual;
    private final int iterations;
    private final int acceptedSteps;
    private final int rejectedSteps;
    private final double finalTrustRadius;
    private final double jacobianConditionNumber;
    private final boolean residualConverged;
    private final boolean distinct;
    private final boolean finiteFractions;
    private final boolean physicalIdentity;
    private final String failureMessage;

    private Result(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne, CandidatePhase bifurcatingPhase,
        CandidatePhase incipientPhase, double temperatureK, double pressureBara, double[] phaseFractions,
        double[][] phaseCompositions, double[] incipientComposition, double[] residual, double initialMaximumResidual,
        double maximumResidual, double equilibriumMaximumResidual, double incipientMaximumResidual,
        double separationResidual, double materialBalanceResidual, int iterations, int acceptedSteps, int rejectedSteps,
        double finalTrustRadius, double jacobianConditionNumber, boolean residualConverged, boolean distinct,
        boolean finiteFractions, boolean physicalIdentity, String failureMessage) {
      this.retainedPhaseZero = retainedPhaseZero;
      this.retainedPhaseOne = retainedPhaseOne;
      this.bifurcatingPhase = bifurcatingPhase;
      this.incipientPhase = incipientPhase;
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.phaseFractions = phaseFractions.clone();
      this.phaseCompositions = copy(phaseCompositions);
      this.incipientComposition = incipientComposition.clone();
      this.residual = residual.clone();
      this.initialMaximumResidual = initialMaximumResidual;
      this.maximumResidual = maximumResidual;
      this.equilibriumMaximumResidual = equilibriumMaximumResidual;
      this.incipientMaximumResidual = incipientMaximumResidual;
      this.separationResidual = separationResidual;
      this.materialBalanceResidual = materialBalanceResidual;
      this.iterations = iterations;
      this.acceptedSteps = acceptedSteps;
      this.rejectedSteps = rejectedSteps;
      this.finalTrustRadius = finalTrustRadius;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.residualConverged = residualConverged;
      this.distinct = distinct;
      this.finiteFractions = finiteFractions;
      this.physicalIdentity = physicalIdentity;
      this.failureMessage = failureMessage;
    }

    private static Result failure(CandidatePhase retainedPhaseZero, CandidatePhase retainedPhaseOne,
        CandidatePhase bifurcatingPhase, CandidatePhase incipientPhase, int componentCount, String failureMessage) {
      return new Result(retainedPhaseZero, retainedPhaseOne, bifurcatingPhase, incipientPhase, Double.NaN, Double.NaN,
          new double[3], new double[3][componentCount], new double[componentCount], new double[4 * componentCount],
          Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0, 0, 0, Double.NaN, Double.NaN,
          false, false, false, false, failureMessage);
    }

    public boolean isConverged() {
      return residualConverged && distinct && finiteFractions && physicalIdentity && failureMessage == null;
    }

    public CandidatePhase getRetainedPhaseZero() {
      return retainedPhaseZero;
    }

    public CandidatePhase getRetainedPhaseOne() {
      return retainedPhaseOne;
    }

    public CandidatePhase getBifurcatingPhase() {
      return bifurcatingPhase;
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

    public double[] getPhaseFractions() {
      return phaseFractions.clone();
    }

    public double getPhaseFraction(int phaseIndex) {
      return phaseFractions[phaseIndex];
    }

    public double[] getPhaseComposition(int phaseIndex) {
      return phaseCompositions[phaseIndex].clone();
    }

    public double[] getIncipientComposition() {
      return incipientComposition.clone();
    }

    public double[] getResidual() {
      return residual.clone();
    }

    public double getInitialMaximumResidual() {
      return initialMaximumResidual;
    }

    public double getMaximumResidual() {
      return maximumResidual;
    }

    public double getEquilibriumMaximumResidual() {
      return equilibriumMaximumResidual;
    }

    public double getIncipientMaximumResidual() {
      return incipientMaximumResidual;
    }

    public double getSeparationResidual() {
      return separationResidual;
    }

    public double getMaterialBalanceResidual() {
      return materialBalanceResidual;
    }

    public int getIterations() {
      return iterations;
    }

    public int getAcceptedSteps() {
      return acceptedSteps;
    }

    public int getRejectedSteps() {
      return rejectedSteps;
    }

    public double getFinalTrustRadius() {
      return finalTrustRadius;
    }

    public double getJacobianConditionNumber() {
      return jacobianConditionNumber;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
