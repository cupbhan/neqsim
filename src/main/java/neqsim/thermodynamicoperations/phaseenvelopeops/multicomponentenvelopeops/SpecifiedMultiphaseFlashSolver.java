package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import neqsim.thermo.phase.LiquidPhaseClassification;

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Trust-region TP flash for two to four caller-specified phase slots.
 *
 * <p>
 * Phase families may repeat, allowing vapor/liquid/liquid/aqueous states to be represented as
 * {@code GAS/OIL/OIL/AQUEOUS}. Component balances are eliminated with a generalized Rachford-Rice denominator. The
 * remaining equations enforce fugacity equality and composition normalization relative to phase slot zero.
 *
 *
 * <p>
 * This solver does not decide whether the requested topology is globally stable. Callers must apply tangent-plane,
 * phase-identity, and topology gates before admitting a result to an engineering phase envelope.
 *
 *
 * @author NeqSim contributors
 * @version 1.0
 */
public final class SpecifiedMultiphaseFlashSolver {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;
  private static final double MINIMUM_PHASE_FRACTION = 1.0e-14;

  private final SystemInterface template;
  private final Object resultOwner = new Object();
  private final CandidatePhase[] phaseSlots;
  private int maximumIterations = 160;
  private double residualTolerance = 1.0e-9;
  private double finiteDifferenceStep = 2.0e-5;
  private double initialTrustRadius = 1.0;
  private double maximumTrustRadius = 100.0;

  /**
   * Creates a specified multiphase flash for two to four phase slots.
   *
   * @param template configured fluid copied privately by this solver
   * @param phaseSlots ordered gas, oil or aqueous phase families
   */
  public SpecifiedMultiphaseFlashSolver(SystemInterface template, CandidatePhase... phaseSlots) {
    if (template == null || template.getNumberOfComponents() < 2 || phaseSlots == null || phaseSlots.length < 2
        || phaseSlots.length > 4) {
      throw new IllegalArgumentException("template and two to four phase slots are required");
    }
    for (CandidatePhase phaseSlot : phaseSlots) {
      if (phaseSlot == null) {
        throw new IllegalArgumentException("every phase slot must specify a phase family");
      }
    }
    this.template = template.clone();
    this.template.init(0);
    this.phaseSlots = phaseSlots.clone();
  }

  /**
   * Sets trust-region and numerical-Jacobian controls.
   *
   * @param maximumIterations positive nonlinear iteration limit
   * @param residualTolerance positive maximum equilibrium residual
   * @param finiteDifferenceStep positive relative Jacobian perturbation
   * @param initialTrustRadius positive initial trust-region radius
   * @param maximumTrustRadius maximum trust-region radius
   * @return computed set numerical controls result
   */
  public SpecifiedMultiphaseFlashSolver setNumericalControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep, double initialTrustRadius, double maximumTrustRadius) {
    if (maximumIterations < 1 || !positive(residualTolerance) || !positive(finiteDifferenceStep)
        || !positive(initialTrustRadius) || !positive(maximumTrustRadius) || initialTrustRadius > maximumTrustRadius) {
      throw new IllegalArgumentException("invalid specified multiphase flash numerical controls");
    }
    this.maximumIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    this.initialTrustRadius = initialTrustRadius;
    this.maximumTrustRadius = maximumTrustRadius;
    return this;
  }

  /**
   * Solves one fixed-temperature, fixed-pressure multiphase equilibrium state.
   *
   * @param temperatureK temperature in kelvin
   * @param pressureBara pressure in bara
   * @param phaseFractions positive phase-fraction seeds
   * @param phaseCompositions one composition seed per phase slot
   * @return immutable equilibrium result and numerical diagnostics
   */
  public Result solve(final double temperatureK, final double pressureBara, double[] phaseFractions,
      double[][] phaseCompositions) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    validateInputs(temperatureK, pressureBara, phaseFractions, phaseCompositions, componentCount);
    double[] normalizedFractions = normalizedFractions(phaseFractions);
    double[][] normalizedCompositions = copy(phaseCompositions);
    for (int phaseIndex = 0; phaseIndex < normalizedCompositions.length; phaseIndex++) {
      normalizedCompositions[phaseIndex] = normalized(normalizedCompositions[phaseIndex]);
    }
    double[] variables = variables(normalizedFractions, normalizedCompositions, componentCount);
    double[] lowerBounds = new double[variables.length];
    double[] upperBounds = new double[variables.length];
    double[] scales = new double[variables.length];
    int logKCount = (phaseSlots.length - 1) * componentCount;
    for (int index = 0; index < variables.length; index++) {
      lowerBounds[index] = index < logKCount ? -80.0 : -32.0;
      upperBounds[index] = index < logKCount ? 80.0 : 32.0;
      scales[index] = 1.0;
    }
    // A private workspace per solve avoids repeated HV/database reconstruction in every Jacobian sample.
    final SystemInterface working = createWorkingSystem(temperatureK, pressureBara, normalizedFractions,
        normalizedCompositions);
    ScaledTrustRegionLeastSquaresSolver trustRegion = new ScaledTrustRegionLeastSquaresSolver(
        trial -> evaluate(working, temperatureK, pressureBara, trial).residual)
        .setNumericalControls(maximumIterations, residualTolerance, finiteDifferenceStep, initialTrustRadius,
            maximumTrustRadius)
        .setVariableScales(scales).setBounds(lowerBounds, upperBounds);
    ScaledTrustRegionLeastSquaresSolver.Result numerical = trustRegion.solve(variables);
    Evaluation evaluation;
    try {
      evaluation = evaluate(working, temperatureK, pressureBara, numerical.getVariables());
    } catch (RuntimeException error) {
      return Result.failure(resultOwner, phaseSlots, temperatureK, pressureBara, componentCount,
          "final residual evaluation failed: " + error.getMessage());
    }
    for (int phaseIndex = 0; phaseIndex < phaseSlots.length; phaseIndex++) {
      evaluation.phaseCompositions[phaseIndex] = normalized(evaluation.phaseCompositions[phaseIndex]);
    }
    boolean physicalIdentity = hasPhysicalPhaseIdentity(evaluation.phaseCompositions);
    String identityDiagnostic = physicalIdentity ? null : phaseIdentityDiagnostic(evaluation.phaseCompositions);
    return new Result(resultOwner, phaseSlots, temperatureK, pressureBara, evaluation.phaseFractions,
        evaluation.phaseCompositions, numerical.getInitialMaximumResidual(), evaluation.maximumResidual,
        materialBalanceResidual(evaluation.phaseFractions, evaluation.phaseCompositions), numerical.getIterations(),
        numerical.getAcceptedSteps(), numerical.getRejectedSteps(), numerical.getFinalTrustRadius(),
        numerical.getJacobianConditionNumber(), numerical.isConverged(), physicalIdentity, identityDiagnostic,
        numerical.getFailureMessage());
  }

  /**
   * Reconstructs a private thermodynamic system from a converged result.
   *
   * @param result result produced by this solver instance
   * @return computed to thermodynamic system result
   */
  public SystemInterface toThermodynamicSystem(Result result) {
    if (result == null || result.owner != resultOwner || !result.isConverged()
        || result.getPhaseCount() != phaseSlots.length) {
      throw new IllegalArgumentException("a converged result with matching phase slots is required");
    }
    CandidatePhase[] resultSlots = result.getPhaseSlots();
    for (int phaseIndex = 0; phaseIndex < phaseSlots.length; phaseIndex++) {
      if (resultSlots[phaseIndex] != phaseSlots[phaseIndex]) {
        throw new IllegalArgumentException("result phase slots do not match this solver");
      }
    }
    return createWorkingSystem(result.getTemperatureK(), result.getPressureBara(), result.phaseFractions,
        result.phaseCompositions);
  }

  /**
   * Replays a converged numerical result and checks finite phases and tangent-plane stability.
   *
   * @param result converged result produced by this solver instance
   * @return scoped physical acceptance diagnostics
   * @throws IllegalArgumentException when the result did not converge or belongs to another solver
   */
  public SpecifiedPhaseEquilibriumValidator.Result validateEquilibrium(Result result) {
    return SpecifiedPhaseEquilibriumValidator.validate(template, toThermodynamicSystem(result), phaseSlots);
  }

  /**
   * Exports an equilibrium only after the independent replay and stability checks pass.
   *
   * @param result converged result produced by this solver instance
   * @return private accepted thermodynamic system
   * @throws IllegalArgumentException when the result cannot be reconstructed
   * @throws IllegalStateException when any physical acceptance check fails
   */
  public SystemInterface toValidatedThermodynamicSystem(Result result) {
    SpecifiedPhaseEquilibriumValidator.Result validation = validateEquilibrium(result);
    if (!validation.isAccepted()) {
      throw new IllegalStateException("specified equilibrium rejected: " + validation.getViolations());
    }
    return toThermodynamicSystem(result);
  }

  /**
   * Constructs logarithmic solver coordinates from the phase seeds.
   *
   * @param fractions fractions
   * @param compositions compositions
   * @param componentCount number of components
   * @return computed variables result
   */
  private double[] variables(double[] fractions, double[][] compositions, int componentCount) {
    int nonReferencePhaseCount = phaseSlots.length - 1;
    double[] variables = new double[nonReferencePhaseCount * (componentCount + 1)];
    for (int phaseIndex = 1; phaseIndex < phaseSlots.length; phaseIndex++) {
      int blockStart = (phaseIndex - 1) * componentCount;
      for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
        double reference = Math.max(compositions[0][componentIndex], MINIMUM_COMPOSITION);
        variables[blockStart + componentIndex] = Math
            .log(Math.max(compositions[phaseIndex][componentIndex], MINIMUM_COMPOSITION) / reference);
      }
      variables[nonReferencePhaseCount * componentCount + phaseIndex - 1] = Math
          .log(fractions[phaseIndex] / fractions[0]);
    }
    return variables;
  }

  /**
   * Numerically safe {@code log(sum(exp(values)))}, used so the inventory never has to be formed directly.
   *
   * @param values values
   * @return computed log sum of exponentials result
   */
  private static double logSumOfExponentials(double[] values) {
    double maximum = Double.NEGATIVE_INFINITY;
    for (int index = 0; index < values.length; index++) {
      maximum = Math.max(maximum, values[index]);
    }
    if (!Double.isFinite(maximum)) {
      return maximum;
    }
    double total = 0.0;
    for (int index = 0; index < values.length; index++) {
      total += Math.exp(values[index] - maximum);
    }
    return maximum + Math.log(total);
  }

  /**
   * Evaluates the nonlinear equilibrium equations at the supplied solver coordinates.
   *
   * @param temperatureK temperature in kelvin
   * @param pressureBara absolute pressure in bara
   * @param variables variables
   * @return computed evaluate result
   *
   * @param working working
   */
  private Evaluation evaluate(SystemInterface working, double temperatureK, double pressureBara, double[] variables) {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    int nonReferencePhaseCount = phaseSlots.length - 1;
    int logKCount = nonReferencePhaseCount * componentCount;
    if (variables.length != logKCount + nonReferencePhaseCount) {
      throw new IllegalStateException("specified multiphase variable count changed");
    }
    double[] fractions = softmaxFractions(variables, logKCount, phaseSlots.length);
    double[] overall = overallComposition();
    double[][] compositions = new double[phaseSlots.length][componentCount];
    double[] sums = new double[phaseSlots.length];
    // Evaluated in log space so the composition floor never enters the generalized Rachford-Rice residual, and the
    // equilibrium ratio is no longer clamped at exp(+-80). A heavy pseudo-component against water reaches K ~ 1e99,
    // i.e. ln K ~ 228, so the old bound truncated the exponent by 148; the true reference-phase fraction near
    // 1e-102 was then floored away as well. Either effect pins the residual to an artefact instead of the
    // equations. The floor belongs only at the EOS input, and the log bound now matches the sibling solvers.
    double[][] logCompositions = new double[phaseSlots.length][componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double[] logKValues = new double[phaseSlots.length];
      double[] denominatorTerms = new double[phaseSlots.length];
      denominatorTerms[0] = Math.log(fractions[0]);
      for (int phaseIndex = 1; phaseIndex < phaseSlots.length; phaseIndex++) {
        double logK = variables[(phaseIndex - 1) * componentCount + componentIndex];
        logKValues[phaseIndex] = Math.max(-700.0, Math.min(700.0, logK));
        denominatorTerms[phaseIndex] = Math.log(fractions[phaseIndex]) + logKValues[phaseIndex];
      }
      double logDenominator = logSumOfExponentials(denominatorTerms);
      if (!Double.isFinite(logDenominator)) {
        throw new IllegalStateException("invalid generalized Rachford-Rice denominator");
      }
      logCompositions[0][componentIndex] = Math.log(overall[componentIndex]) - logDenominator;
      compositions[0][componentIndex] = Math.exp(logCompositions[0][componentIndex]);
      for (int phaseIndex = 1; phaseIndex < phaseSlots.length; phaseIndex++) {
        logCompositions[phaseIndex][componentIndex] = logKValues[phaseIndex] + logCompositions[0][componentIndex];
        compositions[phaseIndex][componentIndex] = Math.exp(logCompositions[phaseIndex][componentIndex]);
      }
    }
    for (int phaseIndex = 0; phaseIndex < phaseSlots.length; phaseIndex++) {
      sums[phaseIndex] = Math.exp(logSumOfExponentials(logCompositions[phaseIndex]));
      if (!Double.isFinite(sums[phaseIndex]) || !(sums[phaseIndex] > 0.0)) {
        throw new IllegalStateException("specified multiphase inventory left the physical domain");
      }
    }
    updateWorkingSystem(working, temperatureK, pressureBara, fractions, compositions);
    double[] residual = new double[variables.length];
    double maximumResidual = 0.0;
    for (int phaseIndex = 1; phaseIndex < phaseSlots.length; phaseIndex++) {
      int blockStart = (phaseIndex - 1) * componentCount;
      for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
        residual[blockStart + componentIndex] = variables[blockStart + componentIndex]
            + working.getPhase(phaseIndex).getComponent(componentIndex).getLogFugacityCoefficient()
            - working.getPhase(0).getComponent(componentIndex).getLogFugacityCoefficient();
        maximumResidual = Math.max(maximumResidual, Math.abs(residual[blockStart + componentIndex]));
      }
      residual[logKCount + phaseIndex - 1] = sums[phaseIndex] - sums[0];
      maximumResidual = Math.max(maximumResidual, Math.abs(residual[logKCount + phaseIndex - 1]));
    }
    if (!Double.isFinite(maximumResidual)) {
      throw new IllegalStateException("specified multiphase residual is non-finite");
    }
    return new Evaluation(residual, maximumResidual, fractions, compositions);
  }

  /**
   * Builds a private EOS state from the requested phase compositions and fractions.
   *
   * @param temperatureK temperature in kelvin
   * @param pressureBara absolute pressure in bara
   * @param fractions fractions
   * @param compositions compositions
   * @return computed create working system result
   */
  private SystemInterface createWorkingSystem(double temperatureK, double pressureBara, double[] fractions,
      double[][] compositions) {
    SystemInterface working = template.clone();
    working.setMultiPhaseCheck(false);
    working.setMaxNumberOfPhases(phaseSlots.length);
    working.setNumberOfPhases(phaseSlots.length);
    PhaseInterface phaseTemplate = template.getPhase(0);
    for (int phaseIndex = 0; phaseIndex < phaseSlots.length; phaseIndex++) {
      // SystemThermo#setPhase reads the old slot temperature before replacement. The fourth
      // slot of ordinary cubic systems is null, so initialize all requested private slots
      // directly before setting the system-wide temperature and pressure.
      working.getPhases()[phaseIndex] = phaseTemplate.clone();
    }
    updateWorkingSystem(working, temperatureK, pressureBara, fractions, compositions);
    return working;
  }

  /**
   * Replaces every trial-dependent value in a solve-local workspace before EOS evaluation.
   *
   * @param working private workspace for this solve
   * @param temperatureK temperature in kelvin
   * @param pressureBara absolute pressure in bara
   * @param fractions phase mole fractions
   * @param compositions phase-by-component mole fractions
   */
  private void updateWorkingSystem(SystemInterface working, double temperatureK, double pressureBara,
      double[] fractions, double[][] compositions) {
    working.setTemperature(temperatureK);
    working.setPressure(pressureBara);
    for (int phaseIndex = 0; phaseIndex < phaseSlots.length; phaseIndex++) {
      working.setPhaseType(phaseIndex, toPhaseType(phaseSlots[phaseIndex]));
      working.setBeta(phaseIndex, fractions[phaseIndex]);
      for (int componentIndex = 0; componentIndex < compositions[phaseIndex].length; componentIndex++) {
        // The EOS requires a strictly positive mole fraction; this is the only place the floor belongs.
        working.getPhase(phaseIndex).getComponent(componentIndex)
            .setx(Math.max(compositions[phaseIndex][componentIndex], MINIMUM_COMPOSITION));
      }
      working.getPhase(phaseIndex).normalize();
      working.setPhaseType(phaseIndex, toPhaseType(phaseSlots[phaseIndex]));
      working.init(1, phaseIndex);
    }
  }

  /**
   * Checks the oil/water composition convention independently of convergence.
   *
   * @param compositions compositions
   * @return true when the documented condition holds
   */
  private boolean hasPhysicalPhaseIdentity(double[][] compositions) {
    for (int phaseIndex = 0; phaseIndex < phaseSlots.length; phaseIndex++) {
      if (phaseSlots[phaseIndex] != CandidatePhase.GAS && LiquidPhaseClassification.classify(template.getPhase(0),
          normalized(compositions[phaseIndex])) != toPhaseType(phaseSlots[phaseIndex])) {
        return false;
      }
    }
    return true;
  }

  /**
   * Describes a mismatch between requested and computed phase families.
   *
   * @param compositions compositions
   * @return computed phase identity diagnostic result
   */
  private String phaseIdentityDiagnostic(double[][] compositions) {
    int waterIndex = waterComponentIndex();
    if (waterIndex < 0) {
      return "aqueous phase requested for a fluid without water";
    }
    StringBuilder diagnostic = new StringBuilder("phase-family mismatch:");
    for (int phaseIndex = 0; phaseIndex < phaseSlots.length; phaseIndex++) {
      diagnostic.append(' ').append(phaseIndex).append('=').append(phaseSlots[phaseIndex]).append("(xWater=")
          .append(normalized(compositions[phaseIndex])[waterIndex]).append(')');
    }
    return diagnostic.toString();
  }

  /**
   * Finds water in the private template component order.
   *
   * @return computed water component index result
   */
  private int waterComponentIndex() {
    for (int componentIndex = 0; componentIndex < template.getPhase(0).getNumberOfComponents(); componentIndex++) {
      if (template.getPhase(0).getComponent(componentIndex).getComponentName().equalsIgnoreCase("water")) {
        return componentIndex;
      }
    }
    return -1;
  }

  /**
   * Returns the normalized overall component inventory of the private template.
   *
   * @return computed overall composition result
   */
  private double[] overallComposition() {
    int componentCount = template.getPhase(0).getNumberOfComponents();
    double[] composition = new double[componentCount];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = Math.max(template.getPhase(0).getComponent(componentIndex).getz(),
          MINIMUM_COMPOSITION);
      total += composition[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }

  /**
   * Measures the largest deviation from the original overall component inventory.
   *
   * @param fractions fractions
   * @param compositions compositions
   * @return computed material balance residual result
   */
  private double materialBalanceResidual(double[] fractions, double[][] compositions) {
    double[] overall = overallComposition();
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < overall.length; componentIndex++) {
      double reconstructed = 0.0;
      for (int phaseIndex = 0; phaseIndex < phaseSlots.length; phaseIndex++) {
        reconstructed += fractions[phaseIndex] * compositions[phaseIndex][componentIndex];
      }
      maximum = Math.max(maximum, Math.abs(reconstructed - overall[componentIndex]));
    }
    return maximum;
  }

  /**
   * Transforms log-ratio coordinates to positive normalized phase fractions.
   *
   * @param variables variables
   * @param logKCount log kcount
   * @param phaseCount phase count
   * @return computed softmax fractions result
   */
  private static double[] softmaxFractions(double[] variables, int logKCount, int phaseCount) {
    double maximum = 0.0;
    for (int phaseIndex = 1; phaseIndex < phaseCount; phaseIndex++) {
      maximum = Math.max(maximum, variables[logKCount + phaseIndex - 1]);
    }
    double[] fractions = new double[phaseCount];
    fractions[0] = Math.exp(-maximum);
    double total = fractions[0];
    for (int phaseIndex = 1; phaseIndex < phaseCount; phaseIndex++) {
      fractions[phaseIndex] = Math.exp(variables[logKCount + phaseIndex - 1] - maximum);
      total += fractions[phaseIndex];
    }
    for (int phaseIndex = 0; phaseIndex < phaseCount; phaseIndex++) {
      fractions[phaseIndex] = Math.max(fractions[phaseIndex] / total, MINIMUM_PHASE_FRACTION);
    }
    return normalizedFractions(fractions);
  }

  /**
   * Returns normalized, finite positive phase fractions.
   *
   * @param phaseFractions phase fractions
   * @return computed normalized fractions result
   */
  private static double[] normalizedFractions(double[] phaseFractions) {
    double[] result = phaseFractions.clone();
    double total = 0.0;
    for (double fraction : result) {
      if (!positive(fraction)) {
        throw new IllegalArgumentException("all specified phase fractions must be positive");
      }
      total += fraction;
    }
    if (!positive(total)) {
      throw new IllegalArgumentException("specified phase-fraction sum must be finite and positive");
    }
    for (int phaseIndex = 0; phaseIndex < result.length; phaseIndex++) {
      result[phaseIndex] /= total;
    }
    return result;
  }

  /**
   * Returns a normalized defensive copy of a composition.
   *
   * @param composition composition
   * @return computed normalized result
   */
  private static double[] normalized(double[] composition) {
    double[] result = composition.clone();
    double total = 0.0;
    for (double value : result) {
      total += value;
    }
    for (int index = 0; index < result.length; index++) {
      result[index] /= total;
    }
    return result;
  }

  /**
   * Checks dimensions and physical domains before starting a numerical solve.
   *
   * @param temperatureK temperature in kelvin
   * @param pressureBara absolute pressure in bara
   * @param fractions fractions
   * @param compositions compositions
   * @param componentCount number of components
   */
  private void validateInputs(double temperatureK, double pressureBara, double[] fractions, double[][] compositions,
      int componentCount) {
    if (!Double.isFinite(temperatureK) || temperatureK < 50.0 || !positive(pressureBara) || fractions == null
        || fractions.length != phaseSlots.length || compositions == null || compositions.length != phaseSlots.length) {
      throw new IllegalArgumentException("invalid specified multiphase flash input");
    }
    normalizedFractions(fractions);
    for (double[] composition : compositions) {
      if (composition == null || composition.length != componentCount) {
        throw new IllegalArgumentException("one composition per phase slot is required");
      }
      validateComposition(composition);
    }
  }

  /**
   * Checks finite nonnegative component values and a valid normalization sum.
   *
   * @param composition composition
   */
  private static void validateComposition(double[] composition) {
    double total = 0.0;
    for (double value : composition) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException("specified phase composition contains an invalid value");
      }
      total += value;
    }
    if (!positive(total)) {
      throw new IllegalArgumentException("specified phase composition cannot be normalized");
    }
  }

  /**
   * Maps the requested physical family to its EOS phase type.
   *
   * @param phase phase
   * @return computed to phase type result
   */
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

  /**
   * Checks that a scalar is finite and strictly positive.
   *
   * @param value value
   * @return true when the documented condition holds
   */
  private static boolean positive(double value) {
    return Double.isFinite(value) && value > 0.0;
  }

  /**
   * Copies every row of a matrix defensively.
   *
   * @param matrix matrix
   * @return computed copy result
   */
  private static double[][] copy(double[][] matrix) {
    double[][] result = new double[matrix.length][];
    for (int row = 0; row < matrix.length; row++) {
      result[row] = matrix[row].clone();
    }
    return result;
  }

  private static final class Evaluation {
    private final double[] residual;
    private final double maximumResidual;
    private final double[] phaseFractions;
    private final double[][] phaseCompositions;

    /**
     * Creates an immutable evaluation.
     *
     * @param residual residual
     * @param maximumResidual maximum residual
     * @param phaseFractions phase fractions
     * @param phaseCompositions phase compositions
     */
    private Evaluation(double[] residual, double maximumResidual, double[] phaseFractions,
        double[][] phaseCompositions) {
      this.residual = residual;
      this.maximumResidual = maximumResidual;
      this.phaseFractions = phaseFractions;
      this.phaseCompositions = phaseCompositions;
    }
  }

  /**
   * Immutable specified multiphase result and trust-region diagnostics.
   *
   * @author NeqSim contributors
   * @version 1.0
   */
  public static final class Result {
    private final Object owner;
    private final CandidatePhase[] phaseSlots;
    private final double temperatureK;
    private final double pressureBara;
    private final double[] phaseFractions;
    private final double[][] phaseCompositions;
    private final double initialMaximumResidual;
    private final double maximumResidual;
    private final double materialBalanceResidual;
    private final int iterations;
    private final int acceptedSteps;
    private final int rejectedSteps;
    private final double finalTrustRadius;
    private final double jacobianConditionNumber;
    private final boolean converged;
    private final boolean physicalPhaseIdentity;
    private final String phaseIdentityDiagnostic;
    private final String failureMessage;

    /**
     * Creates an immutable result.
     *
     * @param owner private identity of the solver that produced this result
     * @param phaseSlots ordered gas, oil or aqueous phase families
     * @param temperatureK temperature in kelvin
     * @param pressureBara absolute pressure in bara
     * @param phaseFractions phase fractions
     * @param phaseCompositions phase compositions
     * @param initialMaximumResidual initial maximum residual
     * @param maximumResidual maximum residual
     * @param materialBalanceResidual material balance residual
     * @param iterations iterations
     * @param acceptedSteps accepted steps
     * @param rejectedSteps rejected steps
     * @param finalTrustRadius final trust radius
     * @param jacobianConditionNumber jacobian condition number
     * @param converged converged
     * @param physicalPhaseIdentity physical phase identity
     * @param phaseIdentityDiagnostic phase identity diagnostic
     * @param failureMessage failure diagnostic, or null for success
     */
    private Result(Object owner, CandidatePhase[] phaseSlots, double temperatureK, double pressureBara,
        double[] phaseFractions, double[][] phaseCompositions, double initialMaximumResidual, double maximumResidual,
        double materialBalanceResidual, int iterations, int acceptedSteps, int rejectedSteps, double finalTrustRadius,
        double jacobianConditionNumber, boolean converged, boolean physicalPhaseIdentity,
        String phaseIdentityDiagnostic, String failureMessage) {
      this.owner = owner;
      this.phaseSlots = phaseSlots.clone();
      this.temperatureK = temperatureK;
      this.pressureBara = pressureBara;
      this.phaseFractions = phaseFractions.clone();
      this.phaseCompositions = copy(phaseCompositions);
      this.initialMaximumResidual = initialMaximumResidual;
      this.maximumResidual = maximumResidual;
      this.materialBalanceResidual = materialBalanceResidual;
      this.iterations = iterations;
      this.acceptedSteps = acceptedSteps;
      this.rejectedSteps = rejectedSteps;
      this.finalTrustRadius = finalTrustRadius;
      this.jacobianConditionNumber = jacobianConditionNumber;
      this.converged = converged;
      this.physicalPhaseIdentity = physicalPhaseIdentity;
      this.phaseIdentityDiagnostic = phaseIdentityDiagnostic;
      this.failureMessage = failureMessage;
    }

    /**
     * Constructs an immutable unsuccessful result with no accepted equilibrium.
     *
     * @param owner private identity of the solver that produced this result
     * @param phaseSlots ordered gas, oil or aqueous phase families
     * @param temperatureK temperature in kelvin
     * @param pressureBara absolute pressure in bara
     * @param componentCount number of components
     * @param failureMessage failure diagnostic, or null for success
     * @return computed failure result
     */
    private static Result failure(Object owner, CandidatePhase[] phaseSlots, double temperatureK, double pressureBara,
        int componentCount, String failureMessage) {
      return new Result(owner, phaseSlots, temperatureK, pressureBara, new double[phaseSlots.length],
          new double[phaseSlots.length][componentCount], Double.NaN, Double.NaN, Double.NaN, 0, 0, 0, Double.NaN,
          Double.NaN, false, false, null, failureMessage);
    }

    /**
     * Reports numerical convergence only; use validateEquilibrium before accepting a physical state.
     *
     * @return true when the documented condition holds
     */
    public boolean isConverged() {
      return converged && failureMessage == null;
    }

    /**
     * Returns the phase slots.
     *
     * @return phase slots
     */
    public CandidatePhase[] getPhaseSlots() {
      return phaseSlots.clone();
    }

    /**
     * Returns the phase count.
     *
     * @return phase count
     */
    public int getPhaseCount() {
      return phaseSlots.length;
    }

    /**
     * Returns the temperature k.
     *
     * @return temperature k
     */
    public double getTemperatureK() {
      return temperatureK;
    }

    /**
     * Returns the pressure bara.
     *
     * @return pressure bara
     */
    public double getPressureBara() {
      return pressureBara;
    }

    /**
     * Returns the phase fractions.
     *
     * @return phase fractions
     */
    public double[] getPhaseFractions() {
      return phaseFractions.clone();
    }

    /**
     * Returns the phase fraction.
     *
     * @param phaseIndex zero-based phase index
     * @return phase fraction
     */
    public double getPhaseFraction(int phaseIndex) {
      return phaseFractions[phaseIndex];
    }

    /**
     * Returns the phase composition.
     *
     * @param phaseIndex zero-based phase index
     * @return phase composition
     */
    public double[] getPhaseComposition(int phaseIndex) {
      return phaseCompositions[phaseIndex].clone();
    }

    /**
     * Returns the phase composition distance.
     *
     * @param firstPhaseIndex zero-based first phase index
     * @param secondPhaseIndex zero-based second phase index
     * @return phase composition distance
     */
    public double getPhaseCompositionDistance(int firstPhaseIndex, int secondPhaseIndex) {
      double distance = 0.0;
      for (int componentIndex = 0; componentIndex < phaseCompositions[firstPhaseIndex].length; componentIndex++) {
        distance += Math.abs(
            phaseCompositions[firstPhaseIndex][componentIndex] - phaseCompositions[secondPhaseIndex][componentIndex]);
      }
      return distance;
    }

    /**
     * Returns the initial maximum residual.
     *
     * @return initial maximum residual
     */
    public double getInitialMaximumResidual() {
      return initialMaximumResidual;
    }

    /**
     * Returns the maximum residual.
     *
     * @return maximum residual
     */
    public double getMaximumResidual() {
      return maximumResidual;
    }

    /**
     * Returns the material balance residual.
     *
     * @return material balance residual
     */
    public double getMaterialBalanceResidual() {
      return materialBalanceResidual;
    }

    /**
     * Returns the iterations.
     *
     * @return iterations
     */
    public int getIterations() {
      return iterations;
    }

    /**
     * Returns the accepted steps.
     *
     * @return accepted steps
     */
    public int getAcceptedSteps() {
      return acceptedSteps;
    }

    /**
     * Returns the rejected steps.
     *
     * @return rejected steps
     */
    public int getRejectedSteps() {
      return rejectedSteps;
    }

    /**
     * Returns the final trust radius.
     *
     * @return final trust radius
     */
    public double getFinalTrustRadius() {
      return finalTrustRadius;
    }

    /**
     * Returns the jacobian condition number.
     *
     * @return jacobian condition number
     */
    public double getJacobianConditionNumber() {
      return jacobianConditionNumber;
    }

    /**
     * Checks the oil/water composition convention independently of convergence.
     *
     * @return true when the documented condition holds
     */
    public boolean hasPhysicalPhaseIdentity() {
      return physicalPhaseIdentity;
    }

    /**
     * Returns the phase identity diagnostic.
     *
     * @return phase identity diagnostic
     */
    public String getPhaseIdentityDiagnostic() {
      return phaseIdentityDiagnostic;
    }

    /**
     * Returns the failure message.
     *
     * @return failure message
     */
    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
