package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import neqsim.thermo.phase.LiquidPhaseClassification;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Replays specified-phase roots and applies conservative finite-phase acceptance checks.
 *
 * <p>
 * Acceptance covers material balance, fugacity equality, phase identity, distinct compositions and the existing
 * gas/oil/water tangent-plane searches. A finite set of stability searches is not a proof of global stability and does
 * not establish experimental accuracy. Critical endpoints and vanishing phases require their dedicated solvers.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
public final class SpecifiedPhaseEquilibriumValidator {
  private static final double BALANCE_TOLERANCE = 1.0e-8;
  private static final double FUGACITY_TOLERANCE = 1.0e-7;
  private static final double MINIMUM_FRACTION = 1.0e-10;
  private static final double MINIMUM_DISTANCE = 1.0e-5;

  /** Prevents instantiation of this stateless validator. */
  private SpecifiedPhaseEquilibriumValidator() {
  }

  /**
   * Checks a private reconstructed equilibrium against the original overall component inventory.
   *
   * @param template original fluid, supplying overall composition and component order
   * @param reconstructed private system carrying the result's normalized phase compositions
   * @param slots requested physical phase families
   * @return immutable acceptance diagnostics; failures and incomplete searches are rejected
   * @throws IllegalArgumentException when systems or phase dimensions are incompatible
   */
  public static Result validate(SystemInterface template, SystemInterface reconstructed, CandidatePhase[] slots) {
    if (template == null || reconstructed == null || slots == null || slots.length < 2 || slots.length > 4
        || reconstructed.getNumberOfPhases() != slots.length
        || template.getNumberOfComponents() != reconstructed.getNumberOfComponents()) {
      throw new IllegalArgumentException("matching template, equilibrium and two to four phase families are required");
    }
    SystemInterface state = reconstructed.clone();
    int count = template.getNumberOfComponents();
    List<String> violations = new ArrayList<String>();
    double balance = 0.0;
    double normalization = 0.0;
    double fractionSum = 0.0;
    for (int component = 0; component < count; component++) {
      String name = template.getPhase(0).getComponent(component).getComponentName();
      if (!name.equals(state.getPhase(0).getComponent(component).getComponentName())) {
        throw new IllegalArgumentException("component order differs from the original fluid");
      }
    }
    for (int phase = 0; phase < slots.length; phase++) {
      if (slots[phase] == null) {
        throw new IllegalArgumentException("every phase family is required");
      }
      double beta = state.getBeta(phase);
      fractionSum += beta;
      if (!Double.isFinite(beta) || beta <= MINIMUM_FRACTION) {
        violations.add("VANISHING_OR_INVALID_PHASE_" + phase);
      }
      double sum = 0.0;
      for (int component = 0; component < count; component++) {
        double x = state.getPhase(phase).getComponent(component).getx();
        if (!Double.isFinite(x) || x < 0.0) {
          violations.add("INVALID_COMPOSITION_" + phase);
        }
        sum += x;
      }
      normalization = Math.max(normalization, Math.abs(sum - 1.0));
      if (state.getPhase(phase).getType() != phaseType(slots[phase])) {
        violations.add("PHASE_TYPE_MISMATCH_" + phase);
      }
      if (slots[phase] != CandidatePhase.GAS
          && LiquidPhaseClassification.classify(state.getPhase(phase)) != phaseType(slots[phase])) {
        violations.add("PHASE_IDENTITY_MISMATCH_" + phase);
      }
      for (int previous = 0; previous < phase; previous++) {
        double distance = 0.0;
        for (int component = 0; component < count; component++) {
          distance += Math.abs(state.getPhase(phase).getComponent(component).getx()
              - state.getPhase(previous).getComponent(component).getx());
        }
        if (!Double.isFinite(distance) || distance <= MINIMUM_DISTANCE) {
          violations.add("COINCIDENT_PHASES_" + previous + "_" + phase);
        }
      }
    }
    normalization = Math.max(normalization, Math.abs(fractionSum - 1.0));
    for (int component = 0; component < count; component++) {
      double inventory = 0.0;
      for (int phase = 0; phase < slots.length; phase++) {
        inventory += state.getBeta(phase) * state.getPhase(phase).getComponent(component).getx();
      }
      balance = Math.max(balance, Math.abs(inventory - template.getPhase(0).getComponent(component).getz()));
    }
    if (!Double.isFinite(normalization) || normalization > BALANCE_TOLERANCE) {
      violations.add("NORMALIZATION_FAILED");
    }
    if (!Double.isFinite(balance) || balance > BALANCE_TOLERANCE) {
      violations.add("MATERIAL_BALANCE_FAILED");
    }
    double fugacity = Double.NaN;
    double minimumTpd = Double.NaN;
    boolean stabilityChecked = false;
    if (violations.isEmpty()) {
      try {
        for (int phase = 0; phase < slots.length; phase++) {
          state.init(1, phase);
        }
        fugacity = 0.0;
        for (int phase = 1; phase < slots.length; phase++) {
          for (int component = 0; component < count; component++) {
            if (template.getPhase(0).getComponent(component).getz() == 0.0) {
              continue;
            }
            double first = Math.log(state.getPhase(0).getComponent(component).getx())
                + state.getPhase(0).getComponent(component).getLogFugacityCoefficient();
            double other = Math.log(state.getPhase(phase).getComponent(component).getx())
                + state.getPhase(phase).getComponent(component).getLogFugacityCoefficient();
            fugacity = Math.max(fugacity, Math.abs(first - other));
          }
        }
        if (!Double.isFinite(fugacity) || fugacity > FUGACITY_TOLERANCE) {
          violations.add("FUGACITY_REPLAY_FAILED");
        } else {
          IncipientPhaseStabilityAnalyzer.Result stability = new IncipientPhaseStabilityAnalyzer(state)
              .setMaximumIterations(800).setDampingFactor(0.2).setTolerances(1.0e-9, -1.0e-8).analyze();
          stabilityChecked = true;
          minimumTpd = Double.POSITIVE_INFINITY;
          for (Candidate candidate : stability.getTrials()) {
            if (!candidate.isConverged() || !Double.isFinite(candidate.getTangentPlaneDistance())) {
              stabilityChecked = false;
              violations.add("STABILITY_SEARCH_INCOMPLETE_" + candidate.getSeedPhase());
            } else {
              minimumTpd = Math.min(minimumTpd, candidate.getTangentPlaneDistance());
            }
          }
          if (!Double.isFinite(minimumTpd)) {
            stabilityChecked = false;
            violations.add("NO_FINITE_STABILITY_RESULT");
          }
          if (stability.hasNonTrivialInstability()) {
            violations.add("ADDITIONAL_PHASE_INSTABILITY");
          }
        }
      } catch (RuntimeException error) {
        violations.add("REPLAY_OR_STABILITY_FAILED: " + error.getMessage());
      }
    }
    return new Result(violations, balance, normalization, fugacity, minimumTpd, stabilityChecked);
  }

  /**
   * Maps a requested physical family to the EOS phase type.
   *
   * @param phase physical family
   * @return corresponding EOS phase type
   */
  private static PhaseType phaseType(CandidatePhase phase) {
    return phase == CandidatePhase.GAS ? PhaseType.GAS
        : phase == CandidatePhase.OIL ? PhaseType.OIL : PhaseType.AQUEOUS;
  }

  /**
   * Immutable replay diagnostics, with stability scoped to the executed trial searches.
   *
   * @author NeqSim contributors
   * @version 1.0
   */
  public static final class Result {
    private final List<String> violations;
    private final double materialBalanceResidual;
    private final double normalizationResidual;
    private final double fugacityResidual;
    private final double minimumTangentPlaneDistance;
    private final boolean stabilityChecked;

    /**
     * Stores a defensive copy of the validation outcome.
     *
     * @param violations failed checks
     * @param balance maximum component inventory error
     * @param normalization maximum mole-fraction normalization error
     * @param fugacity maximum replayed log-fugacity difference
     * @param minimumTpd lowest converged trial tangent-plane distance
     * @param stabilityChecked whether all required trial searches converged
     */
    private Result(List<String> violations, double balance, double normalization, double fugacity, double minimumTpd,
        boolean stabilityChecked) {
      this.violations = Collections.unmodifiableList(new ArrayList<String>(violations));
      this.materialBalanceResidual = balance;
      this.normalizationResidual = normalization;
      this.fugacityResidual = fugacity;
      this.minimumTangentPlaneDistance = minimumTpd;
      this.stabilityChecked = stabilityChecked;
    }

    /**
     * Reports acceptance under the documented finite-phase and stability-search checks.
     *
     * @return true only when every check passed and stability searches completed
     */
    public boolean isAccepted() {
      return violations.isEmpty() && stabilityChecked;
    }

    /**
     * Returns failed checks without exposing mutable state.
     *
     * @return immutable diagnostic list
     */
    public List<String> getViolations() {
      return violations;
    }

    /**
     * Returns the maximum component inventory error.
     *
     * @return dimensionless mole-fraction residual
     */
    public double getMaterialBalanceResidual() {
      return materialBalanceResidual;
    }

    /**
     * Returns the maximum phase or fraction normalization error.
     *
     * @return dimensionless normalization residual
     */
    public double getNormalizationResidual() {
      return normalizationResidual;
    }

    /**
     * Returns the maximum independently replayed log-fugacity difference.
     *
     * @return dimensionless fugacity residual, or NaN if replay was not attempted
     */
    public double getFugacityResidual() {
      return fugacityResidual;
    }

    /**
     * Returns the lowest converged tangent-plane trial value.
     *
     * @return dimensionless tangent-plane distance, or NaN if searches were not attempted
     */
    public double getMinimumTangentPlaneDistance() {
      return minimumTangentPlaneDistance;
    }

    /**
     * Reports completion of every required stability trial.
     *
     * @return true when all gas/oil/water searches applicable to the fluid converged
     */
    public boolean isStabilityChecked() {
      return stabilityChecked;
    }
  }
}
