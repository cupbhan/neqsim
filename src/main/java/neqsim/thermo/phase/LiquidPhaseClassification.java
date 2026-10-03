package neqsim.thermo.phase;

import neqsim.thermo.component.ComponentInterface;

/**
 * Shared mass-based liquid-family convention for cubic EOS states and research acceptance.
 *
 * <p>
 * This labels an already identified liquid root. It does not decide whether a root is liquid, prove equilibrium, or
 * replace stability analysis. Hydrocarbon, inert and TBP components form one family; the remaining components form the
 * aqueous family, matching the personal PhaseEos convention. Equal family masses retain the aqueous label.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
public final class LiquidPhaseClassification {
  /** Prevents construction of the stateless classifier. */
  private LiquidPhaseClassification() {
  }

  /**
   * Classifies the current composition of a liquid phase.
   *
   * @param phase component properties and current mole fractions
   * @return oil or aqueous family
   * @throws IllegalArgumentException when the phase is null
   */
  public static PhaseType classify(PhaseInterface phase) {
    return classify(phase, null);
  }

  /**
   * Classifies a trial liquid composition using the template's component properties.
   *
   * <p>
   * Callers validate finite, nonnegative compositions separately. A null vector selects the phase's current mole
   * fractions. No EOS parameters, component amounts or phase labels are changed.
   *
   * @param phase template supplying component properties in the same order
   * @param composition trial mole fractions, or null for the current phase composition
   * @return oil or aqueous family
   * @throws IllegalArgumentException when the phase or vector dimensions are invalid
   */
  public static PhaseType classify(PhaseInterface phase, double[] composition) {
    if (phase == null || composition != null && composition.length != phase.getNumberOfComponents()) {
      throw new IllegalArgumentException("matching phase and component composition are required");
    }
    double hydrocarbons = 0.0;
    double aqueous = 0.0;
    for (int index = 0; index < phase.getNumberOfComponents(); index++) {
      ComponentInterface component = phase.getComponent(index);
      double amount = (composition == null ? component.getx() : composition[index]) * component.getMolarMass();
      if ((component.isHydrocarbon() || component.isInert() || component.isIsTBPfraction())
          && !component.getName().equals("water") && !component.getName().equals("water_PC")) {
        hydrocarbons += amount;
      } else {
        aqueous += amount;
      }
    }
    return hydrocarbons > aqueous ? PhaseType.OIL : PhaseType.AQUEOUS;
  }
}
