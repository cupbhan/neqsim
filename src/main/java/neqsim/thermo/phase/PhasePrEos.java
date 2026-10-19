package neqsim.thermo.phase;

import neqsim.thermo.component.ComponentPR;

/**
 * PhasePrEos class.
 *
 * @author Even Solbraa
 * @version $Id: $Id
 */
public class PhasePrEos extends PhaseEos {
  /** Serialization version UID. */
  private static final long serialVersionUID = 1000;
  /** Accepted continuation volume; NaN keeps ordinary phase-labelled volume selection. */
  private double volumeRootReference = Double.NaN;

  /**
   * Selects the cubic outer root nearest an accepted continuation volume.
   *
   * @param volume accepted molar volume in internal units (m3/mol times 1e5), or NaN to disable
   * @throws IllegalArgumentException if a finite reference is not positive
   */
  public void setVolumeRootReference(double volume) {
    if (!Double.isNaN(volume) && (!Double.isFinite(volume) || volume <= 0.0)) {
      throw new IllegalArgumentException("Volume reference must be positive and finite, or NaN");
    }
    volumeRootReference = volume;
  }

  /** {@inheritDoc} */
  @Override
  public double molarVolume(double pressure, double temperature, double a, double b, PhaseType pt)
      throws neqsim.util.exception.IsNaNException, neqsim.util.exception.TooManyIterationsException {
    // Old serialized phases initialize this added field to zero: keep their ordinary root policy.
    if (!Double.isFinite(volumeRootReference) || volumeRootReference <= 0.0) {
      return super.molarVolume(pressure, temperature, a, b, pt);
    }
    double vapor = molarVolumeAnalytical(pressure, temperature, PhaseType.GAS);
    double liquid = molarVolumeAnalytical(pressure, temperature, PhaseType.LIQUID);
    boolean validVapor = Double.isFinite(vapor) && vapor > getb();
    boolean validLiquid = Double.isFinite(liquid) && liquid > getb();
    if (!validVapor && !validLiquid) {
      throw new neqsim.util.exception.IsNaNException(this, "molarVolume", "physical cubic volume root");
    }
    double selected = !validLiquid || validVapor
        && Math.abs(Math.log(vapor / volumeRootReference)) < Math.abs(Math.log(liquid / volumeRootReference)) ? vapor
            : liquid;
    setMolarVolume(selected);
    Z = pressure * selected / (R * temperature);
    return selected;
  }

  /**
   * Constructor for PhasePrEos.
   */
  public PhasePrEos() {
    thermoPropertyModelName = "PR-EoS";
    uEOS = 2;
    wEOS = -1;
    delta1 = 1.0 + Math.sqrt(2.0);
    delta2 = 1.0 - Math.sqrt(2.0);
  }

  /** {@inheritDoc} */
  @Override
  public PhasePrEos clone() {
    PhasePrEos clonedPhase = null;
    try {
      clonedPhase = (PhasePrEos) super.clone();
    } catch (Exception ex) {
      logger.error("Cloning failed.", ex);
    }

    return clonedPhase;
  }

  /** {@inheritDoc} */
  @Override
  public void addComponent(String name, double moles, double molesInPhase, int compNumber) {
    super.addComponent(name, molesInPhase, compNumber);
    componentArray[compNumber] = new ComponentPR(name, moles, molesInPhase, compNumber);
  }
}
