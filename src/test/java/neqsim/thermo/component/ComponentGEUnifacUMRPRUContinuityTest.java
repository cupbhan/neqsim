package neqsim.thermo.component;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemUMRPRUMCEos;

/**
 * Checks UMR-PRU activity coefficients and their derivatives on the same logarithmic scale.
 *
 * @author NeqSim developers
 * @version 1.0
 */
class ComponentGEUnifacUMRPRUContinuityTest {
  /** Trace toluene fraction representative of the rich-gas envelope regression. */
  private static final double TRACE_TOLUENE_FRACTION = 0.00025;

  /** Representable small activity coefficients must retain the logarithm used by their derivatives. */
  @Test
  void representableSmallActivityCoefficientRemainsConsistent() {
    for (double temperature : new double[] {115.0, 120.0, 127.0}) {
      ComponentGEInterface component = tolueneActivity(temperature, TRACE_TOLUENE_FRACTION);
      double expected = Math.exp(component.getLnGamma());
      assertTrue(expected > 0.0 && expected < 1.0e-10,
          "The fixture must exercise a representable activity coefficient below the former cutoff");
      assertEquals(expected, component.getGamma(), expected * 1.0e-12,
          "A finite small activity coefficient must not be replaced by an ideal coefficient");
      assertEquals(component.getLnGamma(), Math.log(component.getGamma()), 1.0e-12,
          "The mixing-rule logarithm must agree with the stored log activity coefficient");
    }
  }

  /** Temperature differences of the returned activity coefficient must agree with its analytic derivative. */
  @Test
  void smallActivityCoefficientHasConsistentTemperatureDerivative() {
    double temperature = 120.0;
    ComponentGEInterface component = tolueneActivity(temperature, TRACE_TOLUENE_FRACTION);
    assertTrue(Math.abs(component.getLnGammadt()) > 1.0,
        "The fixture must have a nonzero temperature derivative hidden by the former ideal reset");
    for (double step : new double[] {1.0e-4, 2.0e-4}) {
      double plus = Math.log(tolueneActivity(temperature + step, TRACE_TOLUENE_FRACTION).getGamma());
      double minus = Math.log(tolueneActivity(temperature - step, TRACE_TOLUENE_FRACTION).getGamma());
      assertEquals(component.getLnGammadt(), (plus - minus) / (2.0 * step), 1.0e-5,
          "The returned gamma must follow the same smooth temperature dependence as its derivative");
    }
  }

  /** Ordinary compositions and warmer temperatures must retain finite, consistent activity coefficients. */
  @Test
  void ordinaryActivityCoefficientsRemainConsistent() {
    for (double[] state : new double[][] {{200.0, TRACE_TOLUENE_FRACTION}, {115.0, 0.5}}) {
      ComponentGEInterface component = tolueneActivity(state[0], state[1]);
      assertTrue(Double.isFinite(component.getGamma()) && component.getGamma() > 1.0e-10,
          "The controls must remain outside the small-coefficient regression regime");
      assertEquals(component.getLnGamma(), Math.log(component.getGamma()), 1.0e-12);
    }
  }

  /**
   * Evaluates the activity model directly through a fresh binary UMR-PRU system without a flash.
   *
   * @param temperature temperature in kelvin
   * @param tolueneFraction normalized toluene mole fraction
   * @return initialized toluene activity component
   */
  private static ComponentGEInterface tolueneActivity(double temperature, double tolueneFraction) {
    SystemInterface fluid = new SystemUMRPRUMCEos(temperature, 1.5);
    fluid.addComponent("methane", 1.0 - tolueneFraction);
    fluid.addComponent("Toluene", tolueneFraction);
    fluid.setMixingRule("HV", "UNIFAC_UMRPRU");
    fluid.init(0);
    fluid.init(3);
    PhaseInterface activityPhase = ((PhaseEos) fluid.getPhase(0)).getMixingRule().getGEPhase();
    return (ComponentGEInterface) activityPhase.getComponent("Toluene");
  }
}
