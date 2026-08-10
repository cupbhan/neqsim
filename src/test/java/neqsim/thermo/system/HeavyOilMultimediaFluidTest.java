package neqsim.thermo.system;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Locale;
import neqsim.NeqSimTest;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import org.junit.jupiter.api.Test;

/** Tests for the versioned heavy-oil multimedia delivery fluid. */
class HeavyOilMultimediaFluidTest extends NeqSimTest {
  private static final double TEMPERATURE_K = 273.15 + 220.0;
  private static final double PRESSURE_BARA = 55.0;

  @Test
  void createsAllReleasedDoseAnchors() {
    for (HeavyOilMultimediaFluid.Nh3Dose dose : HeavyOilMultimediaFluid.Nh3Dose.values()) {
      SystemInterface fluid = HeavyOilMultimediaFluid.create(dose, TEMPERATURE_K, PRESSURE_BARA);
      assertEquals(16, fluid.getNumberOfComponents());
      assertEquals(dose.getMolPercent() / 100.0, componentFraction(fluid, "NH3"), 1.0e-12);
      assertEquals(1.0, totalComposition(fluid), 1.0e-12);
      assertTrue(fluid.getMixingRuleName().toLowerCase(Locale.ROOT).contains("tabulated"));
    }
  }

  @Test
  void appliesDoseAndTemperatureDependentInteractionContract() {
    assertDoseInteractions(HeavyOilMultimediaFluid.Nh3Dose.NH3_0P1, 0.0, 0.0);
    assertDoseInteractions(HeavyOilMultimediaFluid.Nh3Dose.NH3_1, -0.11, -0.18);
    assertDoseInteractions(HeavyOilMultimediaFluid.Nh3Dose.NH3_5, -0.125, -0.16);

    SystemInterface fluid = HeavyOilMultimediaFluid.create(5.0, TEMPERATURE_K, PRESSURE_BARA);
    EosMixingRulesInterface rule = mixingRule(fluid);
    assertEquals(0.143, tableValues(rule, fluid, "CO2", "H2O")[2], 1.0e-12);
    assertEquals(-1.90, tableValues(rule, fluid, "N2", "H2O")[4], 1.0e-12);
    assertEquals(-0.40, tableValues(rule, fluid, "N2", "C10-C12")[4], 1.0e-12);
  }

  @Test
  void publicBuilderPresetRunsThreePhaseFlash() {
    SystemInterface fluid = FluidBuilder.heavyOilMultimedia(TEMPERATURE_K, PRESSURE_BARA, 5.0);
    ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
    operations.TPflash();
    assertEquals(3, fluid.getNumberOfPhases());
    assertEquals(1.0, totalPhaseBeta(fluid), 1.0e-10);
  }

  @Test
  void exposesReleaseBoundaryAndLimitations() {
    assertEquals("HeavyOil-DEV-H2O65-CO27-N23-NH3-SRK-v1", HeavyOilMultimediaFluid.CONTRACT_ID);
    assertTrue(HeavyOilMultimediaFluid.isInsideReleasedEnvelope(20.0, 1.0));
    assertTrue(HeavyOilMultimediaFluid.isInsideReleasedEnvelope(250.0, 100.0));
    assertFalse(HeavyOilMultimediaFluid.isInsideReleasedEnvelope(19.9, 55.0));
    assertFalse(HeavyOilMultimediaFluid.isFieldQualified());
    assertFalse(HeavyOilMultimediaFluid.isReactive());
    assertThrows(IllegalArgumentException.class,
        () -> HeavyOilMultimediaFluid.create(2.0, TEMPERATURE_K, PRESSURE_BARA));
  }

  private static void assertDoseInteractions(HeavyOilMultimediaFluid.Nh3Dose dose, double expectedWaterKij,
      double expectedHeavyKij) {
    SystemInterface fluid = HeavyOilMultimediaFluid.create(dose, TEMPERATURE_K, PRESSURE_BARA);
    EosMixingRulesInterface rule = mixingRule(fluid);
    int nh3 = indexOf(fluid, "NH3");
    assertEquals(expectedWaterKij, rule.getBinaryInteractionParameter(nh3, indexOf(fluid, "H2O")), 1.0e-12);
    assertEquals(expectedHeavyKij, rule.getBinaryInteractionParameter(nh3, indexOf(fluid, "C10-C12")), 1.0e-12);
    assertEquals(0.0, rule.getBinaryInteractionParameter(nh3, indexOf(fluid, "CO2")), 1.0e-12);
    assertEquals(0.0, rule.getBinaryInteractionParameter(nh3, indexOf(fluid, "N2")), 1.0e-12);
  }

  private static EosMixingRulesInterface mixingRule(SystemInterface fluid) {
    return ((PhaseEosInterface) fluid.getPhase(0)).getEosMixingRule();
  }

  private static double[] tableValues(EosMixingRulesInterface rule, SystemInterface fluid, String first,
      String second) {
    return rule.getBinaryInteractionParameterTemperatureValues(indexOf(fluid, first), indexOf(fluid, second));
  }

  private static int indexOf(SystemInterface fluid, String canonicalName) {
    for (int index = 0; index < fluid.getPhase(0).getNumberOfComponents(); index++) {
      String name = fluid.getPhase(0).getComponent(index).getComponentName();
      String normalized = name.endsWith("_PC") ? name.substring(0, name.length() - 3) : name;
      if (canonicalName.equals(normalized) || ("H2O".equals(canonicalName) && "water".equals(name))
          || ("CO2".equals(canonicalName) && "carbon dioxide".equals(name))
          || ("N2".equals(canonicalName) && "nitrogen".equals(name))
          || ("NH3".equals(canonicalName) && "ammonia".equals(name))) {
        return index;
      }
    }
    throw new IllegalArgumentException("Component not found: " + canonicalName);
  }

  private static double componentFraction(SystemInterface fluid, String canonicalName) {
    return fluid.getPhase(0).getComponent(indexOf(fluid, canonicalName)).getz();
  }

  private static double totalComposition(SystemInterface fluid) {
    double total = 0.0;
    for (int index = 0; index < fluid.getPhase(0).getNumberOfComponents(); index++) {
      total += fluid.getPhase(0).getComponent(index).getz();
    }
    return total;
  }

  private static double totalPhaseBeta(SystemInterface fluid) {
    double total = 0.0;
    for (int phase = 0; phase < fluid.getNumberOfPhases(); phase++) {
      total += fluid.getPhase(phase).getBeta();
    }
    return total;
  }
}
