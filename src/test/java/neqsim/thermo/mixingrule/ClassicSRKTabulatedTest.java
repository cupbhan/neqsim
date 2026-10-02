package neqsim.thermo.mixingrule;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import org.junit.jupiter.api.Test;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

/** Tests the tabulated temperature-dependent classic EOS mixing rule. */
class ClassicSRKTabulatedTest {

  @Test
  void interpolatesAndClampsSymmetricTemperatureTable() {
    SystemInterface system = createSystem();
    int co2 = system.getPhase(0).getComponent("CO2").getComponentNumber();
    int methane = system.getPhase(0).getComponent("methane").getComponentNumber();
    EosMixingRulesInterface rule = ((PhaseEos) system.getPhase(0)).getMixingRule();
    double[] knots = {300.0, 400.0, 500.0};
    double[] values = {0.05, 0.10, 0.02};

    rule.setBinaryInteractionParameterTemperatureTable(co2, methane, knots, values);
    EosMixingRuleHandler.ClassicSRKTabulated tabulated = (EosMixingRuleHandler.ClassicSRKTabulated) rule;

    assertEquals(0.05, tabulated.getkij(250.0, co2, methane), 1.0e-12);
    assertEquals(0.075, tabulated.getkij(350.0, co2, methane), 1.0e-12);
    assertEquals(0.06, tabulated.getkij(450.0, co2, methane), 1.0e-12);
    assertEquals(0.02, tabulated.getkij(550.0, co2, methane), 1.0e-12);
    assertEquals(tabulated.getkij(350.0, co2, methane), tabulated.getkij(350.0, methane, co2), 0.0);
    assertEquals(0.0005, tabulated.getkijdT(350.0, co2, methane), 1.0e-12);
    assertEquals(-0.0008, tabulated.getkijdT(450.0, co2, methane), 1.0e-12);
    assertEquals(0.0, tabulated.getkijdT(250.0, co2, methane), 0.0);
    assertEquals(0.0, tabulated.getkijdTdT(450.0, co2, methane), 0.0);
    assertArrayEquals(knots, rule.getBinaryInteractionParameterTemperatureKnots(methane, co2), 0.0);
    assertArrayEquals(values, rule.getBinaryInteractionParameterTemperatureValues(methane, co2), 0.0);
  }

  @Test
  void preservesTemperatureTablesWhenSystemIsCloned() {
    SystemInterface system = createSystem();
    int co2 = system.getPhase(0).getComponent("CO2").getComponentNumber();
    int methane = system.getPhase(0).getComponent("methane").getComponentNumber();
    ((PhaseEos) system.getPhase(0)).getMixingRule().setBinaryInteractionParameterTemperatureTable(co2, methane,
        new double[] {300.0, 400.0}, new double[] {0.05, 0.10});

    SystemInterface clonedSystem = system.clone();
    EosMixingRuleHandler.ClassicSRKTabulated clonedRule = (EosMixingRuleHandler.ClassicSRKTabulated) ((PhaseEos) clonedSystem
        .getPhase(0)).getMixingRule();

    assertEquals(0.075, clonedRule.getkij(350.0, co2, methane), 1.0e-12);
    double[] clonedValues = clonedRule.getBinaryInteractionParameterTemperatureValues(co2, methane);
    clonedValues[0] = 0.90;
    assertEquals(0.05, clonedRule.getkij(300.0, co2, methane), 1.0e-12);
  }

  @Test
  void rejectsInvalidTemperatureTables() {
    SystemInterface system = createSystem();
    int co2 = system.getPhase(0).getComponent("CO2").getComponentNumber();
    int methane = system.getPhase(0).getComponent("methane").getComponentNumber();
    EosMixingRulesInterface rule = ((PhaseEos) system.getPhase(0)).getMixingRule();

    assertThrows(IllegalArgumentException.class, () -> rule.setBinaryInteractionParameterTemperatureTable(co2, methane,
        new double[] {300.0}, new double[] {0.05}));
    assertThrows(IllegalArgumentException.class, () -> rule.setBinaryInteractionParameterTemperatureTable(co2, methane,
        new double[] {400.0, 300.0}, new double[] {0.05, 0.10}));
    assertThrows(IllegalArgumentException.class, () -> rule.setBinaryInteractionParameterTemperatureTable(co2, co2,
        new double[] {300.0, 400.0}, new double[] {0.05, 0.10}));
  }

  private static SystemInterface createSystem() {
    SystemInterface system = new SystemSrkEos(350.0, 50.0);
    system.addComponent("CO2", 0.2);
    system.addComponent("methane", 0.8);
    system.setMixingRule(EosMixingRuleType.CLASSIC_T_TABULATED);
    return system;
  }
}
