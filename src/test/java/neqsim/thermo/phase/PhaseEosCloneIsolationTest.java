package neqsim.thermo.phase;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkPenelouxEos;

class PhaseEosCloneIsolationTest extends NeqSimTest {

  @Test
  void huronVidalCloneHasIndependentMutableMixingRuleState() {
    SystemInterface fluid = new SystemSrkPenelouxEos(330.0, 20.0);
    fluid.addComponent("methane", 0.55);
    fluid.addComponent("CO2", 0.15);
    fluid.addComponent("water", 0.30);
    fluid.createDatabase(true);
    fluid.setMixingRule("HV", "NRTL");
    fluid.init(0);
    fluid.init(1);

    PhaseEos original = (PhaseEos) fluid.getPhase(0);
    PhaseEos first = original.clone();
    PhaseEos second = original.clone();

    assertNotSame(original.mixSelect, first.mixSelect);
    assertNotSame(first.mixSelect, second.mixSelect);
    assertNotSame(original.getEosMixingRule(), first.getEosMixingRule());
    assertNotSame(first.getEosMixingRule(), second.getEosMixingRule());
    assertEquals(original.getMixingRuleName(), first.getMixingRuleName());
    assertEquals(original.getEosMixingRule().getBinaryInteractionParameter(0, 1),
        first.getEosMixingRule().getBinaryInteractionParameter(0, 1), 0.0);
  }
}
