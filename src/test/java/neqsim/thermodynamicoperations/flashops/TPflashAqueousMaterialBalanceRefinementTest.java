package neqsim.thermodynamicoperations.flashops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

class TPflashAqueousMaterialBalanceRefinementTest {
  private static final String[] COMPONENTS = { "methane", "nC10", "water" };
  private static final double[] FEED = { 0.20, 0.60, 0.20 };

  @Test
  void ordinaryFlashRejectsNonConservativeAqueousEndpoint() {
    SystemInterface ordinary = createAndFlash(false);
    SystemInterface multiphase = createAndFlash(true);

    assertEquivalentEquilibrium(multiphase, ordinary);
    double firstGibbsEnergy = ordinary.getGibbsEnergy();

    new ThermodynamicOperations(ordinary).TPflash();
    ordinary.init(1);

    assertEquals(firstGibbsEnergy, ordinary.getGibbsEnergy(), 1.0e-8);
    assertEquivalentEquilibrium(multiphase, ordinary);
  }

  private SystemInterface createAndFlash(boolean multiphaseCheck) {
    SystemInterface system = new SystemSrkEos(230.0, 50.0);
    for (int componentIndex = 0; componentIndex < COMPONENTS.length; componentIndex++) {
      system.addComponent(COMPONENTS[componentIndex], FEED[componentIndex]);
    }
    system.setMixingRule(2);
    system.setMultiPhaseCheck(multiphaseCheck);
    new ThermodynamicOperations(system).TPflash();
    system.init(1);
    return system;
  }

  private void assertEquivalentEquilibrium(SystemInterface expected, SystemInterface actual) {
    assertEquals(2, expected.getNumberOfPhases());
    assertEquals(expected.getNumberOfPhases(), actual.getNumberOfPhases());
    assertEquals(PhaseType.OIL, actual.getPhase(0).getType());
    assertEquals(PhaseType.AQUEOUS, actual.getPhase(1).getType());
    assertEquals(expected.getGibbsEnergy(), actual.getGibbsEnergy(), 1.0e-8);

    double betaTotal = 0.0;
    for (int phaseIndex = 0; phaseIndex < expected.getNumberOfPhases(); phaseIndex++) {
      assertEquals(expected.getPhase(phaseIndex).getType(), actual.getPhase(phaseIndex).getType());
      assertEquals(expected.getBeta(phaseIndex), actual.getBeta(phaseIndex), 1.0e-12);
      assertTrue(actual.getBeta(phaseIndex) > 0.0);
      assertTrue(actual.getBeta(phaseIndex) < 1.0);
      assertTrue(Double.isFinite(actual.getPhase(phaseIndex).getZ()));
      assertTrue(actual.getPhase(phaseIndex).getZ() > 0.0);
      betaTotal += actual.getBeta(phaseIndex);

      double compositionTotal = 0.0;
      for (int componentIndex = 0; componentIndex < COMPONENTS.length; componentIndex++) {
        assertEquals(expected.getPhase(phaseIndex).getComponent(componentIndex).getx(),
            actual.getPhase(phaseIndex).getComponent(componentIndex).getx(), 1.0e-12);
        compositionTotal += actual.getPhase(phaseIndex).getComponent(componentIndex).getx();
      }
      assertEquals(1.0, compositionTotal, 1.0e-12);
    }
    assertEquals(1.0, betaTotal, 1.0e-12);
    assertTrue(maximumComponentMaterialBalanceResidual(actual) < 1.0e-8);
    assertTrue(maximumLogFugacityResidual(actual) < 1.0e-8);
  }

  private double maximumComponentMaterialBalanceResidual(SystemInterface system) {
    double maximumResidual = 0.0;
    for (int componentIndex = 0; componentIndex < COMPONENTS.length; componentIndex++) {
      double recoveredFeed = 0.0;
      for (int phaseIndex = 0; phaseIndex < system.getNumberOfPhases(); phaseIndex++) {
        recoveredFeed += system.getBeta(phaseIndex)
            * system.getPhase(phaseIndex).getComponent(componentIndex).getx();
      }
      maximumResidual = Math.max(maximumResidual, Math.abs(FEED[componentIndex] - recoveredFeed));
    }
    return maximumResidual;
  }

  private double maximumLogFugacityResidual(SystemInterface system) {
    double maximumResidual = 0.0;
    for (int componentIndex = 0; componentIndex < COMPONENTS.length; componentIndex++) {
      double firstLogFugacity = Math
          .log(Math.max(system.getPhase(0).getComponent(componentIndex).getx(), Double.MIN_NORMAL))
          + Math.log(system.getPhase(0).getComponent(componentIndex).getFugacityCoefficient());
      double secondLogFugacity = Math
          .log(Math.max(system.getPhase(1).getComponent(componentIndex).getx(), Double.MIN_NORMAL))
          + Math.log(system.getPhase(1).getComponent(componentIndex).getFugacityCoefficient());
      maximumResidual = Math.max(maximumResidual, Math.abs(firstLogFugacity - secondLogFugacity));
    }
    return maximumResidual;
  }
}
