package neqsim.thermo.mixingrule;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryModelProfile;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;

/**
 * Test class for verifying the behavior of EOS (Equation of State) mixing rules in the NeqSim library.
 */
public class EosMixingRulesTest {
  @Test
  void testSetMixingRuleName() {
    neqsim.thermo.system.SystemPrEos testSystem = new neqsim.thermo.system.SystemPrEos(298.0, 10.0);
    testSystem.addComponent("nitrogen", 0.01);
    testSystem.addComponent("CO2", 0.01);
    testSystem.changeComponentName(testSystem.getComponent(0).getName(),
        (testSystem.getComponent(0).getName() + "__well1"));
    testSystem.changeComponentName(testSystem.getComponent(1).getName(),
        (testSystem.getComponent(1).getName() + "__well1"));

    testSystem.addComponent("nitrogen", 0.01);
    testSystem.addComponent("CO2", 0.01);
    testSystem.changeComponentName(testSystem.getComponent(2).getName(),
        (testSystem.getComponent(2).getName() + "__well2"));
    testSystem.changeComponentName(testSystem.getComponent(3).getName(),
        (testSystem.getComponent(3).getName() + "__well2"));

    testSystem.setMixingRule("classic");

    double kij = ((PhaseEos) testSystem.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1);
    double kij2 = ((PhaseEos) testSystem.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(3, 0);

    // Print kij
    assertEquals(-0.019997, kij, 1e-5);
    assertTrue(kij == kij2);
  }

  @Test
  void testMEGOil() {
    neqsim.thermo.system.SystemSrkCPAstatoil testSystem = new neqsim.thermo.system.SystemSrkCPAstatoil(298.0, 10.0);
    testSystem.addTBPfraction("C8", 0.01, 90.9 / 1000.0, 0.9);
    testSystem.addComponent("ethanol", 0.01);

    testSystem.setMixingRule("classic");

    double kij = ((PhaseEos) testSystem.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1);

    // Print kij
    assertEquals(-0.05, kij, 1e-5);
  }

  @Test
  void testHCoilInter() {
    neqsim.thermo.system.SystemPrEos testSystem = new neqsim.thermo.system.SystemPrEos(298.0, 10.0);
    testSystem.addTBPfraction("C8", 0.01, 90.9 / 1000.0, 0.9);
    testSystem.addComponent("CO2", 0.01);
    testSystem.changeComponentName(testSystem.getComponent(0).getName(),
        (testSystem.getComponent(0).getName() + "__well1"));
    testSystem.changeComponentName(testSystem.getComponent(1).getName(),
        (testSystem.getComponent(1).getName() + "__well1"));

    testSystem.addTBPfraction("C8", 0.01, 90.9 / 1000.0, 0.9);
    testSystem.addComponent("CO2", 0.01);
    testSystem.changeComponentName(testSystem.getComponent(2).getName(),
        (testSystem.getComponent(2).getName() + "__well2"));
    testSystem.changeComponentName(testSystem.getComponent(3).getName(),
        (testSystem.getComponent(3).getName() + "__well2"));

    testSystem.setMixingRule("classic");

    double kij = ((PhaseEos) testSystem.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1);
    double kij2 = ((PhaseEos) testSystem.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(3, 0);

    // Print kij
    assertEquals(0.1, kij, 1e-5);
    assertTrue(kij == kij2);
  }

  /** Keeps the official default and applies historical light-pseudo data only on an explicit clone. */
  @Test
  void waterAgainstALightPseudoUsesExplicitCompatibilityProfile() {
    neqsim.thermo.system.SystemSrkEos named = new neqsim.thermo.system.SystemSrkEos(298.0, 10.0);
    named.addComponent("water", 0.5);
    named.addComponent("n-heptane", 0.5);
    named.setMixingRule("classic");
    double namedKij = ((PhaseEos) named.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1);

    neqsim.thermo.system.SystemSrkEos pseudo = new neqsim.thermo.system.SystemSrkEos(298.0, 10.0);
    pseudo.addComponent("water", 0.5);
    pseudo.addTBPfraction("C7", 0.5, 96.0 / 1000.0, 0.75);
    pseudo.setMixingRule("classic");
    double pseudoKij = ((PhaseEos) pseudo.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1);

    assertEquals(0.5, namedKij, 1e-9, "database value for the named molecule");
    assertEquals(0.2, pseudoKij, 1e-9, "ordinary pseudo components retain the official model definition");
    HydrocarbonWaterBoundaryModelProfile profile = new HydrocarbonWaterBoundaryModelProfile("historical-light-pseudo")
        .setLightPseudoWaterCompatibility(true);
    BoundaryFamily family = BoundaryFamily.GO_TO_GOW;
    SystemInterface configured = profile.createTemplate(pseudo, family);
    double configuredKij = ((PhaseEos) configured.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1);
    assertEquals(namedKij, configuredKij, 1e-9, "explicit historical compatibility value");
    assertEquals(pseudoKij, ((PhaseEos) pseudo.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1), 0.0,
        "profile must preserve the official baseline");
    assertTrue(profile.toJson(family).get("lightPseudoWaterCompatibility").getAsBoolean());
  }

  @Test
  void waterAgainstAPseudoBeyondTheDatabaseRangeKeepsTheHistoricalDefault() {
    neqsim.thermo.system.SystemSrkEos testSystem = new neqsim.thermo.system.SystemSrkEos(298.0, 10.0);
    testSystem.addComponent("water", 0.5);
    testSystem.addTBPfraction("C30", 0.5, 420.0 / 1000.0, 0.95);
    testSystem.setMixingRule("classic");

    double kij = ((PhaseEos) testSystem.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1);

    // The database holds no water pair above n-nonane, so heavy fractions are left where they were rather than
    // extrapolated: the one heavy-fluid comparison available gets worse when 0.5 is pushed onto them.
    assertEquals(0.2, kij, 1e-9);
  }

  @Test
  void testCalculatedInteractionParametersUseCriticalVolumes() {
    neqsim.thermo.system.SystemPrEos testSystem = new neqsim.thermo.system.SystemPrEos(298.0, 10.0);
    testSystem.addComponent("methane", 1.0);
    testSystem.addComponent("n-heptane", 1.0);

    try {
      testSystem.calcKIJ(true);
      testSystem.setMixingRule("classic");

      double kij = ((PhaseEos) testSystem.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(0, 1);
      double expected = BIPEstimator.estimateChuehPrausnitz(testSystem.getComponent("methane"),
          testSystem.getComponent("n-heptane"));

      assertEquals(expected, kij, 1e-12);
      assertTrue(kij > 0.0, "Calculated critical-volume kij should not collapse to zero");
    } finally {
      testSystem.calcKIJ(false);
    }
  }
}
