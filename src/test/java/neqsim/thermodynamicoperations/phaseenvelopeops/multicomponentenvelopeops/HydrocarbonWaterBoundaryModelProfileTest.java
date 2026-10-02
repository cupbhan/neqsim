package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.google.gson.JsonObject;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import org.junit.jupiter.api.Test;

/**
 * Tests opt-in family-specific model profiles without changing the baseline fluid.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
class HydrocarbonWaterBoundaryModelProfileTest {

  /**
   * Checks or computes applies independent family multipliers to clones only.
   */
  @Test
  void appliesIndependentFamilyMultipliersToClonesOnly() {
    SystemInterface baseline = createFluid();
    double baselineWaterMethane = kij(baseline, "water", "methane");
    double baselineWaterNc10 = kij(baseline, "water", "nC10");
    double baselineHydrocarbon = kij(baseline, "methane", "nC10");

    HydrocarbonWaterBoundaryModelProfile profile = new HydrocarbonWaterBoundaryModelProfile("software-parity")
        .setWaterKijMultiplier(BoundaryFamily.GO_TO_GOW, 0.98325)
        .setWaterKijMultiplier(BoundaryFamily.OW_TO_GOW, 0.932);

    SystemInterface mainCandidate = profile.createTemplate(baseline, BoundaryFamily.GO_TO_GOW);
    SystemInterface secondaryCandidate = profile.createTemplate(baseline, BoundaryFamily.OW_TO_GOW);

    assertNotSame(baseline, mainCandidate);
    assertNotSame(baseline, secondaryCandidate);
    assertEquals(baselineWaterMethane * 0.98325, kij(mainCandidate, "water", "methane"), 1.0e-12);
    assertEquals(baselineWaterNc10 * 0.98325, kij(mainCandidate, "water", "nC10"), 1.0e-12);
    assertEquals(baselineWaterMethane * 0.932, kij(secondaryCandidate, "water", "methane"), 1.0e-12);
    assertEquals(baselineWaterNc10 * 0.932, kij(secondaryCandidate, "water", "nC10"), 1.0e-12);
    assertEquals(baselineHydrocarbon, kij(mainCandidate, "methane", "nC10"), 0.0);
    assertEquals(baselineHydrocarbon, kij(secondaryCandidate, "methane", "nC10"), 0.0);
    assertEquals(baselineWaterMethane, kij(baseline, "water", "methane"), 0.0);
    assertEquals(baselineWaterNc10, kij(baseline, "water", "nC10"), 0.0);

    SystemInterface identityCandidate = profile.createTemplate(baseline, BoundaryFamily.GW_TO_GOW);
    assertNotSame(baseline, identityCandidate);
    assertEquals(baselineWaterMethane, kij(identityCandidate, "water", "methane"), 0.0);
  }

  /**
   * Checks or computes reports governance and rejects invalid configuration.
   */
  @Test
  void reportsGovernanceAndRejectsInvalidConfiguration() {
    HydrocarbonWaterBoundaryModelProfile profile = new HydrocarbonWaterBoundaryModelProfile("software-parity")
        .setWaterKijMultiplier(BoundaryFamily.GO_TO_GOW, 0.98325);
    JsonObject description = profile.toJson(BoundaryFamily.GO_TO_GOW);

    assertEquals("software-parity", description.get("name").getAsString());
    assertEquals("GO_TO_GOW", description.get("boundaryFamily").getAsString());
    assertEquals(0.98325, description.get("waterKijMultiplier").getAsDouble(), 0.0);
    assertTrue(!description.get("defaultModelModified").getAsBoolean());
    assertTrue(description.get("experimentalCompatibilityProfile").getAsBoolean());
    assertTrue(description.get("requiresIndependentPhysicalValidation").getAsBoolean());
    assertTrue(!description.get("productionPromotionAllowed").getAsBoolean());

    assertThrows(IllegalArgumentException.class, () -> new HydrocarbonWaterBoundaryModelProfile(" "));
    assertThrows(IllegalArgumentException.class, () -> profile.setWaterKijMultiplier(BoundaryFamily.OW_TO_GOW, 0.0));
    assertThrows(IllegalArgumentException.class,
        () -> profile.setWaterKijMultiplier(BoundaryFamily.OW_TO_GOW, Double.NaN));

    SystemInterface dryFluid = new SystemSrkEos(300.0, 10.0);
    dryFluid.addComponent("methane", 1.0);
    dryFluid.setMixingRule("classic");
    assertThrows(IllegalArgumentException.class, () -> profile.createTemplate(dryFluid, BoundaryFamily.GO_TO_GOW));
  }

  /**
   * Requires serializable, explicit research settings and prevents invalid scaled interaction data.
   *
   * @throws Exception when the operation fails
   */
  @Test
  void preservesBaselineAcrossSerializationAndRejectsScalingOverflow() throws Exception {
    HydrocarbonWaterBoundaryModelProfile profile = new HydrocarbonWaterBoundaryModelProfile("research")
        .setWaterKijMultiplier(BoundaryFamily.GO_TO_GOW, 0.98);
    java.io.ByteArrayOutputStream bytes = new java.io.ByteArrayOutputStream();
    try (java.io.ObjectOutputStream output = new java.io.ObjectOutputStream(bytes)) {
      output.writeObject(profile);
    }
    HydrocarbonWaterBoundaryModelProfile recovered;
    try (java.io.ObjectInputStream input = new java.io.ObjectInputStream(
        new java.io.ByteArrayInputStream(bytes.toByteArray()))) {
      recovered = (HydrocarbonWaterBoundaryModelProfile) input.readObject();
    }
    assertEquals(0.98, recovered.getWaterKijMultiplier(BoundaryFamily.GO_TO_GOW), 0.0);
    assertEquals(1.0, recovered.getWaterKijMultiplier(BoundaryFamily.OW_TO_GOW), 0.0);
    SystemInterface baseline = createFluid();
    baseline.setBinaryInteractionParameter("water", "methane", 2.0);
    HydrocarbonWaterBoundaryModelProfile overflowing = new HydrocarbonWaterBoundaryModelProfile("overflow")
        .setWaterKijMultiplier(BoundaryFamily.GO_TO_GOW, Double.MAX_VALUE);
    assertThrows(IllegalArgumentException.class, () -> overflowing.createTemplate(baseline, BoundaryFamily.GO_TO_GOW));
    assertEquals(2.0, kij(baseline, "water", "methane"), 0.0);
  }

  /**
   * Checks or computes create fluid.
   *
   * @return computed create fluid result
   */
  private static SystemInterface createFluid() {
    SystemInterface fluid = new SystemSrkEos(350.0, 100.0);
    fluid.addComponent("methane", 0.20);
    fluid.addComponent("nC10", 0.15);
    fluid.addComponent("water", 0.65);
    fluid.setMixingRule("classic");
    fluid.setBinaryInteractionParameter("water", "methane", 0.45);
    fluid.setBinaryInteractionParameter("water", "nC10", 0.12);
    fluid.setBinaryInteractionParameter("methane", "nC10", 0.03);
    fluid.setMultiPhaseCheck(true);
    fluid.setMaxNumberOfPhases(3);
    fluid.init(0);
    return fluid;
  }

  /**
   * Checks or computes kij.
   *
   * @param fluid fluid
   * @param first first
   * @param second second
   * @return computed kij result
   */
  private static double kij(SystemInterface fluid, String first, String second) {
    int firstIndex = fluid.getPhase(0).getComponent(first).getComponentNumber();
    int secondIndex = fluid.getPhase(0).getComponent(second).getComponentNumber();
    return ((PhaseEosInterface) fluid.getPhase(0)).getEosMixingRule().getBinaryInteractionParameter(firstIndex,
        secondIndex);
  }
}
