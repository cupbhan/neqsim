package neqsim.thermo.util.explicit;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import org.junit.jupiter.api.Test;
import neqsim.thermo.system.SystemElectrolyteCPAstatoil;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Molecular electrolyte accounting and independent public-engine contracts.
 *
 * @author NeqSim development team
 * @version 1.0
 */
class ElectrolyteCompositionTest {
  private static final ObjectMapper JSON = new ObjectMapper();

  /**
   * Create a synthetic gas/brine feed with no simulator parameters.
   *
   * @param temperatureC temperature in degrees Celsius
   * @param pressureMPa absolute pressure in MPa
   * @return independent request
   */
  private ObjectNode request(double temperatureC, double pressureMPa) {
    ObjectNode input = JSON.createObjectNode().put("eos", "ELECTROLYTE-CPA").put("temperatureC", temperatureC)
        .put("pressureMPa", pressureMPa).put("volumeCorrection", true);
    ArrayNode rows = input.putArray("components");
    rows.addObject().put("id", "water-feed").put("name", "H2O").put("kind", "pure").put("z", 0.464);
    rows.addObject().put("id", "salt-feed").put("name", "NaCl").put("kind", "pure").put("z", 0.036);
    rows.addObject().put("id", "gas-feed").put("name", "C1").put("kind", "pure").put("z", 0.49);
    rows.addObject().put("id", "ethane-feed").put("name", "C2").put("kind", "pure").put("z", 0.01);
    return input;
  }

  /** Check formula-unit balances and native physical results at two independent states. */
  @Test
  void conservesFormulaUnitsAndMatchesIndependentNativeSystem() {
    for (double[] state : new double[][] {{80, 20}, {25, 1}}) {
      ObjectNode input = request(state[0], state[1]);
      JsonNode original = input.deepCopy();
      ObjectNode result = ThermodynamicEngine.compositionFlash(input);
      assertTrue(result.path("ok").asBoolean(), result.toString());
      assertEquals(original, input);
      SystemInterface nativeFluid = new SystemElectrolyteCPAstatoil(state[0] + 273.15, state[1] * 10);
      nativeFluid.addComponent("water", 0.464);
      nativeFluid.addComponent("Na+", 0.036);
      nativeFluid.addComponent("Cl-", 0.036);
      nativeFluid.addComponent("methane", 0.49);
      nativeFluid.addComponent("ethane", 0.01);
      nativeFluid.setMixingRule(10);
      nativeFluid.useVolumeCorrection(true);
      nativeFluid.setMultiPhaseCheck(true);
      nativeFluid.setUseIdealGasEnthalpyOfFormation(false);
      new ThermodynamicOperations(nativeFluid).TPflash();
      nativeFluid.init(3);
      nativeFluid.initProperties();
      assertEquals(nativeFluid.getNumberOfPhases(), result.path("phases").size());
      double betaTotal = 0;
      double[] recovered = new double[4];
      for (int p = 0; p < result.path("phases").size(); p++) {
        JsonNode phase = result.path("phases").get(p);
        double beta = phase.path("moleFraction").asDouble();
        betaTotal += beta;
        assertEquals(nativeFluid.getPhase(p).getDensity("kg/m3"), phase.path("density").asDouble(), 1e-5);
        assertEquals(nativeFluid.getBeta(p), phase.path("ionicMoleFraction").asDouble(), 1e-8);
        assertEquals(nativeFluid.getPhase(p).getMolarMass(),
            phase.path("molarMassKgMol").asDouble() * phase.path("formulaUnitFactor").asDouble(), 1e-10);
        double xTotal = 0;
        for (int i = 0; i < recovered.length; i++) {
          JsonNode row = phase.path("composition").get(i);
          assertEquals(input.path("components").get(i).path("id").asText(), row.path("id").asText());
          xTotal += row.path("moleFraction").asDouble();
          recovered[i] += beta * row.path("moleFraction").asDouble();
        }
        assertEquals(1, xTotal, 1e-8);
      }
      assertEquals(1, betaTotal, 1e-9);
      for (int i = 0; i < recovered.length; i++) {
        assertEquals(input.path("components").get(i).path("z").asDouble(), recovered[i], 1e-7);
      }
      assertTrue(result.path("diagnostics").path("maxPhaseChargeResidual").asDouble() < 1e-8);
      assertTrue(result.path("diagnostics").path("maxLogFugacityResidual").asDouble() < 1e-6);
      assertFalse(result.path("components").get(1).has("criticalTemperatureK"));
    }
  }

  /** Reject external salt parameters, missing water and mixed ion/formula-unit feeds. */
  @Test
  void rejectsIncompatibleInputWithoutBorrowingParameters() {
    ObjectNode input = request(80, 20);
    ((ObjectNode) input.path("components").get(1)).put("criticalTemperatureK", 700);
    assertEquals("UNSUPPORTED_PARAMETER", ThermodynamicEngine.compositionFlash(input).path("code").asText());
    input = request(80, 20);
    ((ObjectNode) input.path("components").get(0)).put("name", "C3");
    assertEquals("ELECTROLYTE_WATER_REQUIRED", ThermodynamicEngine.compositionFlash(input).path("code").asText());
    input = request(80, 20);
    ((ObjectNode) input.path("components").get(1)).put("name", "Na+");
    assertEquals("ELECTROLYTE_COMPONENT_SCOPE", ThermodynamicEngine.compositionFlash(input).path("code").asText());
    input = request(80, 20);
    input.put("eos", "SRK");
    assertFalse(ThermodynamicEngine.compositionFlash(input).path("ok").asBoolean());
  }
}
