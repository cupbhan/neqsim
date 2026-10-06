package neqsim.thermo.util.explicit;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.io.IOException;
import org.junit.jupiter.api.Test;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Independent composition modeling, source isolation and physical envelope contracts.
 *
 * @author NeqSim development team
 * @version 1.0
 */
class CompositionThermodynamicEngineTest {
  private static final ObjectMapper JSON = new ObjectMapper();

  /**
   * Build an independent synthetic gas/condensate composition without reference parameters.
   *
   * @return new request
   */
  private ObjectNode request() {
    ObjectNode request = JSON.createObjectNode().put("eos", "PR").put("temperatureC", 30).put("pressureMPa", 5);
    ArrayNode components = request.putArray("components");
    components.addObject().put("id", "gas").put("name", "C1").put("z", 0.7);
    components.addObject().put("id", "middle").put("name", "nC4").put("z", 0.2);
    components.addObject().put("id", "heavy").put("name", "C10").put("z", 0.1);
    return request;
  }

  /**
   * Both EOS families reproduce independently assembled public database systems and their normal properties.
   *
   * @throws IOException when JSON serialization fails
   */
  @Test
  void publicFlashMatchesIndependentDatabaseSystems() throws IOException {
    for (String eos : new String[] {"PR", "SRK"}) {
      ObjectNode input = request().put("eos", eos).put("volumeCorrection", true);
      String before = input.toString();
      JsonNode result = JSON.readTree(ThermodynamicEngine.compositionFlash(before));
      assertTrue(result.path("ok").asBoolean(), result.toString());
      SystemInterface expected = "PR".equals(eos) ? new SystemPrEos(303.15, 50) : new SystemSrkEos(303.15, 50);
      expected.addComponent("methane", 0.7);
      expected.addComponent("n-butane", 0.2);
      expected.addComponent("nC10", 0.1);
      expected.setMixingRule(2);
      expected.useVolumeCorrection(true);
      expected.setMultiPhaseCheck(true);
      expected.setUseIdealGasEnthalpyOfFormation(false);
      new ThermodynamicOperations(expected).TPflash();
      expected.init(3);
      expected.initProperties();
      assertEquals(expected.getNumberOfPhases(), result.path("phases").size());
      for (int p = 0; p < expected.getNumberOfPhases(); p++) {
        PhaseInterface phase = expected.getPhase(p);
        JsonNode actual = result.path("phases").get(p);
        assertEquals(expected.getBeta(p), actual.path("moleFraction").asDouble(), 1e-10);
        assertEquals(phase.getDensity("kg/m3"), actual.path("density").asDouble(), 1e-6);
        assertEquals(phase.getViscosity("cP"), actual.path("viscosity").asDouble(), 1e-9);
        assertEquals(phase.getEnthalpy("kJ/kg"), actual.path("enthalpy").asDouble(), 1e-6);
        assertEquals(phase.getCp("J/kgK"), actual.path("heatCapacity").asDouble(), 1e-6);
        assertEquals(phase.getZvolcorr(), actual.path("compressibilityFactor").asDouble(), 1e-10);
        assertEquals(0, actual.path("unavailableProperties").size());
      }
      assertEquals("independent-neqsim-model", result.path("diagnostics").path("parameterSource").asText());
      assertFalse(result.path("diagnostics").path("formationEnthalpyIncluded").asBoolean());
      assertEquals(before, input.toString());
    }
  }

  /** Verify all advertised short-name bindings use NeqSim's database values, including zero-feed rows. */
  @Test
  void aliasesAndSourceIndicesUseOwnDatabase() {
    ObjectNode input = request();
    ArrayNode components = input.putArray("components");
    String[] names = {"N2", "CO2", "H2O", "MeOH", "C1", "C2", "C3", "iC4", "nC4", "iC5", "nC5", "C6", "C7", "C8", "C9",
        "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20"};
    for (int i = 0; i < names.length; i++) {
      components.addObject().put("id", "component-" + i).put("name", names[i]).put("z", i == 4 ? 1 : 0);
    }
    CompositionModel model = new CompositionModel(input, 303.15, 50);
    assertEquals(names.length, model.bindings.size());
    for (int i = 0; i < names.length; i++) {
      JsonNode binding = model.bindings.get(i);
      SystemInterface independent = new SystemPrEos(303.15, 50);
      independent.addComponent(binding.path("nativeName").asText(), 1);
      assertEquals(i, binding.path("sourceIndex").asInt());
      assertEquals(independent.getComponent(0).getMolarMass(), binding.path("molarMassKgMol").asDouble(), 1e-12);
      assertEquals(independent.getComponent(0).getTC(), binding.path("criticalTemperatureK").asDouble(), 1e-12);
      assertEquals("neqsim-component-database", binding.path("parameterSource").asText());
    }
  }

  /** Unknown petroleum groups require explicit user characterization and never borrow reference parameters. */
  @Test
  void pseudoComponentsRequireUserTbpAndPureComponentsRejectOverrides() {
    ObjectNode input = request();
    ObjectNode heavy = (ObjectNode) input.path("components").get(2);
    heavy.put("name", "C10-C13");
    assertEquals("UNKNOWN_COMPONENT_REQUIRES_TBP", ThermodynamicEngine.compositionFlash(input).path("code").asText());
    heavy.putObject("tbp").put("molarMassKgMol", 0.17).put("densityRelative", 0.79);
    String before = input.toString();
    JsonNode result = ThermodynamicEngine.compositionFlash(input);
    assertTrue(result.path("ok").asBoolean(), result.toString());
    JsonNode binding = result.path("components").get(2);
    assertEquals("neqsim-tbp-characterization", binding.path("parameterSource").asText());
    SystemInterface independent = new SystemPrEos(303.15, 50);
    independent.addTBPfraction("reference", 1, 0.17, 0.79);
    assertEquals(independent.getComponent(0).getTC(), binding.path("criticalTemperatureK").asDouble(), 1e-9);
    assertEquals(before, input.toString());
    heavy.put("name", "C10");
    assertEquals("PURE_COMPONENT_TBP_CONFLICT", ThermodynamicEngine.compositionFlash(input).path("code").asText());
  }

  /** Declared petroleum fractions named C7 cannot be rebound to pure heptane; zero fractions need no estimate. */
  @Test
  void declaredPseudoIdentityOverridesAliasesAndZeroFractionsRemainAudited() {
    for (String kind : new String[] {"pseudo", "plus"}) {
      ObjectNode input = request();
      ObjectNode heavy = (ObjectNode) input.path("components").get(2);
      heavy.put("name", "C7").put("kind", kind);
      assertEquals("UNKNOWN_COMPONENT_REQUIRES_TBP", ThermodynamicEngine.compositionFlash(input).path("code").asText());
      assertEquals("UNKNOWN_COMPONENT_REQUIRES_TBP",
          ThermodynamicEngine.compositionEnvelope(input).path("code").asText());
      heavy.putObject("tbp").put("molarMassKgMol", 0.12).put("densityRelative", 0.77);
      JsonNode result = ThermodynamicEngine.compositionFlash(input);
      assertTrue(result.path("ok").asBoolean(), result.toString());
      assertEquals("neqsim-tbp-characterization", result.path("components").get(2).path("parameterSource").asText());
      assertEquals(0.12, result.path("components").get(2).path("molarMassKgMol").asDouble(), 1e-12);
    }
    ObjectNode input = request();
    ArrayNode original = input.withArray("components").deepCopy();
    ArrayNode components = input.putArray("components");
    components.addObject().put("id", "zero-pseudo").put("name", "C7").put("kind", "pseudo").put("z", 0);
    components.addAll(original);
    input.putArray("interactions").addObject().put("i", 1).put("j", 2).put("value", 0.02);
    String before = input.toString();
    JsonNode result = ThermodynamicEngine.compositionFlash(input);
    assertTrue(result.path("ok").asBoolean(), result.toString());
    assertEquals(1, result.path("omittedComponents").size());
    assertEquals(1, result.path("components").get(0).path("sourceIndex").asInt());
    for (JsonNode phase : result.path("phases")) {
      assertEquals(4, phase.path("composition").size());
      assertEquals("zero-pseudo", phase.path("composition").get(0).path("id").asText());
      assertEquals(0, phase.path("composition").get(0).path("moleFraction").asDouble(), 0);
    }
    assertEquals(before, input.toString());
  }

  /** Critical constants, source molar masses and reference BIP tables cannot cross this API boundary. */
  @Test
  void rejectsSourceParametersAndMalformedInputs() {
    for (String field : new String[] {"molarMass", "mwKgMol", "tcK", "pcBar", "cpJMolK", "penelouxLmol"}) {
      ObjectNode input = request();
      ((ObjectNode) input.path("components").get(0)).put(field, 123);
      assertEquals("UNSUPPORTED_PARAMETER", ThermodynamicEngine.compositionFlash(input).path("code").asText());
    }
    ObjectNode input = request();
    input.putArray("interactions").addObject().put("i", 0).put("j", 1).put("kij", 0.1);
    assertEquals("UNSUPPORTED_PARAMETER", ThermodynamicEngine.compositionFlash(input).path("code").asText());
    input = request().put("pressureMPa", 0);
    assertEquals("INVALID_STATE", ThermodynamicEngine.compositionFlash(input).path("code").asText());
    input = request();
    ((ObjectNode) input.path("components").get(0)).put("z", Double.NaN);
    assertEquals("INVALID_INPUT", ThermodynamicEngine.compositionFlash(input).path("code").asText());
  }

  /** The independent database BIP table is retained except for explicitly supplied local constants. */
  @Test
  void databaseInteractionsRemainIndependentAndOverridesSurviveFlash() {
    ObjectNode input = request();
    ((ObjectNode) input.path("components").get(2)).put("name", "CO2");
    CompositionModel defaultModel = new CompositionModel(input, 303.15, 50);
    double defaultPair = ((PhaseEosInterface) defaultModel.fluid.getPhase(0)).getEosMixingRule()
        .getBinaryInteractionParameter(0, 2);
    assertTrue(Math.abs(defaultPair) > 1e-6);
    input.putArray("interactions").addObject().put("i", 0).put("j", 1).put("value", 0.035);
    JsonNode result = ThermodynamicEngine.compositionFlash(input);
    assertTrue(result.path("ok").asBoolean(), result.toString());
    for (JsonNode pair : result.path("interactions")) {
      if (pair.path("i").asInt() == 1 && pair.path("j").asInt() == 0) {
        assertEquals(0.035, pair.path("value").asDouble(), 1e-12);
        assertEquals("user-constant-override", pair.path("source").asText());
      }
      if (pair.path("i").asInt() == 2 && pair.path("j").asInt() == 0) {
        assertEquals(defaultPair, pair.path("value").asDouble(), 1e-12);
        assertEquals("neqsim-mixing-rule", pair.path("source").asText());
      }
    }
  }

  /** Water and inhibitor phases use normal independently configured NeqSim properties. */
  @Test
  void aqueousPropertiesHaveIndependentAvailability() {
    ObjectNode input = request();
    ArrayNode components = input.putArray("components");
    components.addObject().put("id", "methane").put("name", "C1").put("z", 0.5);
    components.addObject().put("id", "water").put("name", "H2O").put("z", 0.45);
    components.addObject().put("id", "methanol").put("name", "MeOH").put("z", 0.05);
    JsonNode result = ThermodynamicEngine.compositionFlash(input);
    assertTrue(result.path("ok").asBoolean(), result.toString());
    boolean aqueous = false;
    for (JsonNode phase : result.path("phases")) {
      aqueous |= "aqueous".equals(phase.path("phase").asText());
      assertTrue(phase.path("density").asDouble() > 0, phase.toString());
      assertTrue(phase.path("heatCapacity").asDouble() > 0, phase.toString());
      for (String property : new String[] {"density", "viscosity", "enthalpy", "heatCapacity",
          "compressibilityFactor"}) {
        assertTrue(phase.path("propertyAvailability").has(property));
        if (phase.path(property).isNull()) {
          assertEquals("unavailable", phase.path("propertyAvailability").path(property).path("status").asText());
        }
      }
    }
    assertTrue(aqueous, result.toString());
  }

  /** A fresh database model traces physical sides while auditing zero and polar omissions by original index. */
  @Test
  void envelopeCalculatesPhysicalBranchesAndAuditsProjection() {
    ObjectNode input = request();
    ArrayNode components = input.putArray("components");
    components.addObject().put("id", "zero").put("name", "N2").put("z", 0);
    components.addObject().put("id", "wet").put("name", "H2O").put("z", 0.1);
    components.addObject().put("id", "gas").put("name", "C1").put("z", 0.63);
    components.addObject().put("id", "middle").put("name", "nC4").put("z", 0.18);
    components.addObject().put("id", "heavy").put("name", "C10").put("z", 0.09);
    assertEquals("POLAR_FEED_REQUIRES_EXPLICIT_DRY_PROJECTION",
        ThermodynamicEngine.compositionEnvelope(input).path("code").asText());
    input.putObject("envelope").put("hydrocarbonOnly", true);
    String before = input.toString();
    JsonNode result = ThermodynamicEngine.compositionEnvelope(input);
    assertTrue(result.path("ok").asBoolean(), result.toString());
    assertTrue(result.path("segments").size() >= 2, result.toString());
    boolean dew = false, bubble = false;
    for (JsonNode segment : result.path("segments")) {
      dew |= "dew".equals(segment.path("type").asText());
      bubble |= "bubble".equals(segment.path("type").asText());
      for (JsonNode point : segment.path("points")) {
        assertTrue(Double.isFinite(point.path("temperatureC").asDouble()));
        assertTrue(point.path("pressureMPa").asDouble() > 0);
      }
    }
    assertTrue(dew && bubble);
    assertEquals(5, result.path("sourceFeed").size());
    assertEquals(2, result.path("omittedComponents").size());
    assertEquals(2, result.path("feed").get(0).path("sourceIndex").asInt());
    assertEquals(4, result.path("components").get(2).path("sourceIndex").asInt());
    assertEquals(0.7, result.path("feed").get(0).path("moleFraction").asDouble(), 1e-12);
    assertEquals("independent-neqsim-model", result.path("provenance").path("parameterSource").asText());
    assertEquals("unavailable", result.path("availability").path("waterBoundary").path("status").asText());
    assertEquals(before, input.toString());
  }
}
