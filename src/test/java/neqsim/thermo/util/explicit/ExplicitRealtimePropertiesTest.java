package neqsim.thermo.util.explicit;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.io.IOException;
import java.io.InputStream;
import org.junit.jupiter.api.Test;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;

/**
 * Physical derivative, explicit-parameter and immutable-envelope contracts for realtime consumers.
 *
 * @author NeqSim development team
 * @version 1.0
 */
class ExplicitRealtimePropertiesTest {
  private static final ObjectMapper JSON = new ObjectMapper();

  /**
   * Load the independently parameterized public synthetic binary mixture.
   *
   * @return mutable request fixture
   * @throws IOException if the fixture cannot be read
   */
  private ObjectNode model() throws IOException {
    try (InputStream stream = getClass().getResourceAsStream("/neqsim/thermo/explicit/example-request.json")) {
      return (ObjectNode) JSON.readTree(stream);
    }
  }

  /**
   * Fixed-composition heat capacity must reproduce the analytic constant-kij EOS derivative.
   *
   * @throws IOException if the fixture cannot be read
   */
  @Test
  void cpMatchesAnalyticEosAtInteriorAndThermalDomainEdges() throws IOException {
    ObjectNode request = model();
    request.put("mixingRule", "Classic");
    ((ObjectNode) request.path("interactions").path(0)).put("mixingRule", "Classical").put("kijT", 0);
    ArrayNode grid = request.putArray("grid");
    for (double t : new double[] {25, 80, 150}) {
      grid.addObject().put("temperatureC", t).put("pressureMPa", 5);
    }
    JsonNode result = ThermodynamicEngine.calculate(request, true);
    for (JsonNode point : result.path("points")) {
      assertEquals("ok", point.path("status").asText());
      for (JsonNode phase : point.path("phases")) {
        assertEquals("available", phase.path("properties").path("cp").path("status").asText(), phase.toString());
        SystemInterface system = ExplicitProperties.state(request, phase.path("composition"),
            phase.path("type").asText(), point.path("temperatureC").asDouble() + 273.15, 50);
        PhaseInterface analytic = system.getPhase(0);
        double ideal = analytic.getCp0();
        double expected = ideal + (analytic.getCp("J/molK") - ideal) * 8.3147295 / ExplicitFlash.R;
        assertEquals(expected, phase.path("cpJMolK").asDouble(), 0.04);
      }
    }
  }

  /**
   * Huron-Vidal and temperature-dependent interactions must retain available derivative checks.
   *
   * @throws IOException if the fixture cannot be read
   */
  @Test
  void cpRetainsHvAndTemperatureDependentInteractionsWithoutChangingInput() throws IOException {
    ObjectNode request = model();
    String before = request.toString();
    for (JsonNode point : ThermodynamicEngine.calculate(request, true).path("points")) {
      for (JsonNode phase : point.path("phases")) {
        assertEquals("available", phase.path("properties").path("cp").path("status").asText(), phase.toString());
        assertTrue(phase.path("cpJMolK").asDouble() > 0);
      }
    }
    assertEquals(before, request.toString());
    request.putArray("grid").addObject().put("temperatureC", 20).put("pressureMPa", 5);
    JsonNode phase = ThermodynamicEngine.calculate(request, true).path("points").path(0).path("phases").path(0);
    assertTrue(phase.path("cpJMolK").isNull());
    assertEquals("OUTSIDE_PROPERTY_RANGE", phase.path("properties").path("cp").path("code").asText());
  }

  /**
   * A binary hydrocarbon envelope has physical sides, finite extrema and exactly preserved caller data.
   *
   * @throws IOException if the fixture cannot be read
   */
  @Test
  void envelopeUsesRealContinuationAndPreservesExplicitFeed() throws IOException {
    ObjectNode request = model();
    request.put("mixingRule", "Classic");
    ((ObjectNode) request.path("interactions").path(0)).put("mixingRule", "Classical").put("kijT", 0);
    String before = request.toString();
    JsonNode result = ThermodynamicEngine.envelope(request);
    assertTrue("ok".equals(result.path("status").asText()) || "partial".equals(result.path("status").asText()),
        result.toString());
    assertTrue(result.path("segments").size() >= 2, result.toString());
    double highestDew = -Double.MAX_VALUE, highestBubble = -Double.MAX_VALUE;
    double maxT = -Double.MAX_VALUE, maxP = 0;
    for (JsonNode segment : result.path("segments")) {
      for (JsonNode point : segment.path("points")) {
        double t = point.path("temperatureC").asDouble(), p = point.path("pressureMPa").asDouble();
        assertTrue(Double.isFinite(t) && t > -223.15 && Double.isFinite(p) && p > 0);
        maxT = Math.max(maxT, t);
        maxP = Math.max(maxP, p);
        if ("dew".equals(segment.path("type").asText())) {
          highestDew = Math.max(highestDew, t);
        } else {
          highestBubble = Math.max(highestBubble, t);
        }
      }
    }
    assertTrue(highestDew > highestBubble);
    assertEquals(maxT, result.path("cricondentherm").path("temperatureC").asDouble(), 0.1);
    assertEquals(maxP, result.path("cricondenbar").path("pressureMPa").asDouble(), 0.01);
    assertEquals(before, request.toString());
    assertEquals(2, result.path("feed").size());
    assertEquals(64, result.path("provenance").path("modelSha256").asText().length());
  }

  /**
   * A requested projection must disclose every removed polar component and retain its original feed.
   *
   * @throws IOException if the fixture cannot be read
   */
  @Test
  void polarProjectionIsExplicitAndRejectsIncompleteInteractionTables() throws IOException {
    ObjectNode request = model();
    ((ObjectNode) request.path("interactions").path(0)).put("mixingRule", "Classical").put("kijT", 0);
    ObjectNode water = ((ObjectNode) request.path("components").path(0)).deepCopy();
    water.put("name", "H2O").put("engineName", "water").put("templateName", "water").put("z", 0.2);
    ((ObjectNode) request.path("components").path(0)).put("z", 0.4);
    ((ObjectNode) request.path("components").path(1)).put("z", 0.4);
    ((ArrayNode) request.path("components")).add(water);
    JsonNode unavailable = ThermodynamicEngine.envelope(request);
    assertEquals("POLAR_FEED_REQUIRES_EXPLICIT_DRY_PROJECTION", unavailable.path("code").asText());
    assertEquals(0, unavailable.path("segments").size());
    request.putObject("envelope").put("hydrocarbonOnly", true);
    String before = request.toString();
    JsonNode result = ThermodynamicEngine.envelope(request);
    assertEquals("dry-hydrocarbon-projection", result.path("scope").asText());
    assertEquals(3, result.path("sourceFeed").size());
    assertEquals(2, result.path("feed").size());
    assertEquals(0.5, result.path("feed").path(0).path("moleFraction").asDouble(), 1e-12);
    assertEquals("H2O", result.path("omittedComponents").path(0).path("name").asText());
    assertEquals("unavailable", result.path("availability").path("waterBoundary").path("status").asText());
    assertEquals(before, request.toString());
    request.putArray("interactions");
    assertThrows(IllegalArgumentException.class, () -> ThermodynamicEngine.envelope(request));
  }

  /**
   * Unimplemented HV temperature binding is an explicit gap, never a frozen-interaction envelope.
   *
   * @throws IOException if the fixture cannot be read
   */
  @Test
  void envelopeDoesNotSilentlyFreezeHvTemperatureInteractions() throws IOException {
    JsonNode result = ThermodynamicEngine.envelope(model());
    assertEquals("unavailable", result.path("status").asText());
    assertEquals("HV_TEMPERATURE_KIJ_ENVELOPE_UNSUPPORTED", result.path("code").asText());
    assertFalse(result.path("segments").elements().hasNext());
  }

  /**
   * A Classic model must remain Classic even when an inactive pair retains source HV metadata.
   *
   * @throws IOException if the fixture cannot be read
   */
  @Test
  void envelopeRespectsTheSelectedGlobalMixingRule() throws IOException {
    ObjectNode request = model();
    request.put("mixingRule", "Classic");
    JsonNode result = ThermodynamicEngine.envelope(request);
    assertEquals("Classic", result.path("effectiveMixingRule").asText(), result.toString());
    assertTrue(result.path("segments").size() >= 2, result.toString());
  }

  /**
   * A temperature table must replay the same linear kij used by independent fixed-temperature flashes.
   *
   * @throws IOException if the fixture cannot be read
   */
  @Test
  void continuationTemperatureBindingMatchesFreshExplicitSystems() throws IOException {
    ObjectNode request = model();
    request.put("mixingRule", "Classic");
    ((ObjectNode) request.path("interactions").path(0)).put("mixingRule", "Classical");
    SystemInterface continuation = ExplicitFlash.system(request, 300, 50);
    ExplicitEnvelope.bindTemperatureInteractions(continuation, request);
    for (double t : new double[] {280, 350, 450}) {
      continuation.setTemperature(t);
      SystemInterface fixed = ExplicitFlash.system(request, t, 50);
      for (SystemInterface system : new SystemInterface[] {continuation, fixed}) {
        system.setNumberOfPhases(1);
        system.setPhaseType(0, PhaseType.GAS);
        system.setBeta(0, 1);
        system.init(3);
      }
      assertEquals(fixed.getPhase(0).getZ(), continuation.getPhase(0).getZ(), 1e-10);
      for (int i = 0; i < 2; i++) {
        assertEquals(fixed.getPhase(0).getComponent(i).getLogFugacityCoefficient(),
            continuation.getPhase(0).getComponent(i).getLogFugacityCoefficient(), 1e-10);
      }
    }
  }
}
