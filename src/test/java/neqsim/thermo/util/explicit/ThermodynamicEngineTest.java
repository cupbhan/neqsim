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
import java.util.Iterator;
import org.junit.jupiter.api.Test;

/**
 * Standalone consumer, migration, failure-isolation and published-correlation contracts.
 *
 * @author NeqSim development team
 * @version 1.0
 */
public class ThermodynamicEngineTest {
  private static final ObjectMapper JSON = new ObjectMapper();

  /**
   * Read public synthetic test data without depending on the platform or simulator.
   *
   * @param name resource filename
   * @return parsed JSON
   * @throws IOException if reading fails
   */
  private JsonNode resource(String name) throws IOException {
    try (InputStream stream = getClass().getResourceAsStream("/neqsim/thermo/explicit/" + name)) {
      if (stream == null) {
        throw new IOException("Missing test resource " + name);
      }
      return JSON.readTree(stream);
    }
  }

  /**
   * Compare all numerical and availability fields against pre-migration synthetic calculations.
   *
   * @param expected original result
   * @param actual shared engine result
   * @param path diagnostic JSON path
   */
  private void compare(JsonNode expected, JsonNode actual, String path) {
    if (expected.isNumber()) {
      assertTrue(actual.isNumber(), path);
      assertEquals(expected.asDouble(), actual.asDouble(), 1e-10 + 1e-10 * Math.abs(expected.asDouble()), path);
    } else if (expected.isContainerNode()) {
      assertEquals(expected.size(), actual.size(), path);
      if (expected.isArray()) {
        for (int i = 0; i < expected.size(); i++) {
          compare(expected.get(i), actual.path(i), path + "/" + i);
        }
      } else {
        Iterator<String> fields = expected.fieldNames();
        while (fields.hasNext()) {
          String field = fields.next();
          compare(expected.get(field), actual.path(field), path + "/" + field);
        }
      }
    } else {
      assertEquals(expected, actual, path);
    }
  }

  /**
   * Verify both EOS families and mixing rules against the original implementation.
   *
   * @throws IOException if fixture or response parsing fails
   */
  @Test
  void preservesMigrationBaseline() throws IOException {
    JsonNode cases = resource("migration-baseline.json").path("cases");
    assertEquals(4, cases.size());
    for (JsonNode testCase : cases) {
      JsonNode request = testCase.path("input");
      String snapshot = request.toString();
      JsonNode result = JSON.readTree(ThermodynamicEngine.calculate(snapshot));
      compare(testCase.path("expected"), result, request.path("eos").asText());
      assertEquals(snapshot, request.toString(), "caller model must be unchanged");
      JsonNode repeat = ThermodynamicEngine.calculate(request, true);
      compare(result, repeat, "repeat");
    }
  }

  /**
   * Execute the documented external-consumer entry point using its shipped example request.
   *
   * @throws IOException if example parsing fails
   */
  @Test
  void documentationExampleRunsWithoutPlatform() throws IOException {
    String request = resource("example-request.json").toString();
    String response = ThermodynamicEngine.calculate(request);
    JsonNode result = JSON.readTree(response);
    assertEquals("ok", result.path("points").path(0).path("status").asText());
    JsonNode phase = result.path("points").path(0).path("phases").path(0);
    assertTrue(phase.path("densityKgM3").asDouble() > 0);
    assertTrue(phase.path("viscosityPaS").asDouble() > 0);
    assertTrue(phase.path("enthalpyJMol").isNumber());
  }

  /**
   * Invalid states must not suppress valid neighbors or expose unqualified thermal properties.
   *
   * @throws IOException if example parsing fails
   */
  @Test
  void isolatesFailuresAndPropertyLimits() throws IOException {
    ObjectNode request = (ObjectNode) resource("example-request.json");
    ArrayNode grid = (ArrayNode) request.path("grid");
    grid.addObject().put("temperatureC", 20).put("pressureMPa", 1);
    grid.addObject().put("temperatureC", 25).put("pressureMPa", -1);
    JsonNode result = ThermodynamicEngine.calculate(request, true).path("points");
    assertEquals("ok", result.path(0).path("status").asText());
    assertEquals("failed", result.path(3).path("status").asText());
    JsonNode unavailable = result.path(2).path("phases").path(0);
    assertTrue(unavailable.path("enthalpyJMol").isNull());
    assertEquals("OUTSIDE_PROPERTY_RANGE", unavailable.path("properties").path("enthalpy").path("code").asText());
    JsonNode onlyFlash = ThermodynamicEngine.calculate(request, false);
    assertFalse(onlyFlash.path("points").path(0).path("phases").path(0).has("enthalpyJMol"));
    request.put("eos", "typo");
    assertThrows(IllegalArgumentException.class, () -> ThermodynamicEngine.calculate(request, true));
    request.put("eos", "PengRobinson78Peneloux");
    ((ObjectNode) request.path("components").path(0)).put("z", -1);
    assertThrows(IllegalArgumentException.class, () -> ThermodynamicEngine.calculate(request, true));
  }

  /**
   * Exercise the database-component API used by the platform's standard interface.
   *
   * @throws IOException if request or response parsing fails
   */
  @Test
  void standardConsumerUsesSharedFlash() throws IOException {
    String request = "{\"eos\":\"PR\",\"temperatureC\":25,\"pressureMPa\":5,"
        + "\"components\":[{\"id\":\"methane\",\"nativeName\":\"methane\","
        + "\"moleFraction\":1,\"molarMass\":16.043}],\"interactions\":[]}";
    JsonNode result = JSON.readTree(ThermodynamicEngine.standardFlash(request));
    assertTrue(result.path("ok").asBoolean(), result.toString());
    assertEquals("gas", result.path("phases").path(0).path("phase").asText());
    assertTrue(result.path("phases").path(0).path("viscosity").asDouble() > 0);
  }

  /** Verify independent IAPWS R7-97 tables 5 and 35, and R12-08 table 4. */
  @Test
  void publishedWaterVerificationPoints() {
    assertEquals(.00100215168, 1 / PureWaterViscosity.densityKgM3(300, 3), 5e-12);
    assertEquals(.000971180894, 1 / PureWaterViscosity.densityKgM3(300, 80), 5e-13);
    assertEquals(.00120241800, 1 / PureWaterViscosity.densityKgM3(500, 3), 5e-12);
    assertEquals(.00353658941, PureWaterViscosity.saturationPressureMPa(300), 5e-12);
    assertEquals(2.63889776, PureWaterViscosity.saturationPressureMPa(500), 5e-9);
    assertEquals(12.3443146, PureWaterViscosity.saturationPressureMPa(600), 5e-8);
    double[][] samples = {{298.15, 998, 889.735100}, {298.15, 1200, 1437.649467}, {373.15, 1000, 307.883622},
        {433.15, 1, 14.538324}, {433.15, 1000, 217.685358}, {873.15, 1, 32.619287}, {873.15, 100, 35.802262},
        {873.15, 600, 77.430195}, {1173.15, 1, 44.217245}, {1173.15, 100, 47.640433}, {1173.15, 400, 64.154608}};
    for (double[] sample : samples) {
      assertEquals(sample[2], PureWaterViscosity.viscosityPaS(sample[0], sample[1]) * 1e6, 5e-7);
    }
    assertThrows(IllegalArgumentException.class, () -> PureWaterViscosity.densityKgM3(400, .1));
    assertThrows(IllegalArgumentException.class, () -> PureWaterViscosity.densityKgM3(300, Double.NaN));
  }

  /** Verify published methane pressure and viscosity points and the dilute-gas limit. */
  @Test
  void publishedMethaneVerificationPoints() {
    assertEquals(45.387,
        CspViscosity.pressure(10.23, CspViscosity.polynomial(190.555), CspViscosity.exponential(190.555)), .001);
    assertEquals(112.5e-7, CspViscosity.referenceViscosity(300, 1e5 / 101325), .15e-7);
    assertEquals(169.7e-7, CspViscosity.referenceViscosity(500, 1e5 / 101325), .15e-7);
    assertEquals(1, CspViscosity.density(300, 1e-5) * CspViscosity.R_ATM * 300 / 1e-5, 1e-6);
    assertThrows(IllegalArgumentException.class, () -> CspViscosity.density(300, -1));
  }
}
