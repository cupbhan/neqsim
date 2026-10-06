package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;

/**
 * Shared, stateless JSON boundary for explicitly parameterized thermodynamics.
 *
 * <p>
 * No platform files, reference results or simulator installation are used. The explicit profile is opt-in and does not
 * alter the default NeqSim system property methods. Requests and results use degrees Celsius, MPa, kg/m3, J/mol and Pa
 * s as indicated by their field names.
 *
 * @author NeqSim development team
 * @version 1.0
 */
public final class ThermodynamicEngine {
  private static final ObjectMapper JSON = new ObjectMapper();

  /** Prevent construction of the stateless entry point. */
  private ThermodynamicEngine() {
  }

  /**
   * Calculate an explicit-component grid with validated thermal and transport properties.
   *
   * @param request serialized explicit model and grid
   * @return serialized points, phase properties and availability metadata
   * @throws JsonProcessingException if the request or response cannot be serialized
   */
  public static String calculate(String request) throws JsonProcessingException {
    return JSON.writeValueAsString(calculate(JSON.readTree(request), true));
  }

  /**
   * Calculate a grid without reading or modifying the caller's model state.
   *
   * @param input explicit components, interactions, grid and optional property profile
   * @param thermalProperties whether to evaluate validated enthalpy and viscosity
   * @return points with per-point success or failure and per-property availability
   */
  public static ObjectNode calculate(JsonNode input, boolean thermalProperties) {
    validate(input);
    if (!input.path("grid").isArray() || input.path("grid").size() == 0) {
      throw new IllegalArgumentException("A nonempty grid is required");
    }
    ObjectNode response = JSON.createObjectNode();
    ArrayNode points = response.putArray("points");
    for (JsonNode point : input.path("grid")) {
      if (!point.isObject()) {
        throw new IllegalArgumentException("Each grid point must be an object");
      }
      try {
        requireState(point);
        ObjectNode result = ExplicitFlash.calculate(input, point);
        if (thermalProperties && "ok".equals(result.path("status").asText())) {
          for (JsonNode phase : result.path("phases")) {
            ExplicitPropertyCalculator.properties(input, point, (ObjectNode) phase);
          }
        }
        points.add(result);
      } catch (Exception error) {
        ObjectNode failed = point.deepCopy();
        points.add(failed.put("status", "failed").put("error", error.toString()));
      }
    }
    return response;
  }

  /**
   * Run fixed-composition diagnostics without performing phase equilibrium calculations.
   *
   * @param input explicit model with reference points and phase compositions
   * @return isolated phase diagnostics; failed reference points remain identified
   */
  public static ObjectNode diagnose(JsonNode input) {
    validate(input);
    if (!input.path("points").isArray()) {
      throw new IllegalArgumentException("Diagnostic points must be an array");
    }
    ObjectNode response = JSON.createObjectNode();
    ArrayNode points = response.putArray("points");
    for (JsonNode point : input.path("points")) {
      ObjectNode row = points.addObject();
      row.put("temperatureC", ExplicitFlash.number(point, "temperatureC"));
      row.put("pressureMPa", ExplicitFlash.number(point, "pressureMPa"));
      if (!"ok".equals(point.path("status").asText())) {
        row.put("status", "reference-failed");
        continue;
      }
      try {
        requireState(point);
        ArrayNode phases = row.putArray("phases");
        for (JsonNode phase : point.path("phases")) {
          phases.add(ExplicitProperties.phase(input, point, phase));
        }
        row.put("status", "ok");
      } catch (Exception error) {
        row.put("status", "failed").put("error", error.toString());
      }
    }
    return response;
  }

  /**
   * Run a dry database-component flash with constant symmetric interactions.
   *
   * @param request serialized database-component request
   * @return serialized phase results or explicit calculation failure
   * @throws JsonProcessingException if the request or response cannot be serialized
   */
  public static String standardFlash(String request) throws JsonProcessingException {
    return JSON.writeValueAsString(StandardFlash.calculate(JSON.readTree(request)));
  }

  /**
   * Check an explicit model before allocating EOS component arrays.
   *
   * @param input explicit model
   */
  private static void validate(JsonNode input) {
    if (input == null || !input.isObject()) {
      throw new IllegalArgumentException("Explicit model must be an object");
    }
    String eos = input.path("eos").asText();
    if (!("SrkPeneloux".equals(eos) || "SrkPenelouxTemperatureDependent".equals(eos)
        || "PengRobinsonPeneloux".equals(eos) || "PengRobinson78Peneloux".equals(eos))) {
      throw new IllegalArgumentException("Unsupported explicit EOS: " + eos);
    }
    String rule = input.path("mixingRule").asText("Classic");
    if (!("Classic".equals(rule) || "HV".equals(rule))) {
      throw new IllegalArgumentException("Unsupported explicit mixing rule: " + rule);
    }
    JsonNode components = input.path("components");
    if (!components.isArray() || components.size() == 0) {
      throw new IllegalArgumentException("Nonempty explicit components required");
    }
    double total = 0;
    for (JsonNode c : components) {
      double z = ExplicitFlash.number(c, "z");
      if (z < 0) {
        throw new IllegalArgumentException("Negative feed mole fraction");
      }
      total += z;
      for (String key : new String[] {"tcK", "pcBar", "mwKgMol", "omegaA", "omegaB", "tbK"}) {
        if (!(ExplicitFlash.number(c, key) > 0)) {
          throw new IllegalArgumentException("Expected positive " + key);
        }
      }
      for (String key : new String[] {"acentricFactor", "penelouxLmol", "penelouxTLmolK"}) {
        ExplicitFlash.number(c, key);
      }
    }
    if (Math.abs(total - 1) > 1e-8) {
      throw new IllegalArgumentException("Feed mole fractions must sum to one");
    }
    if (!input.path("interactions").isArray()) {
      throw new IllegalArgumentException("Explicit interactions must be an array");
    }
    for (JsonNode pair : input.path("interactions")) {
      for (String key : new String[] {"i", "j"}) {
        if (!pair.path(key).isIntegralNumber() || pair.path(key).asInt() < 0
            || pair.path(key).asInt() >= components.size()) {
          throw new IllegalArgumentException("Invalid interaction component index");
        }
      }
      ExplicitFlash.number(pair, "kij");
      ExplicitFlash.number(pair, "kijT");
    }
  }

  /**
   * Require finite positive absolute temperature and pressure.
   *
   * @param point temperature in degrees Celsius and pressure in MPa
   */
  private static void requireState(JsonNode point) {
    if (!(ExplicitFlash.number(point, "temperatureC") > -273.15) || !(ExplicitFlash.number(point, "pressureMPa") > 0)) {
      throw new IllegalArgumentException("Positive absolute temperature and pressure required");
    }
  }
}
