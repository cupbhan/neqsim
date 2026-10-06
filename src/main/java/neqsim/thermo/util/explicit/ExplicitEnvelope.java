package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Locale;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoHydrocarbonPhaseEnvelopeSolver;

/**
 * Physically classified phase-envelope boundary for explicit EOS models.
 *
 * <p>
 * Water, inhibitor and salt removal is opt-in and fully reported. This boundary never substitutes a database fluid,
 * invents points, or interprets an incomplete continuation as a complete envelope.
 *
 * @author NeqSim development team
 * @version 1.0
 */
final class ExplicitEnvelope {
  /** Prevent construction of the stateless boundary. */
  private ExplicitEnvelope() {
  }

  /**
   * Trace an explicitly parameterized dry two-phase envelope on a private system.
   *
   * @param input explicit model and optional envelope settings
   * @return physical contiguous segments and diagnostic availability
   */
  static ObjectNode calculate(JsonNode input) {
    JsonNode settings = input.path("envelope");
    if (!settings.isMissingNode() && !settings.isObject()) {
      throw new IllegalArgumentException("Envelope settings must be an object");
    }
    Iterator<String> fields = settings.fieldNames();
    while (fields.hasNext()) {
      String field = fields.next();
      if (!Arrays.asList("temperatureMinC", "temperatureMaxC", "pressureMinMPa", "pressureMaxMPa", "hydrocarbonOnly",
          "maximumIterations").contains(field)) {
        throw new IllegalArgumentException("Unsupported envelope setting: " + field);
      }
    }
    if (settings.has("hydrocarbonOnly") && !settings.path("hydrocarbonOnly").isBoolean()) {
      throw new IllegalArgumentException("hydrocarbonOnly must be a boolean");
    }
    double minT = setting(settings, "temperatureMinC", -120) + 273.15;
    double maxT = setting(settings, "temperatureMaxC", 250) + 273.15;
    double minP = setting(settings, "pressureMinMPa", 0.01) * 10;
    double maxP = setting(settings, "pressureMaxMPa", 40) * 10;
    if (minT < 50 || maxT > 2000 || maxT <= minT || minP <= 0 || maxP <= minP || maxP > 5000) {
      throw new IllegalArgumentException(
          "Envelope temperature and pressure ranges must be finite, positive and increasing");
    }
    if (settings.has("maximumIterations") && (!settings.path("maximumIterations").isIntegralNumber()
        || settings.path("maximumIterations").asInt() < 100 || settings.path("maximumIterations").asInt() > 5000)) {
      throw new IllegalArgumentException("Envelope maximumIterations must be an integer from 100 to 5000");
    }
    ObjectNode response = ExplicitFlash.JSON.createObjectNode();
    response.put("eos", input.path("eos").asText());
    response.put("mixingRule", input.path("mixingRule").asText("Classic"));
    response.put("scope", settings.path("hydrocarbonOnly").asBoolean() ? "dry-hydrocarbon-projection" : "full-feed");
    response.putArray("segments");
    ObjectNode availability = response.putObject("availability");
    availability.putObject("waterBoundary").put("status", "unavailable").put("code", "WATER_BOUNDARY_NOT_IMPLEMENTED");
    availability.putObject("threePhaseBoundary").put("status", "unavailable").put("code",
        "THREE_PHASE_BOUNDARY_NOT_IMPLEMENTED");
    ObjectNode provenance = response.putObject("provenance");
    provenance.put("engine", "NeqSim").put("modelSha256", fingerprint(input)).put("callerModelPreserved", true)
        .put("algorithm", "TwoHydrocarbonPhaseEnvelopeSolver");
    String version = ThermodynamicEngine.class.getPackage().getImplementationVersion();
    provenance.put("engineVersion", version == null ? "development" : version);
    ObjectNode projected = project(input, response, settings.path("hydrocarbonOnly").asBoolean());
    if (projected == null) {
      return response;
    }
    if (projected.path("components").size() < 2) {
      return unavailable(response, "MIXTURE_REQUIRES_TWO_COMPONENTS");
    }
    boolean hasHv = false, hasTemperatureKij = false;
    for (JsonNode pair : projected.path("interactions")) {
      hasHv |= "HV".equals(projected.path("mixingRule").asText())
          && "HuronVidal".equals(pair.path("mixingRule").asText());
      hasTemperatureKij |= ExplicitFlash.number(pair, "kijT") != 0;
    }
    if (hasHv && hasTemperatureKij) {
      return unavailable(response, "HV_TEMPERATURE_KIJ_ENVELOPE_UNSUPPORTED");
    }
    // A Huron-Vidal model with only Classic pairs is exactly its Classic submodel.
    projected.put("mixingRule", hasHv ? "HV" : "Classic");
    response.put("effectiveMixingRule", hasHv ? "HV" : "Classic");
    try {
      SystemInterface fluid = ExplicitFlash.system(projected, 0.5 * (minT + maxT), minP);
      if (hasTemperatureKij) {
        bindTemperatureInteractions(fluid, projected);
      }
      TwoHydrocarbonPhaseEnvelopeSolver.Result result = new TwoHydrocarbonPhaseEnvelopeSolver(fluid)
          .setTemperatureRange(minT, maxT).setPressureRange(minP, maxP)
          .setMaximumContinuationIterations(settings.path("maximumIterations").asInt(1200))
          .setStabilityAnalysisEnabled(true).setTopologyFallbackEnabled(false).solve();
      ObjectNode diagnostics = response.putObject("diagnostics");
      diagnostics.put("method", result.getMethod()).put("termination", result.getTerminationSummary())
          .put("envelopeClosed", result.isEnvelopeClosed()).put("criticalPointCount", result.getCriticalPointCount())
          .put("iterationLimitReached", result.isIterationLimitReached()).put("pointCount", result.getPointCount())
          .put("threePhaseStabilityPointCount", result.getThreePhaseStabilityPointCount());
      if (result.getFailureMessage() != null) {
        diagnostics.put("error", result.getFailureMessage());
      }
      if (!result.hasFinitePhysicalValues() || result.getPointCount() == 0) {
        return unavailable(response, "ENVELOPE_CONTINUATION_FAILED");
      }
      boolean dew = false, bubble = false;
      ArrayNode segments = (ArrayNode) response.path("segments");
      for (TwoHydrocarbonPhaseEnvelopeSolver.Segment segment : result.getSegments()) {
        String type = segment.getPhysicalType().name().toLowerCase(Locale.ROOT);
        dew |= "dew".equals(type);
        bubble |= "bubble".equals(type);
        ObjectNode row = segments.addObject().put("type", type);
        ArrayNode points = row.putArray("points");
        double[] temperatures = segment.getTemperaturesK(), pressures = segment.getPressuresBara();
        for (int i = 0; i < temperatures.length; i++) {
          points.addObject().put("temperatureC", temperatures[i] - 273.15).put("pressureMPa", pressures[i] / 10);
        }
      }
      characteristic(response, "criticalPoint", result.getCriticalPoint());
      characteristic(response, "cricondenbar", result.getCricondenbar());
      characteristic(response, "cricondentherm", result.getCricondentherm());
      boolean complete = dew && bubble && result.isEnvelopeClosed() && result.getCriticalPointCount() > 0
          && !result.isIterationLimitReached() && result.getFailureMessage() == null
          && result.getThreePhaseStabilityPointCount() == 0 && completeTermination(result.getTerminationSummary());
      if (complete) {
        complete = crossCheckMaximumTemperature(projected, result, diagnostics);
      }
      response.put("status", complete ? "ok" : "partial");
      if (!complete) {
        response.put("code", "INCOMPLETE_ENVELOPE");
      }
      return response;
    } catch (Exception error) {
      response.putObject("diagnostics").put("error", error.toString());
      return unavailable(response, "ENVELOPE_CALCULATION_FAILED");
    }
  }

  /**
   * Keep source composition, select a dry submodel when explicitly requested, and reindex all interactions.
   *
   * @param input caller model
   * @param response projection audit destination
   * @param dry whether polar components may be removed
   * @return private projected model, or null when projection is unavailable
   */
  private static ObjectNode project(JsonNode input, ObjectNode response, boolean dry) {
    ObjectNode copy = input.deepCopy();
    ArrayNode source = response.putArray("sourceFeed"), omitted = response.putArray("omittedComponents");
    ArrayNode feed = response.putArray("feed"), components = copy.putArray("components");
    int[] indices = new int[input.path("components").size()];
    Arrays.fill(indices, -1);
    double total = 0;
    boolean unrequestedRemoval = false;
    for (int i = 0; i < indices.length; i++) {
      JsonNode c = input.path("components").get(i);
      double z = ExplicitFlash.number(c, "z");
      source.addObject().put("sourceIndex", i).put("name", c.path("name").asText()).put("moleFraction", z);
      String excluded = exclusion(c);
      if (z < 1e-12 || excluded != null) {
        if (z >= 1e-12 && !dry) {
          unrequestedRemoval = true;
        }
        omitted.addObject().put("sourceIndex", i).put("name", c.path("name").asText()).put("moleFraction", z)
            .put("reason", z < 1e-12 ? "numerical-zero" : excluded);
      } else {
        indices[i] = components.size();
        components.add(c.deepCopy());
        feed.addObject().put("sourceIndex", i).put("name", c.path("name").asText()).put("moleFraction", z);
        total += z;
      }
    }
    if (unrequestedRemoval) {
      unavailable(response, "POLAR_FEED_REQUIRES_EXPLICIT_DRY_PROJECTION");
      return null;
    }
    for (int i = 0; i < components.size(); i++) {
      double z = components.get(i).path("z").asDouble() / total;
      ((ObjectNode) components.get(i)).put("z", z);
      ((ObjectNode) feed.get(i)).put("moleFraction", z);
    }
    ArrayNode pairs = copy.putArray("interactions");
    boolean[][] present = new boolean[components.size()][components.size()];
    for (JsonNode pair : input.path("interactions")) {
      int i = indices[pair.path("i").asInt()], j = indices[pair.path("j").asInt()];
      if (i < 0 || j < 0) {
        continue;
      }
      if (i == j || present[i][j]) {
        throw new IllegalArgumentException("Duplicate or diagonal envelope interaction");
      }
      present[i][j] = present[j][i] = true;
      ObjectNode mapped = pair.deepCopy();
      mapped.put("i", i).put("j", j);
      pairs.add(mapped);
    }
    if (pairs.size() != components.size() * (components.size() - 1) / 2) {
      throw new IllegalArgumentException("Envelope model requires every explicit component pair");
    }
    response.put("retainedSourceMoleFraction", total);
    return copy;
  }

  /**
   * Identify components outside the dry hydrocarbon boundary scope using explicit identities.
   *
   * @param component explicit component
   * @return exclusion reason, or null for a supported nonpolar component
   */
  private static String exclusion(JsonNode component) {
    if (component.has("saltIonNumber") || "ion".equalsIgnoreCase(component.path("componentType").asText())) {
      return "salt-or-ion";
    }
    for (String field : new String[] {"name", "engineName", "templateName"}) {
      String name = component.path(field).asText().toLowerCase(Locale.ROOT);
      if ("h2o".equals(name) || "water".equals(name)) {
        return "water";
      }
      if ("meoh".equals(name) || "methanol".equals(name) || "etoh".equals(name) || "ethanol".equals(name)
          || "meg".equals(name) || "deg".equals(name) || "teg".equals(name) || name.contains("glycol")) {
        return "polar-inhibitor";
      }
    }
    return null;
  }

  /**
   * Bind the same linear temperature-dependent Classic interactions across the full physical solver domain.
   *
   * @param fluid private explicit system
   * @param input reindexed dry model
   */
  static void bindTemperatureInteractions(SystemInterface fluid, JsonNode input) {
    fluid.setMixingRule(13);
    for (int ph = 0; ph < fluid.getMaxNumberOfPhases(); ph++) {
      if (fluid.getPhase(ph) == null) {
        continue;
      }
      EosMixingRulesInterface rule = ((PhaseEosInterface) fluid.getPhase(ph)).getEosMixingRule();
      for (JsonNode pair : input.path("interactions")) {
        double k = ExplicitFlash.number(pair, "kij"), slope = ExplicitFlash.number(pair, "kijT");
        rule.setBinaryInteractionParameterTemperatureTable(pair.path("i").asInt(), pair.path("j").asInt(),
            new double[] {0, 10000}, new double[] {k - slope * 288.15, k + slope * (10000 - 288.15)});
      }
    }
  }

  /**
   * Detect a truncated dew side by checking for two phases above its reported maximum temperature.
   *
   * @param input private explicit model
   * @param result continuation result
   * @param diagnostics verification audit destination
   * @return true only when both independent flashes converge without a phase split
   */
  private static boolean crossCheckMaximumTemperature(JsonNode input, TwoHydrocarbonPhaseEnvelopeSolver.Result result,
      ObjectNode diagnostics) {
    double[] maximum = result.getCricondentherm();
    ArrayNode checks = diagnostics.putArray("aboveCricondenthermChecks");
    boolean valid = true;
    for (double pressure : new double[] {maximum[1], result.getCricondenbar()[1] * 0.5}) {
      ObjectNode point = ExplicitFlash.JSON.createObjectNode().put("temperatureC", maximum[0] - 273.15 + 2)
          .put("pressureMPa", pressure / 10);
      try {
        ObjectNode checked = ExplicitFlash.calculate(input, point);
        int count = checked.path("phases").size();
        checks.add(point.put("phaseCount", count).put("status", "ok"));
        valid &= count == 1;
      } catch (Exception error) {
        checks.add(point.put("status", "failed").put("error", error.toString()));
        valid = false;
      }
    }
    diagnostics.put("maximumTemperatureCrossCheckPassed", valid);
    return valid;
  }

  /**
   * Check that continuation ended without a known truncation condition.
   *
   * @param summary solver termination codes
   * @return true only for normal end conditions
   */
  private static boolean completeTermination(String summary) {
    for (String reason : summary.split("\\+")) {
      if (!("PRESSURE_FLOOR".equals(reason) || "CONTINUATION_END".equals(reason))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Serialize a finite physical characteristic point with output units.
   *
   * @param response response object
   * @param key characteristic point name
   * @param point temperature in kelvin and pressure in bara
   */
  private static void characteristic(ObjectNode response, String key, double[] point) {
    if (point.length >= 2 && Double.isFinite(point[0]) && Double.isFinite(point[1]) && point[0] > 0 && point[1] > 0) {
      response.putObject(key).put("temperatureC", point[0] - 273.15).put("pressureMPa", point[1] / 10);
    } else {
      response.putNull(key);
    }
  }

  /**
   * Read an optional finite scalar setting.
   *
   * @param settings envelope settings
   * @param name setting name
   * @param fallback default value
   * @return requested value or default
   */
  private static double setting(JsonNode settings, String name, double fallback) {
    return settings.has(name) ? ExplicitFlash.number(settings, name) : fallback;
  }

  /**
   * Mark unavailable physical coverage without synthesizing a result.
   *
   * @param response response object
   * @param code machine-readable reason
   * @return same response with unavailable status
   */
  private static ObjectNode unavailable(ObjectNode response, String code) {
    return response.put("status", "unavailable").put("code", code);
  }

  /**
   * Hash the exact explicit JSON supplied by the caller for provenance.
   *
   * @param input caller model
   * @return SHA-256 hexadecimal digest
   */
  private static String fingerprint(JsonNode input) {
    try {
      byte[] bytes = MessageDigest.getInstance("SHA-256").digest(input.toString().getBytes(StandardCharsets.UTF_8));
      StringBuilder value = new StringBuilder();
      for (byte b : bytes) {
        value.append(String.format(Locale.ROOT, "%02x", b & 0xff));
      }
      return value.toString();
    } catch (NoSuchAlgorithmException error) {
      throw new IllegalStateException("SHA-256 unavailable", error);
    }
  }
}
