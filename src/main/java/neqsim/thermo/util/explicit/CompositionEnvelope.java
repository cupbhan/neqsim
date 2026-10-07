package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.util.Arrays;
import java.util.Locale;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoHydrocarbonPhaseEnvelopeSolver;

/**
 * Independently parameterized composition envelope using NeqSim database systems and physical branch tracing.
 *
 * @author NeqSim development team
 * @version 1.0
 */
final class CompositionEnvelope {
  /** Prevent construction of the stateless entry point. */
  private CompositionEnvelope() {
  }

  /**
   * Trace actual dry gas/liquid boundaries with explicit feed projection and incomplete-coverage diagnostics.
   *
   * @param input composition and independent model settings
   * @return physical segments, input mapping and independent parameter provenance
   */
  static ObjectNode calculate(JsonNode input) {
    ObjectNode response = ExplicitFlash.JSON.createObjectNode();
    response.putArray("segments");
    try {
      CompositionModel.validate(input);
      JsonNode settings = input.has("envelope") ? input.path("envelope") : ExplicitFlash.JSON.createObjectNode();
      CompositionModel.only(settings, "temperatureMinC", "temperatureMaxC", "pressureMinMPa", "pressureMaxMPa",
          "maximumIterations", "hydrocarbonOnly");
      double minT = setting(settings, "temperatureMinC", -120) + 273.15;
      double maxT = setting(settings, "temperatureMaxC", 250) + 273.15;
      double minP = setting(settings, "pressureMinMPa", 0.01) * 10;
      double maxP = setting(settings, "pressureMaxMPa", 40) * 10;
      CompositionModel.require(minT >= 50 && maxT <= 2000 && maxT > minT && minP > 0 && maxP <= 5000 && maxP > minP,
          "INVALID_ENVELOPE_RANGE", "Envelope bounds must be physical and increasing");
      CompositionModel.require(!settings.has("hydrocarbonOnly") || settings.path("hydrocarbonOnly").isBoolean(),
          "INVALID_INPUT", "hydrocarbonOnly must be boolean");
      CompositionModel.require(
          !settings.has("maximumIterations") || settings.path("maximumIterations").isIntegralNumber()
              && settings.path("maximumIterations").asLong() >= 100
              && settings.path("maximumIterations").asLong() <= 5000,
          "INVALID_INPUT", "maximumIterations must lie between 100 and 5000");
      boolean dry = settings.path("hydrocarbonOnly").asBoolean(false);
      response.put("scope", dry ? "dry-hydrocarbon-projection" : "full-feed");
      ObjectNode projected = project(input, response, dry);
      CompositionModel.require(projected.path("components").size() >= 2, "MIXTURE_REQUIRES_TWO_COMPONENTS",
          "Mixture envelope requires at least two retained components");
      CompositionModel model = new CompositionModel(projected, (minT + maxT) / 2, minP);
      for (int i = 0; i < model.bindings.size(); i++) {
        ((ObjectNode) model.bindings.get(i)).put("sourceIndex",
            response.path("feed").get(i).path("sourceIndex").asInt());
      }
      response.set("components", model.bindings);
      response.set("interactions", model.interactions());
      String mixingRule = "ELECTROLYTE-CPA".equals(input.path("eos").asText()) ? "CPA (10)" : "Classic";
      response.put("eos", input.path("eos").asText()).put("mixingRule", mixingRule)
          .put("effectiveMixingRule", mixingRule)
          .put("volumeCorrection", input.path("volumeCorrection").asBoolean(false));
      ObjectNode availability = response.putObject("availability");
      availability.putObject("waterBoundary").put("status", "unavailable").put("code",
          "WATER_BOUNDARY_NOT_IMPLEMENTED");
      availability.putObject("threePhaseBoundary").put("status", "unavailable").put("code",
          "THREE_PHASE_BOUNDARY_NOT_IMPLEMENTED");
      String version = ThermodynamicEngine.class.getPackage().getImplementationVersion();
      response.putObject("provenance").put("engine", "NeqSim")
          .put("engineVersion", version == null ? "development" : version)
          .put("parameterSource", "independent-neqsim-model").put("callerModelPreserved", true)
          .put("algorithm", "TwoHydrocarbonPhaseEnvelopeSolver");
      TwoHydrocarbonPhaseEnvelopeSolver.Result result = new TwoHydrocarbonPhaseEnvelopeSolver(model.fluid)
          .setTemperatureRange(minT, maxT).setPressureRange(minP, maxP)
          .setMaximumContinuationIterations(settings.path("maximumIterations").asInt(1200))
          .setStabilityAnalysisEnabled(true).setTopologyFallbackEnabled(false).solve();
      ObjectNode diagnostics = response.putObject("diagnostics");
      diagnostics.put("method", result.getMethod()).put("termination", result.getTerminationSummary())
          .put("envelopeClosed", result.isEnvelopeClosed()).put("criticalPointCount", result.getCriticalPointCount())
          .put("iterationLimitReached", result.isIterationLimitReached()).put("pointCount", result.getPointCount())
          .put("threePhaseStabilityPointCount", result.getThreePhaseStabilityPointCount());
      CompositionModel.require(result.hasFinitePhysicalValues() && result.getPointCount() > 0,
          "ENVELOPE_CONTINUATION_FAILED", "No finite physical envelope was produced");
      boolean dew = false, bubble = false;
      ArrayNode segments = (ArrayNode) response.path("segments");
      for (TwoHydrocarbonPhaseEnvelopeSolver.Segment segment : result.getSegments()) {
        String type = segment.getPhysicalType().name().toLowerCase(Locale.ROOT);
        dew |= "dew".equals(type);
        bubble |= "bubble".equals(type);
        ArrayNode points = segments.addObject().put("type", type).putArray("points");
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
          && result.getThreePhaseStabilityPointCount() == 0;
      for (String reason : result.getTerminationSummary().split("\\+")) {
        complete &= "PRESSURE_FLOOR".equals(reason) || "CONTINUATION_END".equals(reason);
      }
      if (complete) {
        complete = crossCheck(model, result, diagnostics);
      }
      response.put("ok", true).put("status", complete ? "ok" : "partial");
      if (!complete) {
        response.put("code", "INCOMPLETE_ENVELOPE");
      }
    } catch (Exception error) {
      response.setAll(CompositionFlash.failure(error));
      response.put("status", "unavailable");
      response.putArray("segments");
    }
    return response;
  }

  /**
   * Build an audited dry projection without using any physical parameters from a reference engine.
   *
   * @param input source composition
   * @param response feed audit destination
   * @param dry whether removal of water and inhibitors is explicitly allowed
   * @return private composition request with reindexed user interactions
   */
  private static ObjectNode project(JsonNode input, ObjectNode response, boolean dry) {
    ObjectNode copy = input.deepCopy();
    ArrayNode components = copy.putArray("components"), feed = response.putArray("feed");
    ArrayNode source = response.putArray("sourceFeed"), omitted = response.putArray("omittedComponents");
    int[] indices = new int[input.path("components").size()];
    Arrays.fill(indices, -1);
    double total = 0;
    for (int i = 0; i < indices.length; i++) {
      JsonNode component = input.path("components").get(i);
      double z = component.path("z").asDouble();
      String name = CompositionModel.resolve(component);
      boolean salt = "ELECTROLYTE-CPA".equals(input.path("eos").asText()) && ElectrolyteComposition.isSalt(component);
      boolean polar = salt || Arrays.asList("water", "methanol", "ethanol", "MEG", "DEG", "TEG").contains(name);
      source.addObject().put("sourceIndex", i).put("id", component.path("id").asText())
          .put("name", component.path("name").asText()).put("moleFraction", z);
      if (z < 1e-12 || polar) {
        CompositionModel.require(z < 1e-12 || dry, "POLAR_FEED_REQUIRES_EXPLICIT_DRY_PROJECTION",
            "Water or inhibitor removal requires hydrocarbonOnly=true");
        omitted.addObject().put("sourceIndex", i).put("id", component.path("id").asText())
            .put("name", component.path("name").asText()).put("moleFraction", z).put("reason",
                z < 1e-12 ? "numerical-zero" : salt ? "salt" : "water".equals(name) ? "water" : "polar-inhibitor");
      } else {
        indices[i] = components.size();
        components.add(component.deepCopy());
        feed.addObject().put("sourceIndex", i).put("id", component.path("id").asText())
            .put("name", component.path("name").asText()).put("moleFraction", z);
        total += z;
      }
    }
    for (int i = 0; i < components.size(); i++) {
      double z = components.get(i).path("z").asDouble() / total;
      ((ObjectNode) components.get(i)).put("z", z);
      ((ObjectNode) feed.get(i)).put("moleFraction", z);
    }
    ArrayNode pairs = copy.putArray("interactions");
    for (JsonNode pair : input.path("interactions")) {
      int i = indices[pair.path("i").asInt()], j = indices[pair.path("j").asInt()];
      if (i >= 0 && j >= 0) {
        ObjectNode mapped = pair.deepCopy();
        pairs.add(mapped.put("i", i).put("j", j));
      }
    }
    response.put("retainedSourceMoleFraction", total);
    return copy;
  }

  /**
   * Reject apparently complete envelopes when independent flashes reveal a missed hotter phase split.
   *
   * @param model independently built private system
   * @param result continuation result
   * @param diagnostics verification destination
   * @return true when both flashes pass equilibrium checks without a phase split
   */
  private static boolean crossCheck(CompositionModel model, TwoHydrocarbonPhaseEnvelopeSolver.Result result,
      ObjectNode diagnostics) {
    ArrayNode checks = diagnostics.putArray("aboveCricondenthermChecks");
    double[] maximum = result.getCricondentherm();
    boolean passed = true;
    for (double pressure : new double[] {maximum[1], result.getCricondenbar()[1] * 0.5}) {
      ObjectNode row = checks.addObject().put("temperatureC", maximum[0] - 273.15 + 2).put("pressureMPa",
          pressure / 10);
      try {
        SystemInterface clone = model.fluid.clone();
        clone.setTemperature(maximum[0] + 2);
        clone.setPressure(pressure);
        new ThermodynamicOperations(clone).TPflash();
        clone.init(1);
        CompositionFlash.audit(clone, model.input);
        row.put("status", "ok").put("phaseCount", clone.getNumberOfPhases());
        passed &= clone.getNumberOfPhases() == 1;
      } catch (Exception error) {
        row.put("status", "failed").put("error", error.toString());
        passed = false;
      }
    }
    diagnostics.put("maximumTemperatureCrossCheckPassed", passed);
    return passed;
  }

  /**
   * Serialize a finite characteristic point using degrees Celsius and MPa.
   *
   * @param response result destination
   * @param name characteristic point name
   * @param values temperature in kelvin and pressure in bara
   */
  private static void characteristic(ObjectNode response, String name, double[] values) {
    if (values.length >= 2 && Double.isFinite(values[0]) && Double.isFinite(values[1]) && values[0] > 0
        && values[1] > 0) {
      response.putObject(name).put("temperatureC", values[0] - 273.15).put("pressureMPa", values[1] / 10);
    } else {
      response.putNull(name);
    }
  }

  /**
   * Read a finite optional envelope setting.
   *
   * @param settings envelope configuration
   * @param name setting name
   * @param fallback default value
   * @return configured or default value
   */
  private static double setting(JsonNode settings, String name, double fallback) {
    return settings.has(name) ? CompositionModel.number(settings, name) : fallback;
  }
}
