package neqsim.mcp.runners;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonPrimitive;
import neqsim.thermo.component.Component;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.mixingrule.HVMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.system.SystemInterface;

/** Deterministic identity of the live NeqSim phase-equilibrium model. */
final class HydrocarbonWaterModelIdentity {
  private HydrocarbonWaterModelIdentity() {
  }

  static JsonObject describe(SystemInterface system, String equationOfState, String mixingRule,
      JsonObject sourceDefinition) {
    if (!(system.getPhase(0) instanceof PhaseEosInterface)) {
      throw new IllegalArgumentException("A cubic-EOS phase is required for model fingerprinting");
    }
    JsonArray composition = composition(system);
    JsonArray componentParameters = componentParameters(system);
    JsonArray interactionParameters = interactionParameters(system);

    String compositionFingerprint = fingerprint(contract("neqsim-composition-v1", composition));
    String componentParameterFingerprint = fingerprint(contract("neqsim-component-parameters-v1", componentParameters));
    String interactionParameterFingerprint = fingerprint(
        contract("neqsim-interaction-parameters-v1", interactionParameters));
    JsonObject full = new JsonObject();
    full.addProperty("schema", "neqsim-live-phase-equilibrium-parameter-identity-v1");
    full.addProperty("equationOfState", equationOfState);
    full.addProperty("mixingRule", mixingRule);
    full.addProperty("componentParameterFingerprint", componentParameterFingerprint);
    full.addProperty("interactionParameterFingerprint", interactionParameterFingerprint);

    JsonObject output = new JsonObject();
    output.addProperty("status", "resolved");
    output.addProperty("schema", "neqsim-live-model-identity-v1");
    output.addProperty("equationOfState", equationOfState);
    output.addProperty("mixingRule", mixingRule);
    output.addProperty("componentCount", system.getNumberOfComponents());
    output.addProperty("compositionFingerprint", compositionFingerprint);
    output.addProperty("compositionFingerprintStatus", "resolved-live-runtime-composition");
    output.addProperty("componentParameterFingerprint", componentParameterFingerprint);
    output.addProperty("interactionParameterFingerprint", interactionParameterFingerprint);
    output.addProperty("fullParameterFingerprint", fingerprint(full));
    output.addProperty("fullParameterFingerprintStatus", "resolved");
    output.add("composition", composition);
    output.add("componentParameters", componentParameters);
    output.add("interactionParameters", interactionParameters);
    if (sourceDefinition != null) {
      output.addProperty("sourceDefinitionFingerprint", fingerprint(canonicalize(sourceDefinition)));
    }
    return output;
  }

  private static JsonObject contract(String schema, JsonArray rows) {
    JsonObject output = new JsonObject();
    output.addProperty("schema", schema);
    output.add("rows", rows);
    return output;
  }

  private static JsonArray composition(SystemInterface system) {
    List<JsonObject> rows = new ArrayList<JsonObject>();
    for (int index = 0; index < system.getNumberOfComponents(); index++) {
      ComponentInterface component = system.getPhase(0).getComponent(index);
      JsonObject row = new JsonObject();
      row.addProperty("component", canonicalComponentName(component.getComponentName()));
      row.addProperty("moleFraction", finite(component.getz(), "component mole fraction"));
      rows.add(row);
    }
    rows.sort(Comparator.comparing(row -> row.get("component").getAsString()));
    JsonArray output = new JsonArray();
    rows.forEach(output::add);
    return output;
  }

  private static JsonArray componentParameters(SystemInterface system) {
    List<JsonObject> rows = new ArrayList<JsonObject>();
    for (int index = 0; index < system.getNumberOfComponents(); index++) {
      ComponentInterface component = system.getPhase(0).getComponent(index);
      JsonObject row = new JsonObject();
      row.addProperty("component", canonicalComponentName(component.getComponentName()));
      row.addProperty("runtimeType", component.getClass().getName());
      row.addProperty("criticalTemperatureK", finite(component.getTC(), "critical temperature"));
      row.addProperty("criticalPressureBara", finite(component.getPC(), "critical pressure"));
      row.addProperty("acentricFactor", finite(component.getAcentricFactor(), "acentric factor"));
      row.addProperty("molarMassKgMol", finite(component.getMolarMass(), "molar mass"));
      row.addProperty("volumeCorrectionConst", finite(component.getVolumeCorrectionConst(), "volume correction"));
      row.addProperty("attractiveTermNumber", component.getAttractiveTermNumber());
      row.add("mathiasCopemanSrk", numbers(component.getMatiascopemanParams()));
      if (component instanceof Component) {
        Component concrete = (Component) component;
        row.add("mathiasCopemanPr", numbers(concrete.getMatiascopemanParamsPR()));
      }
      rows.add(row);
    }
    rows.sort(Comparator.comparing(row -> row.get("component").getAsString()));
    JsonArray output = new JsonArray();
    rows.forEach(output::add);
    return output;
  }

  private static JsonArray interactionParameters(SystemInterface system) {
    PhaseEosInterface phase = (PhaseEosInterface) system.getPhase(0);
    EosMixingRulesInterface mixing = phase.getEosMixingRule();
    JsonArray output = new JsonArray();
    for (int first = 0; first < system.getNumberOfComponents(); first++) {
      for (int second = 0; second < system.getNumberOfComponents(); second++) {
        JsonObject row = new JsonObject();
        row.addProperty("first", canonicalComponentName(system.getPhase(0).getComponent(first).getComponentName()));
        row.addProperty("second", canonicalComponentName(system.getPhase(0).getComponent(second).getComponentName()));
        row.addProperty("kij", finite(mixing.getBinaryInteractionParameter(first, second), "binary interaction"));
        row.addProperty("kijT1",
            finite(mixing.getBinaryInteractionParameterT1(first, second), "temperature binary interaction"));
        if (mixing instanceof HVMixingRulesInterface) {
          HVMixingRulesInterface huronVidal = (HVMixingRulesInterface) mixing;
          row.addProperty("hvDij", finite(huronVidal.getHVDijParameter(first, second), "HV Dij"));
          row.addProperty("hvDijT", finite(huronVidal.getHVDijTParameter(first, second), "HV DijT"));
          row.addProperty("hvAlpha", finite(huronVidal.getHValphaParameter(first, second), "HV alpha"));
        }
        output.add(row);
      }
    }
    return output;
  }

  private static JsonArray numbers(double[] values) {
    JsonArray output = new JsonArray();
    if (values != null) {
      for (double value : values) {
        output.add(finite(value, "component alpha parameter"));
      }
    }
    return output;
  }

  private static double finite(double value, String label) {
    if (!Double.isFinite(value)) {
      throw new IllegalArgumentException(label + " is not finite");
    }
    return value;
  }

  private static String canonicalComponentName(String name) {
    String normalized = String.valueOf(name).trim();
    if (normalized.endsWith("_PC")) {
      normalized = normalized.substring(0, normalized.length() - 3);
    }
    Map<String, String> aliases = Map.ofEntries(Map.entry("water", "H2O"), Map.entry("nitrogen", "N2"),
        Map.entry("methane", "C1"), Map.entry("ethane", "C2"), Map.entry("propane", "C3"), Map.entry("i-butane", "iC4"),
        Map.entry("n-butane", "nC4"), Map.entry("i-pentane", "iC5"), Map.entry("n-pentane", "nC5"),
        Map.entry("n-hexane", "C6"));
    return aliases.getOrDefault(normalized.toLowerCase(Locale.ROOT), normalized);
  }

  private static JsonElement canonicalize(JsonElement element) {
    if (element == null || element.isJsonNull() || element.isJsonPrimitive()) {
      return element == null ? new JsonPrimitive("null") : element.deepCopy();
    }
    if (element.isJsonArray()) {
      JsonArray output = new JsonArray();
      for (JsonElement child : element.getAsJsonArray()) {
        output.add(canonicalize(child));
      }
      return output;
    }
    JsonObject output = new JsonObject();
    element.getAsJsonObject().entrySet().stream().sorted(Map.Entry.comparingByKey())
        .forEach(entry -> output.add(entry.getKey(), canonicalize(entry.getValue())));
    return output;
  }

  private static String fingerprint(JsonElement value) {
    try {
      MessageDigest digest = MessageDigest.getInstance("SHA-256");
      byte[] hash = digest.digest(value.toString().getBytes(StandardCharsets.UTF_8));
      StringBuilder output = new StringBuilder(hash.length * 2);
      for (byte item : hash) {
        output.append(String.format(Locale.ROOT, "%02x", item & 0xff));
      }
      return output.toString();
    } catch (NoSuchAlgorithmException error) {
      throw new IllegalStateException("SHA-256 is unavailable", error);
    }
  }
}
