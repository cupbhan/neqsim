package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.util.database.NeqSimDataBase;

/**
 * Independent NeqSim database and TBP characterization boundary accepting composition only.
 *
 * @author NeqSim development team
 * @version 1.0
 */
final class CompositionModel {
  /** Common composition names mapped to NeqSim catalog identities. */
  private static final Map<String, String> ALIASES = aliases();
  /** Private independently constructed thermodynamic system. */
  final SystemInterface fluid;
  /** Actual modeled component bindings with original source indices. */
  final ArrayNode bindings;
  /** Exactly zero declared fractions omitted without estimating parameters. */
  final ArrayNode omitted;
  /** Private active composition with remapped local interactions. */
  final JsonNode input;

  /**
   * Build a private system without accepting external thermodynamic or transport parameters.
   *
   * @param input composition request
   * @param temperatureK initialization temperature in kelvin
   * @param pressureBar initialization pressure in bara
   */
  CompositionModel(JsonNode input, double temperatureK, double pressureBar) {
    validate(input);
    omitted = ExplicitFlash.JSON.createArrayNode();
    this.input = activeInput(input, omitted);
    fluid = "PR".equals(input.path("eos").asText()) ? new SystemPrEos(temperatureK, pressureBar)
        : new SystemSrkEos(temperatureK, pressureBar);
    fluid.setUseIdealGasEnthalpyOfFormation(false);
    bindings = ExplicitFlash.JSON.createArrayNode();
    Set<String> nativeNames = new HashSet<String>();
    boolean hasTbp = false;
    int index = 0;
    for (JsonNode component : this.input.path("components")) {
      String name = resolve(component);
      double z = number(component, "z");
      boolean tbp = component.has("tbp");
      if (name != null && tbp) {
        throw error("PURE_COMPONENT_TBP_CONFLICT",
            "Pure database components cannot receive TBP parameters: " + component.path("name").asText());
      }
      if (name == null && !tbp) {
        throw error("UNKNOWN_COMPONENT_REQUIRES_TBP",
            "Unknown component requires user TBP molar mass and relative density: " + component.path("name").asText());
      }
      if (tbp) {
        JsonNode values = component.path("tbp");
        only(values, "molarMassKgMol", "densityRelative");
        double mass = number(values, "molarMassKgMol"), density = number(values, "densityRelative");
        require(mass >= 0.05 && mass <= 2 && density >= 0.3 && density <= 1.5, "INVALID_TBP_INPUT",
            "TBP requires molarMassKgMol in [0.05, 2] and densityRelative in [0.3, 1.5]");
        fluid.addTBPfraction("user_tbp_" + index, z, mass, density);
        name = fluid.getPhase(0).getComponent(index).getComponentName();
        hasTbp = true;
      } else {
        require(nativeNames.add(name), "DUPLICATE_COMPONENT_BINDING", "Duplicate database component: " + name);
        fluid.addComponent(name, z);
      }
      ObjectNode binding = bindings.addObject().put("id", component.path("id").asText())
          .put("name", component.path("name").asText()).put("sourceIndex", sourceIndex(input, component))
          .put("nativeName", name)
          .put("parameterSource", tbp ? "neqsim-tbp-characterization" : "neqsim-component-database")
          .put("moleFraction", z);
      if (tbp) {
        binding.set("tbp", component.path("tbp").deepCopy());
        binding.put("characterizationModel", fluid.getCharacterization().getTBPModel().getName());
      }
      index++;
    }
    if (hasTbp) {
      fluid.createDatabase(true);
    }
    fluid.setMixingRule(2);
    fluid.useVolumeCorrection(input.path("volumeCorrection").asBoolean(false));
    fluid.setMultiPhaseCheck(true);
    fluid.init(0);
    fluid.setUseIdealGasEnthalpyOfFormation(false);
    for (int phase = 0; phase < fluid.getMaxNumberOfPhases(); phase++) {
      if (fluid.getPhase(phase) == null) {
        continue;
      }
      EosMixingRulesInterface mixing = ((PhaseEosInterface) fluid.getPhase(phase)).getEosMixingRule();
      for (JsonNode pair : this.input.path("interactions")) {
        mixing.setBinaryInteractionParameter(pair.path("i").asInt(), pair.path("j").asInt(), number(pair, "value"));
      }
    }
    for (int i = 0; i < bindings.size(); i++) {
      ComponentInterface component = fluid.getPhase(0).getComponent(i);
      require(!component.isIsIon(), "ELECTROLYTE_MODEL_REQUIRED", "Ionic components require an electrolyte model");
      require(
          Double.isFinite(component.getTC()) && component.getTC() > 0 && Double.isFinite(component.getPC())
              && component.getPC() > 0 && Double.isFinite(component.getMolarMass()) && component.getMolarMass() > 0,
          "COMPONENT_CHARACTERIZATION_FAILED", "Invalid NeqSim component characterization");
      ObjectNode binding = (ObjectNode) bindings.get(i);
      binding.put("engineMolarMass", component.getMolarMass() * 1000).put("molarMassKgMol", component.getMolarMass())
          .put("criticalTemperatureK", component.getTC()).put("criticalPressureBar", component.getPC())
          .put("acentricFactor", component.getAcentricFactor());
    }
  }

  /**
   * Validate a composition request and reject external parameter payloads rather than silently ignoring them.
   *
   * @param input caller request
   */
  static void validate(JsonNode input) {
    only(input, "eos", "volumeCorrection", "temperatureC", "pressureMPa", "components", "interactions", "envelope");
    require("PR".equals(input.path("eos").asText()) || "SRK".equals(input.path("eos").asText()), "UNSUPPORTED_EOS",
        "Composition calculation requires PR or SRK");
    require(!input.has("volumeCorrection") || input.path("volumeCorrection").isBoolean(), "INVALID_INPUT",
        "volumeCorrection must be boolean");
    JsonNode components = input.path("components");
    require(components.isArray() && components.size() > 0 && components.size() <= 100, "INVALID_INPUT",
        "Between one and 100 components are required");
    Set<String> ids = new HashSet<String>();
    double total = 0;
    for (JsonNode component : components) {
      only(component, "id", "name", "formula", "kind", "z", "tbp");
      require(
          !component.has("kind") || component.path("kind").isTextual()
              && Arrays.asList("pure", "pseudo", "plus", "custom").contains(component.path("kind").asText()),
          "INVALID_INPUT", "Component kind must be pure, pseudo, plus or custom");
      require(
          component.path("id").isTextual() && !component.path("id").asText().trim().isEmpty()
              && ids.add(component.path("id").asText()),
          "INVALID_INPUT", "Component IDs must be unique nonempty strings");
      require(component.path("name").isTextual() && !component.path("name").asText().trim().isEmpty(), "INVALID_INPUT",
          "Component names must be nonempty strings");
      require(!component.has("formula") || component.path("formula").isTextual(), "INVALID_INPUT",
          "formula must be text");
      double z = number(component, "z");
      require(z >= 0 && z <= 1, "INVALID_COMPOSITION", "Mole fractions must lie between zero and one");
      total += z;
    }
    require(Math.abs(total - 1) < 1e-8, "INVALID_COMPOSITION", "Mole fractions must sum to one");
    require(!input.has("interactions") || input.path("interactions").isArray(), "INVALID_INPUT",
        "interactions must be an array");
    Set<String> pairs = new HashSet<String>();
    for (JsonNode pair : input.path("interactions")) {
      only(pair, "i", "j", "value");
      for (String key : new String[] {"i", "j"}) {
        require(
            pair.path(key).isIntegralNumber() && pair.path(key).asLong() >= 0
                && pair.path(key).asLong() < components.size(),
            "INVALID_INTERACTION", "Interaction index is outside the component list");
      }
      int i = pair.path("i").asInt(), j = pair.path("j").asInt();
      require(i != j && pairs.add(Math.min(i, j) + ":" + Math.max(i, j)), "INVALID_INTERACTION",
          "Duplicate or diagonal interaction");
      require(Math.abs(number(pair, "value")) <= 1, "INVALID_INTERACTION", "Constant user BIP must lie in [-1, 1]");
    }
  }

  /**
   * Resolve a composition identity using known aliases and NeqSim's own component catalog.
   *
   * @param component composition row
   * @return NeqSim database name, or null when user characterization is required
   */
  static String resolve(JsonNode component) {
    if (isFraction(component)) {
      return null;
    }
    String byName = databaseName(component.path("name").asText()), byId = databaseName(component.path("id").asText());
    require(byName == null || byId == null || byName.equals(byId), "AMBIGUOUS_COMPONENT_IDENTITY",
        "Component id and name resolve to different database components");
    if (byName != null || byId != null) {
      return byName == null ? byId : byName;
    }
    String formula = component.path("formula").asText().toUpperCase(Locale.ROOT);
    if (Arrays.asList("N2", "CO2", "H2O", "CH4", "C2H6", "C3H8", "H2S", "H2", "O2").contains(formula)) {
      return databaseName(formula);
    }
    return null;
  }

  /**
   * Distinguish declared petroleum fractions from pure compounds with similar names.
   *
   * @param component composition row
   * @return whether the row is explicitly a pseudo or plus fraction
   */
  private static boolean isFraction(JsonNode component) {
    return "pseudo".equals(component.path("kind").asText()) || "plus".equals(component.path("kind").asText());
  }

  /**
   * Preserve original row indices even when an exactly zero fraction needs no characterization.
   *
   * @param original original request
   * @param component retained row
   * @return original component index
   */
  private static int sourceIndex(JsonNode original, JsonNode component) {
    for (int i = 0; i < original.path("components").size(); i++) {
      if (original.path("components").get(i).path("id").equals(component.path("id"))) {
        return i;
      }
    }
    throw error("INVALID_INPUT", "Component source index not found");
  }

  /**
   * Omit exactly zero declared petroleum fractions without inventing their physical parameters.
   *
   * @param original original composition and local model settings
   * @param omitted omitted-row audit destination
   * @return private active request with remapped local BIPs
   */
  private static JsonNode activeInput(JsonNode original, ArrayNode omitted) {
    ObjectNode copy = original.deepCopy();
    ArrayNode components = copy.putArray("components");
    int[] indices = new int[original.path("components").size()];
    Arrays.fill(indices, -1);
    for (int i = 0; i < indices.length; i++) {
      JsonNode row = original.path("components").get(i);
      if (row.path("z").asDouble() == 0 && isFraction(row)) {
        omitted.addObject().put("sourceIndex", i).put("id", row.path("id").asText())
            .put("name", row.path("name").asText()).put("moleFraction", 0).put("reason", "zero-fraction");
      } else {
        indices[i] = components.size();
        components.add(row.deepCopy());
      }
    }
    ArrayNode pairs = copy.putArray("interactions");
    for (JsonNode pair : original.path("interactions")) {
      int i = indices[pair.path("i").asInt()], j = indices[pair.path("j").asInt()];
      if (i >= 0 && j >= 0) {
        ObjectNode mapped = pair.deepCopy();
        pairs.add(mapped.put("i", i).put("j", j));
      }
    }
    return copy;
  }

  /**
   * Resolve one name or alias without estimating physical parameters.
   *
   * @param text user name or identifier
   * @return catalog name, or null
   */
  private static String databaseName(String text) {
    String trimmed = text.trim();
    String name = ALIASES.get(trimmed.toLowerCase(Locale.ROOT).replace(" ", ""));
    if (name != null) {
      return name;
    }
    if (trimmed.isEmpty() || trimmed.length() > 100) {
      return null;
    }
    String aliased = ComponentInterface.getComponentNameFromAlias(trimmed);
    return NeqSimDataBase.hasComponent(aliased) ? aliased : null;
  }

  /**
   * Record the applied independent interaction table, distinguishing user overrides from library defaults.
   *
   * @return complete interaction table
   */
  ArrayNode interactions() {
    ArrayNode pairs = ExplicitFlash.JSON.createArrayNode();
    EosMixingRulesInterface mixing = ((PhaseEosInterface) fluid.getPhase(0)).getEosMixingRule();
    for (int i = 0; i < bindings.size(); i++) {
      for (int j = 0; j < i; j++) {
        boolean overridden = false;
        for (JsonNode pair : input.path("interactions")) {
          overridden |= (pair.path("i").asInt() == i && pair.path("j").asInt() == j)
              || (pair.path("i").asInt() == j && pair.path("j").asInt() == i);
        }
        pairs.addObject().put("i", i).put("j", j).put("value", mixing.getBinaryInteractionParameter(i, j))
            .put("sourceI", bindings.get(i).path("sourceIndex").asInt())
            .put("sourceJ", bindings.get(j).path("sourceIndex").asInt())
            .put("source", overridden ? "user-constant-override" : "neqsim-mixing-rule");
      }
    }
    return pairs;
  }

  /**
   * Read a required finite scalar.
   *
   * @param node containing object
   * @param key field name
   * @return finite numeric value
   */
  static double number(JsonNode node, String key) {
    require(node.path(key).isNumber() && Double.isFinite(node.path(key).asDouble()), "INVALID_INPUT",
        "Finite " + key + " required");
    return node.path(key).asDouble();
  }

  /**
   * Reject undeclared fields on a request object.
   *
   * @param object request object
   * @param names allowed keys
   */
  static void only(JsonNode object, String... names) {
    require(object != null && object.isObject(), "INVALID_INPUT", "An object is required");
    Iterator<String> fields = object.fieldNames();
    while (fields.hasNext()) {
      String field = fields.next();
      require(Arrays.asList(names).contains(field), "UNSUPPORTED_PARAMETER",
          "Composition API does not accept parameter: " + field);
    }
  }

  /**
   * Enforce an input or physical invariant with a machine-readable reason.
   *
   * @param condition required condition
   * @param code failure code
   * @param message explanatory text
   */
  static void require(boolean condition, String code, String message) {
    if (!condition) {
      throw error(code, message);
    }
  }

  /**
   * Build an exception retaining the public failure reason.
   *
   * @param code failure code
   * @param message explanatory text
   * @return input exception
   */
  static IllegalArgumentException error(String code, String message) {
    return new IllegalArgumentException(code + ": " + message);
  }

  /**
   * Map common composition labels to independently parameterized catalog components.
   *
   * @return immutable-by-convention alias dictionary
   */
  private static Map<String, String> aliases() {
    Map<String, String> names = new HashMap<String, String>();
    String[][] groups = {{"nitrogen", "n2", "nitrogen"}, {"CO2", "co2", "carbondioxide"}, {"water", "h2o", "water"},
        {"methanol", "meoh", "methanol"}, {"ethanol", "etoh", "ethanol"}, {"MEG", "meg", "monoethyleneglycol"},
        {"DEG", "deg", "diethyleneglycol"}, {"TEG", "teg", "triethyleneglycol"}, {"methane", "c1", "ch4", "methane"},
        {"ethane", "c2", "c2h6", "ethane"}, {"propane", "c3", "c3h8", "propane"},
        {"i-butane", "ic4", "i-c4", "isobutane", "i-butane"}, {"n-butane", "nc4", "n-c4", "c4", "butane", "n-butane"},
        {"i-pentane", "ic5", "i-c5", "isopentane", "i-pentane"},
        {"n-pentane", "nc5", "n-c5", "c5", "pentane", "n-pentane"}, {"n-hexane", "c6", "nc6", "hexane", "n-hexane"},
        {"n-heptane", "c7", "nc7", "heptane", "n-heptane"}, {"n-octane", "c8", "nc8", "octane", "n-octane"},
        {"n-nonane", "c9", "nc9", "nonane", "n-nonane"}, {"H2S", "h2s", "hydrogensulfide"},
        {"hydrogen", "h2", "hydrogen"}, {"oxygen", "o2", "oxygen"}};
    for (String[] group : groups) {
      for (String name : group) {
        names.put(name.toLowerCase(Locale.ROOT), group[0]);
      }
    }
    for (int i = 10; i <= 20; i++) {
      names.put("c" + i, "nC" + i);
      names.put("nc" + i, "nC" + i);
    }
    return names;
  }
}
