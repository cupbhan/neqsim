package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.util.Locale;

/**
 * Independent electrolyte CPA calculations from molecular formula-unit composition.
 *
 * <p>
 * NaCl is dissociated using NeqSim's own sodium and chloride catalog entries. Phase fractions, compositions and molar
 * masses are returned on the caller's formula-unit basis; ionic accounting remains available in the diagnostic
 * provenance.
 *
 * @author NeqSim development team
 * @version 1.0
 */
final class ElectrolyteComposition {
  /** Prevent construction of this stateless boundary. */
  private ElectrolyteComposition() {
  }

  /**
   * Identify a sodium chloride formula unit without using external salt parameters.
   *
   * @param component declared composition row
   * @return whether the row declares NaCl
   */
  static boolean isSalt(JsonNode component) {
    return "nacl".equals(component.path("name").asText().trim().toLowerCase(Locale.ROOT))
        || "NaCl".equalsIgnoreCase(component.path("formula").asText().trim());
  }

  /**
   * Validate the supported molecular feed and prevent mixed ionic or fitted input.
   *
   * @param input formula-unit composition
   */
  static void validate(JsonNode input) {
    CompositionModel.validate(input);
    CompositionModel.require(!input.has("interactions") || input.path("interactions").size() == 0,
        "ELECTROLYTE_INTERACTION_OVERRIDE_UNSUPPORTED", "Electrolyte CPA uses its own mixing parameters");
    boolean salt = false, water = false;
    for (JsonNode row : input.path("components")) {
      CompositionModel.require(
          !row.has("tbp") && !"pseudo".equals(row.path("kind").asText()) && !"plus".equals(row.path("kind").asText()),
          "ELECTROLYTE_COMPONENT_SCOPE", "Only molecular database components and NaCl are supported");
      if (isSalt(row)) {
        String identity = CompositionModel.resolve(row);
        CompositionModel.require(identity == null || "NaCl".equals(identity), "AMBIGUOUS_COMPONENT_IDENTITY",
            "NaCl cannot also identify a different molecular catalog component");
        CompositionModel.require(!salt, "DUPLICATE_COMPONENT_BINDING", "Only one NaCl input row is permitted");
        salt = true;
      } else {
        String name = CompositionModel.resolve(row);
        CompositionModel.require(name != null && !name.endsWith("+") && !name.endsWith("-"),
            "ELECTROLYTE_COMPONENT_SCOPE", "Provide molecular components, not separate ions");
        water |= "water".equals(name) && row.path("z").asDouble() > 0;
      }
    }
    CompositionModel.require(!salt || water, "ELECTROLYTE_WATER_REQUIRED", "NaCl requires a positive water feed");
  }

  /**
   * Flash a raw molecular feed and restore formula-unit accounting after dissociation.
   *
   * @param input original molecular composition and state
   * @return isolated phase properties and audited formula-unit results
   */
  static ObjectNode calculate(JsonNode input) {
    try {
      validate(input);
      double salt = 0;
      for (JsonNode row : input.path("components")) {
        if (isSalt(row)) {
          salt = row.path("z").asDouble();
        }
      }
      ObjectNode expanded = input.deepCopy();
      ArrayNode components = expanded.putArray("components");
      for (int i = 0; i < input.path("components").size(); i++) {
        JsonNode source = input.path("components").get(i);
        if (isSalt(source)) {
          for (String name : new String[] {"Na+", "Cl-"}) {
            components.addObject().put("id", "ion:" + i + ":" + name).put("name", name).put("kind", "pure").put("z",
                source.path("z").asDouble() / (1 + salt));
          }
        } else {
          ObjectNode copy = source.deepCopy();
          copy.put("id", "molecule:" + i).put("z", source.path("z").asDouble() / (1 + salt));
          components.add(copy);
        }
      }
      ObjectNode result = CompositionFlash.calculate(expanded);
      if (!result.path("ok").asBoolean()) {
        return result;
      }
      ArrayNode ionicBindings = result.path("components").deepCopy();
      double formulaTotal = 0;
      for (JsonNode phase : result.path("phases")) {
        double ions = 0;
        for (JsonNode row : phase.path("composition")) {
          if (row.path("id").asText().startsWith("ion:")) {
            ions += row.path("moleFraction").asDouble();
          }
        }
        double factor = 1 - ions / 2;
        ((ObjectNode) phase).put("ionicMoleFraction", phase.path("moleFraction").asDouble()).put("formulaUnitFactor",
            factor);
        formulaTotal += phase.path("moleFraction").asDouble() * factor;
      }
      CompositionModel.require(Math.abs(formulaTotal - 1 / (1 + salt)) < 1e-7, "EQUILIBRIUM_CHECK_FAILED",
          "Formula-unit total changed during dissociation");
      for (JsonNode phase : result.path("phases")) {
        ObjectNode row = (ObjectNode) phase;
        double factor = row.path("formulaUnitFactor").asDouble();
        row.put("moleFraction", row.path("ionicMoleFraction").asDouble() * factor / formulaTotal);
        row.put("molarMassKgMol", row.path("molarMassKgMol").asDouble() / factor);
        if (row.path("compressibilityFactor").isNumber()) {
          row.put("compressibilityFactor", row.path("compressibilityFactor").asDouble() / factor);
        }
        ArrayNode ionic = row.path("composition").deepCopy();
        ArrayNode molecular = row.putArray("composition");
        for (int i = 0; i < input.path("components").size(); i++) {
          JsonNode source = input.path("components").get(i);
          double x = 0;
          for (JsonNode ion : ionic) {
            String id = ion.path("id").asText();
            if (isSalt(source) ? id.startsWith("ion:" + i + ":") : id.equals("molecule:" + i)) {
              x += ion.path("moleFraction").asDouble();
            }
          }
          molecular.addObject().put("id", source.path("id").asText()).put("moleFraction",
              x / (isSalt(source) ? 2 * factor : factor));
        }
        row.set("ionicComposition", ionic);
      }
      ArrayNode bindings = result.putArray("components");
      for (int i = 0; i < input.path("components").size(); i++) {
        JsonNode source = input.path("components").get(i);
        double mass = 0;
        ObjectNode binding = bindings.addObject().put("id", source.path("id").asText())
            .put("name", source.path("name").asText()).put("sourceIndex", i)
            .put("parameterSource", "neqsim-component-database").put("moleFraction", source.path("z").asDouble());
        ArrayNode nativeNames = binding.putArray("nativeNames");
        for (JsonNode ion : ionicBindings) {
          String id = ion.path("id").asText();
          if (isSalt(source) ? id.startsWith("ion:" + i + ":") : id.equals("molecule:" + i)) {
            mass += ion.path("molarMassKgMol").asDouble();
            nativeNames.add(ion.path("nativeName").asText());
            if (!isSalt(source)) {
              binding.setAll((ObjectNode) ion);
              binding.put("id", source.path("id").asText()).put("sourceIndex", i).put("moleFraction",
                  source.path("z").asDouble());
            }
          }
        }
        binding.put("molarMassKgMol", mass).put("engineMolarMass", mass * 1000);
        if ("nitrogen".equals(binding.path("nativeName").asText())) {
          binding.put("solutePermittivityModel", "constant-liquid-reference-approximation")
              .put("solutePermittivityReferenceTemperatureK", 70.15)
              .put("solutePermittivitySource", "https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular514.pdf");
        }
      }
      ObjectNode diagnostics = (ObjectNode) result.path("diagnostics");
      double balance = 0;
      for (int i = 0; i < input.path("components").size(); i++) {
        double recovered = 0;
        for (JsonNode phase : result.path("phases")) {
          recovered += phase.path("moleFraction").asDouble()
              * phase.path("composition").get(i).path("moleFraction").asDouble();
        }
        balance = Math.max(balance, Math.abs(recovered - input.path("components").get(i).path("z").asDouble()));
      }
      CompositionModel.require(balance < 1e-7, "EQUILIBRIUM_CHECK_FAILED", "Molecular formula units are not conserved");
      diagnostics.put("compositionBasis", "molecular-formula-units").put("dissociation", "NaCl -> Na+ + Cl-")
          .put("maxFormulaUnitBalanceResidual", balance).put("fugacityCheckScope", "neutral-molecular-components")
          .put("ionicPhaseEquilibriumSolver", "NeqSim electrolyte TPflash").put("formulaUnitTotal", formulaTotal)
          .set("ionicComponentBindings", ionicBindings);
      diagnostics.set("ionicInteractions", result.path("interactions"));
      result.putArray("interactions");
      return result;
    } catch (Exception error) {
      return CompositionFlash.failure(error);
    }
  }
}
