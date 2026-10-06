package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.util.Locale;
import java.util.function.DoubleSupplier;
import neqsim.physicalproperties.PhysicalPropertyType;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Independent composition flash using NeqSim's component database, characterization and property methods.
 *
 * @author NeqSim development team
 * @version 1.0
 */
final class CompositionFlash {
  /** Prevent construction of the stateless boundary. */
  private CompositionFlash() {
  }

  /**
   * Calculate equilibrium and isolate unavailable phase properties without importing reference-model parameters.
   *
   * @param input composition and local model settings
   * @return equilibrium, independent parameter provenance, and per-property availability
   */
  static ObjectNode calculate(JsonNode input) {
    ObjectNode response = ExplicitFlash.JSON.createObjectNode();
    try {
      CompositionModel.validate(input);
      double temperature = CompositionModel.number(input, "temperatureC") + 273.15;
      double pressure = CompositionModel.number(input, "pressureMPa") * 10;
      CompositionModel.require(temperature >= 50 && temperature <= 2000 && pressure > 0 && pressure <= 5000,
          "INVALID_STATE", "Temperature or pressure is outside the supported cubic-EOS input bounds");
      CompositionModel model = new CompositionModel(input, temperature, pressure);
      SystemInterface fluid = model.fluid;
      new ThermodynamicOperations(fluid).TPflash();
      fluid.init(1);
      ObjectNode diagnostics = audit(fluid, model.input);
      boolean thermal = true;
      try {
        fluid.init(3);
      } catch (Exception error) {
        thermal = false;
        diagnostics.put("thermalInitializationWarning", error.toString());
        fluid.init(1);
      }
      try {
        fluid.initProperties();
      } catch (Exception error) {
        diagnostics.put("physicalInitializationWarning", error.toString());
      }
      ArrayNode phases = response.putArray("phases");
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        String type = phase.getType().toString().toLowerCase(Locale.ROOT);
        if ("liquid".equals(type)) {
          type = "oil";
        }
        CompositionModel.require("gas".equals(type) || "oil".equals(type) || "aqueous".equals(type),
            "UNSUPPORTED_PHASE", "Unsupported phase: " + type);
        ObjectNode row = phases.addObject().put("phase", type).put("moleFraction", fluid.getBeta(p))
            .put("molarMassKgMol", phase.getMolarMass());
        row.putArray("unavailableProperties");
        row.putObject("propertyMethods");
        row.putObject("propertyAvailability");
        property(row, "density", true, () -> {
          phase.initPhysicalProperties(PhysicalPropertyType.MASS_DENSITY);
          double density = phase.getDensity("kg/m3");
          ((ObjectNode) row.path("propertyMethods")).put("density", phase.getPhysicalProperties().getClass().getName());
          return density;
        });
        property(row, "viscosity", true, () -> {
          phase.initPhysicalProperties(PhysicalPropertyType.DYNAMIC_VISCOSITY);
          double viscosity = phase.getViscosity("cP");
          ((ObjectNode) row.path("propertyMethods")).put("viscosity",
              phase.getPhysicalProperties().getViscosityModel().getClass().getName());
          return viscosity;
        });
        if (thermal) {
          property(row, "enthalpy", false, () -> phase.getEnthalpy("kJ/kg"));
          property(row, "heatCapacity", true, () -> phase.getCp("J/kgK"));
        } else {
          missing(row, "enthalpy", "THERMODYNAMIC_INITIALIZATION_FAILED");
          missing(row, "heatCapacity", "THERMODYNAMIC_INITIALIZATION_FAILED");
        }
        property(row, "compressibilityFactor", true,
            () -> input.path("volumeCorrection").asBoolean(false) ? phase.getZvolcorr() : phase.getZ());
        ArrayNode composition = row.putArray("composition");
        for (JsonNode source : input.path("components")) {
          double x = 0;
          for (int i = 0; i < model.bindings.size(); i++) {
            if (model.bindings.get(i).path("id").equals(source.path("id"))) {
              x = phase.getComponent(i).getx();
            }
          }
          composition.addObject().put("id", source.path("id").asText()).put("moleFraction", x);
        }
      }
      response.set("components", model.bindings);
      response.set("omittedComponents", model.omitted);
      response.set("interactions", model.interactions());
      diagnostics.put("parameterSource", "independent-neqsim-model").put("eos", input.path("eos").asText())
          .put("mixingRule", "CLASSIC (2)").put("volumeCorrection", input.path("volumeCorrection").asBoolean(false))
          .put("formationEnthalpyIncluded", false).put("enthalpyReferenceTemperatureK", 273.15)
          .put("javaVersion", System.getProperty("java.version"));
      response.set("diagnostics", diagnostics);
      response.put("ok", true);
      return response;
    } catch (Exception error) {
      return failure(error);
    }
  }

  /**
   * Check conservation, fugacity equality, selected interactions and enthalpy convention.
   *
   * @param fluid flashed system
   * @param input independent composition request
   * @return equilibrium diagnostics
   */
  static ObjectNode audit(SystemInterface fluid, JsonNode input) {
    double betaTotal = 0, balance = 0, fugacity = 0, total = 0;
    for (JsonNode component : input.path("components")) {
      total += component.path("z").asDouble();
    }
    for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
      PhaseInterface phase = fluid.getPhase(p);
      double beta = fluid.getBeta(p), xTotal = 0;
      CompositionModel.require(Double.isFinite(beta) && beta >= 0 && beta <= 1, "EQUILIBRIUM_CHECK_FAILED",
          "Invalid phase fraction");
      betaTotal += beta;
      for (int i = 0; i < input.path("components").size(); i++) {
        double x = phase.getComponent(i).getx();
        CompositionModel.require(Double.isFinite(x) && x >= 0 && x <= 1, "EQUILIBRIUM_CHECK_FAILED",
            "Invalid phase composition");
        CompositionModel.require(!phase.getComponent(i).isUsingIdealGasEnthalpyOfFormation(),
            "ENTHALPY_REFERENCE_CHANGED", "Formation enthalpy must remain disabled");
        xTotal += x;
        if (p > 0 && input.path("components").get(i).path("z").asDouble() >= 1e-12) {
          double first = fluid.getPhase(0).getComponent(i).getx()
              * fluid.getPhase(0).getComponent(i).getFugacityCoefficient();
          double next = x * phase.getComponent(i).getFugacityCoefficient();
          double residual = Math.abs(Math.log(next / first));
          CompositionModel.require(Double.isFinite(residual), "EQUILIBRIUM_CHECK_FAILED",
              "Nonfinite fugacity residual");
          fugacity = Math.max(fugacity, residual);
        }
      }
      CompositionModel.require(Math.abs(xTotal - 1) < 1e-8, "EQUILIBRIUM_CHECK_FAILED",
          "Phase composition does not sum to one");
      EosMixingRulesInterface mixing = ((PhaseEosInterface) phase).getEosMixingRule();
      for (JsonNode pair : input.path("interactions")) {
        int i = pair.path("i").asInt(), j = pair.path("j").asInt();
        CompositionModel.require(
            Math.abs(mixing.getBinaryInteractionParameter(i, j) - pair.path("value").asDouble()) < 1e-12
                && Math.abs(mixing.getBinaryInteractionParameter(j, i) - pair.path("value").asDouble()) < 1e-12,
            "INTERACTION_CHANGED", "User BIP did not survive the equilibrium calculation");
      }
    }
    for (int i = 0; i < input.path("components").size(); i++) {
      double recovered = 0;
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        recovered += fluid.getBeta(p) * fluid.getPhase(p).getComponent(i).getx();
      }
      balance = Math.max(balance, Math.abs(recovered - input.path("components").get(i).path("z").asDouble() / total));
    }
    CompositionModel.require(Math.abs(betaTotal - 1) < 1e-9 && balance < 1e-7 && fugacity < 1e-6,
        "EQUILIBRIUM_CHECK_FAILED", "Flash did not pass material balance and fugacity checks");
    return ExplicitFlash.JSON.createObjectNode().put("maxMaterialBalanceResidual", balance)
        .put("maxLogFugacityResidual", fugacity).put("phaseFractionTotal", betaTotal)
        .put("inputMoleFractionTotal", total).put("fugacityCheckApplicable", fluid.getNumberOfPhases() > 1);
  }

  /**
   * Evaluate one library property without allowing its failure to erase other properties.
   *
   * @param row phase result
   * @param name property name
   * @param positive whether zero or negative values are invalid
   * @param evaluator isolated property getter
   */
  private static void property(ObjectNode row, String name, boolean positive, DoubleSupplier evaluator) {
    try {
      double value = evaluator.getAsDouble();
      if (!Double.isFinite(value) || (positive && value <= 0)) {
        missing(row, name, "NONFINITE_OR_NONPHYSICAL_PROPERTY");
      } else {
        row.put(name, value);
        ((ObjectNode) row.path("propertyAvailability")).putObject(name).put("status", "available");
      }
    } catch (Exception error) {
      missing(row, name, "PROPERTY_CALCULATION_FAILED");
    }
  }

  /**
   * Mark a single unavailable property explicitly.
   *
   * @param row phase result
   * @param name property name
   * @param reason machine-readable failure reason
   */
  static void missing(ObjectNode row, String name, String reason) {
    row.putNull(name);
    if ("density".equals(name) || "viscosity".equals(name) || "enthalpy".equals(name)) {
      ((ArrayNode) row.path("unavailableProperties")).addObject().put("property", name).put("code", reason);
    }
    ((ObjectNode) row.path("propertyAvailability")).putObject(name).put("status", "unavailable").put("code", reason);
  }

  /**
   * Preserve a useful input error code or report an isolated engine failure.
   *
   * @param error failure
   * @return structured failure without partial values
   */
  static ObjectNode failure(Exception error) {
    String message = error.getMessage() == null ? error.toString() : error.getMessage();
    String code = message.matches("^[A-Z_]+:.*") ? message.substring(0, message.indexOf(':'))
        : "ENGINE_CALCULATION_FAILED";
    return ExplicitFlash.JSON.createObjectNode().put("ok", false).put("code", code).put("message", message);
  }
}
