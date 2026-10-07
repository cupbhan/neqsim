package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import neqsim.thermo.characterization.PlusFractionModel;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;

/**
 * Independent SCN assay closure from a measured composition and bulk heavy-end inputs.
 *
 * <p>
 * Measured carbon-number yields are retained, including the last aggregate plus cut. Ordinary SCN molar masses follow
 * NeqSim's existing Pedersen correlation. By default the terminal-cut molar mass closes the bulk number-average molar
 * mass. Explicit scale-unmeasured mode constrains all unknown cut estimates together. A common multiplier on the
 * Pedersen density curve closes bulk mass/volume density. Defined heavy compounds optionally participate using NeqSim's
 * own catalog values. No reference-engine parameters are accepted. This method retains the terminal plus cut as one EOS
 * pseudo-component; it does not extrapolate its distribution.
 * </p>
 *
 * @author NeqSim development team
 * @version 1.0
 */
final class CompositionCharacterization {
  /** Prevent construction of the stateless adapter. */
  private CompositionCharacterization() {
  }

  /**
   * Generate a mole-basis composition request and a bulk-input closure audit without mutating the input.
   *
   * @param raw raw fractions on one declared basis and bulk heavy-end measurements
   * @return independent effective input and characterization evidence
   */
  static ObjectNode calculate(JsonNode raw) {
    try {
      return resolve(raw);
    } catch (Exception error) {
      return CompositionFlash.failure(error);
    }
  }

  /**
   * Resolve raw SCN properties, close the terminal cut, and convert mass fractions with generated molar masses.
   *
   * @param raw immutable raw composition request
   * @return generated request and audit
   */
  private static ObjectNode resolve(JsonNode raw) {
    CompositionModel.only(raw, "eos", "volumeCorrection", "temperatureC", "pressureMPa", "components", "envelope",
        "basis", "characterization");
    String eos = raw.path("eos").asText(), basis = raw.path("basis").asText();
    CompositionModel.require("PR".equals(eos) || "SRK".equals(eos), "UNSUPPORTED_EOS", "SCN assay requires PR or SRK");
    CompositionModel.require("mole".equals(basis) || "mass".equals(basis), "INVALID_ASSAY_BASIS",
        "Raw SCN composition basis must be mole or mass");
    JsonNode config = raw.path("characterization");
    CompositionModel.only(config, "method", "startCarbonNumber", "terminalCarbonNumber", "molarMassKgMol",
        "densityRelative", "includeDefinedHeavy", "molarMassClosure");
    CompositionModel.require("pedersen-scn-bulk-closure".equals(config.path("method").asText()),
        "UNSUPPORTED_CHARACTERIZATION", "Select pedersen-scn-bulk-closure explicitly");
    int first = carbon(config, "startCarbonNumber"), last = carbon(config, "terminalCarbonNumber");
    CompositionModel.require(first >= 7 && first <= last && last <= 79, "INVALID_ASSAY_RANGE",
        "SCN input range must lie between C7 and C79");
    double bulkMass = CompositionModel.number(config, "molarMassKgMol");
    double bulkDensity = CompositionModel.number(config, "densityRelative");
    CompositionModel.require(bulkMass >= 0.05 && bulkMass <= 2 && bulkDensity >= 0.3 && bulkDensity <= 1.5,
        "INVALID_ASSAY_BULK_INPUT", "Supply a physical bulk molar mass in kg/mol and standard density in g/cm3");
    CompositionModel.require(!config.has("includeDefinedHeavy") || config.path("includeDefinedHeavy").isBoolean(),
        "INVALID_INPUT", "includeDefinedHeavy must be boolean");
    CompositionModel.require(!config.has("molarMassClosure") || config.path("molarMassClosure").isTextual(),
        "INVALID_INPUT", "Molar-mass closure must be text");
    String closure = config.path("molarMassClosure").asText("terminal");
    CompositionModel.require("terminal".equals(closure) || "scale-unmeasured".equals(closure),
        "UNSUPPORTED_CHARACTERIZATION", "Molar-mass closure must be terminal or scale-unmeasured");
    JsonNode rows = raw.path("components");
    CompositionModel.require(rows.isArray() && rows.size() > 0 && rows.size() <= 100, "INVALID_INPUT",
        "Between one and 100 raw components are required");
    int n = rows.size(), terminal = -1;
    double[] amounts = new double[n], masses = new double[n], densities = new double[n];
    boolean[] generated = new boolean[n], inBulk = new boolean[n];
    boolean[] seenCarbon = new boolean[80];
    double total = 0, bulkAmount = 0, knownMass = 0, knownMoles = 0;
    ObjectNode effective = raw.deepCopy();
    effective.remove("basis");
    effective.remove("characterization");
    ArrayNode components = effective.putArray("components");
    SystemInterface catalog = "PR".equals(eos) ? new SystemPrEos(288.15, 1.01325) : new SystemSrkEos(288.15, 1.01325);
    for (int i = 0; i < n; i++) {
      JsonNode row = rows.get(i);
      CompositionModel.only(row, "id", "name", "formula", "kind", "z", "tbp", "carbonNumber");
      amounts[i] = CompositionModel.number(row, "z");
      CompositionModel.require(amounts[i] >= 0 && amounts[i] <= 1, "INVALID_COMPOSITION", "Invalid raw fraction");
      total += amounts[i];
      ObjectNode copy = row.deepCopy();
      copy.remove("carbonNumber");
      components.add(copy);
      int cn = row.has("carbonNumber") ? carbon(row, "carbonNumber") : 0;
      boolean fraction = "pseudo".equals(row.path("kind").asText()) || "plus".equals(row.path("kind").asText());
      if (fraction) {
        CompositionModel.require(cn >= 7 && cn <= last && !seenCarbon[cn], "INVALID_ASSAY_RANGE",
            "Petroleum cuts require unique carbon numbers within the declared terminal range");
        seenCarbon[cn] = true;
        generated[i] = cn >= first;
        inBulk[i] = generated[i];
        if (generated[i]) {
          CompositionModel.require(!row.has("tbp"), "ASSAY_BULK_TBP_CONFLICT",
              "Bulk-input rows cannot also supply individually characterized properties");
          masses[i] = PlusFractionModel.scnMolarMassKgMol(cn);
          densities[i] = PlusFractionModel.scnReferenceDensity(cn);
          if (cn == last) {
            terminal = i;
          }
        } else {
          CompositionModel.only(row.path("tbp"), "molarMassKgMol", "densityRelative");
          masses[i] = CompositionModel.number(row.path("tbp"), "molarMassKgMol");
          densities[i] = CompositionModel.number(row.path("tbp"), "densityRelative");
        }
      } else {
        String name = CompositionModel.resolve(row);
        CompositionModel.require(name != null && !row.has("tbp"), "UNKNOWN_COMPONENT_REQUIRES_TBP",
            "Defined raw components must resolve to NeqSim's own catalog");
        catalog.addComponent(name, 1);
        masses[i] = catalog.getComponent(name).getMolarMass();
        densities[i] = catalog.getComponent(name).getNormalLiquidDensity();
        inBulk[i] = config.path("includeDefinedHeavy").asBoolean(false) && cn >= first;
      }
      CompositionModel.require(Double.isFinite(masses[i]) && masses[i] > 0, "INVALID_ASSAY_INPUT",
          "Component molar mass must be finite and positive");
      if (inBulk[i]) {
        bulkAmount += amounts[i];
        if (!(fraction && cn == last)) {
          knownMass += "mass".equals(basis) ? amounts[i] : amounts[i] * masses[i];
          knownMoles += "mass".equals(basis) ? amounts[i] / masses[i] : amounts[i];
        }
      }
    }
    CompositionModel.require(Math.abs(total - 1) < 1e-8, "INVALID_COMPOSITION", "Raw fractions must sum to one");
    CompositionModel.require(terminal >= 0 && amounts[terminal] > 0 && bulkAmount > 0, "ASSAY_TERMINAL_REQUIRED",
        "A positive terminal plus cut is required");
    for (int cn = first; cn <= last; cn++) {
      CompositionModel.require(seenCarbon[cn], "ASSAY_CARBON_GAP", "Missing carbon-number cut C" + cn);
    }
    double massScale = 1;
    if ("scale-unmeasured".equals(closure)) {
      double fixedMass = 0, fixedMoles = 0, seedMass = 0, seedMoles = 0;
      for (int i = 0; i < n; i++) {
        if (!inBulk[i]) {
          continue;
        }
        double mass = "mass".equals(basis) ? amounts[i] : amounts[i] * masses[i];
        double moles = "mass".equals(basis) ? amounts[i] / masses[i] : amounts[i];
        if (generated[i]) {
          seedMass += mass;
          seedMoles += moles;
        } else {
          fixedMass += mass;
          fixedMoles += moles;
        }
      }
      massScale = "mass".equals(basis) ? seedMoles / (bulkAmount / bulkMass - fixedMoles)
          : (bulkAmount * bulkMass - fixedMass) / seedMass;
      CompositionModel.require(Double.isFinite(massScale) && massScale > 0, "ASSAY_BULK_MASS_INCONSISTENT",
          "Bulk molar mass cannot close after retaining defined heavy compounds");
      for (int i = 0; i < n; i++) {
        if (generated[i]) {
          masses[i] *= massScale;
          double minimum = 12.011 * rows.get(i).path("carbonNumber").asInt() / 1000;
          CompositionModel.require(masses[i] >= minimum && masses[i] <= 2, "ASSAY_BULK_MASS_INCONSISTENT",
              "Constrained SCN molar mass violates its carbon-mass lower bound or the model limit");
        }
      }
    } else {
      double terminalMass = "mass".equals(basis) ? amounts[terminal] / (bulkAmount / bulkMass - knownMoles)
          : (bulkAmount * bulkMass - knownMass) / amounts[terminal];
      CompositionModel.require(Double.isFinite(terminalMass) && terminalMass >= masses[terminal] && terminalMass <= 2,
          "ASSAY_BULK_MASS_INCONSISTENT",
          "Bulk molar mass and measured cut yields imply an invalid terminal-cut molar mass: " + terminalMass);
      masses[terminal] = terminalMass;
    }
    double targetMass = 0, targetMoles = 0, fixedVolume = 0, seedVolume = 0;
    for (int i = 0; i < n; i++) {
      if (!inBulk[i] || amounts[i] == 0) {
        continue;
      }
      double mass = "mass".equals(basis) ? amounts[i] : amounts[i] * masses[i];
      double moles = "mass".equals(basis) ? amounts[i] / masses[i] : amounts[i];
      CompositionModel.require(Double.isFinite(densities[i]) && densities[i] > 0, "ASSAY_DENSITY_UNAVAILABLE",
          "Defined heavy component lacks a normal liquid density");
      targetMass += mass;
      targetMoles += moles;
      if (generated[i]) {
        seedVolume += mass / densities[i];
      } else {
        fixedVolume += mass / densities[i];
      }
    }
    double multiplier = seedVolume / (targetMass / bulkDensity - fixedVolume);
    CompositionModel.require(Double.isFinite(multiplier) && multiplier > 0, "ASSAY_BULK_DENSITY_INCONSISTENT",
        "Bulk density cannot close with the measured composition and own catalog values");
    double moleTotal = 0, recoveredVolume = fixedVolume;
    ArrayNode generatedCuts = ExplicitFlash.JSON.createArrayNode();
    for (int i = 0; i < n; i++) {
      ObjectNode row = (ObjectNode) components.get(i);
      if (generated[i]) {
        densities[i] *= multiplier;
        CompositionModel.require(densities[i] >= 0.3 && densities[i] <= 1.5, "ASSAY_DENSITY_OUT_OF_RANGE",
            "Generated cut density is outside [0.3, 1.5] g/cm3");
        row.putObject("tbp").put("molarMassKgMol", masses[i]).put("densityRelative", densities[i]);
        row.put("kind", "pseudo");
        generatedCuts.addObject().put("id", row.path("id").asText()).put("sourceIndex", i)
            .put("molarMassKgMol", masses[i]).put("densityRelative", densities[i])
            .put("terminalAggregate", i == terminal);
        recoveredVolume += ("mass".equals(basis) ? amounts[i] : amounts[i] * masses[i]) / densities[i];
      }
      double moles = "mass".equals(basis) ? amounts[i] / masses[i] : amounts[i];
      row.put("z", moles);
      moleTotal += moles;
    }
    for (JsonNode row : components) {
      ((ObjectNode) row).put("z", row.path("z").asDouble() / moleTotal);
    }
    CompositionModel.validate(effective);
    ObjectNode result = ExplicitFlash.JSON.createObjectNode().put("ok", true);
    result.set("effectiveInput", effective);
    ObjectNode audit = result.putObject("characterization").put("method", "pedersen-scn-bulk-closure")
        .put("inputBasis", basis).put("terminalCutRepresentation", "aggregate-tbp").put("molarMassClosureMode", closure)
        .put("molarMassCurveScale", massScale).put("bulkMolarMassKgMol", targetMass / targetMoles)
        .put("bulkDensityRelative", targetMass / recoveredVolume)
        .put("molarMassClosureResidual", targetMass / targetMoles - bulkMass)
        .put("densityClosureResidual", targetMass / recoveredVolume - bulkDensity).put("densityCurveScale", multiplier)
        .put("sourceParametersUsed", false);
    audit.set("generatedCuts", generatedCuts);
    return result;
  }

  /**
   * Validate an integral carbon-number input without truncating decimals.
   *
   * @param node input object
   * @param key field name
   * @return carbon number
   */
  private static int carbon(JsonNode node, String key) {
    CompositionModel.require(
        node.path(key).isIntegralNumber() && node.path(key).asInt() >= 0 && node.path(key).asInt() <= 200,
        "INVALID_ASSAY_RANGE", "Carbon number must be an integer in [0, 200]");
    return node.path(key).asInt();
  }
}
