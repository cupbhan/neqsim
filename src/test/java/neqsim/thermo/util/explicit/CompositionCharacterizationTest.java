package neqsim.thermo.util.explicit;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import org.junit.jupiter.api.Test;
import neqsim.thermo.characterization.PlusFractionModel;

/**
 * Raw-assay conservation, independence, immutable inputs and structured failure contracts.
 *
 * @author NeqSim development team
 * @version 1.0
 */
class CompositionCharacterizationTest {
  /** JSON factory for synthetic raw inputs. */
  private static final ObjectMapper JSON = new ObjectMapper();

  /**
   * Create an uncharacterized three-cut assay and methane with no reference parameters.
   *
   * @param basis raw composition basis
   * @return raw assay request
   */
  private ObjectNode input(String basis) {
    ObjectNode input = JSON.createObjectNode().put("eos", "PR").put("basis", basis).put("temperatureC", 80)
        .put("pressureMPa", 20).put("volumeCorrection", true);
    ArrayNode rows = input.putArray("components");
    rows.addObject().put("id", "gas").put("name", "methane").put("kind", "pure").put("z", 0.7);
    for (int cn = 7; cn <= 9; cn++) {
      rows.addObject().put("id", "cut-" + cn).put("name", "C" + cn).put("kind", "pseudo").put("carbonNumber", cn)
          .put("z", 0.1);
    }
    input.putObject("characterization").put("method", "pedersen-scn-bulk-closure").put("startCarbonNumber", 7)
        .put("terminalCarbonNumber", 9).put("molarMassKgMol", 0.13).put("densityRelative", 0.8);
    return input;
  }

  /** Validate both basis conversions against independently recomputed mass, mole and volume balances. */
  @Test
  void bothBasesCloseWithoutChangingRawInputs() {
    for (String basis : new String[] {"mole", "mass"}) {
      ObjectNode raw = input(basis);
      String original = raw.toString();
      JsonNode result = ThermodynamicEngine.compositionCharacterize(raw);
      assertTrue(result.path("ok").asBoolean(), result.toString());
      assertEquals(original, raw.toString());
      JsonNode audit = result.path("characterization");
      assertEquals(0, audit.path("molarMassClosureResidual").asDouble(), 1e-12);
      assertEquals(0, audit.path("densityClosureResidual").asDouble(), 1e-12);
      assertFalse(audit.path("sourceParametersUsed").asBoolean());
      double mole = 0, mass = 0, volume = 0, total = 0;
      JsonNode effective = result.path("effectiveInput");
      for (JsonNode row : effective.path("components")) {
        total += row.path("z").asDouble();
        if (row.has("tbp")) {
          double z = row.path("z").asDouble(), m = row.path("tbp").path("molarMassKgMol").asDouble();
          double rho = row.path("tbp").path("densityRelative").asDouble();
          mole += z;
          mass += z * m;
          volume += z * m / rho;
        }
      }
      assertEquals(1, total, 1e-12);
      assertEquals(0.13, mass / mole, 1e-12);
      assertEquals(0.8, mass / volume, 1e-12);
      assertEquals(PlusFractionModel.scnMolarMassKgMol(7),
          effective.path("components").get(1).path("tbp").path("molarMassKgMol").asDouble(), 1e-12);
      assertEquals("aggregate-tbp", audit.path("terminalCutRepresentation").asText());
    }
  }

  /** Unknown GC cut estimates can be constrained together without changing the source bulk measurements. */
  @Test
  void scaledUnknownCutsCloseBothBasesAndRespectPhysicalBounds() {
    for (String basis : new String[] {"mole", "mass"}) {
      ObjectNode raw = input(basis);
      raw.withObject("characterization").put("molarMassClosure", "scale-unmeasured");
      String original = raw.toString();
      JsonNode result = ThermodynamicEngine.compositionCharacterize(raw);
      assertTrue(result.path("ok").asBoolean(), result.toString());
      assertEquals(original, raw.toString());
      JsonNode audit = result.path("characterization");
      assertEquals("scale-unmeasured", audit.path("molarMassClosureMode").asText());
      assertEquals(0.13, audit.path("bulkMolarMassKgMol").asDouble(), 1e-12);
      assertEquals(0.8, audit.path("bulkDensityRelative").asDouble(), 1e-12);
      double scale = audit.path("molarMassCurveScale").asDouble();
      assertTrue(scale > 0);
      for (int cn = 7; cn <= 9; cn++) {
        double mass = result.path("effectiveInput").path("components").get(cn - 6).path("tbp").path("molarMassKgMol")
            .asDouble();
        assertEquals(PlusFractionModel.scnMolarMassKgMol(cn) * scale, mass, 1e-12);
        assertTrue(mass >= 12.011 * cn / 1000);
      }
      raw.withObject("characterization").put("molarMassKgMol", 0.05);
      assertEquals("ASSAY_BULK_MASS_INCONSISTENT",
          ThermodynamicEngine.compositionCharacterize(raw).path("code").asText());
    }
  }

  /** Verify raw inputs can be flashed directly and share the same result as the independently resolved request. */
  @Test
  void rawFlashMatchesResolvedCompositionForBothEosFamilies() {
    for (String eos : new String[] {"PR", "SRK"}) {
      ObjectNode raw = input("mass").put("eos", eos);
      JsonNode characterized = ThermodynamicEngine.compositionCharacterize(raw);
      JsonNode direct = ThermodynamicEngine.compositionFlash(raw);
      JsonNode resolved = ThermodynamicEngine.compositionFlash(characterized.path("effectiveInput"));
      assertTrue(direct.path("ok").asBoolean(), direct.toString());
      assertTrue(resolved.path("ok").asBoolean(), resolved.toString());
      assertEquals(resolved.path("phases"), direct.path("phases"));
      assertTrue(direct.has("characterization"));
    }
  }

  /** Reject ambiguous input scopes, forbidden EOS parameters, carbon gaps and incompatible measured mass. */
  @Test
  void rejectsInconsistentOrBorrowedParameters() {
    ObjectNode raw = input("mole");
    raw.withObject("characterization").put("molarMassKgMol", 0.08);
    assertEquals("ASSAY_BULK_MASS_INCONSISTENT",
        ThermodynamicEngine.compositionCharacterize(raw).path("code").asText());
    raw = input("mole");
    ((ObjectNode) raw.path("components").get(1)).putObject("tbp").put("molarMassKgMol", 0.13).put("densityRelative",
        0.8);
    assertEquals("ASSAY_BULK_TBP_CONFLICT", ThermodynamicEngine.compositionCharacterize(raw).path("code").asText());
    raw = input("mole");
    raw.withObject("characterization").put("criticalTemperatureK", 600);
    assertEquals("UNSUPPORTED_PARAMETER", ThermodynamicEngine.compositionCharacterize(raw).path("code").asText());
    raw = input("mole");
    ((ObjectNode) raw.path("components").get(2)).put("carbonNumber", 7.5);
    assertEquals("INVALID_ASSAY_RANGE", ThermodynamicEngine.compositionCharacterize(raw).path("code").asText());
    raw = input("mole");
    ((ObjectNode) raw.path("components").get(2)).put("carbonNumber", 6);
    assertEquals("INVALID_ASSAY_RANGE", ThermodynamicEngine.compositionCharacterize(raw).path("code").asText());
  }

  /** Defined GC heavy compounds contribute using own catalog properties, rather than reference-engine molar masses. */
  @Test
  void definedHeavyInputsUseOwnCatalogAndPreserveEarlierMeasuredCuts() {
    ObjectNode raw = input("mole");
    ArrayNode rows = raw.withArray("components");
    ((ObjectNode) rows.get(0)).put("z", 0.6);
    rows.addObject().put("id", "defined-heavy").put("name", "nC7").put("kind", "pure").put("carbonNumber", 7).put("z",
        0.1);
    raw.withObject("characterization").put("includeDefinedHeavy", true);
    JsonNode result = ThermodynamicEngine.compositionCharacterize(raw);
    assertTrue(result.path("ok").asBoolean(), result.toString());
    assertEquals(0.13, result.path("characterization").path("bulkMolarMassKgMol").asDouble(), 1e-12);
    assertFalse(result.path("effectiveInput").path("components").get(4).has("tbp"));
    assertEquals(0.1, result.path("effectiveInput").path("components").get(4).path("z").asDouble(), 1e-12);
  }
}
