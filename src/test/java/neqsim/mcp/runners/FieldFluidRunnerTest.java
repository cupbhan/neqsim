package neqsim.mcp.runners;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;

/** Tests for field-fluid electrolyte and characterized-heavy-oil calculations. */
class FieldFluidRunnerTest extends NeqSimTest {

  @Test
  void retainsWaterRichThreePhaseSolutionAfterSupplementaryStabilityRetry() throws Exception {
    InputStream input = FieldFluidRunnerTest.class
        .getResourceAsStream("/neqsim/mcp/runners/field_wet_classic_payload.json");
    assertNotNull(input);
    try (InputStreamReader reader = new InputStreamReader(input, StandardCharsets.UTF_8)) {
      JsonObject request = JsonParser.parseReader(reader).getAsJsonObject();
      request.addProperty("temperatureC", 215.0);
      request.addProperty("pressureBara", 240.0);
      request.addProperty("temperatureMinC", 215.0);
      request.addProperty("temperatureMaxC", 215.0);
      request.addProperty("eos", "SRK");
      request.addProperty("reactive", false);

      JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request.toString())).getAsJsonObject();
      JsonObject flash = result.getAsJsonObject("flash");
      assertNotNull(flash, result.toString());
      String phases = flash.getAsJsonArray("phases").toString();

      assertEquals("success", result.get("status").getAsString(), result.toString());
      assertEquals(3, flash.get("numberOfPhases").getAsInt(), result.toString());
      assertTrue(phases.contains("gas"), result.toString());
      assertTrue(phases.contains("aqueous"), result.toString());
      assertTrue(phases.contains("oil"), result.toString());
    }
  }

  @Test
  void exposesBenchmarkRegressionWithoutPromotingItToEngineeringAcceptance() {
    String fingerprint = "\"modelFingerprint\":{\"equationOfState\":\"SRK\","
        + "\"volumeTranslation\":\"Peneloux\",\"mixingRule\":\"Huron-Vidal\","
        + "\"compositionFingerprint\":\"fluid-sha256\"," + "\"parameterFingerprint\":\"srk-hv-parameter-sha256\"}";
    String points = "\"points\":[{\"temperatureK\":300.0,\"pressureBara\":100.0},"
        + "{\"temperatureK\":350.0,\"pressureBara\":200.0}]";
    String reference = "{\"identifier\":\"pvtsim-3-aq\",\"source\":\"PVTSIM\"," + fingerprint
        + ",\"boundaryTopology\":\"GO_TO_GOW\",\"topologyResolved\":true,"
        + "\"topologyEvidence\":\"PVTsim 3-Aq contract and phase-side flashes\"," + points + "}";
    String candidate = "{\"identifier\":\"neqsim-go-gow\",\"source\":\"NEQSIM\"," + fingerprint
        + ",\"boundaryTopology\":\"GO_TO_GOW\",\"topologyResolved\":true,"
        + "\"topologyEvidence\":\"global stability gate\"," + points + "}";
    String criteria = "{\"minimumReferenceCoverageFraction\":1.0," + "\"minimumCandidateCoverageFraction\":1.0,"
        + "\"maximumRelativeTemperatureDifference\":0.01," + "\"maximumRmsRelativeTemperatureDifference\":0.01,"
        + "\"maximumP95RelativeTemperatureDifference\":0.01,"
        + "\"maximumSymmetricRelativeTemperatureDifference\":0.01," + "\"criticalPointRequired\":false,"
        + "\"maximumCriticalTemperatureDifference\":0.01," + "\"maximumCriticalPressureDifference\":0.01}";
    String request = "{\"components\":{\"water\":1.0},\"temperatureC\":100.0,"
        + "\"pressureBara\":2.0,\"temperatureMinC\":100.0,\"temperatureMaxC\":100.0,"
        + "\"reactive\":false,\"hydrocarbonWaterBoundaryRegression\":{\"reference\":" + reference + ",\"candidate\":"
        + candidate + ",\"criteria\":" + criteria + "}}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();
    JsonObject regression = result.getAsJsonObject("hydrocarbonWaterBoundaryRegression");
    JsonObject quality = result.getAsJsonObject("qualityGate");

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertTrue(regression.get("accepted").getAsBoolean(), result.toString());
    assertTrue(!regression.get("engineeringEligible").getAsBoolean(), result.toString());
    assertTrue(quality.get("hydrocarbonWaterBoundaryRegressionCalculated").getAsBoolean(), result.toString());
    assertTrue(quality.get("hydrocarbonWaterBoundaryRegressionAccepted").getAsBoolean(), result.toString());
    assertEquals(regression, result.getAsJsonObject("data").getAsJsonObject("hydrocarbonWaterBoundaryRegression"));
  }

  @Test
  void exposesStableRegionScanWithoutRunningStrictBoundaryOrNetworkPromotion() {
    String definition = "{\"format\":\"neqsim-fluid\",\"version\":\"1.0\",\"eos\":\"SRK\","
        + "\"components\":[{\"name\":\"heavy-cut\",\"moleFraction\":1.0,"
        + "\"criticalTemperature\":800.0,\"criticalPressure\":12.0,"
        + "\"acentricFactor\":1.1,\"molarMass\":500.0,\"normalBoilingPoint\":700.0,"
        + "\"criticalVolume\":1.2,\"volumeShift\":0.1,\"parachor\":700.0,"
        + "\"density\":0.95,\"isPseudo\":true}],\"binaryInteractionCoefficients\":[]}";
    String request = "{\"components\":{\"water\":0.65,\"heavyOil\":0.25,\"CO2\":0.07,"
        + "\"nitrogen\":0.03},\"fluidDefinition\":" + definition
        + ",\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":200.0,"
        + "\"temperatureMaxC\":200.0,\"eos\":\"SRK\",\"reactive\":false,"
        + "\"hydrocarbonWaterStableRegionTopologyScan\":{"
        + "\"temperaturesC\":[100.0,200.0],\"pressuresBara\":[10.0,20.0]}}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();
    JsonObject scan = result.getAsJsonObject("hydrocarbonWaterStableRegionTopologyScan");

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertEquals("stable-multiphase-tp-region-scan-only", scan.get("method").getAsString());
    assertTrue(scan.get("diagnosticOnly").getAsBoolean(), result.toString());
    assertTrue(!scan.get("engineeringEligible").getAsBoolean(), result.toString());
    assertTrue(!scan.get("strictBoundaryCalculated").getAsBoolean(), result.toString());
    assertEquals(4, scan.get("requestedGridStateCount").getAsInt());
    assertTrue(
        result.getAsJsonObject("qualityGate").get("hydrocarbonWaterStableRegionTopologyScanCalculated").getAsBoolean());
    assertEquals(scan, result.getAsJsonObject("data").getAsJsonObject("hydrocarbonWaterStableRegionTopologyScan"));
  }

  @Test
  void calculatesSelectedCubicModelForWaterAmmoniaCarbonDioxide() {
    for (String eos : new String[] { "SRK", "PR" }) {
      String request = "{\"components\":{\"water\":0.90,\"ammonia\":0.05,\"CO2\":0.05},"
          + "\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":20.0,"
          + "\"temperatureMaxC\":350.0,\"pointCount\":31,\"eos\":\"" + eos + "\",\"reactive\":false}";

      JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();

      assertEquals("success", result.get("status").getAsString(), result.toString());
      assertEquals(eos, result.get("model").getAsString(), result.toString());
      assertTrue(result.getAsJsonObject("flash").get("numberOfPhases").getAsInt() >= 1);
      assertTrue(!result.getAsJsonObject("speciation").get("enabled").getAsBoolean());
      assertTrue(result.getAsJsonArray("envelope").size() >= 24, result.toString());
    }
  }

  @Test
  void calculatesReactiveWaterAmmoniaCarbonDioxideState() {
    String request = "{\"components\":{\"water\":0.90,\"ammonia\":0.05,\"CO2\":0.05},"
        + "\"temperatureC\":80.0,\"pressureBara\":10.0,\"temperatureMinC\":20.0,"
        + "\"temperatureMaxC\":120.0,\"pointCount\":11,\"reactive\":true}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertEquals("Electrolyte-CPA", result.get("model").getAsString());
    assertTrue(result.getAsJsonObject("flash").get("numberOfPhases").getAsInt() >= 1);
    assertTrue(result.getAsJsonObject("flash").get("reactiveEquilibriumConverged").getAsBoolean());
    assertTrue(result.getAsJsonObject("speciation").get("enabled").getAsBoolean());
    assertTrue(result.getAsJsonArray("envelope").size() >= 8, result.toString());
  }

  @Test
  void excludesDiscontinuousReactiveBoundaryPoints() {
    String request = "{\"components\":{\"water\":0.90,\"ammonia\":0.05,\"CO2\":0.05},"
        + "\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":20.0,"
        + "\"temperatureMaxC\":350.0,\"pointCount\":31,\"reactive\":true}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();
    JsonObject quality = result.getAsJsonObject("qualityGate");

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertEquals("passed", quality.get("verdict").getAsString(), result.toString());
    assertTrue(quality.get("boundaryContinuityPassed").getAsBoolean(), result.toString());
    assertTrue(quality.get("rejectedBoundaryPointCount").getAsInt() >= 1, result.toString());
    assertTrue(result.getAsJsonArray("envelope").size() >= 24, result.toString());
  }

  @Test
  void retainsCharacterizedSubfluidInHeavyOilBlend() {
    String definition = "{\"format\":\"neqsim-fluid\",\"version\":\"1.0\",\"eos\":\"SRK\","
        + "\"components\":[{\"name\":\"heavy-cut\",\"moleFraction\":1.0,"
        + "\"criticalTemperature\":800.0,\"criticalPressure\":12.0,"
        + "\"acentricFactor\":1.1,\"molarMass\":500.0,\"normalBoilingPoint\":700.0,"
        + "\"criticalVolume\":1.2,\"volumeShift\":0.1,\"parachor\":700.0,"
        + "\"density\":0.95,\"isPseudo\":true}],\"binaryInteractionCoefficients\":[]}";
    String request = "{\"components\":{\"water\":0.65,\"heavyOil\":0.25,\"CO2\":0.07,"
        + "\"nitrogen\":0.03},\"fluidDefinition\":" + definition
        + ",\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":200.0,"
        + "\"temperatureMaxC\":200.0,\"pointCount\":11,\"reactive\":false}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertEquals("SRK", result.get("model").getAsString());
    assertTrue(!result.getAsJsonObject("qualityGate").get("boundaryComparisonEligible").getAsBoolean(),
        result.toString());
    assertTrue(
        result.getAsJsonObject("fluid").getAsJsonObject("composition").getAsJsonObject("overall").has("heavy-cut"));
  }

  @Test
  void exposesOptInHydrocarbonWaterBoundaryAuditWithoutPromotingAnchorsToEngineeringCurve() {
    String definition = "{\"format\":\"neqsim-fluid\",\"version\":\"1.0\",\"eos\":\"SRK\","
        + "\"components\":[{\"name\":\"heavy-cut\",\"moleFraction\":1.0,"
        + "\"criticalTemperature\":800.0,\"criticalPressure\":12.0,"
        + "\"acentricFactor\":1.1,\"molarMass\":500.0,\"normalBoilingPoint\":700.0,"
        + "\"criticalVolume\":1.2,\"volumeShift\":0.1,\"parachor\":700.0,"
        + "\"density\":0.95,\"isPseudo\":true}],\"binaryInteractionCoefficients\":[]}";
    String request = "{\"components\":{\"water\":0.65,\"heavyOil\":0.25,\"CO2\":0.07,"
        + "\"nitrogen\":0.03},\"fluidDefinition\":" + definition
        + ",\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":200.0,"
        + "\"temperatureMaxC\":200.0,\"eos\":\"SRK\",\"reactive\":false,"
        + "\"hydrocarbonWaterBoundaryAudit\":{\"temperaturesC\":[100.0,200.0]," + "\"pressuresBara\":[10.0,20.0]}}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();

    assertEquals("success", result.get("status").getAsString(), result.toString());
    JsonObject audit = result.getAsJsonObject("hydrocarbonWaterBoundaryAudit");
    assertEquals("success", audit.get("status").getAsString(), result.toString());
    assertEquals("hydrocarbon-water-two-to-three-phase", audit.get("boundaryTopology").getAsString());
    assertTrue(!audit.get("engineeringEligible").getAsBoolean(), result.toString());
    assertEquals(80, audit.get("continuationPointsPerDirection").getAsInt());
    assertTrue(audit.has("network"), result.toString());
    assertTrue(audit.getAsJsonObject("network").get("attempted").getAsBoolean(), result.toString());
    assertEquals(audit.get("engineeringEligible").getAsBoolean(),
        audit.getAsJsonObject("network").get("engineeringEligible").getAsBoolean());
    assertTrue(result.getAsJsonObject("qualityGate").get("hydrocarbonWaterBoundaryAuditCalculated").getAsBoolean());
    assertEquals(audit, result.getAsJsonObject("data").getAsJsonObject("hydrocarbonWaterBoundaryAudit"));
  }

  @Test
  void auditsHeavyFluidEosTransferWithoutRetuning() {
    String definition = "{\"format\":\"neqsim-fluid\",\"version\":\"1.0\",\"sourceEos\":\"SRK\",\"eos\":\"PR\","
        + "\"components\":[{\"name\":\"heavy-cut\",\"moleFraction\":1.0,"
        + "\"criticalTemperature\":800.0,\"criticalPressure\":12.0,"
        + "\"acentricFactor\":1.1,\"molarMass\":500.0,\"normalBoilingPoint\":700.0,"
        + "\"criticalVolume\":1.2,\"volumeShift\":0.1,\"parachor\":700.0,"
        + "\"density\":0.95,\"isPseudo\":true}],\"binaryInteractionCoefficients\":[]}";
    String request = "{\"components\":{\"water\":0.65,\"heavyOil\":0.25,\"CO2\":0.07,"
        + "\"nitrogen\":0.03},\"fluidDefinition\":" + definition
        + ",\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":200.0,"
        + "\"temperatureMaxC\":200.0,\"pointCount\":11,\"eos\":\"PR\",\"reactive\":false}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertEquals("PR", result.get("model").getAsString());
    assertEquals("SRK", result.get("sourceModel").getAsString());
    assertEquals("PR", result.get("candidateModel").getAsString());
    assertTrue(result.getAsJsonObject("qualityGate").get("modelTransferredFromSource").getAsBoolean());
    assertTrue(result.getAsJsonArray("limitations").toString().contains("without retuning"), result.toString());
  }

  @Test
  void tracesHeavyFluidTwoHydrocarbonPhaseTopology() {
    String definition = "{\"format\":\"neqsim-fluid\",\"version\":\"1.0\",\"eos\":\"SRK\","
        + "\"polarModel\":\"HV\",\"huronVidalInteractionParameters\":[{\"i\":\"heavy-cut\","
        + "\"j\":\"CO2\",\"forward\":100.0,\"reverse\":-50.0}],"
        + "\"components\":[{\"name\":\"heavy-cut\",\"moleFraction\":1.0,"
        + "\"criticalTemperature\":800.0,\"criticalPressure\":12.0,"
        + "\"acentricFactor\":1.1,\"molarMass\":500.0,\"normalBoilingPoint\":700.0,"
        + "\"criticalVolume\":1.2,\"volumeShift\":0.1,\"parachor\":700.0,"
        + "\"density\":0.95,\"isPseudo\":true}],\"binaryInteractionCoefficients\":[]}";
    String request = "{\"components\":{\"water\":0.65,\"heavyOil\":0.25,\"CO2\":0.07,"
        + "\"nitrogen\":0.03},\"fluidDefinition\":" + definition
        + ",\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":20.0,"
        + "\"temperatureMaxC\":350.0,\"phaseMapTemperaturePoints\":21,"
        + "\"phaseMapPressurePoints\":13,\"reactive\":false}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();
    JsonObject quality = result.getAsJsonObject("qualityGate");

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertEquals("michelsen-predictor-corrector", quality.get("boundaryMethod").getAsString());
    assertEquals("hydrocarbon-one-to-two-phase", quality.get("boundaryTopology").getAsString());
    assertEquals("hydrocarbon-one-to-two-phase", result.get("boundaryTopology").getAsString());
    assertTrue(!quality.get("phaseMapCalculated").getAsBoolean(), result.toString());
    assertEquals(0, result.getAsJsonArray("phaseMap").size(), result.toString());
    assertEquals(0, quality.get("stabilityFlashCount").getAsInt(), result.toString());
    assertEquals(0, quality.get("stabilityFailedFlashCount").getAsInt(), result.toString());
    assertTrue(quality.get("envelopeSegmentCount").getAsInt() >= 1, result.toString());
    assertTrue(result.getAsJsonArray("envelope").size() >= 16, result.toString());
  }

  @Test
  void tracesTheHydrocarbonEnvelopeOnItsOwnBasisInsteadOfTheWaterBlendedFeed() {
    String definition = "{\"format\":\"neqsim-fluid\",\"version\":\"1.0\",\"eos\":\"SRK\","
        + "\"components\":[{\"name\":\"heavy-cut\",\"moleFraction\":1.0,"
        + "\"criticalTemperature\":800.0,\"criticalPressure\":12.0,"
        + "\"acentricFactor\":1.1,\"molarMass\":500.0,\"normalBoilingPoint\":700.0,"
        + "\"criticalVolume\":1.2,\"volumeShift\":0.1,\"parachor\":700.0,"
        + "\"density\":0.95,\"isPseudo\":true}],\"binaryInteractionCoefficients\":[]}";
    String request = "{\"components\":{\"water\":0.65,\"heavyOil\":0.25,\"CO2\":0.07,"
        + "\"nitrogen\":0.03},\"fluidDefinition\":" + definition
        + ",\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":20.0,"
        + "\"temperatureMaxC\":350.0,\"phaseMapTemperaturePoints\":21,"
        + "\"phaseMapPressurePoints\":13,\"reactive\":false}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();
    JsonObject quality = result.getAsJsonObject("qualityGate");

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertEquals("hydrocarbon-only", quality.get("envelopeBasis").getAsString(), result.toString());
    assertEquals(0.65, quality.get("envelopeBasisExcludedAqueousFraction").getAsDouble(), 1.0e-12, result.toString());
    assertTrue(result.getAsJsonArray("limitations").toString().contains("hydrocarbon sub-composition"),
        result.toString());
    // The state point keeps the aqueous fraction the envelope basis dropped.
    assertTrue(result.getAsJsonObject("fluid").getAsJsonObject("composition").getAsJsonObject("overall").has("water"),
        result.toString());
    // Whatever the boundary turned out to be, closure must not be claimed while a branch is still at a solver limit.
    if (quality.get("envelopeClosed").getAsBoolean()) {
      assertTrue(!quality.get("envelopeTermination").getAsString().contains("PRESSURE_CEILING"), result.toString());
      assertTrue(quality.get("envelopeCriticalPointCount").getAsInt() >= 1, result.toString());
    }
  }

  @Test
  void calculatesReactiveAqueousPhaseWithImportedHeavyOil() {
    String definition = "{\"format\":\"neqsim-fluid\",\"version\":\"1.0\",\"eos\":\"SRK\","
        + "\"components\":[{\"name\":\"heavy-cut\",\"moleFraction\":1.0,"
        + "\"criticalTemperature\":800.0,\"criticalPressure\":12.0,"
        + "\"acentricFactor\":1.1,\"molarMass\":500.0,\"normalBoilingPoint\":700.0,"
        + "\"criticalVolume\":1.2,\"volumeShift\":0.1,\"parachor\":700.0,"
        + "\"density\":0.95,\"isPseudo\":true}],\"binaryInteractionCoefficients\":[]}";
    String request = "{\"components\":{\"water\":0.60,\"heavyOil\":0.25,\"ammonia\":0.08,"
        + "\"CO2\":0.07},\"fluidDefinition\":" + definition
        + ",\"temperatureC\":200.0,\"pressureBara\":15.0,\"temperatureMinC\":200.0,"
        + "\"temperatureMaxC\":200.0,\"pointCount\":11,\"reactive\":true}";

    JsonObject result = JsonParser.parseString(FieldFluidRunner.run(request)).getAsJsonObject();

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertTrue(result.get("model").getAsString().startsWith("Electrolyte-CPA"));
    assertTrue(result.getAsJsonObject("flash").get("reactiveEquilibriumConverged").getAsBoolean());
    assertTrue(result.getAsJsonObject("speciation").get("enabled").getAsBoolean());
    assertTrue(
        result.getAsJsonObject("fluid").getAsJsonObject("composition").getAsJsonObject("overall").has("heavy-cut"));
  }
}
