package neqsim.mcp.runners;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.junit.jupiter.api.Test;

class HydrocarbonWaterBoundaryRegressionRunnerTest {

  @Test
  void serializesSameTopologySamePressureRegressionContract() {
    JsonObject result = HydrocarbonWaterBoundaryRegressionRunner
        .run(JsonParser.parseString(request(true, "SRK", true)).getAsJsonObject());

    assertEquals("success", result.get("status").getAsString());
    assertEquals("OW->GOW:GAS", result.get("boundaryTopology").getAsString());
    assertTrue(result.get("comparisonEligible").getAsBoolean(), result.toString());
    assertTrue(result.get("thresholdsFrozen").getAsBoolean(), result.toString());
    assertTrue(result.get("accepted").getAsBoolean(), result.toString());
    assertFalse(result.get("engineeringEligible").getAsBoolean(), result.toString());
    assertEquals(1.0, result.getAsJsonObject("metrics").getAsJsonObject("referenceToCandidate")
        .get("pressureCoverageFraction").getAsDouble(), 1.0e-12);
  }

  @Test
  void rejectsUnresolvedPvtsimThreeHcTopologyBeforeCalculatingErrors() {
    JsonObject result = HydrocarbonWaterBoundaryRegressionRunner
        .run(JsonParser.parseString(request(false, "SRK", true)).getAsJsonObject());

    assertFalse(result.get("comparisonEligible").getAsBoolean(), result.toString());
    assertFalse(result.has("metrics"), result.toString());
    assertTrue(result.getAsJsonArray("violations").toString().contains("REFERENCE_TOPOLOGY_UNRESOLVED"),
        result.toString());
  }

  @Test
  void rejectsDifferentEquationOfStateBeforeCalculatingErrors() {
    JsonObject result = HydrocarbonWaterBoundaryRegressionRunner
        .run(JsonParser.parseString(request(true, "PR78", true)).getAsJsonObject());

    assertFalse(result.get("comparisonEligible").getAsBoolean(), result.toString());
    assertTrue(result.getAsJsonArray("violations").toString().contains("EQUATION_OF_STATE_MISMATCH"),
        result.toString());
  }

  @Test
  void keepsMetricsPendingUntilAcceptanceThresholdsAreFrozen() {
    JsonObject result = HydrocarbonWaterBoundaryRegressionRunner
        .run(JsonParser.parseString(request(true, "SRK", false)).getAsJsonObject());

    assertTrue(result.get("comparisonEligible").getAsBoolean(), result.toString());
    assertFalse(result.get("thresholdsFrozen").getAsBoolean(), result.toString());
    assertTrue(result.get("benchmarkPending").getAsBoolean(), result.toString());
    assertFalse(result.get("accepted").getAsBoolean(), result.toString());
    assertTrue(result.has("metrics"), result.toString());
  }

  private static String request(boolean referenceTopologyResolved, String candidateEos, boolean includeCriteria) {
    String model = "\"modelFingerprint\":{\"equationOfState\":\"%s\","
        + "\"volumeTranslation\":\"Peneloux\",\"mixingRule\":\"Huron-Vidal\","
        + "\"compositionFingerprint\":\"fluid-22-sha256\"," + "\"parameterFingerprint\":\"srk-hv-parameter-sha256\"}";
    String points = "\"points\":[{\"temperatureK\":300.0,\"pressureBara\":100.0},"
        + "{\"temperatureK\":350.0,\"pressureBara\":200.0}," + "{\"temperatureK\":400.0,\"pressureBara\":100.0},"
        + "{\"temperatureK\":350.0,\"pressureBara\":50.0}," + "{\"temperatureK\":300.0,\"pressureBara\":100.0}]";
    String reference = "{\"identifier\":\"pvtsim-3-hc\",\"source\":\"PVTSIM\"," + String.format(model, "SRK")
        + ",\"boundaryTopology\":\"OW_TO_GOW\",\"topologyResolved\":" + referenceTopologyResolved
        + ",\"topologyEvidence\":\"independent phase-side flashes\"," + points
        + ",\"criticalPoint\":{\"temperatureK\":350.0,\"pressureBara\":200.0}}";
    String candidatePoints = points.replace("300.0", "302.0").replace("350.0", "352.0").replace("400.0", "402.0");
    String candidate = "{\"identifier\":\"neqsim-loop\",\"source\":\"NEQSIM\"," + String.format(model, candidateEos)
        + ",\"boundaryTopology\":\"OW_TO_GOW\",\"topologyResolved\":true,"
        + "\"topologyEvidence\":\"global stability and closed-loop gate\"," + candidatePoints
        + ",\"criticalPoint\":{\"temperatureK\":352.0,\"pressureBara\":200.0}}";
    String criteria = includeCriteria
        ? ",\"criteria\":{\"minimumReferenceCoverageFraction\":1.0," + "\"minimumCandidateCoverageFraction\":1.0,"
            + "\"maximumRelativeTemperatureDifference\":0.01," + "\"maximumRmsRelativeTemperatureDifference\":0.01,"
            + "\"maximumP95RelativeTemperatureDifference\":0.01,"
            + "\"maximumSymmetricRelativeTemperatureDifference\":0.01," + "\"criticalPointRequired\":true,"
            + "\"maximumCriticalTemperatureDifference\":0.01," + "\"maximumCriticalPressureDifference\":0.01}"
        : "";
    return "{\"reference\":" + reference + ",\"candidate\":" + candidate + criteria + "}";
  }
}
