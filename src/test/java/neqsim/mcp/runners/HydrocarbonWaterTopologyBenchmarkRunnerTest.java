package neqsim.mcp.runners;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.nio.file.Path;
import java.nio.file.Paths;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import org.junit.jupiter.api.Test;

class HydrocarbonWaterTopologyBenchmarkRunnerTest {
  private static final Path MANIFEST = Paths.get(
      "task_solve/2026-08-02_重做高含水重油2_hc相包络求解器/step1_scope_and_research/references/ma-2021-acs-omega/benchmark_manifest.json")
      .toAbsolutePath().normalize();

  @Test
  void independentlyVerifiesManifestArtifactsAndAcceptsCompleteTopology() {
    JsonObject result = HydrocarbonWaterTopologyBenchmarkRunner.run(request(true));

    assertEquals("success", result.get("status").getAsString());
    assertTrue(result.get("sourceContractVerified").getAsBoolean(), result.toString());
    assertTrue(result.get("topologyAccepted").getAsBoolean(), result.toString());
    assertFalse(result.get("numericComparisonEligible").getAsBoolean());
    assertTrue(result.get("numericBenchmarkPending").getAsBoolean());
    assertFalse(result.get("engineeringEligible").getAsBoolean());
    assertEquals(2, result.getAsJsonArray("artifactAudit").size());
  }

  @Test
  void reportsMissingGoToGowForCurrentOwToGowOnlyCandidate() {
    JsonObject result = HydrocarbonWaterTopologyBenchmarkRunner.run(request(false));

    assertTrue(result.get("sourceContractVerified").getAsBoolean(), result.toString());
    assertFalse(result.get("topologyAccepted").getAsBoolean());
    assertTrue(result.getAsJsonArray("missingBoundaryTopologies").toString().contains("GO->GOW:AQUEOUS"),
        result.toString());
    assertFalse(result.get("engineeringEligible").getAsBoolean());
  }

  private static Path benchmarkManifest() {
    try {
      return Paths.get(HydrocarbonWaterTopologyBenchmarkRunnerTest.class
          .getResource("/neqsim/mcp/runners/ma-2021-acs-omega/benchmark_manifest.json").toURI());
    } catch (Exception exception) {
      throw new IllegalStateException("benchmark manifest test resource is unavailable", exception);
    }
  }

  private static JsonObject request(boolean complete) {
    JsonObject request = new JsonObject();
    request.addProperty("referenceManifestPath", benchmarkManifest().toString());
    JsonObject candidate = new JsonObject();
    candidate.addProperty("identifier", complete ? "complete-test-network" : "current-srk-ow-gow-loop");
    candidate.addProperty("topologyEvidenceResolved", true);
    candidate.addProperty("internalQualityEligible", true);
    JsonArray regions = new JsonArray();
    regions.add("OW");
    regions.add("GOW");
    if (complete) {
      regions.add("GO");
    }
    candidate.add("observedRegions", regions);
    JsonArray boundaries = new JsonArray();
    boundaries.add("OW_TO_GOW");
    if (complete) {
      boundaries.add("GO_TO_GOW");
    }
    candidate.add("observedBoundaryTopologies", boundaries);
    request.add("candidate", candidate);
    return request;
  }
}
