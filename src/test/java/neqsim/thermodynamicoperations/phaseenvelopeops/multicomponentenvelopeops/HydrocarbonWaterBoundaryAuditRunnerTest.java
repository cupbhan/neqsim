package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.mcp.runners.HydrocarbonWaterBoundaryAuditRunner;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;

/** End-to-end MCP audit regression using the published Fluid 1 topology. */
class HydrocarbonWaterBoundaryAuditRunnerTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void fluidOnePseudoArcCheckpointRestartAddsStrictPointAndReturnsTwoReusableStates() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double[] temperaturesK = new double[] { 180.0, 200.0, 220.0, 240.0, 250.0, 260.0, 270.0, 280.0, 300.0 };
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).setStableScanControls(32, 0.25)
        .discoverFromStableRegionTransitions(temperaturesK, new double[] { 10.0, 20.0 }).getDiscovery();
    HydrocarbonWaterBoundaryAnchorDiscoverer.Branch branch = discovery.getBranches().get(0);

    JsonObject diagnostic = HydrocarbonWaterBoundaryAuditRunner.runPseudoArcRestartDiagnostic(fluid,
        branch.getPoints().get(0).toContinuationState(), branch.getPoints().get(1).toContinuationState(), 1, 0.25,
        1.0e-4, 2.0e-5, 170.0, 310.0, 1.0, 30.0);

    assertTrue(diagnostic.get("completedRequestedPoints").getAsBoolean(), diagnostic.toString());
    assertEquals(1, diagnostic.get("newStrictPointCount").getAsInt(), diagnostic.toString());
    JsonArray restartStates = diagnostic.getAsJsonObject("trace").getAsJsonArray("restartStates");
    assertEquals(2, restartStates.size(), diagnostic.toString());
    assertEquals(fluid.getPhase(0).getNumberOfComponents(),
        restartStates.get(0).getAsJsonObject().getAsJsonArray("phaseZeroComposition").size());

    JsonObject topology = HydrocarbonWaterBoundaryAuditRunner.runRestartStateTopologyDiagnostic(fluid,
        branch.getPoints().get(0).toContinuationState(), BoundaryFamily.GW_TO_GOW);
    assertEquals("REGULAR_BOUNDARY_POINT", topology.get("endpointClassification").getAsString(), topology.toString());
    assertTrue(topology.get("independentBranchSupport").getAsBoolean(), topology.toString());

    JsonObject seeds = HydrocarbonWaterBoundaryAuditRunner.runStrictBranchSeedDiagnostic(fluid, temperaturesK,
        new double[] { 10.0, 20.0 }, BoundaryFamily.GW_TO_GOW);
    assertEquals(1, seeds.get("branchCount").getAsInt(), seeds.toString());
    assertEquals(2, seeds.getAsJsonArray("branches").get(0).getAsJsonObject().getAsJsonArray("restartStates").size(),
        seeds.toString());

    HydrocarbonWaterBoundaryAnchorDiscoverer.Result denseDiscovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .discoverFromStableRegionTransitions(temperaturesK, new double[] { 10.0, 10.1, 10.2 }).getDiscovery();
    HydrocarbonWaterBoundaryAnchorDiscoverer.Branch denseBranch = denseDiscovery.getBranches().get(0);
    assertTrue(denseBranch.getPoints().size() >= 3, denseDiscovery.toString());
    JsonObject switched = HydrocarbonWaterBoundaryAuditRunner.runBranchSwitchDiagnostic(fluid,
        denseBranch.getPoints().get(0).toContinuationState(), denseBranch.getPoints().get(1).toContinuationState(),
        denseBranch.getPoints().get(2).toContinuationState(), 1, 0.25, 1.0e-4, 2.0e-5, 170.0, 310.0, 1.0, 30.0);
    assertTrue(switched.getAsJsonObject("seedQualityGate").get("passed").getAsBoolean(), switched.toString());
    assertTrue(switched.get("switchAccepted").getAsBoolean(), switched.toString());
    assertTrue(switched.getAsJsonObject("continuationQualityGate").get("passed").getAsBoolean(), switched.toString());
    assertEquals(1, switched.get("newStrictPointCount").getAsInt(), switched.toString());
  }

  @Test
  @Tag("slow")
  void fluidOneHybridDiagnosticSelectsRegularPressureForDenseAnchorsInBothDirections() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double[] temperaturesK = new double[] { 180.0, 200.0, 220.0, 240.0, 250.0, 260.0, 270.0, 280.0, 300.0 };

    JsonObject diagnostic = HydrocarbonWaterBoundaryAuditRunner.runHybridBranchDiagnostic(fluid, temperaturesK,
        new double[] { 10.0, 10.5, 11.0 }, BoundaryFamily.GW_TO_GOW, 1, 170.0, 310.0, 1.0, 30.0);

    assertEquals("success", diagnostic.get("status").getAsString(), diagnostic.toString());
    assertTrue(diagnostic.get("assembled").getAsBoolean(), diagnostic.toString());
    assertTrue(diagnostic.get("regularPressurePrimarySelected").getAsBoolean(), diagnostic.toString());
    assertEquals(3, diagnostic.get("anchorCount").getAsInt(), diagnostic.toString());
    assertEquals(2, diagnostic.get("newStrictPointCount").getAsInt(), diagnostic.toString());
    assertEquals(5, diagnostic.getAsJsonArray("points").size(), diagnostic.toString());
    assertTrue(
        diagnostic.getAsJsonObject("regularPressureBackwardTrace").get("completedRequestedPoints").getAsBoolean(),
        diagnostic.toString());
    assertTrue(diagnostic.getAsJsonObject("regularPressureForwardTrace").get("completedRequestedPoints").getAsBoolean(),
        diagnostic.toString());
    assertTrue(!diagnostic.get("engineeringEligible").getAsBoolean(), diagnostic.toString());
  }

  @Test
  @Tag("slow")
  void fluidOneRegularPressureDiagnosticAddsOneStrictPointWithoutRunningTheGlobalNetwork() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double[] temperaturesK = new double[] { 180.0, 200.0, 220.0, 240.0, 250.0, 260.0, 270.0, 280.0, 300.0 };

    JsonObject diagnostic = HydrocarbonWaterBoundaryAuditRunner.runRegularPressureDiagnostic(fluid, temperaturesK,
        new double[] { 10.0, 20.0 }, BoundaryFamily.GW_TO_GOW, 1);

    assertEquals("success", diagnostic.get("status").getAsString(), diagnostic.toString());
    assertTrue(diagnostic.get("completedRequestedPoints").getAsBoolean(), diagnostic.toString());
    assertEquals(1, diagnostic.get("newStrictPointCount").getAsInt(), diagnostic.toString());
    assertEquals(3, diagnostic.getAsJsonArray("points").size(), diagnostic.toString());
    assertEquals(3, diagnostic.getAsJsonArray("continuationStates").size(), diagnostic.toString());
    assertEquals(2, diagnostic.getAsJsonArray("restartStates").size(), diagnostic.toString());
    assertTrue(
        diagnostic.getAsJsonArray("points").get(2).getAsJsonObject().get("globalStabilityAccepted").getAsBoolean(),
        diagnostic.toString());
    assertTrue(diagnostic.getAsJsonArray("attempts").get(0).getAsJsonObject().get("accepted").getAsBoolean(),
        diagnostic.toString());
  }

  @Test
  @Tag("slow")
  void explicitStableBracketCorrectionRecomputesAdjacencyAndReturnsStrictReusableState() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    HydrocarbonWaterStableRegionTransitionScanner.Result scan = new HydrocarbonWaterStableRegionTransitionScanner(fluid)
        .scan(new double[] { 180.0, 200.0, 220.0, 240.0, 250.0, 260.0, 270.0, 280.0, 300.0 }, new double[] { 10.0 });
    HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket = scan.getBrackets().stream()
        .filter(candidate -> candidate.getFamily() == BoundaryFamily.GW_TO_GOW).findFirst().orElseThrow();

    JsonObject diagnostic = HydrocarbonWaterBoundaryAuditRunner.runStableBracketCorrectionDiagnostic(fluid,
        bracket.getLowerTemperatureState().getTemperatureK(), bracket.getUpperTemperatureState().getTemperatureK(),
        bracket.getPressureBara(), BoundaryFamily.GW_TO_GOW);

    assertEquals("success", diagnostic.get("status").getAsString(), diagnostic.toString());
    assertTrue(diagnostic.get("strictCorrectionConverged").getAsBoolean(), diagnostic.toString());
    assertTrue(diagnostic.get("ordinaryRootCount").getAsInt() >= 1, diagnostic.toString());
    assertTrue(diagnostic.getAsJsonArray("roots").get(0).getAsJsonObject().has("continuationState"),
        diagnostic.toString());
    assertTrue(!diagnostic.get("engineeringEligible").getAsBoolean(), diagnostic.toString());
  }

  @Test
  @Tag("slow")
  void modelSeededDiagnosticRecorrectsSerializedStateInsteadOfTrustingItsSourceModel() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double[] temperaturesK = new double[] { 180.0, 200.0, 220.0, 240.0, 250.0, 260.0, 270.0, 280.0, 300.0 };
    HydrocarbonWaterBoundaryAnchorDiscoverer.Branch branch = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .discoverFromStableRegionTransitions(temperaturesK, new double[] { 10.0, 20.0 }).getDiscovery().getBranches()
        .get(0);
    TwoToThreePhaseArcLengthCorrector.State seed = branch.getPoints().get(0).toContinuationState();

    JsonObject diagnostic = HydrocarbonWaterBoundaryAuditRunner.runModelSeededFixedPressureDiagnostic(fluid, seed,
        BoundaryFamily.GW_TO_GOW, seed.getTemperatureK() - 20.0, seed.getTemperatureK() + 20.0);

    assertEquals("success", diagnostic.get("status").getAsString(), diagnostic.toString());
    assertTrue(diagnostic.get("seedOnly").getAsBoolean(), diagnostic.toString());
    assertTrue(diagnostic.get("strictRootCount").getAsInt() >= 1, diagnostic.toString());
    assertTrue(diagnostic.get("ordinaryRootCount").getAsInt() >= 1, diagnostic.toString());
    assertTrue(!diagnostic.get("engineeringEligible").getAsBoolean(), diagnostic.toString());
    assertTrue(diagnostic.getAsJsonArray("roots").get(0).getAsJsonObject().has("continuationState"),
        diagnostic.toString());
  }

  @Test
  @Tag("slow")
  void fluidOneAuditClosesDomainToCepBranchButRejectsSparseDiscoveryCoverage() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double[] temperaturesK = new double[] { 180.0, 200.0, 220.0, 240.0, 250.0, 260.0, 270.0, 280.0, 300.0 };
    double[] pressuresBara = new double[] { 10.0, 20.0 };

    JsonObject audit = HydrocarbonWaterBoundaryAuditRunner.run(fluid, temperaturesK, pressuresBara, 100);
    JsonObject network = audit.getAsJsonObject("network");

    assertEquals("success", audit.get("status").getAsString(), audit.toString());
    assertTrue(network.get("attempted").getAsBoolean(), audit.toString());
    assertTrue(network.get("branchNetworkEligible").getAsBoolean(), audit.toString());
    assertTrue(!network.get("discoveryCoveragePassed").getAsBoolean(), audit.toString());
    assertTrue(!network.get("engineeringEligible").getAsBoolean(), audit.toString());
    assertTrue(!audit.get("engineeringEligible").getAsBoolean(), audit.toString());
    assertTrue(!audit.getAsJsonObject("discoveryCoverage").get("passed").getAsBoolean(), audit.toString());
    assertEquals(1, network.get("branchCount").getAsInt(), audit.toString());
    assertEquals(2, network.get("attachedEndpointCount").getAsInt(), audit.toString());
    JsonArray endpoints = network.getAsJsonArray("endpoints");
    assertTrue(containsEndpoint(endpoints, "DOMAIN_EXIT"), audit.toString());
    assertTrue(containsEndpoint(endpoints, "CRITICAL_END_POINT"), audit.toString());
    assertEquals(0, network.getAsJsonArray("violations").size(), audit.toString());
    assertTrue(network.getAsJsonObject("qualityGate").get("accepted").getAsBoolean(), audit.toString());
    assertTrue(network.getAsJsonArray("branches").get(0).getAsJsonObject().getAsJsonArray("points").size() > 40,
        audit.toString());
  }

  @Test
  @Tag("slow")
  void fluidTwoAuditClosesGoToGowFromDeclaredDomainToOilSpinodal() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperaturesK = new double[] { 196.75, 233.15, 273.15, 323.15, 373.15, 423.15, 473.15, 498.15, 523.15,
        548.15, 573.15, 623.15, 673.15, 723.15, 753.15 };
    double[] pressuresBara = new double[] { 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0, 150.0, 200.0, 250.0,
        260.0, 270.0, 275.0 };

    JsonObject audit = HydrocarbonWaterBoundaryAuditRunner.run(fluid, temperaturesK, pressuresBara, 200);
    JsonObject network = audit.getAsJsonObject("network");

    assertTrue(network.get("branchNetworkEligible").getAsBoolean(), audit.toString());
    assertTrue(network.get("discoveryCoveragePassed").getAsBoolean(), audit.toString());
    assertTrue(network.get("engineeringEligible").getAsBoolean(), audit.toString());
    JsonObject coverage = audit.getAsJsonObject("discoveryCoverage");
    assertTrue(coverage.get("passed").getAsBoolean(), audit.toString());
    assertTrue(coverage.getAsJsonObject("transitionBracketCounts").get("GO_TO_GOW").getAsInt() > 0, audit.toString());
    assertEquals(0, coverage.getAsJsonObject("transitionBracketCounts").get("GW_TO_GOW").getAsInt(), audit.toString());
    assertEquals(0, coverage.getAsJsonObject("transitionBracketCounts").get("OW_TO_GOW").getAsInt(), audit.toString());
    assertTrue(coverage.getAsJsonArray("observedStableRegions").toString().contains("GO"), audit.toString());
    assertTrue(coverage.getAsJsonArray("observedStableRegions").toString().contains("GOW"), audit.toString());
    assertEquals(1, network.get("branchCount").getAsInt(), audit.toString());
    assertEquals(2, network.get("attachedEndpointCount").getAsInt(), audit.toString());
    JsonArray endpoints = network.getAsJsonArray("endpoints");
    assertTrue(containsEndpoint(endpoints, "DOMAIN_EXIT"), audit.toString());
    assertTrue(containsEndpoint(endpoints, "RETAINED_PHASE_SPINODAL"), audit.toString());
    assertEquals(0, network.getAsJsonArray("violations").size(), audit.toString());
    JsonObject spinodal = endpoint(endpoints, "RETAINED_PHASE_SPINODAL");
    assertEquals(504.366971105, spinodal.get("temperatureK").getAsDouble(), 2.0e-3);
    assertEquals(275.169677734, spinodal.get("pressureBara").getAsDouble(), 2.0e-3);
  }

  private static boolean containsEndpoint(JsonArray endpoints, String type) {
    return endpoint(endpoints, type) != null;
  }

  private static JsonObject endpoint(JsonArray endpoints, String type) {
    for (int index = 0; index < endpoints.size(); index++) {
      JsonObject endpoint = endpoints.get(index).getAsJsonObject();
      if (type.equals(endpoint.get("type").getAsString()) && endpoint.get("attachmentCount").getAsInt() == 1) {
        return endpoint;
      }
    }
    return null;
  }
}
