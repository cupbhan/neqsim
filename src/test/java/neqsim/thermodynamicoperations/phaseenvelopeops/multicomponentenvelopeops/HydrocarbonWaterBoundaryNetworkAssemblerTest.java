package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.EndpointEvidence;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.End;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.EvidenceSource;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.NetworkBranch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Real-fluid regression for physical endpoint attachment and global branch-network gates. */
class HydrocarbonWaterBoundaryNetworkAssemblerTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void oneStrictAnchorRemainsIsolatedByDefaultUntilTopologyIsClassified() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    List<AnchorPoint> anchors = Collections
        .singletonList(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 10.0), fluid));
    List<HydrocarbonWaterBoundaryAnchorDiscoverer.Branch> discovered = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).clusterCorrectedAnchors(anchors);

    HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .assembleBranches(discovered, 1).getBranches().get(0);

    assertTrue(branch.isIsolated());
    assertEquals(0, branch.getGeneratedSeedPointCount());
  }

  @Test
  @Tag("slow")
  void oneStrictAnchorGeneratesASecondGloballyStableContinuationSeed() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    List<AnchorPoint> anchors = Collections
        .singletonList(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 10.0), fluid));
    List<HydrocarbonWaterBoundaryAnchorDiscoverer.Branch> discovered = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).clusterCorrectedAnchors(anchors);

    HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .setCorrectorControls(80, 1.0e-8, 2.0e-5).setStepControls(0.25, 1.0e-4, 0.5, 16, 0.35)
        .setSingleAnchorContinuationEnabled(true).assembleBranches(discovered, 1).getBranches().get(0);

    assertFalse(branch.isIsolated());
    assertEquals(1, branch.getAnchorCount());
    assertEquals(1, branch.getGeneratedSeedPointCount());
    assertNotNull(branch.getBackwardTrace());
    assertNotNull(branch.getForwardTrace());
    assertTrue(branch.getEvidencePoints().size() >= 3);
  }

  @Test
  @Tag("slow")
  void sparseStrictAnchorsAreSubdividedWithoutRelaxingTheContinuityGate() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    List<AnchorPoint> anchors = new ArrayList<AnchorPoint>();
    anchors.add(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 10.0), fluid));
    anchors.add(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 50.0), fluid));
    List<HydrocarbonWaterBoundaryAnchorDiscoverer.Branch> discovered = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).setClusteringControls(2.01, 4.0, 0.5).clusterCorrectedAnchors(anchors);

    HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .setCorrectorControls(80, 1.0e-8, 2.0e-5).setStepControls(0.25, 1.0e-4, 0.5, 16, 0.35)
        .setInternalBridgeControls(0.75, 0.20, 12).assembleBranches(discovered, 1).getBranches().get(0);

    assertTrue(branch.getInternalBridgePointCount() > 0, "A pressure ratio of five requires strict bridge points");
    List<TwoToThreePhaseBoundaryQualityGate.EvidencePoint> points = branch.getEvidencePoints();
    for (int index = 1; index < points.size(); index++) {
      double logPressureJump = Math
          .abs(Math.log(points.get(index).getPressureBara() / points.get(index - 1).getPressureBara()));
      assertTrue(logPressureJump <= 0.75 + 1.0e-10,
          "Every assembled adjacent segment must satisfy the configured pressure-continuity limit");
    }
  }

  @Test
  @Tag("slow")
  void fluidOneConnectsDeclaredDomainExitToIndependentlySolvedCriticalEndpoint() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    List<AnchorPoint> anchors = new ArrayList<AnchorPoint>();
    anchors.add(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 10.0), fluid));
    anchors.add(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 20.0), fluid));
    List<HydrocarbonWaterBoundaryAnchorDiscoverer.Branch> discovered = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).setClusteringControls(2.01, 4.0, 0.5).clusterCorrectedAnchors(anchors);
    HydrocarbonWaterBoundaryBranchAssembler.Result assembly = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .setCorrectorControls(80, 1.0e-8, 2.0e-5).setStepControls(1.0, 1.0e-4, 1.0, 16, 0.35)
        .assembleBranches(discovered, 400);
    assertEquals(1, assembly.getBranches().size());
    HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch = assembly.getBranches().get(0);
    assertFalse(branch.getForwardTrace().hasCompletedRequestedPoints());
    assertFalse(branch.getBackwardTrace().hasCompletedRequestedPoints());

    TwoToThreePhaseArcLengthCorrector.State forwardEnd = branch.getForwardTrace().getPoints()
        .get(branch.getForwardTrace().getPoints().size() - 1);
    HydrocarbonWaterCriticalEndpointSolver.Result critical = new HydrocarbonWaterCriticalEndpointSolver(fluid,
        CandidatePhase.GAS, CandidatePhase.AQUEOUS, CandidatePhase.OIL, CandidatePhase.GAS)
        .setNumericalControls(12, 1.0e-6, 1.0e-9, 1.0e-3, 2.0e-4, 1.0e-7, 1.0e-4).solve(forwardEnd);
    assertTrue(critical.isPhysicalEndpoint(), critical.getFailureMessage());

    EndpointEvidence criticalEvidence = EndpointEvidence.criticalEndpoint("fluid1-gw-gow-cep", BoundaryFamily.GW_TO_GOW,
        critical);
    EndpointEvidence domainEvidence = EndpointEvidence.domainExit("fluid1-gw-gow-low-pressure-domain",
        BoundaryFamily.GW_TO_GOW, branch.getBackwardTrace(), 50.0, 2500.0, 1.0e-6, 1.0e6, 1.0e-5);
    HydrocarbonWaterBoundaryNetworkAssembler.Result network = new HydrocarbonWaterBoundaryNetworkAssembler()
        .setAttachmentTolerances(0.05, 0.25, 1.0e-6, 1.0e-6)
        .assemble(assembly, Arrays.asList(criticalEvidence, domainEvidence));

    assertTrue(network.isEngineeringEligible(), network.getViolations().toString());
    assertNotNull(network.getQualityReport());
    assertTrue(network.getQualityReport().isAccepted());
    assertEquals(1, network.getBranches().size());
    NetworkBranch networkBranch = network.getBranches().get(0);
    assertEquals(SpecialPointType.DOMAIN_EXIT, networkBranch.getStartAttachment().getType());
    assertEquals(EvidenceSource.DOMAIN_BOUND_CLASSIFIER, networkBranch.getStartAttachment().getSource());
    assertEquals(SpecialPointType.CRITICAL_END_POINT, networkBranch.getEndAttachment().getType());
    assertEquals(EvidenceSource.CRITICAL_ENDPOINT_SOLVER, networkBranch.getEndAttachment().getSource());
    assertEquals(206.97619256, networkBranch.getEndAttachment().getTemperatureK(), 2.0e-4);
    assertEquals(36.61109242, networkBranch.getEndAttachment().getPressureBara(), 2.0e-4);
    assertEquals(1, network.getEndpoints().get(0).getAttachmentCount());
    assertEquals(1, network.getEndpoints().get(1).getAttachmentCount());
    assertEquals(0, criticalEvidence.getAttachmentCount(),
        "Network assembly must not mutate caller-owned endpoint evidence");
  }

  @Test
  @Tag("slow")
  void fluidTwoAttachesRefinedOilSpinodalButKeepsOtherEndDiagnostic() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperatures = new double[] { 423.15, 448.15, 473.15, 493.15, 503.15, 508.15, 513.15, 523.15, 548.15,
        573.15 };
    double[] pressures = new double[] { 250.0, 260.0, 270.0, 275.0 };
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).discoverFromStableRegionTransitions(temperatures, pressures)
        .getDiscovery();
    HydrocarbonWaterBoundaryBranchAssembler.Result assembly = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .setCorrectorControls(80, 1.0e-8, 2.0e-5).setStepControls(0.25, 1.0e-4, 0.5, 16, 0.35).assemble(discovery, 100);

    HydrocarbonWaterBoundaryNetworkAssembler.Result network = new HydrocarbonWaterBoundaryNetworkAssembler()
        .assemble(assembly, Collections.<EndpointEvidence>emptyList());
    assertFalse(network.isEngineeringEligible(), "The unobserved opposite endpoint must keep the network diagnostic");
    assertEquals(1, network.getBranches().size());
    NetworkBranch branch = network.getBranches().get(0);
    assertNotNull(branch.getEndAttachment());
    assertEquals(SpecialPointType.RETAINED_PHASE_SPINODAL, branch.getEndAttachment().getType());
    assertEquals(EvidenceSource.TERMINATION_CLASSIFIER, branch.getEndAttachment().getSource());
    assertEquals(504.366971105, branch.getEndAttachment().getTemperatureK(), 2.0e-3);
    assertEquals(275.169677734, branch.getEndAttachment().getPressureBara(), 2.0e-3);
    assertEquals(1, network.getEndpoints().stream()
        .filter(endpoint -> endpoint.getType() == SpecialPointType.RETAINED_PHASE_SPINODAL).count());
    assertTrue(network.getViolations().stream()
        .anyMatch(violation -> violation.endsWith(":START") && violation.startsWith("UNATTACHED_BRANCH_ENDPOINT:")));
  }

  @Test
  @Tag("slow")
  void targetBranchMergeTrimsTheExcludedTailBeforeTheQualityGate() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    List<AnchorPoint> anchors = new ArrayList<AnchorPoint>();
    anchors.add(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 10.0), fluid));
    anchors.add(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 20.0), fluid));
    List<HydrocarbonWaterBoundaryAnchorDiscoverer.Branch> discovered = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).setClusteringControls(2.01, 4.0, 0.5).clusterCorrectedAnchors(anchors);
    HydrocarbonWaterBoundaryBranchAssembler.Result assembly = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .setCorrectorControls(80, 1.0e-8, 2.0e-5).setStepControls(0.25, 1.0e-4, 0.5, 16, 0.35)
        .assembleBranches(discovered, 8);
    List<TwoToThreePhaseBoundaryQualityGate.EvidencePoint> original = assembly.getBranches().get(0).getEvidencePoints();
    int mergeIndex = original.size() - 3;
    EndpointEvidence merge = EndpointEvidence.targetBranchMerge("synthetic-proven-target-merge",
        BoundaryFamily.GW_TO_GOW, original.get(mergeIndex), original.get(mergeIndex + 1), End.END);

    HydrocarbonWaterBoundaryNetworkAssembler.Result network = new HydrocarbonWaterBoundaryNetworkAssembler()
        .assemble(assembly, Collections.singletonList(merge));

    NetworkBranch branch = network.getBranches().get(0);
    assertNotNull(branch.getEndAttachment());
    assertEquals(SpecialPointType.TARGET_BRANCH_MERGE, branch.getEndAttachment().getType());
    assertEquals(EvidenceSource.SECONDARY_STATIONARY_BRANCH_TRACKER, branch.getEndAttachment().getSource());
    assertEquals(2, branch.getTrimmedPointCount());
    assertEquals(mergeIndex + 1, branch.getPoints().size());
    assertEquals(original.get(mergeIndex).getTemperatureK(),
        branch.getPoints().get(branch.getPoints().size() - 1).getTemperatureK(), 1.0e-12);
    assertEquals(1, network.getEndpoints().stream().filter(
        endpoint -> endpoint.getType() == SpecialPointType.TARGET_BRANCH_MERGE && endpoint.getAttachmentCount() == 1)
        .count());
  }

  private static TwoToThreePhaseBoundaryPointSolver.Result solveFluidOnePoint(SystemInterface fluid,
      double pressureBara) {
    IncipientPhaseBoundaryPointSolver.Result aqueousDew = new IncipientPhaseBoundaryPointSolver(fluid,
        CandidatePhase.GAS, CandidatePhase.AQUEOUS).solve(pressureBara, 200.0, 450.0);
    assertTrue(aqueousDew.isConverged(), aqueousDew.getFailureMessage());
    TwoToThreePhaseBoundaryPointSolver.Result point = new TwoToThreePhaseBoundaryPointSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.AQUEOUS, CandidatePhase.OIL).setNumericalControls(36, 80, 1.0e-5, 1.0e-8).solve(pressureBara,
            150.0, aqueousDew.getTemperatureK() - 0.25, 1.0 - 1.0e-6, overallComposition(fluid),
            aqueousDew.getIncipientComposition());
    assertTrue(point.isConverged(), point.getFailureMessage());
    return point;
  }

  private static double[] overallComposition(SystemInterface fluid) {
    double[] composition = new double[fluid.getPhase(0).getNumberOfComponents()];
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = fluid.getPhase(0).getComponent(componentIndex).getz();
    }
    return composition;
  }
}
