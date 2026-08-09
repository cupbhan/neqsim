package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;

class HydrocarbonWaterBoundaryBranchAssemblerTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void densePressureAnchorsSelectRegularPressureBeforePseudoArcLength() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double[] temperaturesK = new double[] { 180.0, 200.0, 220.0, 240.0, 250.0, 260.0, 270.0, 280.0, 300.0 };
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).setStableScanControls(32, 0.25)
        .discoverFromStableRegionTransitions(temperaturesK, new double[] { 10.0, 10.5, 11.0 }).getDiscovery();

    HydrocarbonWaterBoundaryBranchAssembler.Result assembled = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .assemble(discovery, 1);

    assertEquals(1, assembled.getBranches().size());
    HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch = assembled.getBranches().get(0);
    assertEquals(3, branch.getAnchorCount());
    assertTrue(branch.getBackwardTrace() == null);
    assertTrue(branch.getForwardTrace() == null);
    assertNotNull(branch.getRegularBackwardTrace());
    assertNotNull(branch.getRegularForwardTrace());
    assertTrue(branch.getRegularBackwardTrace().hasCompletedRequestedPoints(),
        branch.getRegularBackwardTrace().getFailureMessage());
    assertTrue(branch.getRegularForwardTrace().hasCompletedRequestedPoints(),
        branch.getRegularForwardTrace().getFailureMessage());
    assertEquals(5, branch.getEvidencePoints().size());
  }

  @Test
  @Tag("slow")
  void fixedPressureFallbackPreservesARegularBranchWhenPseudoArcCannotStart() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double[] temperaturesK = new double[] { 180.0, 200.0, 220.0, 240.0, 250.0, 260.0, 270.0, 280.0, 300.0 };
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).setStableScanControls(32, 0.25)
        .discoverFromStableRegionTransitions(temperaturesK, new double[] { 10.0, 20.0 }).getDiscovery();

    HydrocarbonWaterBoundaryBranchAssembler.Result assembled = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .setCorrectorControls(1, 1.0e-14, 2.0e-5).setStepControls(0.25, 0.01, 0.5, 2, 0.35).assemble(discovery, 1);

    assertEquals(1, assembled.getBranches().size());
    HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch = assembled.getBranches().get(0);
    assertEquals(BoundaryFamily.GW_TO_GOW, branch.getFamily());
    assertTrue(branch.getBackwardTrace().getAcceptedCorrections().isEmpty());
    assertTrue(branch.getForwardTrace().getAcceptedCorrections().isEmpty());
    assertNotNull(branch.getRegularBackwardTrace());
    assertNotNull(branch.getRegularForwardTrace());
    assertTrue(branch.getRegularBackwardTrace().hasCompletedRequestedPoints(),
        branch.getRegularBackwardTrace().getFailureMessage());
    assertTrue(branch.getRegularForwardTrace().hasCompletedRequestedPoints(),
        branch.getRegularForwardTrace().getFailureMessage());
    assertEquals(4, branch.getEvidencePoints().size());
  }
}
