package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterIsolatedBoundaryRootClassifier.RootTopology;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterIsolatedBoundaryRootClassifierTest extends NeqSimTest {

  @Test
  void rejectsInvalidInputsAndControls() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    HydrocarbonWaterIsolatedBoundaryRootClassifier classifier = new HydrocarbonWaterIsolatedBoundaryRootClassifier(
        fluid);
    assertThrows(IllegalArgumentException.class, () -> classifier.classify(null));
    assertThrows(IllegalArgumentException.class, () -> classifier.setNeighborhoodControls(0, 0.1, 9, 0.1, 5.0));
    assertThrows(IllegalArgumentException.class, () -> classifier.setNeighborhoodControls(2, 0.1, 8, 0.1, 5.0));
    assertThrows(IllegalArgumentException.class, () -> classifier.setQualityControls(0.35, 1.1));
  }

  @Test
  @Tag("slow")
  void fluidOneOrdinaryRootObtainsIndependentLocalBranchSupportWithoutMutatingCaller() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double originalTemperature = fluid.getTemperature();
    double originalPressure = fluid.getPressure();
    AnchorPoint anchor = AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 10.0), fluid);

    HydrocarbonWaterIsolatedBoundaryRootClassifier.Result result = new HydrocarbonWaterIsolatedBoundaryRootClassifier(
        fluid).setNeighborhoodControls(2, 0.08, 17, 0.10, 12.0).classify(anchor);

    assertTrue(result.hasIndependentBranchSupport(), result.getDiagnostic());
    assertTrue(result.getTopology() == RootTopology.REGULAR_BRANCH_SUPPORTED
        || result.getTopology() == RootTopology.ONE_SIDED_BRANCH_SUPPORTED
        || result.getTopology() == RootTopology.PRESSURE_TURNING_POINT_CANDIDATE);
    assertFalse(result.getMatchingRoots().isEmpty());
    assertTrue(result.getSuccessfulGridFraction() >= 0.90);

    HydrocarbonWaterBoundaryAnchorDiscoverer discoverer = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid);
    List<AnchorPoint> independentlyRecovered = new ArrayList<AnchorPoint>();
    for (TwoToThreePhaseBoundaryPointSolver.Result root : result.getMatchingRoots()) {
      independentlyRecovered.add(AnchorPoint.from(BoundaryFamily.GW_TO_GOW, root, fluid));
    }
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result initial = discoverer.setCorrectionControls(24, 80, 1.0e-5, 1.0e-8)
        .discoverFromStableRegionTransitions(new double[] {230.0, 260.0}, new double[] {10.0}).getDiscovery();
    assertEquals(1, initial.getCorrectedAnchorCount());
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result augmented = discoverer.augmentWithIndependentAnchors(initial,
        independentlyRecovered);
    assertTrue(augmented.getCorrectedAnchorCount() > initial.getCorrectedAnchorCount());
    assertTrue(augmented.getBranches().stream().anyMatch(branch -> !branch.isIsolated()));
    assertEquals(originalTemperature, fluid.getTemperature(), 0.0);
    assertEquals(originalPressure, fluid.getPressure(), 0.0);
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
