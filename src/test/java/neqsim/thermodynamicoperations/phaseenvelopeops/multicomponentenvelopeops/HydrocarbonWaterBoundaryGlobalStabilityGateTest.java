package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterBoundaryGlobalStabilityGateTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void acceptsIndependentlyKnownFluidOneBoundaryRoot() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    TwoToThreePhaseBoundaryPointSolver.Result root = fluidOneRootAtTenBara(fluid);

    HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = new HydrocarbonWaterBoundaryGlobalStabilityGate(
        fluid).evaluate(root);

    print("Fluid 1 physical root", stability);
    assertTrue(stability.isTargetMatched(), stability.getFailureMessage());
    assertTrue(stability.isAccepted(), stability.getFailureMessage());
  }

  @Test
  @Tag("slow")
  void rejectsFluidTwoMetastableRootsInsideStableThreePhaseRegion() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .setCorrectionControls(32, 60, 2.0e-5, 1.0e-8)
        .discover(new double[] {180.0, 220.0, 260.0, 300.0, 350.0, 425.0, 500.0}, new double[] {50.0}, 100.0, 650.0);

    int evaluated = 0;
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.EndpointCandidate candidate : discovery.getEndpointCandidates()) {
      HydrocarbonWaterBoundaryEndpointClassifier.Result classification = candidate.getClassification();
      if (classification
          .getClassification() != HydrocarbonWaterBoundaryEndpointClassifier.Classification.METASTABLE_BOUNDARY_ROOT) {
        continue;
      }
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = classification.getGlobalStabilityResult();
      assertTrue(stability != null, "Every metastable classification must retain its global TPD evidence");
      print("Fluid 2 " + candidate.getFamily() + " at " + stability.getBoundaryRoot().getTemperatureK() + " K",
          stability);
      assertFalse(stability.isAccepted(),
          "A zero-TPD root inside the independently stable GOW region must not seed a boundary branch");
      evaluated++;
    }
    assertTrue(discovery.getBranches().isEmpty(), "Rejected metastable roots must never seed ordinary branches");
    assertTrue(evaluated >= 4, "The regression must preserve all observed Fluid 2 metastable-root evidence");
  }

  private static TwoToThreePhaseBoundaryPointSolver.Result fluidOneRootAtTenBara(SystemInterface fluid) {
    IncipientPhaseBoundaryPointSolver.Result aqueousDew = new IncipientPhaseBoundaryPointSolver(fluid,
        CandidatePhase.GAS, CandidatePhase.AQUEOUS).solve(10.0, 225.0, 425.0);
    assertTrue(aqueousDew.isConverged(), aqueousDew.getFailureMessage());
    TwoToThreePhaseBoundaryPointSolver.Result root = new TwoToThreePhaseBoundaryPointSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.AQUEOUS, CandidatePhase.OIL).setNumericalControls(36, 80, 1.0e-5, 1.0e-8).solve(10.0, 190.0,
            aqueousDew.getTemperatureK() - 0.25, 1.0 - 1.0e-6, overallComposition(fluid),
            aqueousDew.getIncipientComposition());
    assertTrue(root.isConverged(), root.getFailureMessage());
    return root;
  }

  private static double[] overallComposition(SystemInterface fluid) {
    double[] composition = new double[fluid.getPhase(0).getNumberOfComponents()];
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = fluid.getPhase(0).getComponent(componentIndex).getz();
    }
    return composition;
  }

  private static void print(String label, HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability) {
    System.out.printf("%s: accepted=%s target=%s targetDistance=%.8g minimumTPD=%.8g trials=%d failure=%s%n", label,
        stability.isAccepted(), stability.isTargetMatched(), stability.getTargetCompositionDistance(),
        stability.getMinimumNonTrivialTangentPlaneDistance(), stability.getTrials().size(),
        stability.getFailureMessage());
  }
}
