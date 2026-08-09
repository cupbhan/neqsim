package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterRegularPressureBoundaryTracer.Attempt;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterRegularPressureBoundaryTracer.TerminationReason;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterRegularPressureBoundaryTracerTest extends NeqSimTest {

  @Test
  void rejectsInvalidSeedsAndControls() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    HydrocarbonWaterRegularPressureBoundaryTracer tracer = new HydrocarbonWaterRegularPressureBoundaryTracer(fluid,
        BoundaryFamily.GW_TO_GOW);

    assertThrows(IllegalArgumentException.class, () -> tracer.setStepControls(0.0, 0.001, 0.1, 4));
    assertThrows(IllegalArgumentException.class, () -> tracer.setStepControls(0.01, 0.02, 0.1, 4));
    assertThrows(IllegalArgumentException.class, () -> tracer.setCorrectionControls(10, 8.0, 0.04, 0.35, 30.0));
    assertThrows(IllegalArgumentException.class, () -> tracer.setDomainBounds(400.0, 300.0, 1.0, 100.0));
    assertThrows(IllegalArgumentException.class, () -> tracer.trace(null, null, 1));
  }

  @Test
  @Tag("slow")
  void fluidOneAddsAStableStrictPointWithPressureAsTheRegularParameterWithoutMutatingCaller() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double originalTemperature = fluid.getTemperature();
    double originalPressure = fluid.getPressure();
    AnchorPoint first = AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 10.0), fluid);
    AnchorPoint second = AnchorPoint.from(BoundaryFamily.GW_TO_GOW, solveFluidOnePoint(fluid, 20.0), fluid);

    HydrocarbonWaterRegularPressureBoundaryTracer.Result result = new HydrocarbonWaterRegularPressureBoundaryTracer(
        fluid, BoundaryFamily.GW_TO_GOW).setStepControls(0.04, 0.0025, 0.08, 6)
        .setCorrectionControls(11, 8.0, 0.04, 0.35, 30.0).trace(first, second, 1);

    assertTrue(result.hasCompletedRequestedPoints(), result.getFailureMessage());
    assertEquals(TerminationReason.REQUESTED_POINT_COUNT, result.getTerminationReason());
    assertEquals(3, result.getPoints().size());
    AnchorPoint accepted = result.getPoints().get(2);
    assertTrue(accepted.getPressureBara() > second.getPressureBara());
    assertTrue(accepted.getBoundaryRoot().isConverged());
    assertTrue(accepted.getGlobalStabilityResult().isAccepted());
    assertTrue(accepted.getRetainedFlashResidual() <= 1.0e-5);
    assertTrue(accepted.getTangentPlaneDistance() <= 1.0e-8);
    assertTrue(accepted.getStationarityResidual() <= 1.0e-5);
    Attempt acceptedAttempt = result.getAttempts().stream().filter(attempt -> attempt.getAcceptedAnchor() != null)
        .findFirst().orElse(null);
    assertNotNull(acceptedAttempt);
    assertTrue(acceptedAttempt.getMatchingBracketCount() >= 1);
    assertTrue(acceptedAttempt.getStrictOrdinaryRootCount() >= 1);
    assertEquals(null, acceptedAttempt.getFailureMessage());
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
