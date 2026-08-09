package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterBoundaryEventDetectorTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void fluidOneDetectsInteriorTemperatureAndPressureMaximaOnPseudoArcBranch() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double originalTemperature = fluid.getTemperature();
    double originalPressure = fluid.getPressure();
    TwoToThreePhaseBoundaryPointSolver.Result pointTen = solvePoint(fluid, 10.0);
    TwoToThreePhaseBoundaryPointSolver.Result pointTwenty = solvePoint(fluid, 20.0);
    TwoToThreePhasePseudoArcLengthTracer.Result trace = new TwoToThreePhasePseudoArcLengthTracer(fluid,
        CandidatePhase.GAS, CandidatePhase.AQUEOUS, CandidatePhase.OIL).setCorrectorControls(60, 1.0e-8, 2.0e-5)
        .setStepControls(0.25, 0.02, 0.25, 8).setMaximumCompositionJump(0.35)
        .trace(TwoToThreePhaseArcLengthCorrector.State.from(pointTen),
            TwoToThreePhaseArcLengthCorrector.State.from(pointTwenty), 50);
    assertTrue(trace.hasCompletedRequestedPoints(), trace.getFailureMessage());
    HydrocarbonWaterBoundaryEventDetector.Result events = new HydrocarbonWaterBoundaryEventDetector()
        .detect(trace.getPoints());
    HydrocarbonWaterBoundaryEventDetector.Event maximumTemperature = events.getEvent(SpecialPointType.CRICONDENTHERM);
    HydrocarbonWaterBoundaryEventDetector.Event maximumPressure = events.getEvent(SpecialPointType.CRICONDENBAR);
    assertNotNull(maximumTemperature);
    assertNotNull(maximumPressure);
    System.out.printf("Fluid 1 cricondentherm: T=%.9f K P=%.9f bara bracket=%d..%d quality=%.5g%n",
        maximumTemperature.getTemperatureK(), maximumTemperature.getPressureBara(),
        maximumTemperature.getLowerPointIndex(), maximumTemperature.getUpperPointIndex(),
        maximumTemperature.getQualityMeasure());
    System.out.printf("Fluid 1 cricondenbar: T=%.9f K P=%.9f bara bracket=%d..%d quality=%.5g%n",
        maximumPressure.getTemperatureK(), maximumPressure.getPressureBara(), maximumPressure.getLowerPointIndex(),
        maximumPressure.getUpperPointIndex(), maximumPressure.getQualityMeasure());
    assertTrue(maximumTemperature.isInterpolated());
    assertTrue(maximumTemperature.getTemperatureK() > 270.895 && maximumTemperature.getTemperatureK() < 270.900);
    assertTrue(maximumTemperature.getPressureBara() > 48.0 && maximumTemperature.getPressureBara() < 51.0);
    assertTrue(maximumPressure.isInterpolated());
    assertTrue(maximumPressure.getPressureBara() > 92.62 && maximumPressure.getPressureBara() < 92.64);
    assertTrue(maximumPressure.getTemperatureK() > 245.0 && maximumPressure.getTemperatureK() < 248.0);
    assertTrue(events.getMaximumTemperaturePointIndex() < events.getMaximumPressurePointIndex());
    assertNull(events.getEvent(SpecialPointType.CRITICAL_END_POINT),
        "A truncated non-coalesced branch must not be mislabeled as a critical endpoint");
    assertEquals(originalTemperature, fluid.getTemperature(), 0.0);
    assertEquals(originalPressure, fluid.getPressure(), 0.0);
  }

  private static TwoToThreePhaseBoundaryPointSolver.Result solvePoint(SystemInterface fluid, double pressureBara) {
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
