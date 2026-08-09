package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.EnumSet;
import java.util.Set;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterThreePhaseTopologyFinderTest extends NeqSimTest {

  @Test
  void rejectsInvalidGlobalControls() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    assertThrows(IllegalArgumentException.class, () -> new HydrocarbonWaterThreePhaseTopologyFinder(fluid)
        .setNumericalControls(2, 12, 8, 20, 1.0e-5, 1.0e-3, 5.0, 1.0e-8, 1.0e-7));
  }

  @Test
  @Tag("slow")
  void fluidOneSearchesEveryMotherTopologyWithoutFabricatingPoint() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double originalTemperature = fluid.getTemperature();
    double originalPressure = fluid.getPressure();
    HydrocarbonWaterThreePhaseTopologyFinder.Result result = new HydrocarbonWaterThreePhaseTopologyFinder(fluid)
        .setNumericalControls(3, 12, 8, 30, 1.0e-4, 1.0e-3, 10.0, 1.0e-8, 1.0e-7).find(1.0, 80.0, 150.0, 500.0);

    assertEquals(3, result.getSearches().size());
    Set<CandidatePhase> mothers = EnumSet.noneOf(CandidatePhase.class);
    for (HydrocarbonWaterThreePhaseTopologyFinder.MotherSearch search : result.getSearches()) {
      mothers.add(search.getMotherPhase());
      assertEquals(3, search.getResult().getSamples().size());
      System.out.printf("global three-phase mother=%s physical=%d trials=%d minGap=%.8g diagnostic=%s%n",
          search.getMotherPhase(), search.getResult().getPhysicalPoints().size(),
          search.getResult().getCoupledTrials().size(), search.getResult().getMinimumTemperatureGapK(),
          search.getResult().getDiagnostic());
    }
    assertEquals(EnumSet.allOf(CandidatePhase.class), mothers);
    for (HydrocarbonWaterThreePhaseTopologyFinder.PhysicalPoint physical : result.getPhysicalPoints()) {
      assertTrue(physical.getPoint().isConverged());
      assertTrue(physical.getPoint().hasDistinctPhases());
      assertTrue(physical.getPoint().getMaximumResidual() <= 1.0e-7);
    }
    assertEquals(originalTemperature, fluid.getTemperature(), 0.0);
    assertEquals(originalPressure, fluid.getPressure(), 0.0);
  }
}
