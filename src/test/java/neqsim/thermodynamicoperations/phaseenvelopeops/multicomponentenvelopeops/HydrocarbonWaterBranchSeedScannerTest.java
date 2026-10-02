package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterBranchSeedScannerTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void scanPreservesTemplateAndReturnsSameStatePhaseInventory() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double originalTemperatureK = fluid.getTemperature();
    double originalPressureBara = fluid.getPressure();

    HydrocarbonWaterBranchSeedScanner.Result result = new HydrocarbonWaterBranchSeedScanner(fluid,
        CandidatePhase.AQUEOUS).scan(new double[] {213.0, 214.0}, new double[] {65.0});

    assertEquals(2, result.getEvaluatedPointCount());
    assertEquals(originalTemperatureK, fluid.getTemperature(), 0.0);
    assertEquals(originalPressureBara, fluid.getPressure(), 0.0);
    assertTrue(result.getFailures().stream().anyMatch(
        failure -> failure.getTemperatureK() == 213.0 && "REQUESTED_PHASE_ABSENT".equals(failure.getReason())));
    HydrocarbonWaterBranchSeedScanner.Seed seed = result.getSeeds().stream()
        .filter(candidate -> candidate.getTemperatureK() == 214.0).findFirst().orElseThrow();
    assertTrue(seed.hasPhase(CandidatePhase.GAS));
    assertTrue(seed.hasPhase(CandidatePhase.OIL));
    assertTrue(seed.hasPhase(CandidatePhase.AQUEOUS));
    assertEquals(1.0, seed.getPhaseFraction(CandidatePhase.GAS) + seed.getPhaseFraction(CandidatePhase.OIL)
        + seed.getPhaseFraction(CandidatePhase.AQUEOUS), 1.0e-10);
    for (int componentIndex = 0; componentIndex < fluid.getPhase(0).getNumberOfComponents(); componentIndex++) {
      double reconstructed = seed.getPhaseFraction(CandidatePhase.GAS)
          * seed.getPhaseComposition(CandidatePhase.GAS)[componentIndex]
          + seed.getPhaseFraction(CandidatePhase.OIL) * seed.getPhaseComposition(CandidatePhase.OIL)[componentIndex]
          + seed.getPhaseFraction(CandidatePhase.AQUEOUS)
              * seed.getPhaseComposition(CandidatePhase.AQUEOUS)[componentIndex];
      assertEquals(fluid.getPhase(0).getComponent(componentIndex).getz(), reconstructed, 1.0e-10);
    }
    double retainedFraction = seed.getPhaseFraction(CandidatePhase.OIL) + seed.getPhaseFraction(CandidatePhase.AQUEOUS);
    TwoToThreePhaseArcLengthCorrector.State oilAqueousOrdering = TwoToThreePhaseArcLengthCorrector.State.create(
        CandidatePhase.OIL, CandidatePhase.AQUEOUS, CandidatePhase.GAS, seed.getTemperatureK(), seed.getPressureBara(),
        seed.getPhaseFraction(CandidatePhase.OIL) / retainedFraction, seed.getPhaseComposition(CandidatePhase.OIL),
        seed.getPhaseComposition(CandidatePhase.AQUEOUS), seed.getPhaseComposition(CandidatePhase.GAS));
    IncipientPhaseCurvatureAnalyzer.Result curvature = new IncipientPhaseCurvatureAnalyzer(fluid, CandidatePhase.OIL,
        CandidatePhase.GAS).analyze(oilAqueousOrdering);
    assertTrue(Double.isFinite(curvature.getMinimumEigenvalue()),
        "OW->GOW ordering must support oil as the retained curvature reference");
  }

  @Test
  @Tag("slow")
  void fluidTwoCoarseScanExposesItsObservedMultiphaseSeedDomain() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperaturesK = new double[] {180.0, 220.0, 260.0, 300.0, 350.0, 425.0, 500.0, 650.0};
    double[] pressuresBara = new double[] {0.5, 2.0, 10.0, 50.0, 100.0, 250.0, 500.0};
    HydrocarbonWaterBranchSeedScanner.Result result = new HydrocarbonWaterBranchSeedScanner(fluid,
        CandidatePhase.AQUEOUS).scan(temperaturesK, pressuresBara);
    long threePhaseSeeds = result.getSeeds().stream().filter(seed -> seed.hasPhase(CandidatePhase.GAS)
        && seed.hasPhase(CandidatePhase.OIL) && seed.hasPhase(CandidatePhase.AQUEOUS)).count();
    for (HydrocarbonWaterBranchSeedScanner.Seed seed : result.getSeeds()) {
      System.out.printf("Fluid 2 grid seed T=%.3f K P=%.3f bara phases=%d G=%s O=%s W=%s%n", seed.getTemperatureK(),
          seed.getPressureBara(), seed.getStablePhaseCount(), seed.hasPhase(CandidatePhase.GAS),
          seed.hasPhase(CandidatePhase.OIL), seed.hasPhase(CandidatePhase.AQUEOUS));
    }
    System.out.printf("Fluid 2 grid inventory: evaluations=%d aqueousSeeds=%d threePhaseSeeds=%d failures=%d%n",
        result.getEvaluatedPointCount(), result.getSeeds().size(), threePhaseSeeds, result.getFailures().size());
    assertEquals(temperaturesK.length * pressuresBara.length, result.getEvaluatedPointCount());
    assertTrue(!result.getSeeds().isEmpty(), "The high-water Fluid 2 grid must expose aqueous composition seeds");
  }
}
