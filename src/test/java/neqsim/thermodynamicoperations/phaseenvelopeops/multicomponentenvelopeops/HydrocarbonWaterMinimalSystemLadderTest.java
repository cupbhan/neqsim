package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkPenelouxEos;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBranchSeedScanner.Seed;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/** Minimal composition ladder for reproducible hydrocarbon-water topology development. */
class HydrocarbonWaterMinimalSystemLadderTest extends NeqSimTest {
  private static final double[] TEMPERATURES_K = new double[] {180.0, 220.0, 260.0, 300.0, 350.0, 425.0};
  private static final double[] PRESSURES_BARA = new double[] {1.0, 5.0, 10.0, 25.0, 50.0, 100.0};

  @Test
  @Tag("slow")
  void stableFlashExposesGasOilWaterSeedsAtEveryCompositionLevel() {
    for (SystemInterface fluid : fluids()) {
      HydrocarbonWaterBranchSeedScanner.Result scan = new HydrocarbonWaterBranchSeedScanner(fluid,
          CandidatePhase.AQUEOUS).setCompositionDeduplicationTolerance(1.0e-8).scan(TEMPERATURES_K, PRESSURES_BARA);
      long gasOilWaterSeeds = scan.getSeeds().stream().filter(seed -> seed.hasPhase(CandidatePhase.GAS))
          .filter(seed -> seed.hasPhase(CandidatePhase.OIL)).filter(seed -> seed.hasPhase(CandidatePhase.AQUEOUS))
          .count();
      assertTrue(gasOilWaterSeeds > 0,
          fluid.getNumberOfComponents() + "-component ladder level returned no GOW seed; seeds="
              + scan.getSeeds().size() + " failures=" + scan.getFailures().size());
    }
  }

  @Test
  @Tag("slow")
  void oneStrictBoundaryRootPerLevelReplaysInAnIndependentCorrector() {
    for (SystemInterface fluid : fluids()) {
      HydrocarbonWaterBranchSeedScanner.Result scan = new HydrocarbonWaterBranchSeedScanner(fluid,
          CandidatePhase.AQUEOUS).setCompositionDeduplicationTolerance(1.0e-8).scan(TEMPERATURES_K, PRESSURES_BARA);
      List<Seed> gasOilWaterSeeds = new ArrayList<Seed>();
      for (Seed seed : scan.getSeeds()) {
        if (seed.hasPhase(CandidatePhase.GAS) && seed.hasPhase(CandidatePhase.OIL)
            && seed.hasPhase(CandidatePhase.AQUEOUS)) {
          gasOilWaterSeeds.add(seed);
        }
      }
      boolean reproducibleRootFound = false;
      String lastFailure = "no GOW seed";
      for (Seed seed : gasOilWaterSeeds) {
        TwoToThreePhaseBoundaryPointSolver.RootSet roots = new HydrocarbonWaterBoundarySeedCorrector(fluid,
            CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS).setNumericalControls(24, 60, 1.0e-5, 1.0e-8)
            .correctAtSeedPressure(seed, TEMPERATURES_K[0], TEMPERATURES_K[TEMPERATURES_K.length - 1]);
        lastFailure = roots.getFailureMessage();
        for (TwoToThreePhaseBoundaryPointSolver.Result root : roots.getRoots()) {
          TwoToThreePhaseArcLengthCorrector.State state = TwoToThreePhaseArcLengthCorrector.State.from(root);
          TwoToThreePhaseArcLengthCorrector.EquationResiduals replay = new TwoToThreePhaseArcLengthCorrector(fluid,
              CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS).replayEquationResiduals(state);
          if (replay.getThermodynamicMaximumResidual() <= 1.0e-8) {
            reproducibleRootFound = true;
            break;
          }
          lastFailure = "independent replay residual=" + replay.getThermodynamicMaximumResidual();
        }
        if (reproducibleRootFound) {
          break;
        }
      }
      assertTrue(reproducibleRootFound, fluid.getNumberOfComponents()
          + "-component ladder level returned no reproducible GO->GOW root; " + lastFailure);
    }
  }

  private static List<SystemInterface> fluids() {
    List<SystemInterface> fluids = new ArrayList<SystemInterface>();
    fluids.add(fluid(false, false));
    fluids.add(fluid(true, false));
    fluids.add(fluid(true, true));
    return fluids;
  }

  private static SystemInterface fluid(boolean carbonDioxide, boolean nitrogen) {
    SystemInterface fluid = new SystemSrkPenelouxEos(298.15, 10.0);
    fluid.addComponent("methane", 0.20);
    fluid.addComponent("n-hexane", 0.20);
    fluid.addComponent("water", carbonDioxide ? nitrogen ? 0.53 : 0.55 : 0.60);
    if (carbonDioxide) {
      fluid.addComponent("CO2", 0.05);
    }
    if (nitrogen) {
      fluid.addComponent("nitrogen", 0.02);
    }
    fluid.createDatabase(true);
    fluid.setMixingRule("HV", "NRTL");
    fluid.setMultiPhaseCheck(true);
    fluid.init(0);
    return fluid;
  }
}
