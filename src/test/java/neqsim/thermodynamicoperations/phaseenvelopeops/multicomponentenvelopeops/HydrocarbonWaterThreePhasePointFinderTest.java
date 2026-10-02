package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterThreePhasePointFinderTest extends NeqSimTest {

  @Test
  void rejectsDuplicatePhaseFamilies() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    assertThrows(IllegalArgumentException.class, () -> new HydrocarbonWaterThreePhasePointFinder(fluid,
        CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.OIL));
  }

  @Test
  @Tag("slow")
  void fluidOneSearchReportsIndependentBoundaryEvidenceWithoutFabricatingPoint() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    double originalTemperature = fluid.getTemperature();
    double originalPressure = fluid.getPressure();
    HydrocarbonWaterThreePhasePointFinder.Result result = new HydrocarbonWaterThreePhasePointFinder(fluid,
        CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS)
        .setNumericalControls(9, 28, 24, 50, 1.0e-5, 1.0e-3, 8.0, 1.0e-8, 1.0e-7).find(0.5, 100.0, 150.0, 900.0);

    for (HydrocarbonWaterThreePhasePointFinder.Sample sample : result.getSamples()) {
      System.out.printf("three-phase scan P=%.8g bara usable=%s oilT=%.8g waterT=%.8g gap=%.8g%n",
          sample.getPressureBara(), sample.isUsable(), sample.getFirstBoundary().getTemperatureK(),
          sample.getSecondBoundary().getTemperatureK(), sample.getTemperatureGapK());
    }
    for (HydrocarbonWaterThreePhasePointFinder.CoupledTrial trial : result.getCoupledTrials()) {
      System.out.printf("three-phase candidate source=%s seedP=%.8g gap=%.8g coupled=%s%n", trial.getSource(),
          trial.getSeed().getPressureBara(), trial.getSeed().getTemperatureGapK(), trial.getPoint());
    }
    System.out.printf("Fluid 1 simultaneous-incidence result: physical=%s points=%d minimumGap=%.8g diagnostic=%s%n",
        result.hasPhysicalPoint(), result.getPhysicalPoints().size(), result.getMinimumTemperatureGapK(),
        result.getDiagnostic());

    assertTrue(result.getSamples().size() >= 9);
    for (ThreePhasePointSolver.Result point : result.getPhysicalPoints()) {
      assertTrue(point.isConverged());
      assertTrue(point.hasDistinctPhases());
      assertTrue(point.getMaximumResidual() <= 1.0e-7);
    }
    if (result.hasPhysicalPoint()) {
      assertEquals(null, result.getDiagnostic());
    } else {
      assertFalse(result.getDiagnostic() == null || result.getDiagnostic().isEmpty());
    }
    assertEquals(originalTemperature, fluid.getTemperature(), 0.0);
    assertEquals(originalPressure, fluid.getPressure(), 0.0);
  }

  @Test
  @Tag("slow")
  void diagnoseFluidOneOilStationaryPointAtFiftyBara() {
    SystemInterface source = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    for (double temperatureK : new double[] {150.0, 225.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0}) {
      SystemInterface grid = source.clone();
      grid.setTemperature(temperatureK);
      grid.setPressure(50.0);
      grid.setMultiPhaseCheck(true);
      grid.setMaxNumberOfPhases(3);
      new ThermodynamicOperations(grid).TPflash();
      grid.init(1);
      double[] oilSeed = null;
      StringBuilder phaseTypes = new StringBuilder();
      for (int phaseIndex = 0; phaseIndex < grid.getNumberOfPhases(); phaseIndex++) {
        phaseTypes.append(grid.getPhase(phaseIndex).getType()).append('/');
        if (grid.getPhase(phaseIndex).getType() == PhaseType.OIL
            || grid.getPhase(phaseIndex).getType() == PhaseType.LIQUID) {
          oilSeed = composition(grid, phaseIndex);
        }
      }
      SystemInterface gas = source.clone();
      gas.setNumberOfPhases(1);
      gas.setTemperature(temperatureK);
      gas.setPressure(50.0);
      gas.setPhaseType(0, PhaseType.GAS);
      for (int componentIndex = 0; componentIndex < gas.getPhase(0).getNumberOfComponents(); componentIndex++) {
        gas.getPhase(0).getComponent(componentIndex).setx(gas.getPhase(0).getComponent(componentIndex).getz());
      }
      gas.getPhase(0).normalize();
      gas.init(1);
      if (oilSeed == null) {
        oilSeed = composition(gas, 0);
      }
      IncipientPhaseStabilityAnalyzer.Candidate successive = new IncipientPhaseStabilityAnalyzer(gas)
          .setMaximumIterations(500).setDampingFactor(0.2).analyzeCandidate(CandidatePhase.OIL, oilSeed);
      IncipientPhaseStationaryPointSolver.Result newton = new IncipientPhaseStationaryPointSolver(gas,
          CandidatePhase.OIL).setNumericalControls(100, 1.0e-9, 2.0e-5).solve(oilSeed);
      int waterIndex = -1;
      for (int componentIndex = 0; componentIndex < gas.getPhase(0).getNumberOfComponents(); componentIndex++) {
        if (gas.getPhase(0).getComponent(componentIndex).getComponentName().equalsIgnoreCase("water")) {
          waterIndex = componentIndex;
          break;
        }
      }
      assertTrue(waterIndex >= 0);
      System.out.printf(
          "oil stationary P=50 T=%.2f grid=%s seedWater=%.6g SS(conv=%s trivial=%s physical=%s TPD=%.8g "
              + "r=%.3g water=%.6g) Newton(conv=%s trivial=%s physical=%s TPD=%.8g r=%.3g water=%.6g "
              + "cond=%.3g failure=%s)%n",
          temperatureK, phaseTypes, oilSeed[waterIndex], successive.isConverged(), successive.isTrivial(),
          successive.getPhase(), successive.getTangentPlaneDistance(), successive.getStationarityResidual(),
          successive.getComposition()[waterIndex], newton.isConverged(), newton.isTrivial(), newton.getPhysicalPhase(),
          newton.getTangentPlaneDistance(), newton.getStationarityResidual(), newton.getComposition()[waterIndex],
          newton.getJacobianConditionNumber(), newton.getFailureMessage());
      double[] dryOilSeed = oilSeed.clone();
      double dryTotal = 0.0;
      for (int componentIndex = 0; componentIndex < dryOilSeed.length; componentIndex++) {
        if (gas.getPhase(0).getComponent(componentIndex).getComponentName().equalsIgnoreCase("water")) {
          dryOilSeed[componentIndex] *= 1.0e-12;
        }
        dryTotal += dryOilSeed[componentIndex];
      }
      for (int componentIndex = 0; componentIndex < dryOilSeed.length; componentIndex++) {
        dryOilSeed[componentIndex] /= dryTotal;
      }
      IncipientPhaseStabilityAnalyzer.Candidate drySuccessive = new IncipientPhaseStabilityAnalyzer(gas)
          .setMaximumIterations(500).setDampingFactor(0.2).analyzeCandidate(CandidatePhase.OIL, dryOilSeed);
      IncipientPhaseStationaryPointSolver.Result dryNewton = new IncipientPhaseStationaryPointSolver(gas,
          CandidatePhase.OIL).setNumericalControls(100, 1.0e-9, 2.0e-5).solve(dryOilSeed);
      System.out.printf(
          "  dry seed SS(conv=%s trivial=%s physical=%s TPD=%.8g r=%.3g) "
              + "Newton(conv=%s trivial=%s physical=%s TPD=%.8g r=%.3g failure=%s)%n",
          drySuccessive.isConverged(), drySuccessive.isTrivial(), drySuccessive.getPhase(),
          drySuccessive.getTangentPlaneDistance(), drySuccessive.getStationarityResidual(), dryNewton.isConverged(),
          dryNewton.isTrivial(), dryNewton.getPhysicalPhase(), dryNewton.getTangentPlaneDistance(),
          dryNewton.getStationarityResidual(), dryNewton.getFailureMessage());
    }
  }

  @Test
  @Tag("slow")
  void diagnoseWaterCompositionHomotopyForOilStationaryPoint() {
    SystemInterface source = LindeloffMichelsenReferenceFluidTest.fluidOne(false);
    SystemInterface grid = source.clone();
    grid.setTemperature(225.0);
    grid.setPressure(50.0);
    grid.setMultiPhaseCheck(true);
    grid.setMaxNumberOfPhases(3);
    new ThermodynamicOperations(grid).TPflash();
    grid.init(1);
    double[] seed = null;
    for (int phaseIndex = 0; phaseIndex < grid.getNumberOfPhases(); phaseIndex++) {
      if (grid.getPhase(phaseIndex).getType() == PhaseType.OIL
          || grid.getPhase(phaseIndex).getType() == PhaseType.LIQUID) {
        seed = composition(grid, phaseIndex);
      }
    }
    assertTrue(seed != null);
    for (double waterScale : new double[] {1.0e-8, 1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 0.03, 0.1, 0.2, 0.4, 0.6, 0.8,
        1.0}) {
      SystemInterface gas = source.clone();
      gas.setNumberOfPhases(1);
      gas.setTemperature(225.0);
      gas.setPressure(50.0);
      gas.setPhaseType(0, PhaseType.GAS);
      int waterIndex = -1;
      for (int componentIndex = 0; componentIndex < gas.getPhase(0).getNumberOfComponents(); componentIndex++) {
        double value = gas.getPhase(0).getComponent(componentIndex).getz();
        if (gas.getPhase(0).getComponent(componentIndex).getComponentName().equalsIgnoreCase("water")) {
          waterIndex = componentIndex;
          value *= waterScale;
        }
        gas.getPhase(0).getComponent(componentIndex).setx(value);
      }
      gas.getPhase(0).normalize();
      gas.init(1);
      IncipientPhaseStationaryPointSolver.Result result = new IncipientPhaseStationaryPointSolver(gas,
          CandidatePhase.OIL).setNumericalControls(120, 1.0e-9, 2.0e-5).solve(seed);
      System.out.printf(
          "water homotopy scale=%.8g conv=%s trivial=%s physical=%s TPD=%.8g r=%.3g water=%.8g "
              + "cond=%.3g failure=%s%n",
          waterScale, result.isConverged(), result.isTrivial(), result.getPhysicalPhase(),
          result.getTangentPlaneDistance(), result.getStationarityResidual(), result.getComposition()[waterIndex],
          result.getJacobianConditionNumber(), result.getFailureMessage());
      if (result.isConverged() && !result.isTrivial() && result.getPhysicalPhase() == CandidatePhase.OIL) {
        seed = result.getComposition();
      }
    }
  }

  private static double[] composition(SystemInterface system, int phaseIndex) {
    double[] composition = new double[system.getPhase(phaseIndex).getNumberOfComponents()];
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = system.getPhase(phaseIndex).getComponent(componentIndex).getx();
    }
    return composition;
  }
}
