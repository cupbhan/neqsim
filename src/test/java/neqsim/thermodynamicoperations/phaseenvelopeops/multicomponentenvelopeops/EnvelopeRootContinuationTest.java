package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import neqsim.thermo.phase.PhasePrEos;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

/**
 * Public synthetic regressions for volume-root continuity in distributed heavy feeds.
 *
 * @author NeqSim developers
 * @version 1.0
 */
class EnvelopeRootContinuationTest {
  /** PR volume references preserve an outer root and can be explicitly cleared. */
  @Test
  void referencePreventsRootRelabellingJump() {
    SystemInterface fluid = new SystemPrEos(150.0, 10.0);
    fluid.addComponent("methane", 1.0);
    fluid.setMixingRule("classic");
    fluid.init(0);
    fluid.setBeta(0.5);
    fluid.init(3);
    PhasePrEos phase = (PhasePrEos) fluid.getPhase(0);
    double gasVolume = phase.getMolarVolume();
    phase.setVolumeRootReference(gasVolume);
    fluid.setPhaseType(0, PhaseType.LIQUID);
    fluid.init(3);
    assertEquals(gasVolume, phase.getMolarVolume(), gasVolume * 1.0e-10);
    assertTrue(phase.getdPdVTn() < 0.0, "Continuation must not select the mechanically unstable middle root");
    phase.setVolumeRootReference(Double.NaN);
    fluid.init(3);
    assertTrue(phase.getMolarVolume() < gasVolume * 0.5, "Clearing the reference restores ordinary root selection");
  }

  /** Requires a complete physical boundary without changing the caller or its composition. */
  @Test
  void distributedHeavyFeedCloses() {
    SystemInterface fluid = new SystemPrEos(300.0, 1.0);
    fluid.addComponent("methane", 0.75);
    fluid.addComponent("ethane", 0.06);
    fluid.addComponent("propane", 0.03);
    double[] masses = {0.1, 0.15, 0.25, 0.4};
    double[] densities = {0.75, 0.8, 0.85, 0.9};
    double[] fractions = {0.5, 0.3, 0.15, 0.05};
    for (int i = 0; i < masses.length; i++) {
      fluid.addTBPfraction("fraction" + i, 0.16 * fractions[i], masses[i], densities[i]);
    }
    fluid.setMixingRule("classic");
    TwoHydrocarbonPhaseEnvelopeSolver.Result result = new TwoHydrocarbonPhaseEnvelopeSolver(fluid)
        .setPressureRange(0.1, 400.0).setMaximumContinuationIterations(1800).setStabilityAnalysisEnabled(false)
        .setTopologyFallbackEnabled(false).solve();
    assertTrue(result.isEnvelopeClosed(), result.getTerminationSummary());
    assertTrue(result.getCriticalPointCount() >= 1);
    assertTrue(result.getPointCount() > 50);
    assertTrue(result.getCricondentherm()[0] > 600.0 && result.getCricondentherm()[0] < 700.0);
    assertTrue(result.getCricondenbar()[1] > 350.0 && result.getCricondenbar()[1] < 400.0);
    assertEquals(300.0, fluid.getTemperature(), 0.0);
    assertEquals(1.0, fluid.getPressure(), 0.0);
    assertEquals(7, fluid.getNumberOfComponents());
    assertEquals(0.75, fluid.getComponent("methane").getNumberOfmoles(), 1.0e-15);
  }
}
