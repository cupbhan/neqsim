package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

/**
 * Checks composition derivatives for extreme partition coefficients at incipient phases.
 *
 * @author NeqSim developers
 * @version 1.0
 */
class EnvelopeDerivativeStabilityTest {
  /** A bounded-domain truncation disappears only when the traced maximum lies below the ceiling. */
  @Test
  void explicitPressureExtensionIsIndependentOfHigherCeiling() {
    SystemInterface fluid = new SystemPrEos(300.0, 1.0);
    fluid.addComponent("methane", 0.8);
    fluid.addComponent("ethane", 0.06);
    fluid.addComponent("propane", 0.03);
    double[] masses = {0.1, 0.15, 0.25, 0.4};
    double[] densities = {0.75, 0.8, 0.85, 0.9};
    double[] fractions = {0.5, 0.3, 0.15, 0.05};
    for (int i = 0; i < masses.length; i++) {
      fluid.addTBPfraction("fraction" + i, 0.11 * fractions[i], masses[i], densities[i]);
    }
    fluid.setMixingRule("classic");
    TwoHydrocarbonPhaseEnvelopeSolver.Result limited = trace(fluid, 400.0);
    TwoHydrocarbonPhaseEnvelopeSolver.Result extended = trace(fluid, 500.0);
    TwoHydrocarbonPhaseEnvelopeSolver.Result control = trace(fluid, 600.0);
    assertFalse(limited.isEnvelopeClosed());
    assertTrue(limited.getTerminationSummary().contains("PRESSURE_CEILING"));
    assertTrue(extended.isEnvelopeClosed(), extended.getTerminationSummary());
    assertTrue(control.isEnvelopeClosed(), control.getTerminationSummary());
    assertTrue(extended.getCricondenbar()[1] > 400.0 && extended.getCricondenbar()[1] < 450.0);
    assertEquals(extended.getCricondenbar()[1], control.getCricondenbar()[1], 0.0);
    assertEquals(extended.getCricondentherm()[0], control.getCricondentherm()[0], 0.0);
    assertEquals(extended.getPointCount(), control.getPointCount());
    assertEquals(300.0, fluid.getTemperature(), 0.0);
    assertEquals(1.0, fluid.getPressure(), 0.0);
  }

  /**
   * Runs a boundary trace without changing continuation tolerances or the caller.
   *
   * @param fluid original public synthetic feed
   * @param ceiling pressure termination limit in bara
   * @return traced boundary and diagnostic status
   */
  private TwoHydrocarbonPhaseEnvelopeSolver.Result trace(SystemInterface fluid, double ceiling) {
    return new TwoHydrocarbonPhaseEnvelopeSolver(fluid).setPressureRange(0.1, ceiling)
        .setMaximumContinuationIterations(1800).setStabilityAnalysisEnabled(false).setTopologyFallbackEnabled(false)
        .solve();
  }

  /** The material-balance Jacobian must retain small nonzero composition derivatives. */
  @Test
  void extremePartitionCoefficientPreservesDerivative() {
    SystemInterface fluid = new SystemPrEos(200.0, 30.0);
    fluid.addComponent("methane", 0.5);
    fluid.addComponent("ethane", 0.5);
    fluid.setMixingRule("classic");
    fluid.init(0);
    fluid.setBeta(1.0 - 1.0e-10);
    for (double partition : new double[] {1.0e-19, 1.0e14, 1.0e20}) {
      for (int phase = 0; phase < 2; phase++) {
        fluid.getPhase(phase).getComponent(0).setK(partition);
        fluid.getPhase(phase).getComponent(1).setK(1.0);
      }
      SysNewtonRhapsonPhaseEnvelope solver = new SysNewtonRhapsonPhaseEnvelope(fluid, 2, 2);
      solver.init();
      solver.setJac();
      double beta = fluid.getBeta();
      double denominator = (1.0 - beta) + beta * partition;
      double expected = 0.5 * partition / (denominator * denominator);
      assertEquals(expected, solver.Jac.get(2, 0), Math.abs(expected) * 1.0e-12,
          "Derivative must be algebraically consistent with the stable Rachford-Rice denominator");
    }
  }
}
