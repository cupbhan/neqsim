package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;
import org.junit.jupiter.api.Test;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Regression for loss of continuation direction through integer rounding.
 *
 * @author NeqSim developers
 * @version 1.0
 */
class EnvelopeSmallStepTest {
  /** Numerical interruption must not be advertised as a closed envelope. */
  @Test
  void numericalInterruptionIsNotClosure() {
    SystemInterface fluid = new SystemPrEos(300.0, 1.0);
    fluid.addComponent("methane", 0.85);
    fluid.addComponent("ethane", 0.07);
    fluid.addComponent("propane", 0.03);
    fluid.addTBPfraction("heavy", 0.05, 0.5, 0.9);
    fluid.setMixingRule("classic");
    TwoHydrocarbonPhaseEnvelopeSolver.Result result = new TwoHydrocarbonPhaseEnvelopeSolver(fluid)
        .setPressureRange(0.1, 2000.0).setMaximumContinuationIterations(1800).setStabilityAnalysisEnabled(false)
        .setTopologyFallbackEnabled(false).solve();
    assertTrue(result.getPointCount() > 0);
    assertTrue(result.getTerminationSummary().contains("NUMERICAL_FAILURE"));
    assertFalse(result.isEnvelopeClosed());
  }

  /** A pressure limit retains the converged portion without claiming a complete envelope. */
  @Test
  void pressureLimitRetainsPoints() {
    SystemInterface fluid = new SystemPrEos(250.0, 1.0);
    fluid.addComponent("methane", 0.8);
    fluid.addComponent("ethane", 0.15);
    fluid.addComponent("propane", 0.05);
    fluid.setMixingRule("classic");
    TwoHydrocarbonPhaseEnvelopeSolver.Result result = new TwoHydrocarbonPhaseEnvelopeSolver(fluid)
        .setPressureRange(0.1, 2.0).setStabilityAnalysisEnabled(false).setTopologyFallbackEnabled(false).solve();
    assertTrue(result.getTerminationSummary().contains("PRESSURE_CEILING"));
    assertTrue(result.getPointCount() > 0, "Retain converged points for bounded-domain diagnostics");
    assertTrue(result.hasFinitePhysicalValues());
    assertFalse(result.isEnvelopeClosed());
    assertTrue(Double.isNaN(result.getCricondenbar()[1]), "A truncated curve has no validated global extremum");
  }

  /**
   * Verifies a small nonzero step still produces a physical equilibrium point.
   *
   * @throws Exception if the initial dew point cannot be calculated
   */
  @Test
  void preservesSmallNonzeroStep() throws Exception {
    SystemInterface fluid = new SystemPrEos(220.0, 1.0);
    fluid.addComponent("methane", 0.8);
    fluid.addComponent("ethane", 0.15);
    fluid.addComponent("propane", 0.05);
    fluid.setMixingRule("classic");
    new ThermodynamicOperations(fluid).dewPointTemperatureFlash();
    SysNewtonRhapsonPhaseEnvelope solver = new SysNewtonRhapsonPhaseEnvelope(fluid, 2, 3);
    for (int point = 1; point <= 4; point++) {
      solver.calcInc(point);
      solver.solve(point);
    }
    solver.ds = 1.0e-9;
    solver.calcInc(5);
    assertTrue(Math.abs(solver.ds) > 0.0, "A nonzero step must not be rounded to zero");
    solver.solve(5);
    assertTrue(solver.fvec.norm2() < 1.0e-8, "Original equilibrium tolerance remains mandatory");
    assertTrue(fluid.getTemperature() > 0.0 && fluid.getPressure() > 0.0);
  }
}
