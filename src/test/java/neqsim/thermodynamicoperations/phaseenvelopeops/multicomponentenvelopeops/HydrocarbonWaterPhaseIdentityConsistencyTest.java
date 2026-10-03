package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import org.junit.jupiter.api.Test;

/**
 * Checks that research acceptance and EOS initialization use the same liquid-family convention.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
class HydrocarbonWaterPhaseIdentityConsistencyTest {
  /** Water-rich vapor must keep its evaluated gas identity in both stability solvers. */
  @Test
  void waterVaporIsNotRelabeledAqueousByItsComposition() {
    SystemInterface steam = new SystemSrkEos(600.0, 1.0);
    steam.addComponent("water", 1.0);
    steam.setMixingRule("classic");
    new neqsim.thermodynamicoperations.ThermodynamicOperations(steam).TPflash();
    steam.init(1);
    assertEquals(PhaseType.GAS, steam.getPhase(0).getType());
    IncipientPhaseStabilityAnalyzer.Candidate trial = new IncipientPhaseStabilityAnalyzer(steam)
        .analyzeCandidate(CandidatePhase.GAS, new double[] {1.0});
    assertTrue(trial.isConverged());
    assertEquals(CandidatePhase.GAS, trial.getPhase());
    IncipientPhaseStationaryPointSolver.Result stationary = new IncipientPhaseStationaryPointSolver(steam,
        CandidatePhase.GAS).solve(new double[] {1.0});
    assertTrue(stationary.isConverged());
    assertEquals(CandidatePhase.GAS, stationary.getPhysicalPhase());
  }

  /** A water-mole-rich oil must retain its family while a duplicated root is still rejected. */
  @Test
  void waterMoleMajorityDoesNotTurnHeavyLiquidIntoAnAqueousIdentityFailure() {
    SystemInterface fluid = fluid(0.6);
    SpecifiedMultiphaseFlashSolver solver = new SpecifiedMultiphaseFlashSolver(fluid, CandidatePhase.OIL,
        CandidatePhase.OIL);
    SpecifiedMultiphaseFlashSolver.Result root = solver.solve(300.0, 100.0, new double[] {0.5, 0.5},
        new double[][] {{0.4, 0.6}, {0.4, 0.6}});
    assertTrue(root.isConverged(), root.getFailureMessage());
    assertTrue(root.hasPhysicalPhaseIdentity(), root.getPhaseIdentityDiagnostic());
    SystemInterface state = solver.toThermodynamicSystem(root);
    assertEquals(PhaseType.OIL, state.getPhase(0).getType());
    assertEquals(PhaseType.OIL, state.getPhase(1).getType());
    SpecifiedPhaseEquilibriumValidator.Result checked = solver.validateEquilibrium(root);
    assertFalse(checked.isAccepted(), "Coincident slots do not establish a stable two-phase equilibrium");
    assertTrue(checked.getViolations().toString().contains("COINCIDENT_PHASES"));
    assertFalse(checked.getViolations().toString().contains("PHASE_IDENTITY_MISMATCH"));
  }

  /** Independent replay must check fugacity instead of rejecting a correctly classified liquid pair. */
  @Test
  void replayUsesMassFamilyBeforeTestingEquilibrium() {
    SystemInterface template = fluid(0.75);
    SystemInterface state = template.clone();
    double[][] compositions = {{0.4, 0.6}, {0.1, 0.9}};
    CandidatePhase[] slots = {CandidatePhase.OIL, CandidatePhase.AQUEOUS};
    for (int phase = 0; phase < 2; phase++) {
      state.setBeta(phase, 0.5);
      state.setPhaseType(phase, phase == 0 ? PhaseType.OIL : PhaseType.AQUEOUS);
      for (int component = 0; component < 2; component++) {
        state.getPhase(phase).getComponent(component).setx(compositions[phase][component]);
      }
      state.init(1, phase);
    }
    assertEquals(PhaseType.OIL, state.getPhase(0).getType());
    assertEquals(PhaseType.AQUEOUS, state.getPhase(1).getType());
    SpecifiedPhaseEquilibriumValidator.Result checked = SpecifiedPhaseEquilibriumValidator.validate(template, state,
        slots);
    assertFalse(checked.isAccepted(), "Initialized liquid compositions have not been equilibrated");
    assertFalse(checked.getViolations().toString().contains("PHASE_IDENTITY_MISMATCH"));
    assertTrue(checked.getViolations().toString().contains("FUGACITY_REPLAY_FAILED"));
  }

  /**
   * Creates a dense hydrocarbon/water mixture without asserting that it is at equilibrium.
   *
   * @param water water mole fraction
   * @return initialized SRK template in component order decane, water
   */
  private static SystemInterface fluid(double water) {
    SystemInterface fluid = new SystemSrkEos(300.0, 100.0);
    fluid.addComponent("nC10", 1.0 - water);
    fluid.addComponent("water", water);
    fluid.setMixingRule("classic");
    fluid.init(0);
    return fluid;
  }
}
