package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import org.junit.jupiter.api.Test;

/**
 * Physical regression tests for the specified two-to-four-phase research solver.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
class SpecifiedMultiphaseFlashAcceptanceTest {
  /** Ensures changing PT between calls leaves no residual state in a reused solver. */
  @Test
  void solveLocalWorkspaceDoesNotLeakBetweenOperatingPoints() {
    SystemInterface first = ResearchEquilibriumFixtures.equilibrium(false);
    SystemInterface second = first.clone();
    second.setTemperature(400.0);
    second.setPressure(30.0);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(second).TPflash();
    assertEquals(2, second.getNumberOfPhases());
    SpecifiedMultiphaseFlashSolver solver = new SpecifiedMultiphaseFlashSolver(first, CandidatePhase.GAS,
        CandidatePhase.OIL);
    for (SystemInterface point : new SystemInterface[] {first, second, first}) {
      double[][] seeds = ResearchEquilibriumFixtures.compositions(first);
      SpecifiedMultiphaseFlashSolver.Result root = solver.solve(point.getTemperature(), point.getPressure(),
          ResearchEquilibriumFixtures.fractions(first), seeds);
      assertTrue(root.isConverged(), root.getFailureMessage());
      for (int phase = 0; phase < 2; phase++) {
        assertArrayEquals(ResearchEquilibriumFixtures.compositions(point)[phase], root.getPhaseComposition(phase),
            1.0e-7);
      }
      assertTrue(solver.validateEquilibrium(root).isAccepted());
    }
    assertEquals(350.0, first.getTemperature(), 0.0);
    assertEquals(10.0, first.getPressure(), 0.0);
  }

  /** Finds a metastable gas/oil root and requires rejection when an aqueous phase lowers Gibbs energy. */
  @Test
  void rejectsConvergedMetastableRootWithAnAdditionalWaterPhase() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double temperature = 500.0;
    double pressure = 250.0;
    double[][] seeds = new double[2][fluid.getNumberOfComponents()];
    for (int component = 0; component < fluid.getNumberOfComponents(); component++) {
      neqsim.thermo.component.ComponentInterface data = fluid.getPhase(0).getComponent(component);
      double k = data.getPC() / pressure
          * Math.exp(5.373 * (1.0 + data.getAcentricFactor()) * (1.0 - data.getTC() / temperature));
      seeds[0][component] = Math.max(data.getz() * k, 1.0e-100);
      seeds[1][component] = Math.max(data.getz() / k, 1.0e-100);
      if ("water".equalsIgnoreCase(data.getComponentName())) {
        seeds[1][component] *= 1.0e-8;
      }
    }
    SpecifiedMultiphaseFlashSolver solver = new SpecifiedMultiphaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL).setNumericalControls(160, 1.0e-10, 2.0e-5, 1.0, 100.0);
    SpecifiedMultiphaseFlashSolver.Result root = solver.solve(temperature, pressure, new double[] {0.65, 0.35}, seeds);
    assertTrue(root.isConverged(), root.getFailureMessage());
    assertTrue(root.hasPhysicalPhaseIdentity(), root.getPhaseIdentityDiagnostic());
    SpecifiedPhaseEquilibriumValidator.Result checked = solver.validateEquilibrium(root);
    assertFalse(checked.isAccepted());
    assertTrue(checked.getViolations().toString().contains("ADDITIONAL_PHASE_INSTABILITY"),
        checked.getViolations().toString());
    assertThrows(IllegalStateException.class, () -> solver.toValidatedThermodynamicSystem(root));
  }

  /** Verifies perturbed two- and three-phase seeds against the independent ordinary TP flash. */
  @Test
  void reproducesTwoAndThreePhaseFlashWithConservedNormalizedResults() {
    for (boolean wet : new boolean[] {false, true}) {
      SystemInterface fluid = ResearchEquilibriumFixtures.equilibrium(wet);
      double[][] original = ResearchEquilibriumFixtures.compositions(fluid);
      double[][] seeds = ResearchEquilibriumFixtures.compositions(fluid);
      for (double[] row : seeds) {
        for (int index = 0; index < row.length; index++) {
          row[index] *= index % 2 == 0 ? 1.03 : 0.97;
        }
      }
      SpecifiedMultiphaseFlashSolver solver = new SpecifiedMultiphaseFlashSolver(fluid,
          ResearchEquilibriumFixtures.slots(fluid)).setNumericalControls(160, 1.0e-10, 2.0e-5, 1.0, 100.0);
      SpecifiedMultiphaseFlashSolver.Result result = solver.solve(350.0, 10.0,
          ResearchEquilibriumFixtures.fractions(fluid), seeds);
      assertTrue(result.isConverged(), result.getFailureMessage());
      assertTrue(result.getMaterialBalanceResidual() < 1.0e-9);
      assertTrue(result.getIterations() > 0, "perturbed seeds must exercise the nonlinear solver");
      SpecifiedPhaseEquilibriumValidator.Result checked = solver.validateEquilibrium(result);
      assertTrue(checked.isAccepted(), checked.getViolations().toString());
      assertTrue(checked.isStabilityChecked());
      assertTrue(checked.getFugacityResidual() < 1.0e-7);
      assertEquals(fluid.getNumberOfPhases(), solver.toValidatedThermodynamicSystem(result).getNumberOfPhases());
      for (int phase = 0; phase < original.length; phase++) {
        assertArrayEquals(original[phase], result.getPhaseComposition(phase), 1.0e-7);
        assertArrayEquals(original[phase], ResearchEquilibriumFixtures.compositions(fluid)[phase], 0.0);
        double sum = 0.0;
        for (double value : result.getPhaseComposition(phase)) {
          sum += value;
        }
        assertEquals(1.0, sum, 1.0e-14);
      }
      double[] exposed = result.getPhaseComposition(0);
      exposed[0] = -1.0;
      assertTrue(result.getPhaseComposition(0)[0] >= 0.0);
      SpecifiedMultiphaseFlashSolver foreign = new SpecifiedMultiphaseFlashSolver(fluid,
          ResearchEquilibriumFixtures.slots(fluid));
      assertThrows(IllegalArgumentException.class, () -> foreign.toThermodynamicSystem(result));
    }
  }

  /** Ensures a numerically exact duplicated liquid is never accepted as four distinct phases. */
  @Test
  void rejectsCoincidentFourPhaseRootAndVanishingPhase() {
    SystemInterface fluid = ResearchEquilibriumFixtures.equilibrium(false);
    double[][] phases = ResearchEquilibriumFixtures.compositions(fluid);
    double gas = fluid.getBeta(0);
    double oil = fluid.getBeta(1);
    SpecifiedMultiphaseFlashSolver solver = new SpecifiedMultiphaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.OIL, CandidatePhase.OIL);
    SpecifiedMultiphaseFlashSolver.Result root = solver.solve(350.0, 10.0,
        new double[] {gas, oil / 3.0, oil / 3.0, oil / 3.0},
        new double[][] {phases[0], phases[1], phases[1], phases[1]});
    assertTrue(root.isConverged(), root.getFailureMessage());
    SpecifiedPhaseEquilibriumValidator.Result checked = solver.validateEquilibrium(root);
    assertFalse(checked.isAccepted());
    assertTrue(checked.getViolations().toString().contains("COINCIDENT_PHASES"));
    assertThrows(IllegalStateException.class, () -> solver.toValidatedThermodynamicSystem(root));
    SystemInterface diagnostic = solver.toThermodynamicSystem(root);
    diagnostic.setBeta(3, 1.0e-14);
    assertTrue(SpecifiedPhaseEquilibriumValidator.validate(fluid, diagnostic, root.getPhaseSlots()).getViolations()
        .toString().contains("VANISHING_OR_INVALID_PHASE"));
  }

  /** Verifies rejected data are reported before numerical iteration. */
  @Test
  void rejectsInvalidSeedsAndOverflowingTotals() {
    SystemInterface fluid = ResearchEquilibriumFixtures.equilibrium(false);
    SpecifiedMultiphaseFlashSolver solver = new SpecifiedMultiphaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL);
    double[][] phases = ResearchEquilibriumFixtures.compositions(fluid);
    assertThrows(IllegalArgumentException.class, () -> solver.solve(Double.NaN, 10.0, new double[] {0.5, 0.5}, phases));
    assertThrows(IllegalArgumentException.class,
        () -> solver.solve(350.0, 10.0, new double[] {Double.MAX_VALUE, Double.MAX_VALUE}, phases));
    assertThrows(IllegalArgumentException.class, () -> solver.solve(350.0, 10.0, new double[] {0.5, 0.5},
        new double[][] {{Double.MAX_VALUE, Double.MAX_VALUE}, phases[1]}));
    assertThrows(IllegalArgumentException.class, () -> solver.solve(350.0, 10.0, new double[] {0.0, 1.0}, phases));
    assertThrows(IllegalArgumentException.class,
        () -> solver.solve(350.0, 10.0, new double[] {0.5, 0.5}, new double[][] {{-0.1, 1.1}, phases[1]}));
  }

  /** Checks the independent replay catches inventory, fugacity and physical-family corruption. */
  @Test
  void rejectsInvalidPhysicalStateOnReplay() {
    SystemInterface fluid = ResearchEquilibriumFixtures.equilibrium(true);
    CandidatePhase[] slots = ResearchEquilibriumFixtures.slots(fluid);
    SystemInterface bad = fluid.clone();
    bad.setBeta(0, bad.getBeta(0) * 0.5);
    assertTrue(SpecifiedPhaseEquilibriumValidator.validate(fluid, bad, slots).getViolations().toString()
        .contains("MATERIAL_BALANCE_FAILED"));
    CandidatePhase[] wrong = slots.clone();
    wrong[2] = CandidatePhase.OIL;
    assertTrue(SpecifiedPhaseEquilibriumValidator.validate(fluid, fluid, wrong).getViolations().toString()
        .contains("PHASE_IDENTITY_MISMATCH"));
    SystemInterface wrongTemperature = fluid.clone();
    wrongTemperature.setTemperature(400.0);
    assertTrue(SpecifiedPhaseEquilibriumValidator.validate(fluid, wrongTemperature, slots).getViolations().toString()
        .contains("FUGACITY_REPLAY_FAILED"));
  }
}
