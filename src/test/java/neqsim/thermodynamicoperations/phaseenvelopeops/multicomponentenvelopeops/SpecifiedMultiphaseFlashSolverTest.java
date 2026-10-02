package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.Arrays;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class SpecifiedMultiphaseFlashSolverTest extends NeqSimTest {
  private static final Logger logger = LogManager.getLogger(SpecifiedMultiphaseFlashSolverTest.class);

  /**
   * Checks or computes four specified slots reproduce one gas oil equilibrium without changing the template.
   */
  @Test
  @Tag("slow")
  void fourSpecifiedSlotsReproduceOneGasOilEquilibriumWithoutChangingTheTemplate() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double originalTemperature = fluid.getTemperature();
    double originalPressure = fluid.getPressure();
    double temperatureK = 500.0;
    double pressureBara = 250.0;
    SpecifiedTwoPhaseFlashSolver.Result twoPhase = new SpecifiedTwoPhaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL).setNumericalControls(100, 1.0e-10, 2.0e-5).solve(temperatureK, pressureBara, 0.65,
            wilsonSeed(fluid, temperatureK, pressureBara, true), wilsonSeed(fluid, temperatureK, pressureBara, false));
    assertTrue(twoPhase.isConverged(), twoPhase.getFailureMessage());

    double oilFraction = 1.0 - twoPhase.getBeta();
    SpecifiedMultiphaseFlashSolver solver = new SpecifiedMultiphaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.OIL, CandidatePhase.OIL).setNumericalControls(80, 1.0e-9, 2.0e-5, 0.5, 20.0);
    SpecifiedMultiphaseFlashSolver.Result result = solver.solve(temperatureK, pressureBara,
        new double[] {twoPhase.getBeta(), 0.5 * oilFraction, 0.3 * oilFraction, 0.2 * oilFraction},
        new double[][] {twoPhase.getPhaseZeroComposition(), twoPhase.getPhaseOneComposition(),
            twoPhase.getPhaseOneComposition(), twoPhase.getPhaseOneComposition()});

    assertTrue(result.isConverged(), result.getFailureMessage());
    assertTrue(result.getMaximumResidual() <= 1.0e-9);
    assertTrue(result.getMaterialBalanceResidual() <= 1.0e-9);
    assertTrue(result.getPhaseCompositionDistance(1, 2) <= 1.0e-8);
    assertTrue(result.getPhaseCompositionDistance(2, 3) <= 1.0e-8);
    assertTrue(result.hasPhysicalPhaseIdentity(), result.getPhaseIdentityDiagnostic());
    assertEquals(4, solver.toThermodynamicSystem(result).getNumberOfPhases());
    assertEquals(originalTemperature, fluid.getTemperature(), 0.0);
    assertEquals(originalPressure, fluid.getPressure(), 0.0);
  }

  /**
   * Checks or computes fluid two high pressure stationary modes are audited as four phase seeds.
   */
  @Test
  @Tag("slow")
  void fluidTwoHighPressureStationaryModesAreAuditedAsFourPhaseSeeds() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double temperatureK = 504.5;
    double pressureBara = 280.0;
    SpecifiedTwoPhaseFlashSolver retainedSolver = new SpecifiedTwoPhaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL).setNumericalControls(120, 1.0e-10, 2.0e-5);
    SpecifiedTwoPhaseFlashSolver.Result retained = retainedSolver.solve(temperatureK, pressureBara, 0.655,
        wilsonSeed(fluid, temperatureK, pressureBara, true), wilsonSeed(fluid, temperatureK, pressureBara, false));
    assertTrue(retained.isConverged(), retained.getFailureMessage());
    IncipientPhaseStabilityAnalyzer.Result stability = new IncipientPhaseStabilityAnalyzer(
        retainedSolver.toThermodynamicSystem(retained)).setMaximumIterations(800).setDampingFactor(0.15).analyze();
    Candidate oil = null;
    Candidate aqueous = null;
    for (Candidate candidate : stability.getCandidates()) {
      logger.info("Fluid 2 stationary seed phase={} TPD={} residual={} trivial={} water={}", candidate.getPhase(),
          candidate.getTangentPlaneDistance(), candidate.getStationarityResidual(), candidate.isTrivial(),
          water(candidate.getComposition(), fluid));
      if (candidate.isConverged() && candidate.getPhase() == CandidatePhase.OIL && !candidate.isTrivial()) {
        oil = candidate;
      }
      if (candidate.isConverged() && candidate.getPhase() == CandidatePhase.AQUEOUS && !candidate.isTrivial()) {
        aqueous = candidate;
      }
    }
    TwoToThreePhaseArcLengthCorrector.State homogeneousOil = TwoToThreePhaseArcLengthCorrector.State.create(
        CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.OIL, temperatureK, pressureBara, retained.getBeta(),
        retained.getPhaseZeroComposition(), retained.getPhaseOneComposition(), retained.getPhaseOneComposition());
    IncipientPhaseStationarityJacobianAnalyzer.Result oilMode = new IncipientPhaseStationarityJacobianAnalyzer(fluid,
        CandidatePhase.OIL, CandidatePhase.OIL).setFiniteDifferenceStep(2.0e-4).analyze(homogeneousOil);
    IncipientPhaseStabilityAnalyzer explicitAnalyzer = new IncipientPhaseStabilityAnalyzer(
        retainedSolver.toThermodynamicSystem(retained)).setMaximumIterations(800).setDampingFactor(0.15);
    for (double amplitude : new double[] {-0.02, 0.02, -0.05, 0.05, -0.1, 0.1, -0.2, 0.2, -0.5, 0.5, -1.0, 1.0}) {
      Candidate trial = explicitAnalyzer.analyzeCandidate(CandidatePhase.OIL,
          perturbAlongMode(retained.getPhaseOneComposition(), oilMode, amplitude));
      logger.info("Fluid 2 explicit oil seed amplitude={} TPD={} residual={} trivial={} water={} failure={}", amplitude,
          trial.getTangentPlaneDistance(), trial.getStationarityResidual(), trial.isTrivial(),
          water(trial.getComposition(), fluid), trial.getFailureMessage());
      if (trial.isConverged() && trial.getPhase() == CandidatePhase.OIL && !trial.isTrivial()) {
        oil = trial;
        break;
      }
    }
    assertNotNull(oil, "a distinct oil stationary mode is required at the post-spinodal audit state");
    assertNotNull(aqueous, "the water-rich stationary mode is required at the post-spinodal audit state");

    SpecifiedMultiphaseFlashSolver solver = new SpecifiedMultiphaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.OIL, CandidatePhase.AQUEOUS)
        .setNumericalControls(240, 1.0e-8, 1.0e-5, 0.25, 40.0);
    double retainedOilFraction = 1.0 - retained.getBeta();
    SpecifiedMultiphaseFlashSolver.Result physicalVlla = null;
    StringBuilder diagnostics = new StringBuilder();
    for (double oilTwoFraction : new double[] {1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2}) {
      for (double aqueousFraction : new double[] {1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2}) {
        double retainedScale = 1.0 - oilTwoFraction - aqueousFraction;
        SpecifiedMultiphaseFlashSolver.Result trial = solver.solve(temperatureK, pressureBara,
            new double[] {retainedScale * retained.getBeta(), retainedScale * retainedOilFraction, oilTwoFraction,
                aqueousFraction},
            new double[][] {retained.getPhaseZeroComposition(), retained.getPhaseOneComposition(), oil.getComposition(),
                aqueous.getComposition()});
        logger.info(
            "Fluid 2 VLLA seed oil2={} water={} converged={} physical={} fractions={} residual={} balance={} oilDistance={} failure={} identity={}",
            oilTwoFraction, aqueousFraction, trial.isConverged(), trial.hasPhysicalPhaseIdentity(),
            Arrays.toString(trial.getPhaseFractions()), trial.getMaximumResidual(), trial.getMaterialBalanceResidual(),
            trial.getPhaseCompositionDistance(1, 2), trial.getFailureMessage(), trial.getPhaseIdentityDiagnostic());
        diagnostics.append(" seed=").append(oilTwoFraction).append('/').append(aqueousFraction).append(" converged=")
            .append(trial.isConverged()).append(" physical=").append(trial.hasPhysicalPhaseIdentity())
            .append(" fractions=").append(Arrays.toString(trial.getPhaseFractions())).append(" residual=")
            .append(trial.getMaximumResidual()).append(" oilDistance=").append(trial.getPhaseCompositionDistance(1, 2))
            .append(" identity=").append(trial.getPhaseIdentityDiagnostic()).append(';');
        if (trial.isConverged() && trial.hasPhysicalPhaseIdentity() && trial.getPhaseCompositionDistance(1, 2) > 1.0e-5
            && trial.getPhaseFraction(2) > 1.0e-10 && trial.getPhaseFraction(3) > 1.0e-10) {
          physicalVlla = trial;
        }
      }
    }
    assertNull(physicalVlla,
        "an arbitrary post-spinodal PT state must not be promoted to a finite VLLA equilibrium;" + diagnostics);
  }

  /**
   * Checks or computes wilson seed.
   *
   * @param fluid fluid
   * @param temperatureK temperature in kelvin
   * @param pressureBara absolute pressure in bara
   * @param gas gas
   * @return computed wilson seed result
   */
  private static double[] wilsonSeed(SystemInterface fluid, double temperatureK, double pressureBara, boolean gas) {
    int componentCount = fluid.getPhase(0).getNumberOfComponents();
    double[] seed = new double[componentCount];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      ComponentInterface component = fluid.getPhase(0).getComponent(componentIndex);
      double wilsonK = component.getPC() / pressureBara
          * Math.exp(5.373 * (1.0 + component.getAcentricFactor()) * (1.0 - component.getTC() / temperatureK));
      wilsonK = Math.max(1.0e-20, Math.min(1.0e20, wilsonK));
      seed[componentIndex] = Math.max(component.getz() * (gas ? wilsonK : 1.0 / wilsonK), 1.0e-100);
      if (!gas && component.getComponentName().equalsIgnoreCase("water")) {
        seed[componentIndex] *= 1.0e-8;
      }
      total += seed[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      seed[componentIndex] /= total;
    }
    return seed;
  }

  /**
   * Checks or computes water.
   *
   * @param composition composition
   * @param fluid fluid
   * @return computed water result
   */
  private static double water(double[] composition, SystemInterface fluid) {
    return composition[fluid.getPhase(0).getComponent("water").getComponentNumber()];
  }

  /**
   * Checks or computes perturb along mode.
   *
   * @param base base
   * @param mode mode
   * @param amplitude amplitude
   * @return computed perturb along mode result
   */
  private static double[] perturbAlongMode(double[] base, IncipientPhaseStationarityJacobianAnalyzer.Result mode,
      double amplitude) {
    int referenceIndex = mode.getReferenceComponentIndex();
    double[] eigenvector = mode.getBifurcationEigenvector();
    double[] logarithms = new double[base.length];
    double reference = Math.max(base[referenceIndex], 1.0e-100);
    int coordinateIndex = 0;
    double maximum = 0.0;
    for (int componentIndex = 0; componentIndex < base.length; componentIndex++) {
      if (componentIndex != referenceIndex) {
        logarithms[componentIndex] = Math.log(Math.max(base[componentIndex], 1.0e-100) / reference)
            + amplitude * eigenvector[coordinateIndex++];
      }
      maximum = Math.max(maximum, logarithms[componentIndex]);
    }
    double[] composition = new double[base.length];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = Math.exp(Math.max(-700.0, logarithms[componentIndex] - maximum));
      total += composition[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] /= total;
    }
    return composition;
  }
}
