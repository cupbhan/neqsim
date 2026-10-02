package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertEquals;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.jupiter.api.Test;

/**
 * Regression of the retained-phase corrector and its physical rejection checks.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
class RetainedPhaseBifurcationCorrectorTest {
  private static final Logger logger = LogManager.getLogger(RetainedPhaseBifurcationCorrectorTest.class);

  /**
   * Replays a known cubic three-phase root but rejects relabelling water as a second oil phase.
   */
  @Test
  void convergedAqueousRootCannotBePromotedToASecondOilPhase() {
    SystemInterface fluid = ResearchEquilibriumFixtures.equilibrium(true);
    double[][] x = ResearchEquilibriumFixtures.compositions(fluid);
    double retainedBeta = fluid.getBeta(0) / (fluid.getBeta(0) + fluid.getBeta(1));
    TwoToThreePhaseArcLengthCorrector.State homogeneous = TwoToThreePhaseArcLengthCorrector.State.create(
        CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.OIL, 350.0, 10.0, retainedBeta, x[0], x[1], x[1]);
    IncipientPhaseStationarityJacobianAnalyzer.Result mode = new IncipientPhaseStationarityJacobianAnalyzer(fluid,
        CandidatePhase.OIL, CandidatePhase.OIL).analyze(homogeneous);
    assertTrue(mode.hasRealBifurcationMode());
    int reference = mode.getReferenceComponentIndex();
    double[] vector = mode.getBifurcationEigenvector();
    double separation = 0.0;
    int coordinate = 0;
    for (int component = 0; component < x[0].length; component++) {
      if (component != reference) {
        separation += vector[coordinate++]
            * (Math.log(x[2][component] / x[2][reference]) - Math.log(x[1][component] / x[1][reference]));
      }
    }
    TwoToThreePhaseArcLengthCorrector.State seed = TwoToThreePhaseArcLengthCorrector.State.create(CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.OIL, 350.0, 10.0, retainedBeta, x[0], x[1], x[2]);
    RetainedPhaseBifurcationCorrector corrector = new RetainedPhaseBifurcationCorrector(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.OIL).setNumericalControls(160, 1.0e-8, 2.0e-5);
    RetainedPhaseBifurcationCorrector.Result root = corrector.correct(seed, mode, 10.0, separation, fluid.getBeta(2));
    assertTrue(root.isConverged(), root.getFailureMessage());
    assertTrue(root.getMaterialBalanceResidual() < 1.0e-8);
    assertEquals(350.0, root.getTemperatureK(), 1.0e-5);
    assertEquals(10.0, root.getPressureBara(), 1.0e-9);
    SpecifiedPhaseEquilibriumValidator.Result checked = corrector.validateEquilibrium(root);
    assertFalse(checked.isAccepted());
    assertTrue(checked.getViolations().toString().contains("PHASE_IDENTITY_MISMATCH"));
    assertThrows(IllegalStateException.class, () -> corrector.toValidatedThermodynamicSystem(root));
    RetainedPhaseBifurcationCorrector foreign = new RetainedPhaseBifurcationCorrector(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.OIL);
    assertThrows(IllegalArgumentException.class, () -> foreign.toThermodynamicSystem(root));
    assertThrows(IllegalArgumentException.class, () -> corrector.correct(seed, null, 10.0, 0.1, 0.01));
    assertThrows(IllegalArgumentException.class, () -> corrector.correct(seed, mode, 10.0, 0.0, 0.01));
    assertThrows(IllegalArgumentException.class, () -> corrector.correct(seed, mode, 10.0, 0.1, Double.NaN));
    TwoToThreePhaseArcLengthCorrector.State invalid = TwoToThreePhaseArcLengthCorrector.State.create(CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.OIL, Double.NaN, 10.0, retainedBeta, x[0], x[1], x[2]);
    assertThrows(IllegalArgumentException.class, () -> corrector.correct(invalid, mode, 10.0, 0.1, 0.01));
    assertEquals(350.0, fluid.getTemperature(), 0.0);
    assertEquals(10.0, fluid.getPressure(), 0.0);
  }

  /**
   * Replays the historical high-pressure spinodal without inventing a finite second oil phase.
   */
  @Test
  void fluidTwoSpinodalNullModeDoesNotInventAFiniteThreePhaseOilSplit() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double temperatureK = 504.366971105;
    double pressureBara = 275.169677734;
    SpecifiedTwoPhaseFlashSolver.Result retained = new SpecifiedTwoPhaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL).setNumericalControls(120, 1.0e-10, 2.0e-5).solve(temperatureK, pressureBara, 0.655,
            wilsonSeed(fluid, temperatureK, pressureBara, true), wilsonSeed(fluid, temperatureK, pressureBara, false));
    assertTrue(retained.isConverged(), retained.getFailureMessage());
    TwoToThreePhaseArcLengthCorrector.State homogeneousOil = TwoToThreePhaseArcLengthCorrector.State.create(
        CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.OIL, temperatureK, pressureBara, retained.getBeta(),
        retained.getPhaseZeroComposition(), retained.getPhaseOneComposition(), retained.getPhaseOneComposition());
    IncipientPhaseStationarityJacobianAnalyzer.Result mode = new IncipientPhaseStationarityJacobianAnalyzer(fluid,
        CandidatePhase.OIL, CandidatePhase.OIL).setFiniteDifferenceStep(2.0e-4).analyze(homogeneousOil);
    logger.info("spinodal predictor eigen={} residual={}", mode.getBifurcationEigenvalue(), mode.getMaximumResidual());
    RetainedPhaseBifurcationCorrector corrector = new RetainedPhaseBifurcationCorrector(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.OIL).setNumericalControls(160, 1.0e-9, 2.0e-5);
    RetainedPhaseBifurcationCorrector.Result accepted = null;
    for (double targetSeparation : new double[] {-0.01, 0.01, -0.02, 0.02, -0.05, 0.05}) {
      double fractionSeed = 1.0e-4;
      double firstAmplitude = -targetSeparation * fractionSeed / (1.0 - retained.getBeta());
      TwoToThreePhaseArcLengthCorrector.State predictor = TwoToThreePhaseArcLengthCorrector.State.create(
          CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.OIL, temperatureK, pressureBara, retained.getBeta(),
          retained.getPhaseZeroComposition(), perturbAlongMode(retained.getPhaseOneComposition(), mode, firstAmplitude),
          perturbAlongMode(retained.getPhaseOneComposition(), mode, targetSeparation));
      RetainedPhaseBifurcationCorrector.Result correction = corrector.correct(predictor, mode, pressureBara,
          targetSeparation, fractionSeed);
      logger.info(
          "bordered switch separation={} converged={} T={} K newFraction={} residual={} equilibrium={} pressure={} separationResidual={} balance={} condition={} oilDistance={} water=({}, {}) failure={}",
          targetSeparation, correction.isConverged(), correction.getTemperatureK(), correction.getBifurcatingFraction(),
          correction.getMaximumResidual(), correction.getEquilibriumMaximumResidual(), correction.getPressureResidual(),
          correction.getSeparationResidual(), correction.getMaterialBalanceResidual(),
          correction.getJacobianConditionNumber(), correction.getBifurcatingPhaseDistance(),
          water(correction.getPhaseOneComposition(), fluid), water(correction.getBifurcatingComposition(), fluid),
          correction.getFailureMessage());
      if (correction.isConverged() && correction.getBifurcatingPhaseDistance() > 1.0e-5
          && water(correction.getPhaseOneComposition(), fluid) < 0.5
          && water(correction.getBifurcatingComposition(), fluid) < 0.5) {
        accepted = correction;
        break;
      }
    }
    assertNull(accepted,
        "the Fluid 2 retained-OIL spinodal is a stationary-branch stability limit, not evidence for a finite "
            + "G/O1/O2 equilibrium");
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

  /**
   * Checks or computes water.
   *
   * @param result result produced by this solver instance
   * @param fluid fluid
   * @param phaseIndex zero-based phase index
   * @return computed water result
   */
  private static double water(SpecifiedThreePhaseFlashSolver.Result result, SystemInterface fluid, int phaseIndex) {
    int waterIndex = fluid.getPhase(0).getComponent("water").getComponentNumber();
    return result.getPhaseComposition(phaseIndex)[waterIndex];
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
   * Checks or computes distance.
   *
   * @param first first
   * @param second second
   * @return computed distance result
   */
  private static double distance(double[] first, double[] second) {
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

}
