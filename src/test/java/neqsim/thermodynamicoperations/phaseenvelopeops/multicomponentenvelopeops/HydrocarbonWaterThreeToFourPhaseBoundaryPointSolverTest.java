package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterThreeToFourPhaseBoundaryPointSolverTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void fluidTwoSpinodalDoesNotInventGasOilOilWithIncipientWaterBoundary() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double temperatureK = 504.366971105;
    double pressureBara = 275.169677734;
    SpecifiedTwoPhaseFlashSolver retainedSolver = new SpecifiedTwoPhaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL).setNumericalControls(160, 1.0e-10, 2.0e-5);
    SpecifiedTwoPhaseFlashSolver.Result retained = retainedSolver.solve(temperatureK, pressureBara, 0.655,
        wilsonSeed(fluid, temperatureK, pressureBara, true), wilsonSeed(fluid, temperatureK, pressureBara, false));
    assertTrue(retained.isConverged(), retained.getFailureMessage());

    IncipientPhaseStabilityAnalyzer analyzer = new IncipientPhaseStabilityAnalyzer(
        retainedSolver.toThermodynamicSystem(retained)).setMaximumIterations(1000).setDampingFactor(0.15);
    Candidate aqueous = analyzer.analyzeCandidate(CandidatePhase.AQUEOUS, aqueousSeed(fluid));
    assertTrue(aqueous.isConverged(), aqueous.getFailureMessage());
    assertTrue(aqueous.getPhase() == CandidatePhase.AQUEOUS);
    assertTrue(water(aqueous.getComposition(), fluid) >= 0.5);

    TwoToThreePhaseArcLengthCorrector.State boundaryState = TwoToThreePhaseArcLengthCorrector.State.create(
        CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS, temperatureK, pressureBara, retained.getBeta(),
        retained.getPhaseZeroComposition(), retained.getPhaseOneComposition(), aqueous.getComposition());
    TwoToThreePhaseArcLengthCorrector.State homogeneousOil = TwoToThreePhaseArcLengthCorrector.State.create(
        CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.OIL, temperatureK, pressureBara, retained.getBeta(),
        retained.getPhaseZeroComposition(), retained.getPhaseOneComposition(), retained.getPhaseOneComposition());
    IncipientPhaseStationarityJacobianAnalyzer.Result mode = new IncipientPhaseStationarityJacobianAnalyzer(fluid,
        CandidatePhase.OIL, CandidatePhase.OIL).setFiniteDifferenceStep(2.0e-4).analyze(homogeneousOil);
    IncipientPhaseCurvatureAnalyzer.Result tpdCurvature = new IncipientPhaseCurvatureAnalyzer(fluid, CandidatePhase.OIL,
        CandidatePhase.OIL).setFiniteDifferenceStep(2.0e-4).analyze(homogeneousOil);

    HydrocarbonWaterThreeToFourPhaseBoundaryPointSolver solver = new HydrocarbonWaterThreeToFourPhaseBoundaryPointSolver(
        fluid, CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.OIL, CandidatePhase.AQUEOUS)
        .setNumericalControls(integerProperty("neqsim.fluid2.maximumIterations", 240), 1.0e-8, 2.0e-5, 0.25, 40.0);
    HydrocarbonWaterThreeToFourPhaseBoundaryPointSolver.Result accepted = null;
    StringBuilder diagnostics = new StringBuilder(" modeEigen=").append(mode.getBifurcationEigenvalue())
        .append(" symmetricMinimum=").append(mode.getMinimumEigenvalue()).append(" minimumSingular=")
        .append(mode.getMinimumSingularValue()).append(" antisymmetry=").append(mode.getMaximumAntisymmetry())
        .append(" tpdMinimum=").append(tpdCurvature.getMinimumEigenvalue()).append(" tpdGradient=")
        .append(tpdCurvature.getMaximumGradient());
    for (double separation : new double[] { doubleProperty("neqsim.fluid2.separation", 0.01) }) {
      for (double fraction : new double[] { doubleProperty("neqsim.fluid2.fraction", 1.0e-4) }) {
        HydrocarbonWaterThreeToFourPhaseBoundaryPointSolver.Result trial = solver.solve(boundaryState, mode, separation,
            fraction);
        diagnostics.append(" sep=").append(separation).append(" fraction=").append(fraction).append(" converged=")
            .append(trial.isConverged()).append(" T=").append(trial.getTemperatureK()).append(" P=")
            .append(trial.getPressureBara()).append(" betas=")
            .append(java.util.Arrays.toString(trial.getPhaseFractions())).append(" residual=")
            .append(trial.getMaximumResidual()).append(" equilibrium=").append(trial.getEquilibriumMaximumResidual())
            .append(" incipient=").append(trial.getIncipientMaximumResidual()).append(" separationResidual=")
            .append(trial.getSeparationResidual()).append(" balance=").append(trial.getMaterialBalanceResidual())
            .append(" accepted/rejected=").append(trial.getAcceptedSteps()).append('/').append(trial.getRejectedSteps())
            .append(" radius=").append(trial.getFinalTrustRadius()).append(" condition=")
            .append(trial.getJacobianConditionNumber()).append(" maxResidualIndex=")
            .append(maximumResidualIndex(trial.getResidual())).append(" initial=")
            .append(trial.getInitialMaximumResidual()).append(" failure=").append(trial.getFailureMessage())
            .append(';');
        if (trial.isConverged()) {
          accepted = trial;
          break;
        }
      }
      if (accepted != null) {
        break;
      }
    }
    assertNull(accepted,
        "Fluid 2 uses a continuous retained-OIL stationary-branch switch at this limit; a finite G/O1/O2 + "
            + "incipient-water root would be a false topology: " + diagnostics);
  }

  private static double[] aqueousSeed(SystemInterface fluid) {
    int componentCount = fluid.getPhase(0).getNumberOfComponents();
    double[] seed = new double[componentCount];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      ComponentInterface component = fluid.getPhase(0).getComponent(componentIndex);
      seed[componentIndex] = component.getComponentName().equalsIgnoreCase("water") ? 0.999
          : Math.max(component.getz(), 1.0e-100) * 1.0e-4;
      total += seed[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      seed[componentIndex] /= total;
    }
    return seed;
  }

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

  private static double water(double[] composition, SystemInterface fluid) {
    return composition[fluid.getPhase(0).getComponent("water").getComponentNumber()];
  }

  private static int maximumResidualIndex(double[] residual) {
    int result = -1;
    double maximum = -1.0;
    for (int index = 0; index < residual.length; index++) {
      if (Math.abs(residual[index]) > maximum) {
        maximum = Math.abs(residual[index]);
        result = index;
      }
    }
    return result;
  }

  private static double doubleProperty(String name, double defaultValue) {
    String value = System.getProperty(name);
    return value == null ? defaultValue : Double.parseDouble(value);
  }

  private static int integerProperty(String name, int defaultValue) {
    String value = System.getProperty(name);
    return value == null ? defaultValue : Integer.parseInt(value);
  }
}
