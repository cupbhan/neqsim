package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterSecondaryOilTransitionTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void fluidTwoTracksTheDistinctOilStationaryBranchBackToItsZeroTpdLimit() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double temperatureK = 504.366971105;
    double pressureBara = 275.169677734;
    SpecifiedTwoPhaseFlashSolver retainedSolver = new SpecifiedTwoPhaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL).setNumericalControls(160, 1.0e-10, 2.0e-5);
    SpecifiedTwoPhaseFlashSolver.Result retained = retainedSolver.solve(temperatureK, pressureBara, 0.655,
        wilsonSeed(fluid, temperatureK, pressureBara, true), wilsonSeed(fluid, temperatureK, pressureBara, false));
    assertTrue(retained.isConverged(), retained.getFailureMessage());
    Candidate aqueous = new IncipientPhaseStabilityAnalyzer(retainedSolver.toThermodynamicSystem(retained))
        .setMaximumIterations(1000).setDampingFactor(0.15).analyzeCandidate(CandidatePhase.AQUEOUS, aqueousSeed(fluid));
    assertTrue(aqueous.isConverged(), aqueous.getFailureMessage());

    TwoToThreePhaseArcLengthCorrector.State state = TwoToThreePhaseArcLengthCorrector.State.create(CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.AQUEOUS, temperatureK, pressureBara, retained.getBeta(),
        retained.getPhaseZeroComposition(), retained.getPhaseOneComposition(), aqueous.getComposition());
    List<TwoToThreePhaseBoundaryPointSolver.Result> roots = new ArrayList<TwoToThreePhaseBoundaryPointSolver.Result>();
    roots.add(TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(state, 1.0e-9));
    TwoToThreePhaseArcLengthCorrector corrector = new TwoToThreePhaseArcLengthCorrector(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL, CandidatePhase.AQUEOUS).setNumericalControls(200, 1.0e-8, 2.0e-5);
    for (double targetPressure : new double[] { 275.2, 275.3, 275.5, 276.0, 277.0, 278.0, 279.0, 280.0, 282.0,
        284.0 }) {
      TwoToThreePhaseArcLengthCorrector.Result corrected = corrector.correctAtPressure(state, targetPressure);
      assertTrue(corrected.isConverged(), "P=" + targetPressure + " residual=" + corrected.getMaximumResidual()
          + " failure=" + corrected.getFailureMessage());
      state = corrected.getState();
      roots.add(TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(state,
          corrected.getThermodynamicMaximumResidual()));
    }

    Candidate alternateOil = distinctOilCandidate(fluid, state);
    assertNotNull(alternateOil, "a distinct terminal oil stationary point is required");
    HydrocarbonWaterSecondaryStationaryBranchTracker.Result tracked = new HydrocarbonWaterSecondaryStationaryBranchTracker(
        fluid).setNumericalControls(160, 1000, 1.0e-10, 2.0e-5, 0.15, 0.25, 1.0e-8)
        .trackBackward(roots, CandidatePhase.OIL, alternateOil.getComposition());

    StringBuilder diagnostics = new StringBuilder(" terminalTPD=").append(alternateOil.getTangentPlaneDistance())
        .append(" zeroBrackets=").append(tracked.getZeroTpdBrackets().size()).append(" branchLimits=")
        .append(tracked.getBranchLimitBrackets().size()).append(" termination=")
        .append(tracked.getTerminationMessage());
    for (HydrocarbonWaterSecondaryStationaryBranchTracker.Sample sample : tracked.getSamples()) {
      diagnostics.append(" [P=").append(sample.getBoundaryRoot().getPressureBara()).append(" T=")
          .append(sample.getBoundaryRoot().getTemperatureK()).append(" tracked=").append(sample.isTracked())
          .append(" tpd=").append(sample.getTangentPlaneDistance()).append(" jump=").append(sample.getCompositionJump())
          .append(" distance=")
          .append(sample.getCandidate() == null ? Double.NaN
              : distance(sample.getCandidate().getComposition(), sample.getRetainedFlash().getPhaseOneComposition()))
          .append(" failure=").append(sample.getFailureMessage()).append(']');
    }
    assertTrue(tracked.getZeroTpdBrackets().isEmpty(),
        "the Fluid 2 secondary-oil branch is expected to leave zero TPD continuously at coalescence;" + diagnostics);
    HydrocarbonWaterSecondaryStationaryBranchTracker.Sample first = tracked.getSamples().get(0);
    assertTrue(first.isTracked(), first.getFailureMessage());
    assertTrue(Math.abs(first.getTangentPlaneDistance()) <= 1.0e-6, diagnostics.toString());
    assertTrue(
        distance(first.getCandidate().getComposition(), first.getRetainedFlash().getPhaseOneComposition()) <= 1.0e-4,
        diagnostics.toString());
    assertTrue(tracked.getSamples().stream().allMatch(sample -> sample.getTangentPlaneDistance() <= 1.0e-8),
        diagnostics.toString());

    TwoToThreePhaseBoundaryPointSolver.Result terminalRoot = roots.get(roots.size() - 1);
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejected = new HydrocarbonWaterBoundaryGlobalStabilityGate(fluid)
        .evaluate(terminalRoot);
    assertTrue(!rejected.isAccepted() && rejected.isTargetMatched(), rejected.getFailureMessage());
    HydrocarbonWaterRetainedPhaseBranchSwitcher.Result switched = new HydrocarbonWaterRetainedPhaseBranchSwitcher(fluid)
        .setFlashControls(200, 1.0e-10, 2.0e-5).setBoundaryControls(24, 80, 1.0e-5, 1.0e-8, 1.0e-4)
        .switchAndCorrect(rejected, alternateOil, state.getTemperatureK() - 40.0, state.getTemperatureK() + 40.0);
    StringBuilder switchDiagnostics = new StringBuilder(" failure=").append(switched.getFailureMessage())
        .append(" distinctBranches=").append(switched.getDistinctRetainedBranchCount()).append(" corrections=")
        .append(switched.getCorrections().size());
    for (HydrocarbonWaterRetainedPhaseBranchSwitcher.Correction correction : switched.getCorrections()) {
      SpecifiedTwoPhaseFlashSolver.Result branch = correction.getRetainedBranch();
      switchDiagnostics.append(" [beta=").append(branch.getBeta()).append(" GOdistance=")
          .append(distance(branch.getPhaseZeroComposition(), branch.getPhaseOneComposition())).append(" gasWater=")
          .append(water(branch.getPhaseZeroComposition(), fluid)).append(" oilWater=")
          .append(water(branch.getPhaseOneComposition(), fluid)).append(" roots=")
          .append(correction.getRootSet().getRoots().size()).append(" scanFailure=")
          .append(correction.getRootSet().getFailureMessage());
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : correction.getClassifications()) {
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = classification.getGlobalStabilityResult();
        switchDiagnostics.append(" {T=").append(classification.getBoundaryRoot().getTemperatureK()).append(" class=")
            .append(classification.getClassification()).append(" minTPD=")
            .append(stability == null ? Double.NaN : stability.getMinimumNonTrivialTangentPlaneDistance())
            .append(" failure=").append(classification.getFailureMessage()).append('}');
      }
      switchDiagnostics.append(']');
    }
    assertEquals(1, switched.getDistinctRetainedBranchCount(),
        "the explicitly tracked lower-TPD oil mode did not seed a distinct retained GO branch: " + switchDiagnostics);
    assertFalse(switched.hasGloballyStableBoundaryRoot(), switchDiagnostics.toString());
    assertTrue(switched.getCorrections().stream().flatMap(correction -> correction.getClassifications().stream())
        .allMatch(classification -> classification
            .getClassification() == HydrocarbonWaterBoundaryEndpointClassifier.Classification.PHASE_COALESCENCE_CANDIDATE),
        switchDiagnostics.toString());
  }

  private static Candidate distinctOilCandidate(SystemInterface fluid, TwoToThreePhaseArcLengthCorrector.State state) {
    SpecifiedTwoPhaseFlashSolver retainedSolver = new SpecifiedTwoPhaseFlashSolver(fluid, CandidatePhase.GAS,
        CandidatePhase.OIL).setNumericalControls(160, 1.0e-10, 2.0e-5);
    SpecifiedTwoPhaseFlashSolver.Result retained = retainedSolver.solve(state.getTemperatureK(),
        state.getPressureBara(), state.getBeta(), state.getPhaseZeroComposition(), state.getPhaseOneComposition());
    if (!retained.isConverged()) {
      return null;
    }
    TwoToThreePhaseArcLengthCorrector.State homogeneous = TwoToThreePhaseArcLengthCorrector.State.create(
        CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.OIL, state.getTemperatureK(), state.getPressureBara(),
        retained.getBeta(), retained.getPhaseZeroComposition(), retained.getPhaseOneComposition(),
        retained.getPhaseOneComposition());
    IncipientPhaseStationarityJacobianAnalyzer.Result mode = new IncipientPhaseStationarityJacobianAnalyzer(fluid,
        CandidatePhase.OIL, CandidatePhase.OIL).setFiniteDifferenceStep(2.0e-4).analyze(homogeneous);
    IncipientPhaseStabilityAnalyzer analyzer = new IncipientPhaseStabilityAnalyzer(
        retainedSolver.toThermodynamicSystem(retained)).setMaximumIterations(1000).setDampingFactor(0.15);
    Candidate best = null;
    for (double amplitude : new double[] { -0.02, 0.02, -0.05, 0.05, -0.1, 0.1, -0.2, 0.2, -0.5, 0.5, -1.0, 1.0 }) {
      Candidate candidate = analyzer.analyzeCandidate(CandidatePhase.OIL,
          perturbAlongMode(retained.getPhaseOneComposition(), mode, amplitude));
      if (candidate.isConverged() && candidate.getPhase() == CandidatePhase.OIL && !candidate.isTrivial()
          && water(candidate.getComposition(), fluid) < 0.5
          && (best == null || candidate.getTangentPlaneDistance() < best.getTangentPlaneDistance())) {
        best = candidate;
      }
    }
    return best;
  }

  private static double[] aqueousSeed(SystemInterface fluid) {
    double[] seed = new double[fluid.getPhase(0).getNumberOfComponents()];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
      ComponentInterface component = fluid.getPhase(0).getComponent(componentIndex);
      seed[componentIndex] = component.getComponentName().equalsIgnoreCase("water") ? 0.999
          : Math.max(component.getz(), 1.0e-100) * 1.0e-4;
      total += seed[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
      seed[componentIndex] /= total;
    }
    return seed;
  }

  private static double[] wilsonSeed(SystemInterface fluid, double temperatureK, double pressureBara, boolean gas) {
    double[] seed = new double[fluid.getPhase(0).getNumberOfComponents()];
    double total = 0.0;
    for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
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
    for (int componentIndex = 0; componentIndex < seed.length; componentIndex++) {
      seed[componentIndex] /= total;
    }
    return seed;
  }

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

  private static double water(double[] composition, SystemInterface fluid) {
    return composition[fluid.getPhase(0).getComponent("water").getComponentNumber()];
  }

  private static double distance(double[] first, double[] second) {
    double result = 0.0;
    for (int index = 0; index < first.length; index++) {
      result += Math.abs(first[index] - second[index]);
    }
    return result;
  }
}
