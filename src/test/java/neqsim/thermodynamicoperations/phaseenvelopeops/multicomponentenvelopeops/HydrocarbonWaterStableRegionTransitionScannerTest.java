package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterStableRegionTransitionScannerTest extends NeqSimTest {

  @Test
  void rejectsInvalidControls() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    HydrocarbonWaterStableRegionTransitionScanner scanner = new HydrocarbonWaterStableRegionTransitionScanner(fluid);
    assertThrows(IllegalArgumentException.class, () -> scanner.setNumericalControls(0, 1.0e-4, 1.0e-9));
    assertThrows(IllegalArgumentException.class, () -> scanner.scan(new double[] {300.0}, new double[] {50.0}));
  }

  @Test
  @Tag("slow")
  void fluidTwoInventoriesStableRegionTransitionsAtFiftyBaraWithoutMutatingCaller() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double originalTemperature = fluid.getTemperature();
    double originalPressure = fluid.getPressure();
    double[] temperatures = new double[] {196.75, 213.15, 233.15, 253.15, 273.15, 298.15, 323.15, 348.15, 373.15,
        398.15, 423.15, 448.15, 473.15, 498.15, 523.15, 548.15, 573.15, 623.15, 673.15, 723.15, 753.15};

    HydrocarbonWaterStableRegionTransitionScanner.Result result = new HydrocarbonWaterStableRegionTransitionScanner(
        fluid).setNumericalControls(80, 1.0e-4, 1.0e-9).scan(temperatures, new double[] {50.0});

    for (HydrocarbonWaterStableRegionTransitionScanner.StableState state : result.getStates()) {
      System.out.printf("Fluid 2 stable state T=%.6f K P=%.3f bara region=%s rawPhases=%d%n", state.getTemperatureK(),
          state.getPressureBara(), state.getRegion(), state.getRawPhaseCount());
    }
    for (HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket : result.getBrackets()) {
      System.out.printf(
          "Fluid 2 stable transition family=%s T=[%.9f, %.9f] K width=%.3g twoPhase=%s threePhase=%s bisections=%d%n",
          bracket.getFamily(), bracket.getLowerTemperatureState().getTemperatureK(),
          bracket.getUpperTemperatureState().getTemperatureK(), bracket.getTemperatureWidthK(),
          bracket.getTwoPhaseState().getRegion(), bracket.getThreePhaseState().getRegion(), bracket.getBisections());
      assertTrue(bracket.getIncipientPhaseFractionOnThreePhaseSide() > 0.0);
      assertTrue(bracket.getRetainedPhaseZeroSplitOnTwoPhaseSide() > 0.0);
      assertTrue(bracket.getRetainedPhaseZeroSplitOnTwoPhaseSide() < 1.0);
      assertTrue(bracket.getRetainedPhaseZeroSplitOnThreePhaseSide() > 0.0);
      assertTrue(bracket.getRetainedPhaseZeroSplitOnThreePhaseSide() < 1.0);
      assertTrue(bracket.getRetainedPhaseSplitJump() >= 0.0);
      assertTrue(bracket.getRetainedPhaseZeroCompositionJump() >= 0.0);
      assertTrue(bracket.getRetainedPhaseOneCompositionJump() >= 0.0);
      HydrocarbonWaterStableRegionBoundaryCorrector.Result corrected = new HydrocarbonWaterStableRegionBoundaryCorrector(
          fluid).setNumericalControls(8, 80, 1.0e-7, 1.0e-8).correct(bracket);
      System.out.printf("  strict correction converged=%s roots=%d accepted=%d failure=%s%n", corrected.isConverged(),
          corrected.getRootSet() == null ? 0 : corrected.getRootSet().getRoots().size(),
          corrected.getAcceptedOrdinaryRootCount(), corrected.getFailureMessage());
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : corrected.getClassifications()) {
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = classification.getGlobalStabilityResult();
        System.out.printf("    T=%.9f K class=%s globalTPD=%s failure=%s%n",
            classification.getBoundaryRoot().getTemperatureK(), classification.getClassification(),
            stability == null ? "n/a" : Double.toString(stability.getMinimumNonTrivialTangentPlaneDistance()),
            classification.getFailureMessage());
      }
      assertTrue(corrected.isConverged(), corrected.getFailureMessage());
      assertFalse(corrected.getIntervalAttempts().isEmpty());
      assertEquals("TWO_PHASE_SIDE", corrected.getIntervalAttempts().get(0).getSeedSource());
    }
    for (HydrocarbonWaterStableRegionTransitionScanner.Failure failure : result.getFailures()) {
      System.out.printf("Fluid 2 stable transition diagnostic T=%.6f K P=%.3f bara reason=%s%n",
          failure.getTemperatureK(), failure.getPressureBara(), failure.getReason());
    }
    assertEquals(temperatures.length, result.getStates().size());
    assertTrue(result.getFlashEvaluations() >= temperatures.length);
    assertEquals(1, result.getBrackets().size());
    assertEquals(HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily.GO_TO_GOW,
        result.getBrackets().get(0).getFamily());
    assertEquals(originalTemperature, fluid.getTemperature(), 0.0);
    assertEquals(originalPressure, fluid.getPressure(), 0.0);
  }

  @Test
  @Tag("slow")
  void fluidTwoBuildsGloballyStableGoToGowAnchorsAcrossPressureLevels() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperatures = new double[] {196.75, 233.15, 273.15, 323.15, 373.15, 423.15, 473.15, 498.15, 523.15,
        548.15, 573.15, 623.15, 673.15, 723.15, 753.15};

    HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).setCorrectionControls(32, 80, 1.0e-5, 1.0e-8)
        .discoverFromStableRegionTransitions(temperatures, new double[] {40.0, 50.0, 60.0});
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = stable.getDiscovery();

    for (HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket : stable.getStableRegionScan()
        .getBrackets()) {
      System.out.printf("Fluid 2 stable discovery bracket family=%s P=%.9f T=[%.9f, %.9f] width=%.8g%n",
          bracket.getFamily(), bracket.getPressureBara(), bracket.getLowerTemperatureState().getTemperatureK(),
          bracket.getUpperTemperatureState().getTemperatureK(), bracket.getTemperatureWidthK());
    }
    for (HydrocarbonWaterStableRegionBoundaryCorrector.Result correction : stable.getCorrections()) {
      System.out.printf("  correction P=%.9f converged=%s roots=%d accepted=%d failure=%s%n",
          correction.getBracket().getPressureBara(), correction.isConverged(),
          correction.getRootSet() == null ? 0 : correction.getRootSet().getRoots().size(),
          correction.getAcceptedOrdinaryRootCount(), correction.getFailureMessage());
      if (correction.getRootSet() != null) {
        System.out.printf("    scan TPD=[%.8g, %.8g] usable=%d flash=%d stability=%d%n",
            correction.getRootSet().getMinimumTangentPlaneDistance(),
            correction.getRootSet().getMaximumTangentPlaneDistance(),
            correction.getRootSet().getUsableStabilityStates(), correction.getRootSet().getFlashEvaluations(),
            correction.getRootSet().getStabilityEvaluations());
      }
    }
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.Branch branch : discovery.getBranches()) {
      System.out.printf("Fluid 2 stable branch=%s family=%s anchors=%d continuation=%s%n", branch.getIdentifier(),
          branch.getFamily(), branch.getPoints().size(), branch.canSeedContinuation());
      for (HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint anchor : branch.getPoints()) {
        System.out.printf("  stable anchor T=%.9f K P=%.9f bara beta=%.9g minimumTPD=%.8g%n", anchor.getTemperatureK(),
            anchor.getPressureBara(), anchor.getBeta(), new HydrocarbonWaterBoundaryGlobalStabilityGate(fluid)
                .evaluate(anchor.getBoundaryRoot()).getMinimumNonTrivialTangentPlaneDistance());
      }
    }
    assertTrue(stable.getStableRegionScan().getBrackets().size() >= 3,
        "Each pressure level must expose at least one stable adjacent-region bracket");
    assertTrue(discovery.getCorrectedAnchorCount() >= 3,
        "Stable topology brackets must produce globally stable strict anchors");
    assertTrue(
        discovery.getBranches().stream().anyMatch(HydrocarbonWaterBoundaryAnchorDiscoverer.Branch::canSeedContinuation),
        "Pressure-adjacent stable anchors must seed continuation");
    assertTrue(
        discovery.getEndpointCandidates().stream().noneMatch(candidate -> candidate.getClassification()
            .getClassification() == HydrocarbonWaterBoundaryEndpointClassifier.Classification.METASTABLE_BOUNDARY_ROOT),
        "Stable-bracket discovery must not reintroduce the rejected interior metastable roots");

    HydrocarbonWaterBoundaryBranchAssembler.Result assembled = new HydrocarbonWaterBoundaryBranchAssembler(fluid)
        .setCorrectorControls(60, 1.0e-8, 2.0e-5).setStepControls(0.05, 0.005, 0.1, 8, 0.35).assemble(discovery, 2);
    assertEquals(1, assembled.getBranches().size());
    HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch = assembled.getBranches().get(0);
    assertTrue(branch.getBackwardTrace().hasCompletedRequestedPoints(), branch.getBackwardTrace().getFailureMessage());
    assertTrue(branch.getForwardTrace().hasCompletedRequestedPoints(), branch.getForwardTrace().getFailureMessage());
    assertEquals(7, branch.getEvidencePoints().size());
    assertEquals(2, branch.getBackwardTrace().getGlobalStabilityEvidence().stream()
        .filter(HydrocarbonWaterBoundaryGlobalStabilityGate.Result::isAccepted).count());
    assertEquals(2, branch.getForwardTrace().getGlobalStabilityEvidence().stream()
        .filter(HydrocarbonWaterBoundaryGlobalStabilityGate.Result::isAccepted).count());
  }

  @Test
  @Tag("slow")
  void fluidTwoPreservesBroadPressureStableTopologyInventory() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperatures = new double[] {196.75, 233.15, 273.15, 323.15, 373.15, 423.15, 473.15, 498.15, 523.15,
        548.15, 573.15, 623.15, 673.15, 723.15, 753.15};
    double[] pressures = new double[] {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0, 150.0, 200.0, 250.0, 300.0,
        400.0, 500.0};

    HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).setCorrectionControls(32, 80, 1.0e-5, 1.0e-8)
        .discoverFromStableRegionTransitions(temperatures, pressures);
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = stable.getDiscovery();

    System.out.printf(
        "Fluid 2 broad stable topology: states=%d brackets=%d failures=%d corrections=%d anchors=%d branches=%d isolated=%d withheld=%d%n",
        stable.getStableRegionScan().getStates().size(), stable.getStableRegionScan().getBrackets().size(),
        stable.getStableRegionScan().getFailures().size(), stable.getCorrections().size(),
        discovery.getCorrectedAnchorCount(), discovery.getBranches().size(), discovery.getIsolatedBranchCount(),
        discovery.getEndpointCandidateCount());
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.Branch branch : discovery.getBranches()) {
      System.out.printf("  broad branch=%s family=%s anchors=%d continuation=%s%n", branch.getIdentifier(),
          branch.getFamily(), branch.getPoints().size(), branch.canSeedContinuation());
      for (HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint anchor : branch.getPoints()) {
        assertTrue(anchor.getGlobalStabilityResult().isAccepted());
        System.out.printf("    P=%.9f bara T=%.9f K beta=%.9g%n", anchor.getPressureBara(), anchor.getTemperatureK(),
            anchor.getBeta());
      }
    }
    for (HydrocarbonWaterStableRegionBoundaryCorrector.Result correction : stable.getCorrections()) {
      System.out.printf("  broad correction family=%s P=%.9f accepted=%d intervals=%d failure=%s%n",
          correction.getBracket().getFamily(), correction.getBracket().getPressureBara(),
          correction.getAcceptedOrdinaryRootCount(), correction.getIntervalAttempts().size(),
          correction.getFailureMessage());
    }
    assertEquals(stable.getStableRegionScan().getBrackets().size(), stable.getCorrections().size());
    assertTrue(discovery.getCorrectedAnchorCount() > 0);
    assertTrue(discovery.getBranches().stream()
        .anyMatch(HydrocarbonWaterBoundaryAnchorDiscoverer.Branch::canSeedContinuation));
    assertEquals(1, discovery.getBranches().size(),
        "A sufficiently sampled pressure grid must join the same composition-continuous physical branch");
  }

  @Test
  @Tag("slow")
  void fluidTwoLongStableGoToGowContinuationReportsBothTerminations() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperatures = new double[] {196.75, 233.15, 273.15, 323.15, 373.15, 423.15, 473.15, 498.15, 523.15,
        548.15, 573.15, 623.15, 673.15, 723.15, 753.15};
    double[] pressures = new double[] {0.5, 1.0, 2.0, 5.0, 100.0, 150.0, 200.0, 250.0};
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).discoverFromStableRegionTransitions(temperatures, pressures)
        .getDiscovery();
    assertEquals(2, discovery.getBranches().size(),
        "The intentionally sparse middle pressure gap must preserve two seedable partial segments");

    HydrocarbonWaterBoundaryAnchorDiscoverer.Branch lowPressure = discovery.getBranches().get(0);
    HydrocarbonWaterBoundaryAnchorDiscoverer.Branch highPressure = discovery.getBranches().get(1);
    TwoToThreePhasePseudoArcLengthTracer.Result reverse = new TwoToThreePhasePseudoArcLengthTracer(fluid,
        lowPressure.getFamily().getRetainedPhaseZero(), lowPressure.getFamily().getRetainedPhaseOne(),
        lowPressure.getFamily().getIncipientPhase()).setCorrectorControls(80, 1.0e-8, 2.0e-5)
        .setStepControls(0.5, 1.0e-4, 1.0, 16).setMaximumCompositionJump(0.35)
        .trace(lowPressure.getPoints().get(1).toContinuationState(),
            lowPressure.getPoints().get(0).toContinuationState(), 200);
    int highLast = highPressure.getPoints().size() - 1;
    TwoToThreePhasePseudoArcLengthTracer.Result forward = new TwoToThreePhasePseudoArcLengthTracer(fluid,
        highPressure.getFamily().getRetainedPhaseZero(), highPressure.getFamily().getRetainedPhaseOne(),
        highPressure.getFamily().getIncipientPhase()).setCorrectorControls(80, 1.0e-8, 2.0e-5)
        .setStepControls(0.5, 1.0e-4, 1.0, 16).setMaximumCompositionJump(0.35)
        .trace(highPressure.getPoints().get(highLast - 1).toContinuationState(),
            highPressure.getPoints().get(highLast).toContinuationState(), 200);

    printTermination("reverse", reverse);
    printTermination("forward", forward);
    assertTrue(reverse.getPoints().size() >= 3);
    assertTrue(forward.getPoints().size() >= 3);
    assertTrue(!reverse.hasCompletedRequestedPoints(), "The reverse trace must reach a classified numerical boundary");
    assertTrue(!forward.hasCompletedRequestedPoints(), "The forward trace must reach a classified numerical boundary");
    assertEquals(reverse.getAcceptedCorrections().size(), reverse.getGlobalStabilityEvidence().stream()
        .filter(HydrocarbonWaterBoundaryGlobalStabilityGate.Result::isAccepted).count());
    assertEquals(forward.getAcceptedCorrections().size(), forward.getGlobalStabilityEvidence().stream()
        .filter(HydrocarbonWaterBoundaryGlobalStabilityGate.Result::isAccepted).count());
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejected = forward.getGlobalStabilityEvidence()
        .get(forward.getGlobalStabilityEvidence().size() - 1);
    assertTrue(!rejected.isAccepted());
    System.out.printf("  forward rejected root=(%.9f K, %.9f bara) target=%s minimumTPD=%.9g unstablePhase=%s%n",
        rejected.getBoundaryRoot().getTemperatureK(), rejected.getBoundaryRoot().getPressureBara(),
        rejected.getBoundaryRoot().getIncipientPhase(), rejected.getMinimumNonTrivialTangentPlaneDistance(),
        rejected.getMostUnstableCandidate() == null ? "none" : rejected.getMostUnstableCandidate().getPhase());
    TwoToThreePhaseArcLengthCorrector.State last = forward.getPoints().get(forward.getPoints().size() - 1);
    System.out.printf("  forward last phase distances GO=%.9g GW=%.9g OW=%.9g%n",
        distance(last.getPhaseZeroComposition(), last.getPhaseOneComposition()),
        distance(last.getPhaseZeroComposition(), last.getIncipientComposition()),
        distance(last.getPhaseOneComposition(), last.getIncipientComposition()));
    int physicalCriticalEndpoints = 0;
    for (CandidatePhase criticalPhase : new CandidatePhase[] {CandidatePhase.GAS, CandidatePhase.OIL}) {
      HydrocarbonWaterCriticalEndpointSolver.Result critical = new HydrocarbonWaterCriticalEndpointSolver(fluid,
          last.getRetainedPhaseZero(), last.getRetainedPhaseOne(), last.getIncipientPhase(), criticalPhase)
          .setNumericalControls(24, 1.0e-6, 1.0e-9, 1.0e-3, 2.0e-4, 1.0e-7, 1.0e-4).solve(last);
      System.out.printf(
          "  Fluid 2 GO->GOW CEP criticalPhase=%s physical=%s T=%.9f K P=%.9f bara eigen=%.9g third=%.9g gradient=%.9g distance=%.9g failure=%s%n",
          criticalPhase, critical.isPhysicalEndpoint(), critical.getTemperatureK(), critical.getPressureBara(),
          critical.getMinimumEigenvalue(), critical.getThirdDirectionalDerivative(),
          critical.getHomogeneousTpdGradient(), critical.getPhaseDistance(), critical.getFailureMessage());
      physicalCriticalEndpoints += critical.isPhysicalEndpoint() ? 1 : 0;
    }
    assertEquals(0, physicalCriticalEndpoints,
        "A retained-phase spinodal with distinct gas/oil phases must not be mislabeled as a CEP");
    assertTrue(distance(last.getPhaseZeroComposition(), last.getPhaseOneComposition()) > 1.0e-3);
    HydrocarbonWaterBoundaryGlobalStabilityGate.RetainedPhaseLocalStability retainedOil = rejected
        .getRetainedPhaseStability().stream().filter(local -> local.getPhase() == CandidatePhase.OIL).findFirst()
        .orElseThrow();
    assertTrue(retainedOil.getMinimumEigenvalue() < 0.0);
    HydrocarbonWaterBoundaryEventDetector detector = new HydrocarbonWaterBoundaryEventDetector();
    HydrocarbonWaterBoundaryEventDetector.Result events = detector
        .attachGlobalStabilityLimit(detector.detect(forward.getPoints()), forward);
    HydrocarbonWaterBoundaryEventDetector.Event spinodal = events
        .getEvent(HydrocarbonWaterPhaseTopology.SpecialPointType.RETAINED_PHASE_SPINODAL);
    assertNotNull(spinodal);
    assertEquals(CandidatePhase.OIL, spinodal.getDestabilizingPhase());
    assertEquals(rejected.getBoundaryRoot().getTemperatureK(), spinodal.getTemperatureK(), 0.0);
    assertEquals(rejected.getBoundaryRoot().getPressureBara(), spinodal.getPressureBara(), 0.0);
    assertEquals(Math.abs(retainedOil.getMinimumEigenvalue()), spinodal.getQualityMeasure(), 0.0);
    assertTrue(spinodal.getDiagnostic().contains("stationary Jacobian"));
  }

  @Test
  @Tag("slow")
  void fluidTwoRefinesHighPressureRetainedOilSpinodal() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperatures = new double[] {423.15, 448.15, 473.15, 493.15, 503.15, 508.15, 513.15, 523.15, 548.15,
        573.15};
    double[] pressures = new double[] {50.0, 100.0, 150.0, 200.0, 225.0, 240.0, 250.0, 260.0, 270.0, 275.0, 280.0,
        282.0, 284.0, 285.0, 287.0, 290.0, 295.0, 300.0};
    HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).setCorrectionControls(32, 80, 1.0e-5, 1.0e-8)
        .discoverFromStableRegionTransitions(temperatures, pressures);
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = stable.getDiscovery();
    System.out.printf("Fluid 2 high-P corrected topology brackets=%d anchors=%d branches=%d withheld=%d%n",
        stable.getStableRegionScan().getBrackets().size(), discovery.getCorrectedAnchorCount(),
        discovery.getBranches().size(), discovery.getEndpointCandidateCount());
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.Branch branch : discovery.getBranches()) {
      System.out.printf("  high-P corrected branch=%s anchors=%d%n", branch.getIdentifier(), branch.getPoints().size());
      for (HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint anchor : branch.getPoints()) {
        System.out.printf("    P=%.9f bara T=%.9f K beta=%.9g%n", anchor.getPressureBara(), anchor.getTemperatureK(),
            anchor.getBeta());
      }
    }
    for (HydrocarbonWaterStableRegionBoundaryCorrector.Result correction : stable.getCorrections()) {
      System.out.printf("  high-P corrected attempt P=%.9f accepted=%d intervals=%d failure=%s%n",
          correction.getBracket().getPressureBara(), correction.getAcceptedOrdinaryRootCount(),
          correction.getIntervalAttempts().size(), correction.getFailureMessage());
    }
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.EndpointCandidate candidate : discovery.getEndpointCandidates()) {
      HydrocarbonWaterBoundaryEndpointClassifier.Result classification = candidate.getClassification();
      System.out.printf("  high-P withheld P=%.9f T=%.9f class=%s minimumTPD=%s failure=%s%n",
          classification.getBoundaryRoot().getPressureBara(), classification.getBoundaryRoot().getTemperatureK(),
          classification.getClassification(),
          classification.getGlobalStabilityResult() == null ? "n/a"
              : Double.toString(classification.getGlobalStabilityResult().getMinimumNonTrivialTangentPlaneDistance()),
          classification.getFailureMessage());
      if (classification.getGlobalStabilityResult() != null) {
        for (HydrocarbonWaterBoundaryGlobalStabilityGate.RetainedPhaseLocalStability local : classification
            .getGlobalStabilityResult().getRetainedPhaseStability()) {
          System.out.printf("    retained=%s minimumEigenvalue=%.9g localAccepted=%s%n", local.getPhase(),
              local.getMinimumEigenvalue(), local.isAccepted());
        }
      }
    }
    List<TwoToThreePhaseBoundaryPointSolver.Result> curvatureRoots = new ArrayList<>();
    for (HydrocarbonWaterStableRegionBoundaryCorrector.Result correction : stable.getCorrections()) {
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : correction.getAllClassifications()) {
        if (classification.getGlobalStabilityResult() != null) {
          curvatureRoots.add(classification.getBoundaryRoot());
        }
      }
    }
    curvatureRoots.sort((first, second) -> Double.compare(first.getPressureBara(), second.getPressureBara()));
    assertTrue(curvatureRoots.size() >= 2);
    TwoToThreePhaseBoundaryPointSolver.Result curvatureTemplate = curvatureRoots.get(0);
    HydrocarbonWaterRetainedPhaseSpinodalSolver spinodalSolver = new HydrocarbonWaterRetainedPhaseSpinodalSolver(fluid,
        curvatureTemplate.getRetainedPhaseZero(), curvatureTemplate.getRetainedPhaseOne(),
        curvatureTemplate.getIncipientPhase(), CandidatePhase.OIL)
        .setNumericalControls(40, 80, 1.0e-8, 2.0e-5, 2.0e-4, 1.0e-4, 1.0e-6);
    TwoToThreePhaseBoundaryPointSolver.Result spinodalFirst = null;
    TwoToThreePhaseBoundaryPointSolver.Result spinodalSecond = null;
    TwoToThreePhaseBoundaryPointSolver.Result previousCurvatureRoot = null;
    double previousEigenvalue = Double.NaN;
    for (TwoToThreePhaseBoundaryPointSolver.Result root : curvatureRoots) {
      HydrocarbonWaterRetainedPhaseSpinodalSolver.CurvatureState curvature = spinodalSolver.analyze(root);
      System.out.printf("  retained-OIL Jacobian P=%.9f bara T=%.9f K eigen=%.9g residual=%.9g antisym=%.9g%n",
          root.getPressureBara(), root.getTemperatureK(), curvature.getMinimumEigenvalue(),
          curvature.getCurvature().getMaximumResidual(), curvature.getCurvature().getMaximumAntisymmetry());
      if (previousCurvatureRoot != null
          && Math.copySign(1.0, previousEigenvalue) != Math.copySign(1.0, curvature.getMinimumEigenvalue())) {
        spinodalFirst = previousCurvatureRoot;
        spinodalSecond = root;
        break;
      }
      previousCurvatureRoot = root;
      previousEigenvalue = curvature.getMinimumEigenvalue();
    }
    assertNotNull(spinodalFirst, "the pressure inventory must bracket the retained-oil spinodal");
    HydrocarbonWaterRetainedPhaseSpinodalSolver.Result spinodal = spinodalSolver.solve(spinodalFirst, spinodalSecond);
    System.out.printf(
        "Fluid 2 retained-OIL spinodal converged=%s P=%.9f bara T=%.9f K eigen=%.9g bracket=[%.9f, %.9f] width=%.9g iterations=%d failure=%s%n",
        spinodal.isConverged(), spinodal.getSpinodal().getState().getPressureBara(),
        spinodal.getSpinodal().getState().getTemperatureK(), spinodal.getSpinodal().getMinimumEigenvalue(),
        spinodal.getLowerPressure().getState().getPressureBara(),
        spinodal.getUpperPressure().getState().getPressureBara(), spinodal.getPressureBracketWidthBara(),
        spinodal.getIterations(), spinodal.getFailureMessage());
    assertTrue(spinodal.isConverged(), spinodal.getFailureMessage());
    HydrocarbonWaterRetainedPhaseSpinodalBranchInitializer.Result postSpinodal = new HydrocarbonWaterRetainedPhaseSpinodalBranchInitializer(
        fluid, CandidatePhase.OIL).setNumericalControls(5.0, 100, 1.0e-8, 2.0e-5, 200, 1.0e-9, 2.0e-5, 1.0e-4)
        .initialize(spinodal);
    System.out.printf("Fluid 2 post-spinodal G/O1/O2 initialized=%s attempts=%d failure=%s%n",
        postSpinodal.isInitialized(), postSpinodal.getAttempts().size(), postSpinodal.getFailureMessage());
    for (HydrocarbonWaterRetainedPhaseSpinodalBranchInitializer.Attempt attempt : postSpinodal.getAttempts()) {
      SpecifiedThreePhaseFlashSolver.Result flash = attempt.getFlash();
      System.out.printf(
          "  split=%.6g modes=(%+.6g,%+.6g) preconditioned=%s gibbsError=%.9g flash=%s thermo=%s balance=%s oilDistance=%s usable=%s failure=%s%n",
          attempt.getSplitFraction(), attempt.getFirstModeAmplitude(), attempt.getSecondModeAmplitude(),
          attempt.getPreconditioned().isConverged(), attempt.getPreconditioned().getFinalError(),
          flash != null && flash.isConverged(),
          flash == null ? "n/a" : Double.toString(flash.getThermodynamicMaximumResidual()),
          flash == null ? "n/a" : Double.toString(flash.getMaterialBalanceResidual()),
          flash == null ? "n/a" : Double.toString(flash.getPhaseCompositionDistance(1, 2)), attempt.isUsable(),
          flash == null ? attempt.getPreconditioned().getFailureMessage() : flash.getFailureMessage());
    }
    assertFalse(postSpinodal.isInitialized(),
        "Fluid 2 must not promote the water-rich third slot returned by the post-spinodal flash to a second oil "
            + "phase: " + postSpinodal.getFailureMessage());
    assertTrue(discovery.getCorrectedAnchorCount() > 0);
    assertTrue(discovery.getBranches().stream()
        .anyMatch(HydrocarbonWaterBoundaryAnchorDiscoverer.Branch::canSeedContinuation));
    assertTrue(discovery.getEndpointCandidateCount() >= 1);
    HydrocarbonWaterBoundaryAnchorDiscoverer.Branch branch = discovery.getBranches().get(0);
    int lastIndex = branch.getPoints().size() - 1;
    TwoToThreePhasePseudoArcLengthTracer.Result forward = new TwoToThreePhasePseudoArcLengthTracer(fluid,
        branch.getFamily().getRetainedPhaseZero(), branch.getFamily().getRetainedPhaseOne(),
        branch.getFamily().getIncipientPhase()).setCorrectorControls(80, 1.0e-8, 2.0e-5)
        .setStepControls(0.1, 1.0e-5, 0.25, 20).setMaximumCompositionJump(0.35)
        .trace(branch.getPoints().get(lastIndex - 1).toContinuationState(),
            branch.getPoints().get(lastIndex).toContinuationState(), 100);
    printTermination("high-pressure refined forward", forward);
    assertTrue(!forward.hasCompletedRequestedPoints());
    assertTrue(forward.getPoints().get(forward.getPoints().size() - 1)
        .getPressureBara() >= spinodal.getLowerPressure().getState().getPressureBara() - 5.0);
    assertTrue(forward.getPoints().get(forward.getPoints().size() - 1)
        .getPressureBara() <= spinodal.getUpperPressure().getState().getPressureBara() + 5.0);
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejected = forward.getGlobalStabilityEvidence()
        .get(forward.getGlobalStabilityEvidence().size() - 1);
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result precedingAccepted = null;
    for (HydrocarbonWaterBoundaryGlobalStabilityGate.Result evidence : forward.getGlobalStabilityEvidence()) {
      if (evidence.isAccepted()) {
        precedingAccepted = evidence;
      }
    }
    if (precedingAccepted == null) {
      precedingAccepted = branch.getPoints().get(lastIndex).getGlobalStabilityResult();
    }
    assertNotNull(precedingAccepted);
    System.out.printf(
        "Fluid 2 stability bracket accepted=(%.9f K, %.9f bara, minTPD=%.9g) rejected=(%.9f K, %.9f bara, minTPD=%.9g)%n",
        precedingAccepted.getBoundaryRoot().getTemperatureK(), precedingAccepted.getBoundaryRoot().getPressureBara(),
        precedingAccepted.getMinimumNonTrivialTangentPlaneDistance(), rejected.getBoundaryRoot().getTemperatureK(),
        rejected.getBoundaryRoot().getPressureBara(), rejected.getMinimumNonTrivialTangentPlaneDistance());
    if (rejected.getMostUnstableCandidate() == null
        || rejected.getMostUnstableCandidate().getPhase() != CandidatePhase.OIL
        || rejected.getMostUnstableCandidate().getTangentPlaneDistance() >= -1.0e-3) {
      HydrocarbonWaterBoundaryGlobalStabilityGate.RetainedPhaseLocalStability retainedOil = rejected
          .getRetainedPhaseStability().stream().filter(local -> local.getPhase() == CandidatePhase.OIL).findFirst()
          .orElseThrow();
      System.out.printf("Fluid 2 continuation stopped first at retained-OIL local stability eigenvalue=%.9g%n",
          retainedOil.getMinimumEigenvalue());
      assertTrue(!retainedOil.isAccepted());
      assertTrue(retainedOil.getMinimumEigenvalue() < -1.0e-6);
      return;
    }
    double[] rejectedGas = rejected.getBoundaryRoot().getPhaseZeroComposition();
    double[] rejectedOil = rejected.getBoundaryRoot().getPhaseOneComposition();
    double[] rejectedWater = rejected.getBoundaryRoot().getIncipientComposition();
    double[] alternateOil = rejected.getMostUnstableCandidate().getComposition();
    for (int componentIndex = 0; componentIndex < alternateOil.length; componentIndex++) {
      System.out.printf("  component=%s gas=%.9g oil=%.9g water=%.9g alternateOil=%.9g%n",
          fluid.getPhase(0).getComponent(componentIndex).getComponentName(), rejectedGas[componentIndex],
          rejectedOil[componentIndex], rejectedWater[componentIndex], alternateOil[componentIndex]);
    }
    List<TwoToThreePhaseBoundaryPointSolver.Result> trackedBoundaryRoots = new ArrayList<>();
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint anchor : branch.getPoints()) {
      trackedBoundaryRoots.add(anchor.getBoundaryRoot());
    }
    for (TwoToThreePhaseArcLengthCorrector.Result correction : forward.getAcceptedCorrections()) {
      trackedBoundaryRoots.add(TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(correction.getState(),
          correction.getThermodynamicMaximumResidual()));
    }
    trackedBoundaryRoots.add(rejected.getBoundaryRoot());
    HydrocarbonWaterSecondaryStationaryBranchTracker.Result tracked = new HydrocarbonWaterSecondaryStationaryBranchTracker(
        fluid).setNumericalControls(100, 500, 1.0e-9, 2.0e-5, 0.2, 0.25, 1.0e-7)
        .trackBackward(trackedBoundaryRoots, CandidatePhase.OIL, alternateOil);
    System.out.printf(
        "Fluid 2 explicit secondary-OIL tracking samples=%d zeroBrackets=%d branchLimits=%d reachesFirst=%s failure=%s%n",
        tracked.getSamples().size(), tracked.getZeroTpdBrackets().size(), tracked.getBranchLimitBrackets().size(),
        tracked.reachesFirstBoundaryRoot(), tracked.getTerminationMessage());
    for (HydrocarbonWaterSecondaryStationaryBranchTracker.Sample sample : tracked.getSamples()) {
      System.out.printf("  tracked P=%.9f bara T=%.9f K TPD=%.9g jump=%.9g accepted=%s failure=%s%n",
          sample.getBoundaryRoot().getPressureBara(), sample.getBoundaryRoot().getTemperatureK(),
          sample.getTangentPlaneDistance(), sample.getCompositionJump(), sample.isTracked(),
          sample.getFailureMessage());
    }
    HydrocarbonWaterSecondaryStationaryBranchTracker.Sample terminalTracked = tracked.getSamples()
        .get(tracked.getSamples().size() - 1);
    assertTrue(terminalTracked.isTracked(), terminalTracked.getFailureMessage());
    assertTrue(terminalTracked.getTangentPlaneDistance() < -1.0e-3);
    assertTrue(distance(terminalTracked.getCandidate().getComposition(), alternateOil) < 1.0e-5);
    assertEquals(1, tracked.getBranchLimitBrackets().size());
    HydrocarbonWaterRetainedPhaseBranchSwitcher.Result switched = new HydrocarbonWaterRetainedPhaseBranchSwitcher(fluid)
        .setBoundaryControls(32, 80, 1.0e-5, 1.0e-8, 1.0e-4).switchAndCorrect(rejected,
            rejected.getBoundaryRoot().getTemperatureK() - 40.0, rejected.getBoundaryRoot().getTemperatureK() + 40.0);
    long convergedFlashes = switched.getAttempts().stream().filter(attempt -> attempt.getFlash().isConverged()).count();
    long distinctFlashes = switched.getAttempts().stream()
        .filter(HydrocarbonWaterRetainedPhaseBranchSwitcher.Attempt::isDistinctBranch).count();
    int correctedRoots = switched.getCorrections().stream()
        .mapToInt(correction -> correction.getRootSet().getRoots().size()).sum();
    System.out.printf(
        "Fluid 2 retained-phase switch replacement=%s attempts=%d converged=%d distinctAttempts=%d distinctBranches=%d correctedRoots=%d accepted=%d failure=%s%n",
        switched.getReplacementPhase(), switched.getAttempts().size(), convergedFlashes, distinctFlashes,
        switched.getDistinctRetainedBranchCount(), correctedRoots, switched.getAcceptedBoundaryRoots().size(),
        switched.getFailureMessage());
    for (int correctionIndex = 0; correctionIndex < switched.getCorrections().size(); correctionIndex++) {
      HydrocarbonWaterRetainedPhaseBranchSwitcher.Correction correction = switched.getCorrections()
          .get(correctionIndex);
      System.out.printf("  switched branch=%d beta=%.9g roots=%d scanFailure=%s%n", correctionIndex,
          correction.getRetainedBranch().getBeta(), correction.getRootSet().getRoots().size(),
          correction.getRootSet().getFailureMessage());
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : correction.getClassifications()) {
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = classification.getGlobalStabilityResult();
        System.out.printf("    root T=%.9f K class=%s minTPD=%s unstable=%s failure=%s%n",
            classification.getBoundaryRoot().getTemperatureK(), classification.getClassification(),
            stability == null ? "n/a" : Double.toString(stability.getMinimumNonTrivialTangentPlaneDistance()),
            stability == null || stability.getMostUnstableCandidate() == null ? "none"
                : stability.getMostUnstableCandidate().getPhase(),
            classification.getFailureMessage());
      }
    }
    assertEquals(CandidatePhase.OIL, switched.getReplacementPhase());
    assertTrue(switched.getAttempts().size() >= 28);
    assertTrue(convergedFlashes > 0);
    TwoToThreePhaseBoundaryPointSolver.Result coalescenceCandidate = switched.getCorrections().get(0)
        .getClassifications().get(0).getBoundaryRoot();
    int retainedPairEndpoints = 0;
    for (CandidatePhase criticalPhase : new CandidatePhase[] {CandidatePhase.GAS, CandidatePhase.OIL}) {
      HydrocarbonWaterRetainedPairCriticalEndpointSolver.Result endpoint = new HydrocarbonWaterRetainedPairCriticalEndpointSolver(
          fluid, coalescenceCandidate.getRetainedPhaseZero(), coalescenceCandidate.getRetainedPhaseOne(),
          coalescenceCandidate.getIncipientPhase(), criticalPhase)
          .setNumericalControls(24, 1.0e-7, 1.0e-3, 2.0e-4, 1.0e-7, 1.0e-4).solve(coalescenceCandidate);
      System.out.printf(
          "  retained-pair CEP root=%s physical=%s T=%.9f K P=%.9f bara eigen=%.9g third=%.9g gradient=%.9g aqueousTPD=%.9g aqueousDistance=%.9g failure=%s%n",
          criticalPhase, endpoint.isPhysicalEndpoint(), endpoint.getTemperatureK(), endpoint.getPressureBara(),
          endpoint.getMinimumEigenvalue(), endpoint.getThirdDirectionalDerivative(),
          endpoint.getHomogeneousTpdGradient(), endpoint.getIncipientTangentPlaneDistance(),
          endpoint.getIncipientPhaseDistance(), endpoint.getFailureMessage());
      retainedPairEndpoints += endpoint.isPhysicalEndpoint() ? 1 : 0;
    }
    assertTrue(retainedPairEndpoints <= 1,
        "gas- and oil-root evaluations may not create duplicate retained-pair critical endpoints");
  }

  @Test
  @Tag("slow")
  void fluidTwoInventoriesHighPressureStableRegionsForBranchSwitching() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperatures = new double[] {423.15, 448.15, 473.15, 493.15, 503.15, 508.15, 513.15, 523.15, 548.15,
        573.15, 623.15, 673.15, 723.15, 753.15};
    double[] pressures = new double[] {250.0, 275.0, 285.0, 300.0, 325.0, 350.0, 400.0, 500.0};
    HydrocarbonWaterStableRegionTransitionScanner.Result result = new HydrocarbonWaterStableRegionTransitionScanner(
        fluid).setNumericalControls(80, 1.0e-4, 1.0e-9).scan(temperatures, pressures);

    for (HydrocarbonWaterStableRegionTransitionScanner.StableState state : result.getStates()) {
      System.out.printf("Fluid 2 high-P stable state P=%.3f bara T=%.6f K region=%s fractions=(G=%.8g O=%.8g W=%.8g)%n",
          state.getPressureBara(), state.getTemperatureK(), state.getRegion(),
          state.getPhaseFraction(CandidatePhase.GAS), state.getPhaseFraction(CandidatePhase.OIL),
          state.getPhaseFraction(CandidatePhase.AQUEOUS));
    }
    for (HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket : result.getBrackets()) {
      System.out.printf("  high-P transition P=%.3f family=%s T=%.9f K%n", bracket.getPressureBara(),
          bracket.getFamily(), bracket.getMidpointTemperatureK());
    }
    for (HydrocarbonWaterStableRegionTransitionScanner.Failure failure : result.getFailures()) {
      System.out.printf("  high-P diagnostic P=%.3f T=%.6f K reason=%s%n", failure.getPressureBara(),
          failure.getTemperatureK(), failure.getReason());
    }
    assertTrue(result.getStates().size() > 0);
  }

  @Test
  @Tag("slow")
  void fluidTwoChecksFourPhaseFlashAtHighPressureGlobalStabilityLimit() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    fluid.setTemperature(503.985596882);
    fluid.setPressure(284.124866867);
    fluid.setMultiPhaseCheck(true);
    fluid.setMaxNumberOfPhases(4);
    RuntimeException failure = assertThrows(RuntimeException.class, () -> new ThermodynamicOperations(fluid).TPflash());
    System.out.printf("Fluid 2 native four-phase flash unavailable at T=%.9f K P=%.9f bara: %s%n",
        fluid.getTemperature(), fluid.getPressure(), failure.toString());
  }

  private static void printTermination(String direction, TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    TwoToThreePhaseArcLengthCorrector.State last = trace.getPoints().get(trace.getPoints().size() - 1);
    System.out.printf(
        "Fluid 2 %s long GO->GOW trace points=%d attempts=%d rejected=%d last=(%.9f K, %.9f bara) beta=%.9g reason=%s failure=%s globalTrials=%d%n",
        direction, trace.getPoints().size(), trace.getAttemptedCorrections(), trace.getRejectedCorrections(),
        last.getTemperatureK(), last.getPressureBara(), last.getBeta(), trace.getTerminationReason(),
        trace.getFailureMessage(), trace.getGlobalStabilityEvidence().size());
  }

  private static double distance(double[] first, double[] second) {
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }
}
