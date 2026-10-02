package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

class HydrocarbonWaterBoundaryAnchorDiscovererTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void fluidTwoPreservesEveryCorrectedRootAcrossAllBoundaryFamiliesAtFiftyBara() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result result = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .setCorrectionControls(32, 60, 2.0e-5, 1.0e-8)
        .discover(new double[] {180.0, 220.0, 260.0, 300.0, 350.0, 425.0, 500.0}, new double[] {50.0}, 100.0, 650.0);

    int correctedRoots = 0;
    int ordinaryAnchors = 0;
    int endpointCandidates = 0;
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.CorrectionAttempt attempt : result.getAttempts()) {
      correctedRoots += attempt.getCorrectedRootCount();
      ordinaryAnchors += attempt.getOrdinaryAnchorCount();
      endpointCandidates += attempt.getEndpointCandidateCount();
      System.out.printf(
          "Fluid 2 P=50 family=%s seeds=%d roots=%d ordinary=%d endpoints=%d flash=%d stability=%d failure=%s%n",
          attempt.getFamily(), attempt.getThreePhaseSeedCount(), attempt.getCorrectedRootCount(),
          attempt.getOrdinaryAnchorCount(), attempt.getEndpointCandidateCount(), attempt.getFlashEvaluations(),
          attempt.getStabilityEvaluations(), attempt.getFailureMessage());
    }
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.EndpointCandidate candidate : result.getEndpointCandidates()) {
      HydrocarbonWaterBoundaryEndpointClassifier.Result classification = candidate.getClassification();
      System.out.printf("  endpoint family=%s class=%s T=%.9f K P=%.9f bara minDistance=%.8g failure=%s%n",
          candidate.getFamily(), classification.getClassification(), classification.getBoundaryRoot().getTemperatureK(),
          classification.getBoundaryRoot().getPressureBara(), classification.getMinimumPhaseDistance(),
          classification.getFailureMessage());
    }
    List<HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint> allAnchors = new ArrayList<>();
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.Branch branch : result.getBranches()) {
      for (HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint anchor : branch.getPoints()) {
        allAnchors.add(anchor);
        System.out.printf("  anchor branch=%s family=%s T=%.9f K P=%.9f bara beta=%.9g%n", branch.getIdentifier(),
            branch.getFamily(), anchor.getTemperatureK(), anchor.getPressureBara(), anchor.getBeta());
        printStableFlashComparison(fluid, anchor);
      }
    }
    printCrossFamilyCoincidences(allAnchors);
    assertEquals(3, result.getAttempts().size());
    assertEquals(21, result.getGridEvaluationCount());
    assertEquals(correctedRoots, ordinaryAnchors + endpointCandidates,
        "Every corrected root must be routed to an ordinary anchor or an endpoint candidate");
    assertEquals(ordinaryAnchors, result.getCorrectedAnchorCount());
    assertEquals(endpointCandidates, result.getEndpointCandidateCount());
    assertTrue(correctedRoots > 0,
        "The observed Fluid 2 three-phase seed domain must yield corrected boundary evidence");
    assertEquals(0, ordinaryAnchors,
        "Every finite-fraction Fluid 2 root at 50 bara lies inside the stable GOW region and must be withheld");
    assertEquals(7, endpointCandidates);
    assertTrue(result.getBranches().isEmpty(), "Metastable roots must not seed an ordinary continuation branch");
    assertEquals(4,
        result.getEndpointCandidates().stream().filter(candidate -> candidate.getClassification()
            .getClassification() == HydrocarbonWaterBoundaryEndpointClassifier.Classification.METASTABLE_BOUNDARY_ROOT)
            .count());
    assertEquals(1,
        result.getEndpointCandidates().stream()
            .filter(candidate -> candidate.getClassification()
                .getClassification() == HydrocarbonWaterBoundaryEndpointClassifier.Classification.CRITICAL_END_POINT)
            .count());
    assertEquals(2, result.getEndpointCandidates().stream().filter(candidate -> candidate.getClassification()
        .getClassification() == HydrocarbonWaterBoundaryEndpointClassifier.Classification.PHASE_COALESCENCE_CANDIDATE)
        .count());
    TwoToThreePhaseBoundaryQualityGate.EnvelopeReport quality = new TwoToThreePhaseBoundaryQualityGate().validate(
        java.util.Collections.<TwoToThreePhaseBoundaryQualityGate.Branch>emptyList(), result.getEndpointCandidates());
    assertTrue(!quality.isAccepted(), "Withheld roots and unattached special points must keep the envelope blocked");
    assertEquals(endpointCandidates, quality.getEndpointCandidateCount());
  }

  @Test
  @Tag("slow")
  void fluidTwoLinksCompositionContinuousRootsAcrossPressureAnchors() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result result = new HydrocarbonWaterBoundaryAnchorDiscoverer(fluid)
        .setCorrectionControls(32, 60, 2.0e-5, 1.0e-8)
        .discover(new double[] {180.0, 220.0, 260.0, 300.0, 350.0, 425.0, 500.0}, new double[] {40.0, 50.0, 60.0},
            100.0, 650.0);

    int preserved = result.getEndpointCandidateCount();
    int continuationBranches = 0;
    for (HydrocarbonWaterBoundaryAnchorDiscoverer.Branch branch : result.getBranches()) {
      preserved += branch.getPoints().size();
      continuationBranches += branch.canSeedContinuation() ? 1 : 0;
      System.out.printf("Fluid 2 linked branch=%s family=%s anchors=%d continuation=%s%n", branch.getIdentifier(),
          branch.getFamily(), branch.getPoints().size(), branch.canSeedContinuation());
      for (HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint anchor : branch.getPoints()) {
        System.out.printf("  T=%.9f K P=%.9f bara beta=%.9g%n", anchor.getTemperatureK(), anchor.getPressureBara(),
            anchor.getBeta());
      }
    }
    int correctedRoots = result.getAttempts().stream()
        .mapToInt(HydrocarbonWaterBoundaryAnchorDiscoverer.CorrectionAttempt::getCorrectedRootCount).sum();
    assertEquals(correctedRoots, preserved);
    assertEquals(0, continuationBranches,
        "Composition-continuous but globally metastable roots must not seed pseudo-arclength continuation");
    assertTrue(result.getEndpointCandidates().stream().anyMatch(candidate -> candidate.getClassification()
        .getClassification() == HydrocarbonWaterBoundaryEndpointClassifier.Classification.METASTABLE_BOUNDARY_ROOT));
  }

  private static void printCrossFamilyCoincidences(List<HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint> anchors) {
    for (int firstIndex = 0; firstIndex < anchors.size(); firstIndex++) {
      for (int secondIndex = firstIndex + 1; secondIndex < anchors.size(); secondIndex++) {
        HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint first = anchors.get(firstIndex);
        HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint second = anchors.get(secondIndex);
        if (first.getFamily() == second.getFamily()
            || Math.abs(first.getTemperatureK() - second.getTemperatureK()) > 1.0e-5
            || Math.abs(first.getPressureBara() - second.getPressureBara()) > 1.0e-5) {
          continue;
        }
        System.out.printf(
            "  cross-family coincidence T=%.9f K families=%s/%s mapped distances G=%.8g O=%.8g W=%.8g "
                + "swappedGO=(%.8g, %.8g)%n",
            first.getTemperatureK(), first.getFamily(), second.getFamily(),
            distance(composition(first, CandidatePhase.GAS), composition(second, CandidatePhase.GAS)),
            distance(composition(first, CandidatePhase.OIL), composition(second, CandidatePhase.OIL)),
            distance(composition(first, CandidatePhase.AQUEOUS), composition(second, CandidatePhase.AQUEOUS)),
            distance(composition(first, CandidatePhase.GAS), composition(second, CandidatePhase.OIL)),
            distance(composition(first, CandidatePhase.OIL), composition(second, CandidatePhase.GAS)));
      }
    }
  }

  private static void printStableFlashComparison(SystemInterface fluid,
      HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint anchor) {
    HydrocarbonWaterBranchSeedScanner.Result stable = new HydrocarbonWaterBranchSeedScanner(fluid,
        CandidatePhase.AQUEOUS).setCompositionDeduplicationTolerance(1.0e-10)
        .scan(new double[] {anchor.getTemperatureK()}, new double[] {anchor.getPressureBara()});
    if (stable.getSeeds().isEmpty()) {
      System.out.printf("    stable TP comparison unavailable failures=%d%n", stable.getFailures().size());
      return;
    }
    HydrocarbonWaterBranchSeedScanner.Seed seed = stable.getSeeds().get(0);
    System.out.printf("    stable TP phases=%d fractions G=%.9g O=%.9g W=%.9g", seed.getStablePhaseCount(),
        seed.getPhaseFraction(CandidatePhase.GAS), seed.getPhaseFraction(CandidatePhase.OIL),
        seed.getPhaseFraction(CandidatePhase.AQUEOUS));
    for (CandidatePhase phase : CandidatePhase.values()) {
      if (seed.hasPhase(phase)) {
        System.out.printf(" d%s=%.8g", phase.name().substring(0, 1),
            distance(composition(anchor, phase), seed.getPhaseComposition(phase)));
      }
    }
    System.out.println();
  }

  private static double[] composition(HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint anchor,
      CandidatePhase phase) {
    if (anchor.getFamily().getRetainedPhaseZero() == phase) {
      return anchor.getPhaseZeroComposition();
    }
    if (anchor.getFamily().getRetainedPhaseOne() == phase) {
      return anchor.getPhaseOneComposition();
    }
    if (anchor.getFamily().getIncipientPhase() == phase) {
      return anchor.getIncipientComposition();
    }
    throw new IllegalArgumentException("phase is not present in anchor topology");
  }

  private static double distance(double[] first, double[] second) {
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }
}
