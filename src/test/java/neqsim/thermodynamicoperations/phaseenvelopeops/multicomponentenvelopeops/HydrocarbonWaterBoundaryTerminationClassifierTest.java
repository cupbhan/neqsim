package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryGlobalStabilityGate.RetainedPhaseLocalStability;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryTerminationClassifier.Type;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.Branch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EnvelopeReport;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EvidencePoint;

/** Regression tests for physical versus unresolved continuation termination. */
class HydrocarbonWaterBoundaryTerminationClassifierTest extends NeqSimTest {

  @Test
  @Tag("slow")
  void refinesFluidTwoRetainedOilSpinodalInsteadOfReportingCoverageFailure() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    double[] temperatures = new double[] { 423.15, 448.15, 473.15, 493.15, 503.15, 508.15, 513.15, 523.15, 548.15,
        573.15 };
    double[] pressures = new double[] { 275.0, 280.0 };
    HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult discovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        fluid).setCorrectionControls(32, 80, 1.0e-5, 1.0e-8)
        .discoverFromStableRegionTransitions(temperatures, pressures);

    TwoToThreePhaseBoundaryPointSolver.Result acceptedRoot = null;
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result rejectedEvidence = null;
    for (HydrocarbonWaterStableRegionBoundaryCorrector.Result correction : discovery.getCorrections()) {
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : correction.getAllClassifications()) {
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = classification.getGlobalStabilityResult();
        if (stability == null) {
          continue;
        }
        if (Math.abs(classification.getBoundaryRoot().getPressureBara() - 275.0) < 0.1 && stability.isAccepted()) {
          acceptedRoot = classification.getBoundaryRoot();
        }
        if (Math.abs(classification.getBoundaryRoot().getPressureBara() - 280.0) < 0.1 && !stability.isAccepted()
            && hasUnstableRetainedOil(stability)) {
          rejectedEvidence = stability;
        }
      }
    }

    assertNotNull(acceptedRoot, "275 bara must supply the stable side of the spinodal bracket");
    assertNotNull(rejectedEvidence, "280 bara must supply the retained-oil-unstable side of the spinodal bracket");
    HydrocarbonWaterBoundaryTerminationClassifier.Result termination = new HydrocarbonWaterBoundaryTerminationClassifier(
        fluid).setNumericalControls(40, 80, 1.0e-8, 2.0e-5, 2.0e-4, 1.0e-4, 1.0e-6)
        .classify(acceptedRoot, rejectedEvidence);

    assertEquals(Type.RETAINED_PHASE_SPINODAL, termination.getType(), termination.getDiagnostic());
    assertTrue(termination.isPhysicalEndpoint(), termination.getDiagnostic());
    assertEquals(CandidatePhase.OIL, termination.getDestabilizingPhase());
    assertEquals(504.366971105, termination.getTemperatureK(), 2.0e-3);
    assertEquals(275.169677734, termination.getPressureBara(), 2.0e-3);
    assertTrue(termination.getQualityMeasure() <= 1.1e-6, termination.getDiagnostic());
    assertNotNull(termination.getSpinodalResult());
    assertTrue(termination.getSpinodalResult().getPressureBracketWidthBara() <= 1.0e-3, termination.getDiagnostic());

    Branch physicalSegment = new Branch("fluid2-go-to-gow-physical-termination",
        BoundaryFamily.GO_TO_GOW.getDefinition(), CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS,
        Arrays.asList(EvidencePoint.from(acceptedRoot), EvidencePoint.from(termination.getBoundaryRoot())),
        SpecialPointType.DOMAIN_EXIT, SpecialPointType.RETAINED_PHASE_SPINODAL, null, termination);
    EnvelopeReport report = new TwoToThreePhaseBoundaryQualityGate().validate(Arrays.asList(physicalSegment));
    assertTrue(report.isAccepted(),
        report.getViolations().toString() + " " + report.getBranchReports().get(0).getViolations());
  }

  private static boolean hasUnstableRetainedOil(HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability) {
    for (RetainedPhaseLocalStability retained : stability.getRetainedPhaseStability()) {
      if (retained.getPhase() == CandidatePhase.OIL && !retained.isAccepted()
          && retained.getMinimumEigenvalue() < 0.0) {
        return true;
      }
    }
    return false;
  }
}
