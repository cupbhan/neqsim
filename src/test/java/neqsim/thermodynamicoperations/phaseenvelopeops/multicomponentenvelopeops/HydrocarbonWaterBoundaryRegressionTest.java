package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import org.junit.jupiter.api.Test;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.AcceptanceCriteria;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.Curve;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.ModelFingerprint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.Point;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.ReferenceSource;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.BoundaryDefinition;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;

class HydrocarbonWaterBoundaryRegressionTest {
  private static final BoundaryDefinition OW_TO_GOW = new BoundaryDefinition(Region.OIL_AQUEOUS,
      Region.GAS_OIL_AQUEOUS);
  private static final ModelFingerprint MODEL = new ModelFingerprint("SRK", "Peneloux", "Huron-Vidal",
      "fluid-22-sha256", "srk-hv-parameter-sha256");

  @Test
  void comparesClosedMultivaluedCurvesAtTheSamePressureAndPassesFrozenThresholds() {
    Curve reference = curve("pvtsim", ReferenceSource.PVTSIM, true, 0.0);
    Curve candidate = curve("neqsim", ReferenceSource.NEQSIM, true, 2.0);
    AcceptanceCriteria criteria = new AcceptanceCriteria(1.0, 1.0, 0.01, 0.01, 0.01, 0.01, true, 0.01, 0.01);

    HydrocarbonWaterBoundaryRegression.Result result = new HydrocarbonWaterBoundaryRegression().compare(reference,
        candidate, criteria);

    assertTrue(result.isComparisonEligible());
    assertTrue(result.isThresholdsFrozen());
    assertFalse(result.isBenchmarkPending());
    assertTrue(result.isAccepted(), result.getViolations().toString());
    assertEquals(1.0, result.getMetrics().getReferenceCoverageFraction(), 1.0e-12);
    assertEquals(1.0, result.getMetrics().getCandidateCoverageFraction(), 1.0e-12);
    assertEquals(2.0 / 300.0, result.getMetrics().getMaximumRelativeTemperatureDifference(), 1.0e-12);
    assertEquals(2.0, result.getMetrics().getCriticalPointMetrics().getAbsoluteTemperatureDifferenceK(), 1.0e-12);
  }

  @Test
  void computesMetricsButCannotPassBeforeThresholdsAreFrozen() {
    HydrocarbonWaterBoundaryRegression.Result result = new HydrocarbonWaterBoundaryRegression().compare(
        curve("paper", ReferenceSource.PAPER, true, 0.0), curve("neqsim", ReferenceSource.NEQSIM, true, 0.0), null);

    assertTrue(result.isComparisonEligible());
    assertFalse(result.isThresholdsFrozen());
    assertTrue(result.isBenchmarkPending());
    assertFalse(result.isAccepted());
    assertTrue(result.getViolations().isEmpty());
  }

  @Test
  void rejectsAmbiguousPvtsimThreeHcWithoutSideFlashTopologyEvidence() {
    HydrocarbonWaterBoundaryRegression.Result result = new HydrocarbonWaterBoundaryRegression().compare(
        curve("pvtsim-3-hc", ReferenceSource.PVTSIM, false, 0.0), curve("neqsim", ReferenceSource.NEQSIM, true, 0.0),
        criteria());

    assertFalse(result.isComparisonEligible());
    assertFalse(result.isAccepted());
    assertTrue(result.getViolations().contains("REFERENCE_TOPOLOGY_UNRESOLVED"));
  }

  @Test
  void rejectsDifferentModelFingerprintBeforeCalculatingCurveError() {
    ModelFingerprint wrongModel = new ModelFingerprint("PR78", "Peneloux", "Huron-Vidal", "fluid-22-sha256",
        "srk-hv-parameter-sha256");
    Curve candidate = new Curve("neqsim-pr", ReferenceSource.NEQSIM, wrongModel, OW_TO_GOW, true,
        "independent global stability gate", points(0.0), new Point(352.0, 202.0));

    HydrocarbonWaterBoundaryRegression.Result result = new HydrocarbonWaterBoundaryRegression()
        .compare(curve("pvtsim", ReferenceSource.PVTSIM, true, 0.0), candidate, criteria());

    assertFalse(result.isComparisonEligible());
    assertTrue(result.getViolations().contains("EQUATION_OF_STATE_MISMATCH"));
    assertEquals(null, result.getMetrics());
  }

  @Test
  void failsWhenSamePressureTemperatureDifferenceExceedsThreshold() {
    HydrocarbonWaterBoundaryRegression.Result result = new HydrocarbonWaterBoundaryRegression().compare(
        curve("pvtsim", ReferenceSource.PVTSIM, true, 0.0), curve("neqsim", ReferenceSource.NEQSIM, true, 20.0),
        criteria());

    assertTrue(result.isComparisonEligible());
    assertFalse(result.isAccepted());
    assertTrue(result.getViolations().contains("MAXIMUM_SAME_PRESSURE_TEMPERATURE_DIFFERENCE_ABOVE_THRESHOLD"));
  }

  @Test
  void rejectsUnresolvedParameterIdentityBeforeCalculatingErrors() {
    ModelFingerprint unresolved = new ModelFingerprint("SRK", "Peneloux", "Huron-Vidal", "fluid-22-sha256",
        "not-proven");
    Curve candidate = new Curve("neqsim", ReferenceSource.NEQSIM, unresolved, OW_TO_GOW, true, "global stability gate",
        points(0.0), null);

    HydrocarbonWaterBoundaryRegression.Result result = new HydrocarbonWaterBoundaryRegression()
        .compare(curve("pvtsim", ReferenceSource.PVTSIM, true, 0.0), candidate, criteria());

    assertFalse(result.isComparisonEligible());
    assertTrue(result.getViolations().contains("PARAMETER_FINGERPRINT_UNRESOLVED"));
    assertEquals(null, result.getMetrics());
  }

  @Test
  void rejectsUnresolvedCompositionIdentityBeforeCalculatingErrors() {
    ModelFingerprint unresolved = new ModelFingerprint("SRK", "Peneloux", "Huron-Vidal", "not-proven",
        "srk-hv-parameter-sha256");
    Curve candidate = new Curve("neqsim", ReferenceSource.NEQSIM, unresolved, OW_TO_GOW, true, "global stability gate",
        points(0.0), null);

    HydrocarbonWaterBoundaryRegression.Result result = new HydrocarbonWaterBoundaryRegression()
        .compare(curve("pvtsim", ReferenceSource.PVTSIM, true, 0.0), candidate, criteria());

    assertFalse(result.isComparisonEligible());
    assertTrue(result.getViolations().contains("COMPOSITION_FINGERPRINT_UNRESOLVED"));
    assertEquals(null, result.getMetrics());
  }

  private static Curve curve(String identifier, ReferenceSource source, boolean topologyResolved,
      double temperatureShiftK) {
    return new Curve(identifier, source, MODEL, OW_TO_GOW, topologyResolved,
        topologyResolved ? "independent phase-side flashes" : "PVTsim 3-HC label only", points(temperatureShiftK),
        new Point(350.0 + temperatureShiftK, 200.0));
  }

  private static java.util.List<Point> points(double shiftK) {
    return Arrays.asList(new Point(300.0 + shiftK, 100.0), new Point(350.0 + shiftK, 200.0),
        new Point(400.0 + shiftK, 100.0), new Point(350.0 + shiftK, 50.0), new Point(300.0 + shiftK, 100.0));
  }

  private static AcceptanceCriteria criteria() {
    return new AcceptanceCriteria(1.0, 1.0, 0.01, 0.01, 0.01, 0.01, false, 0.01, 0.01);
  }
}
