package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.Arrays;
import java.util.LinkedHashSet;
import org.junit.jupiter.api.Test;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.BoundaryDefinition;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;

class HydrocarbonWaterTopologyBenchmarkTest {
  private static final String OW_TO_GOW = new BoundaryDefinition(Region.OIL_AQUEOUS, Region.GAS_OIL_AQUEOUS).getCode();
  private static final String GO_TO_GOW = new BoundaryDefinition(Region.GAS_OIL, Region.GAS_OIL_AQUEOUS).getCode();

  @Test
  void acceptsRequiredPaperAdjacencyButNeverClaimsNumericOrEngineeringEligibility() {
    HydrocarbonWaterTopologyBenchmark.Reference reference = reference(true);
    HydrocarbonWaterTopologyBenchmark.Candidate candidate = new HydrocarbonWaterTopologyBenchmark.Candidate("candidate",
        true, true, regions(Region.GAS_OIL, Region.OIL_AQUEOUS, Region.GAS_OIL_AQUEOUS), strings(OW_TO_GOW, GO_TO_GOW));

    HydrocarbonWaterTopologyBenchmark.Result result = new HydrocarbonWaterTopologyBenchmark().compare(reference,
        candidate);

    assertTrue(result.isComparisonEligible());
    assertTrue(result.isTopologyAccepted(), result.getViolations().toString());
    assertFalse(result.isNumericComparisonEligible());
    assertTrue(result.isNumericBenchmarkPending());
    assertFalse(result.isEngineeringEligible());
  }

  @Test
  void rejectsCandidateThatOnlyContainsTheCurrentOwToGowLoop() {
    HydrocarbonWaterTopologyBenchmark.Candidate candidate = new HydrocarbonWaterTopologyBenchmark.Candidate(
        "current-loop", true, true, regions(Region.OIL_AQUEOUS, Region.GAS_OIL_AQUEOUS), strings(OW_TO_GOW));

    HydrocarbonWaterTopologyBenchmark.Result result = new HydrocarbonWaterTopologyBenchmark().compare(reference(true),
        candidate);

    assertFalse(result.isTopologyAccepted());
    assertTrue(result.getMissingRegions().contains(Region.GAS_OIL));
    assertTrue(result.getMissingBoundaryCodes().contains(GO_TO_GOW));
  }

  @Test
  void rejectsUnverifiedPaperSourceContract() {
    HydrocarbonWaterTopologyBenchmark.Candidate candidate = new HydrocarbonWaterTopologyBenchmark.Candidate("candidate",
        true, true, regions(Region.GAS, Region.GAS_OIL, Region.OIL_AQUEOUS, Region.GAS_OIL_AQUEOUS),
        strings(OW_TO_GOW, GO_TO_GOW));

    HydrocarbonWaterTopologyBenchmark.Result result = new HydrocarbonWaterTopologyBenchmark().compare(reference(false),
        candidate);

    assertFalse(result.isComparisonEligible());
    assertFalse(result.isTopologyAccepted());
    assertTrue(result.getViolations().contains("REFERENCE_SOURCE_CONTRACT_UNVERIFIED"));
  }

  private static HydrocarbonWaterTopologyBenchmark.Reference reference(boolean verified) {
    return new HydrocarbonWaterTopologyBenchmark.Reference("ma-2021", verified, true,
        regions(Region.GAS_OIL, Region.OIL_AQUEOUS, Region.GAS_OIL_AQUEOUS), strings(OW_TO_GOW, GO_TO_GOW));
  }

  @SafeVarargs
  private static <T> java.util.Set<T> strings(T... values) {
    return new LinkedHashSet<T>(Arrays.asList(values));
  }

  private static java.util.Set<Region> regions(Region... values) {
    return new LinkedHashSet<Region>(Arrays.asList(values));
  }
}
