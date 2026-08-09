package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import org.junit.jupiter.api.Test;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;

class HydrocarbonWaterClosedBoundaryLoopAssemblerTest {

  @Test
  void closesTwoIdentityCompatiblePathsWithoutDuplicatingTheMergePoint() {
    State common = state(300.0, 250.0, 0.0);
    State firstMiddle = state(280.0, 150.0, 0.02);
    State secondMiddle = state(320.0, 200.0, -0.02);
    State firstMerge = state(350.0, 300.0, 0.01);
    State secondMerge = state(350.0001, 300.0002, 0.0101);

    HydrocarbonWaterClosedBoundaryLoopAssembler.Result result = new HydrocarbonWaterClosedBoundaryLoopAssembler()
        .assemble(Arrays.asList(common, firstMiddle, firstMerge), Arrays.asList(common, secondMiddle, secondMerge));

    assertTrue(result.isAccepted(), result.getFailureMessage());
    assertEquals(5, result.getStates().size());
    assertEquals(2, result.getMergePointIndex());
    assertEquals(SpecialPointType.TARGET_BRANCH_MERGE, result.getMergePointType());
    assertEquals(result.getStates().get(0).getTemperatureK(),
        result.getStates().get(result.getStates().size() - 1).getTemperatureK(), 0.0);
    assertEquals(1, result.getStates().stream().filter(point -> point.getTemperatureK() == firstMerge.getTemperatureK()
        && point.getPressureBara() == firstMerge.getPressureBara()).count());
  }

  @Test
  void rejectsPathsWhoseMergeCompositionsBelongToDifferentIdentities() {
    State common = state(300.0, 250.0, 0.0);
    State firstMerge = state(350.0, 300.0, 0.01);
    State unrelatedMerge = state(350.0, 300.0, 0.3);

    HydrocarbonWaterClosedBoundaryLoopAssembler.Result result = new HydrocarbonWaterClosedBoundaryLoopAssembler()
        .assemble(Arrays.asList(common, state(320.0, 270.0, 0.02), firstMerge),
            Arrays.asList(common, state(310.0, 260.0, -0.02), unrelatedMerge));

    assertFalse(result.isAccepted());
    assertEquals("the two paths do not share the same target-branch merge identity", result.getFailureMessage());
  }

  @Test
  void ordinaryClosureUsesBookkeepingSeamWithoutInventingTargetBranchMerge() {
    State common = state(300.0, 250.0, 0.0);
    State firstSeam = state(350.0, 300.0, 0.01);
    State secondSeam = state(350.0001, 300.0002, 0.0101);

    HydrocarbonWaterClosedBoundaryLoopAssembler.Result result = new HydrocarbonWaterClosedBoundaryLoopAssembler()
        .assembleOrdinaryClosure(Arrays.asList(common, state(320.0, 270.0, 0.02), firstSeam),
            Arrays.asList(common, state(310.0, 260.0, -0.02), secondSeam));

    assertTrue(result.isAccepted(), result.getFailureMessage());
    assertEquals(SpecialPointType.CLOSED_LOOP_SEAM, result.getMergePointType());
  }

  private static State state(double temperatureK, double pressureBara, double shift) {
    return State.create(CandidatePhase.OIL, CandidatePhase.AQUEOUS, CandidatePhase.GAS, temperatureK, pressureBara, 0.4,
        new double[] { 0.8 - shift, 0.1 + shift, 0.1 }, new double[] { 0.01, 0.01, 0.98 },
        new double[] { 0.7 - shift, 0.2 + shift, 0.1 });
  }
}
