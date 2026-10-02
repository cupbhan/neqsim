package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.Collections;
import org.junit.jupiter.api.Test;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterContinuationChainAssembler.Result;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterContinuationChainAssembler.Segment;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;

class HydrocarbonWaterContinuationChainAssemblerTest {

  @Test
  void joinsTwoCheckpointSegmentsThroughTheirTwoStateIdentityOverlap() {
    State first = state(250.0, 100.0, 0.0);
    State overlapOne = state(260.0, 110.0, 0.01);
    State overlapTwo = state(270.0, 120.0, 0.02);
    State last = state(280.0, 130.0, 0.03);

    Result result = new HydrocarbonWaterContinuationChainAssembler()
        .assemble(Arrays.asList(new Segment("low", Arrays.asList(first, overlapOne, overlapTwo)),
            new Segment("high", Arrays.asList(overlapOne, overlapTwo, last))));

    assertTrue(result.isSingleConnectedChain(), result.getViolations().toString());
    assertEquals(1, result.getChains().size());
    assertEquals(4, result.getChains().get(0).getStates().size());
    assertEquals(1, result.getJoins().size());
    assertEquals(2, result.getJoins().get(0).getOverlapStateCount());
    assertEquals(0.0, result.getJoins().get(0).getMaximumPtDistance(), 0.0);
    assertEquals(0.0, result.getJoins().get(0).getMaximumCompositionDistance(), 0.0);
  }

  @Test
  void reversesOneSegmentWhenTheIdentityOverlapIsAtTwoEnds() {
    State first = state(250.0, 100.0, 0.0);
    State overlapOne = state(260.0, 110.0, 0.01);
    State overlapTwo = state(270.0, 120.0, 0.02);
    State last = state(280.0, 130.0, 0.03);

    Result result = new HydrocarbonWaterContinuationChainAssembler()
        .assemble(Arrays.asList(new Segment("first", Arrays.asList(first, overlapOne, overlapTwo)),
            new Segment("reversed", Arrays.asList(last, overlapTwo, overlapOne))));

    assertTrue(result.isSingleConnectedChain(), result.getViolations().toString());
    assertEquals(4, result.getChains().get(0).getStates().size());
    assertTrue(result.getJoins().get(0).isReversedSecond());
    assertEquals(2, result.getJoins().get(0).getOverlapStateCount());
  }

  @Test
  void preservesDisconnectedSegmentsInsteadOfInventingAGapConnection() {
    Result result = new HydrocarbonWaterContinuationChainAssembler()
        .assemble(Arrays.asList(new Segment("first", Collections.singletonList(state(250.0, 100.0, 0.0))),
            new Segment("second", Collections.singletonList(state(350.0, 300.0, 0.5)))));

    assertFalse(result.isSingleConnectedChain());
    assertEquals(2, result.getChains().size());
    assertEquals(Collections.singletonList("DISCONNECTED_CHECKPOINT_CHAINS:2"), result.getViolations());
  }

  private static State state(double temperatureK, double pressureBara, double shift) {
    return State.create(CandidatePhase.OIL, CandidatePhase.AQUEOUS, CandidatePhase.GAS, temperatureK, pressureBara, 0.4,
        new double[] {0.8 - shift, 0.1 + shift, 0.1}, new double[] {0.01, 0.01, 0.98},
        new double[] {0.7 - shift, 0.2 + shift, 0.1});
  }
}
