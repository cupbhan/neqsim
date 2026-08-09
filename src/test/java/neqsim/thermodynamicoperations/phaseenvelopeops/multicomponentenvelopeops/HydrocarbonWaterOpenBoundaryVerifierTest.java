package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.Arrays;
import java.util.List;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;

/** Contract regressions that prevent incomplete open/network evidence from entering thermodynamic verification. */
class HydrocarbonWaterOpenBoundaryVerifierTest extends NeqSimTest {

  @Test
  void rejectsAnOpenBranchWithoutADeclaredDomainCrossing() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    List<State> states = Arrays.asList(
        state(fluid, CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS, 300.0, 1.0),
        state(fluid, CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS, 310.0, 2.0),
        state(fluid, CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS, 320.0, 3.0));

    assertThrows(IllegalArgumentException.class,
        () -> new HydrocarbonWaterOpenBoundaryVerifier(fluid).verify(states, 0.05, states.get(2)));
  }

  @Test
  void fullNetworkRejectsTwoCopiesOfTheSameBoundaryFamily() {
    SystemInterface fluid = LindeloffMichelsenReferenceFluidTest.fluidTwo(false);
    State first = state(fluid, CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS, 300.0, 0.04);
    State second = state(fluid, CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS, 310.0, 0.1);
    State third = state(fluid, CandidatePhase.GAS, CandidatePhase.OIL, CandidatePhase.AQUEOUS, 320.0, 1.0);
    List<State> duplicateFamilyLoop = Arrays.asList(first, second, third, first);
    List<State> open = Arrays.asList(first, second, third);

    assertThrows(IllegalArgumentException.class,
        () -> new HydrocarbonWaterFullNetworkVerifier(fluid).verify(duplicateFamilyLoop, open, 0.05, third));
  }

  private static State state(SystemInterface fluid, CandidatePhase phaseZero, CandidatePhase phaseOne,
      CandidatePhase incipient, double temperatureK, double pressureBara) {
    int componentCount = fluid.getPhase(0).getNumberOfComponents();
    double[] phaseZeroComposition = new double[componentCount];
    double[] phaseOneComposition = new double[componentCount];
    double[] incipientComposition = new double[componentCount];
    phaseZeroComposition[0] = 1.0;
    phaseOneComposition[Math.min(1, componentCount - 1)] = 1.0;
    incipientComposition[componentCount - 1] = 1.0;
    return State.create(phaseZero, phaseOne, incipient, temperatureK, pressureBara, 0.5, phaseZeroComposition,
        phaseOneComposition, incipientComposition);
  }
}
