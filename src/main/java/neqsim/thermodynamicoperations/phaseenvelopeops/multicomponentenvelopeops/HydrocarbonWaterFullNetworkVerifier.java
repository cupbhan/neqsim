package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.Branch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EnvelopeReport;

/** Verifies the first-scope OW-to-GOW closed loop and GO-to-GOW open branch as one physical network. */
public final class HydrocarbonWaterFullNetworkVerifier {
  private final SystemInterface template;

  /** Creates a verifier from the exact model and fluid definition under test. */
  public HydrocarbonWaterFullNetworkVerifier(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Independently re-solves both branches and applies one combined topology/continuity quality gate. */
  public Result verify(List<State> closedLoopStates, List<State> openBranchStates, double minimumPressureBara,
      State firstRejectedHighState) {
    BoundaryFamily closedFamily = boundaryFamily(
        closedLoopStates == null || closedLoopStates.isEmpty() ? null : closedLoopStates.get(0));
    BoundaryFamily openFamily = boundaryFamily(
        openBranchStates == null || openBranchStates.isEmpty() ? null : openBranchStates.get(0));
    EnumSet<BoundaryFamily> observed = EnumSet.of(closedFamily, openFamily);
    EnumSet<BoundaryFamily> required = EnumSet.of(BoundaryFamily.OW_TO_GOW, BoundaryFamily.GO_TO_GOW);
    if (!observed.equals(required)) {
      throw new IllegalArgumentException("full first-scope network requires exactly OW_TO_GOW and GO_TO_GOW");
    }

    HydrocarbonWaterClosedBoundaryLoopVerifier.Result closed = new HydrocarbonWaterClosedBoundaryLoopVerifier(template)
        .verify(closedLoopStates);
    HydrocarbonWaterOpenBoundaryVerifier.Result open = new HydrocarbonWaterOpenBoundaryVerifier(template)
        .verify(openBranchStates, minimumPressureBara, firstRejectedHighState);
    Branch closedBranch = new Branch("verified-network-ow-to-gow-closed-loop", closedFamily.getDefinition(),
        closedFamily.getRetainedPhaseZero(), closedFamily.getRetainedPhaseOne(), closedFamily.getIncipientPhase(),
        closed.getEvidence(), null, null);
    Branch openBranch = new Branch("verified-network-go-to-gow-domain-to-spinodal", openFamily.getDefinition(),
        openFamily.getRetainedPhaseZero(), openFamily.getRetainedPhaseOne(), openFamily.getIncipientPhase(),
        open.getEvidence(), SpecialPointType.DOMAIN_EXIT, SpecialPointType.RETAINED_PHASE_SPINODAL, null,
        open.getHighTermination());
    EnvelopeReport combined = new TwoToThreePhaseBoundaryQualityGate()
        .validate(Arrays.asList(closedBranch, openBranch));
    boolean internalQualityEligible = closed.isInternalQualityEligible() && open.isInternalQualityEligible()
        && combined.isAccepted();
    return new Result(closed, open, combined, internalQualityEligible);
  }

  private static BoundaryFamily boundaryFamily(State state) {
    if (state != null) {
      for (BoundaryFamily family : BoundaryFamily.values()) {
        if (family.getRetainedPhaseZero() == state.getRetainedPhaseZero()
            && family.getRetainedPhaseOne() == state.getRetainedPhaseOne()
            && family.getIncipientPhase() == state.getIncipientPhase()) {
          return family;
        }
      }
    }
    throw new IllegalArgumentException("a supported hydrocarbon-water boundary state is required");
  }

  /** Branch-level independent results plus the combined network-quality verdict. */
  public static final class Result {
    private final HydrocarbonWaterClosedBoundaryLoopVerifier.Result closedLoop;
    private final HydrocarbonWaterOpenBoundaryVerifier.Result openBranch;
    private final EnvelopeReport combinedQualityReport;
    private final boolean internalQualityEligible;

    private Result(HydrocarbonWaterClosedBoundaryLoopVerifier.Result closedLoop,
        HydrocarbonWaterOpenBoundaryVerifier.Result openBranch, EnvelopeReport combinedQualityReport,
        boolean internalQualityEligible) {
      this.closedLoop = closedLoop;
      this.openBranch = openBranch;
      this.combinedQualityReport = combinedQualityReport;
      this.internalQualityEligible = internalQualityEligible;
    }

    public HydrocarbonWaterClosedBoundaryLoopVerifier.Result getClosedLoop() {
      return closedLoop;
    }

    public HydrocarbonWaterOpenBoundaryVerifier.Result getOpenBranch() {
      return openBranch;
    }

    public EnvelopeReport getCombinedQualityReport() {
      return combinedQualityReport;
    }

    public boolean isInternalQualityEligible() {
      return internalQualityEligible;
    }
  }
}
