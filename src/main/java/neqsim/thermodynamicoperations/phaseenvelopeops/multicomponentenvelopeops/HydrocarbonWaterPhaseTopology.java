package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Locale;
import java.util.Set;

/**
 * Physical phase-region topology for non-reactive gas-oil-aqueous PT envelopes.
 *
 * <p>
 * The topology contract deliberately separates thermodynamic phase classification from numerical continuation. A valid
 * phase boundary connects two stable regions that differ by exactly one phase. This prevents diagnostic TP-grid states
 * from being joined across a three-phase point and presented as one continuous engineering boundary.
 * </p>
 */
public final class HydrocarbonWaterPhaseTopology {
  private HydrocarbonWaterPhaseTopology() {
  }

  /** Fluid phase identities used by the non-reactive hydrocarbon-water envelope. */
  public enum Phase {
    GAS, OIL, AQUEOUS
  }

  /** All physically addressable stable phase regions for a gas-oil-water system. */
  public enum Region {
    GAS("G", Phase.GAS), OIL("O", Phase.OIL), AQUEOUS("W", Phase.AQUEOUS), GAS_OIL("GO", Phase.GAS, Phase.OIL),
    GAS_AQUEOUS("GW", Phase.GAS, Phase.AQUEOUS), OIL_AQUEOUS("OW", Phase.OIL, Phase.AQUEOUS),
    GAS_OIL_AQUEOUS("GOW", Phase.GAS, Phase.OIL, Phase.AQUEOUS);

    private final String code;
    private final Set<Phase> phases;

    Region(String code, Phase... phases) {
      this.code = code;
      EnumSet<Phase> phaseSet = EnumSet.noneOf(Phase.class);
      Collections.addAll(phaseSet, phases);
      this.phases = Collections.unmodifiableSet(phaseSet);
    }

    /** @return compact stable-region code used in result serialization */
    public String getCode() {
      return code;
    }

    /** @return immutable phase set */
    public Set<Phase> getPhases() {
      return phases;
    }

    /**
     * Tests whether two regions may share one physical saturation boundary.
     *
     * @param other candidate adjacent region
     * @return true when one region contains every phase of the other plus exactly one phase
     */
    public boolean isAdjacentTo(Region other) {
      if (other == null || other == this || Math.abs(phases.size() - other.phases.size()) != 1) {
        return false;
      }
      Set<Phase> smaller = phases.size() < other.phases.size() ? phases : other.phases;
      Set<Phase> larger = phases.size() < other.phases.size() ? other.phases : phases;
      return larger.containsAll(smaller);
    }

    /**
     * Maps converged NeqSim phase type names to one stable topology region.
     *
     * @param phaseTypeNames NeqSim phase names such as gas, oil, and aqueous
     * @return classified region
     */
    public static Region fromPhaseTypeNames(Collection<String> phaseTypeNames) {
      if (phaseTypeNames == null || phaseTypeNames.isEmpty()) {
        throw new IllegalArgumentException("at least one phase type is required");
      }
      EnumSet<Phase> classified = EnumSet.noneOf(Phase.class);
      for (String phaseTypeName : phaseTypeNames) {
        if (phaseTypeName == null) {
          continue;
        }
        String normalized = phaseTypeName.trim().toLowerCase(Locale.ROOT);
        if (normalized.equals("gas")) {
          classified.add(Phase.GAS);
        } else if (normalized.equals("oil") || normalized.equals("liquid")) {
          classified.add(Phase.OIL);
        } else if (normalized.equals("aqueous") || normalized.equals("water")) {
          classified.add(Phase.AQUEOUS);
        }
      }
      for (Region region : values()) {
        if (region.phases.equals(classified)) {
          return region;
        }
      }
      throw new IllegalArgumentException("unsupported phase set " + classified);
    }
  }

  /** Special points that may terminate, join, or characterize physical boundary segments. */
  public enum SpecialPointType {
    THREE_PHASE_POINT, CRITICAL_POINT, CRITICAL_END_POINT, PHASE_COALESCENCE_CANDIDATE, STABILITY_STATIONARY_FOLD,
    RETAINED_PHASE_SPINODAL, ALTERNATE_STATIONARY_BRANCH_ONSET, TARGET_BRANCH_MERGE, CLOSED_LOOP_SEAM, CRICONDENBAR,
    CRICONDENTHERM, DOMAIN_EXIT
  }

  /** Immutable declaration of one oriented physical boundary. */
  public static final class BoundaryDefinition {
    private final Region lowerPhaseRegion;
    private final Region higherPhaseRegion;
    private final Phase incipientPhase;

    /**
     * Creates a boundary from the region with fewer phases to the region with one additional phase.
     *
     * @param lowerPhaseRegion stable region on the lower-phase-count side
     * @param higherPhaseRegion stable region on the higher-phase-count side
     */
    public BoundaryDefinition(Region lowerPhaseRegion, Region higherPhaseRegion) {
      if (lowerPhaseRegion == null || higherPhaseRegion == null || !lowerPhaseRegion.isAdjacentTo(higherPhaseRegion)
          || higherPhaseRegion.getPhases().size() != lowerPhaseRegion.getPhases().size() + 1
          || !higherPhaseRegion.getPhases().containsAll(lowerPhaseRegion.getPhases())) {
        throw new IllegalArgumentException("a boundary must add exactly one phase to an adjacent stable region");
      }
      this.lowerPhaseRegion = lowerPhaseRegion;
      this.higherPhaseRegion = higherPhaseRegion;
      EnumSet<Phase> difference = EnumSet.copyOf(higherPhaseRegion.getPhases());
      difference.removeAll(lowerPhaseRegion.getPhases());
      this.incipientPhase = difference.iterator().next();
    }

    /** @return lower-phase-count side */
    public Region getLowerPhaseRegion() {
      return lowerPhaseRegion;
    }

    /** @return higher-phase-count side */
    public Region getHigherPhaseRegion() {
      return higherPhaseRegion;
    }

    /** @return phase that appears when crossing toward the higher-phase-count side */
    public Phase getIncipientPhase() {
      return incipientPhase;
    }

    /** @return stable serialization label such as {@code G->GO:OIL} */
    public String getCode() {
      return lowerPhaseRegion.getCode() + "->" + higherPhaseRegion.getCode() + ":" + incipientPhase.name();
    }
  }
}
