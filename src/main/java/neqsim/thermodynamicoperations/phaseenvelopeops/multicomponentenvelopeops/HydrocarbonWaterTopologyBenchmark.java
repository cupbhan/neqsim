package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;

/**
 * Topology-only paper benchmark for non-reactive gas-oil-aqueous PT envelopes.
 *
 * <p>
 * This contract is intentionally separate from {@link HydrocarbonWaterBoundaryRegression}. A paper figure produced with
 * another thermodynamic model may prove that required phase-region adjacencies exist, but it cannot establish
 * same-model coordinate accuracy or engineering release eligibility.
 * </p>
 */
public final class HydrocarbonWaterTopologyBenchmark {

  /** Compare an audited paper topology contract with independently classified candidate evidence. */
  public Result compare(Reference reference, Candidate candidate) {
    if (reference == null || candidate == null) {
      throw new IllegalArgumentException("topology benchmark requires reference and candidate");
    }
    List<String> violations = new ArrayList<String>();
    if (!reference.isSourceContractVerified()) {
      violations.add("REFERENCE_SOURCE_CONTRACT_UNVERIFIED");
    }
    if (!reference.isTopologyOnly()) {
      violations.add("REFERENCE_SCOPE_IS_NOT_TOPOLOGY_ONLY");
    }
    if (!candidate.isTopologyEvidenceResolved()) {
      violations.add("CANDIDATE_TOPOLOGY_EVIDENCE_UNRESOLVED");
    }

    Set<Region> missingRegions = new LinkedHashSet<Region>(reference.getRequiredRegions());
    missingRegions.removeAll(candidate.getObservedRegions());
    if (!missingRegions.isEmpty()) {
      violations.add("REQUIRED_PHASE_REGIONS_MISSING");
    }
    Set<String> missingBoundaries = new LinkedHashSet<String>(reference.getRequiredBoundaryCodes());
    missingBoundaries.removeAll(candidate.getObservedBoundaryCodes());
    if (!missingBoundaries.isEmpty()) {
      violations.add("REQUIRED_TWO_TO_THREE_PHASE_BOUNDARIES_MISSING");
    }

    boolean comparisonEligible = reference.isSourceContractVerified() && reference.isTopologyOnly()
        && candidate.isTopologyEvidenceResolved();
    boolean topologyAccepted = comparisonEligible && missingRegions.isEmpty() && missingBoundaries.isEmpty();
    return new Result(comparisonEligible, topologyAccepted, missingRegions, missingBoundaries, violations);
  }

  /** Audited paper topology declaration. */
  public static final class Reference {
    private final String identifier;
    private final boolean sourceContractVerified;
    private final boolean topologyOnly;
    private final Set<Region> requiredRegions;
    private final Set<String> requiredBoundaryCodes;

    public Reference(String identifier, boolean sourceContractVerified, boolean topologyOnly,
        Set<Region> requiredRegions, Set<String> requiredBoundaryCodes) {
      requireText(identifier, "reference identifier");
      if (requiredRegions == null || requiredRegions.isEmpty() || requiredBoundaryCodes == null
          || requiredBoundaryCodes.isEmpty()) {
        throw new IllegalArgumentException("paper topology reference requires regions and boundary families");
      }
      this.identifier = identifier;
      this.sourceContractVerified = sourceContractVerified;
      this.topologyOnly = topologyOnly;
      this.requiredRegions = Collections.unmodifiableSet(new LinkedHashSet<Region>(requiredRegions));
      this.requiredBoundaryCodes = Collections.unmodifiableSet(new LinkedHashSet<String>(requiredBoundaryCodes));
    }

    public String getIdentifier() {
      return identifier;
    }

    public boolean isSourceContractVerified() {
      return sourceContractVerified;
    }

    public boolean isTopologyOnly() {
      return topologyOnly;
    }

    public Set<Region> getRequiredRegions() {
      return requiredRegions;
    }

    public Set<String> getRequiredBoundaryCodes() {
      return requiredBoundaryCodes;
    }
  }

  /** Independently classified NeqSim topology evidence. */
  public static final class Candidate {
    private final String identifier;
    private final boolean topologyEvidenceResolved;
    private final boolean internalQualityEligible;
    private final Set<Region> observedRegions;
    private final Set<String> observedBoundaryCodes;

    public Candidate(String identifier, boolean topologyEvidenceResolved, boolean internalQualityEligible,
        Set<Region> observedRegions, Set<String> observedBoundaryCodes) {
      requireText(identifier, "candidate identifier");
      this.identifier = identifier;
      this.topologyEvidenceResolved = topologyEvidenceResolved;
      this.internalQualityEligible = internalQualityEligible;
      this.observedRegions = observedRegions == null ? Collections.<Region>emptySet()
          : Collections.unmodifiableSet(new LinkedHashSet<Region>(observedRegions));
      this.observedBoundaryCodes = observedBoundaryCodes == null ? Collections.<String>emptySet()
          : Collections.unmodifiableSet(new LinkedHashSet<String>(observedBoundaryCodes));
    }

    public String getIdentifier() {
      return identifier;
    }

    public boolean isTopologyEvidenceResolved() {
      return topologyEvidenceResolved;
    }

    public boolean isInternalQualityEligible() {
      return internalQualityEligible;
    }

    public Set<Region> getObservedRegions() {
      return observedRegions;
    }

    public Set<String> getObservedBoundaryCodes() {
      return observedBoundaryCodes;
    }
  }

  /** Topology-only result; numerical and engineering eligibility remain false by design. */
  public static final class Result {
    private final boolean comparisonEligible;
    private final boolean topologyAccepted;
    private final Set<Region> missingRegions;
    private final Set<String> missingBoundaryCodes;
    private final List<String> violations;

    private Result(boolean comparisonEligible, boolean topologyAccepted, Set<Region> missingRegions,
        Set<String> missingBoundaryCodes, List<String> violations) {
      this.comparisonEligible = comparisonEligible;
      this.topologyAccepted = topologyAccepted;
      this.missingRegions = Collections.unmodifiableSet(new LinkedHashSet<Region>(missingRegions));
      this.missingBoundaryCodes = Collections.unmodifiableSet(new LinkedHashSet<String>(missingBoundaryCodes));
      this.violations = Collections.unmodifiableList(new ArrayList<String>(violations));
    }

    public boolean isComparisonEligible() {
      return comparisonEligible;
    }

    public boolean isTopologyAccepted() {
      return topologyAccepted;
    }

    public boolean isNumericComparisonEligible() {
      return false;
    }

    public boolean isNumericBenchmarkPending() {
      return true;
    }

    public boolean isEngineeringEligible() {
      return false;
    }

    public Set<Region> getMissingRegions() {
      return missingRegions;
    }

    public Set<String> getMissingBoundaryCodes() {
      return missingBoundaryCodes;
    }

    public List<String> getViolations() {
      return violations;
    }
  }

  private static void requireText(String value, String label) {
    if (value == null || value.trim().isEmpty()) {
      throw new IllegalArgumentException(label + " is required");
    }
  }
}
