package neqsim.mcp.runners;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.Branch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryEndpointClassifier;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryGlobalStabilityGate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryBranchAssembler;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterIsolatedBoundaryRootClassifier;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.EndpointAttachment;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.EndpointEvidence;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.End;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.NetworkBranch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryNetworkAssembler.TargetBranchMergeTrimResult;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterCriticalEndpointSolver;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryTerminationClassifier;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterClosedBoundaryLoopAssembler;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterClosedBoundaryLoopVerifier;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterOpenBoundaryVerifier;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterFullNetworkVerifier;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterRegularPressureBoundaryTracer;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterSecondaryStationaryBranchTracker;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterStableRegionBoundaryCorrector;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterStableRegionTransitionScanner;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.SeedEvaluation;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStationarityJacobianAnalyzer;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EvidencePoint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryPointSolver;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhasePseudoArcLengthTracer;

/**
 * Produces an opt-in inventory and physical network audit of hydrocarbon-water two-to-three-phase boundaries.
 *
 * <p>
 * Corrected anchors alone are never promoted. The runner exposes stable-region brackets, strict roots, disconnected
 * branches, withheld endpoint roots, and independently refined endpoints, then performs bidirectional continuation and
 * network assembly. Only a complete network that passes both endpoint topology and per-point quality gates is marked
 * engineering eligible.
 * </p>
 */
public final class HydrocarbonWaterBoundaryAuditRunner {
  private static final double MAXIMUM_MATCHED_COMPOSITION_JUMP = 0.35;
  private static final double MAXIMUM_BRANCH_SWITCH_COMPOSITION_JUMP = 0.15;
  private static final double MAXIMUM_BRANCH_SWITCH_PT_DISTANCE = 0.02;
  private static final double MINIMUM_BRANCH_SWITCH_TARGET_SEED_PT_DISTANCE = 1.0e-5;
  private static final double MINIMUM_BRANCH_SWITCH_FORWARD_ADVANCE_RATIO = 0.02;
  private static final int DEFAULT_CONTINUATION_POINTS_PER_DIRECTION = 80;
  // A stable TP phase-count switch can be displaced slightly from the specified-two-phase zero-TPD
  // root because its retained phase split and compositions come from one side of the bracket. One
  // attempt examines only the already-refined (often sub-kelvin) bracket and defeats the corrector's
  // controlled interval-expansion algorithm. Six attempts cover the corrector's complete bounded
  // sequence up to its 25 K half-width without permitting an unbounded search.
  private static final int STABLE_BRACKET_CORRECTION_INTERVAL_ATTEMPTS = 6;

  private HydrocarbonWaterBoundaryAuditRunner() {
  }

  /**
   * Runs stable-region discovery and retained-phase termination refinement.
   *
   * @param template fully configured non-reactive EOS system
   * @param temperaturesK stable-region temperature probes in kelvin
   * @param pressuresBara stable-region pressure probes in bara
   * @return JSON audit, marked engineering eligible only when the complete physical network passes
   */
  public static JsonObject run(SystemInterface template, double[] temperaturesK, double[] pressuresBara) {
    return run(template, temperaturesK, pressuresBara, DEFAULT_CONTINUATION_POINTS_PER_DIRECTION);
  }

  /**
   * Runs the audit with an explicit pseudo-arclength extension budget per discovered branch end.
   *
   * @param template fully configured non-reactive EOS system
   * @param temperaturesK stable-region temperature probes in kelvin
   * @param pressuresBara stable-region pressure probes in bara
   * @param continuationPointsPerDirection maximum accepted corrections requested from each anchor end
   * @return JSON audit with strict anchors and the best available physical branch network
   */
  public static JsonObject run(SystemInterface template, double[] temperaturesK, double[] pressuresBara,
      int continuationPointsPerDirection) {
    if (template == null || temperaturesK == null || temperaturesK.length < 2 || pressuresBara == null
        || pressuresBara.length < 2 || continuationPointsPerDirection < 1 || continuationPointsPerDirection > 400) {
      throw new IllegalArgumentException("a thermodynamic template, two temperatures, and two pressures are required");
    }
    long started = System.currentTimeMillis();
    long discoveryStarted = System.currentTimeMillis();
    HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        template).setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).setStableScanControls(32, 0.25)
        .setMaximumCorrectionIntervalAttempts(STABLE_BRACKET_CORRECTION_INTERVAL_ATTEMPTS)
        .setClusteringControls(MAXIMUM_MATCHED_COMPOSITION_JUMP, 1.75, 0.5)
        .discoverFromStableRegionTransitions(temperaturesK, pressuresBara);
    long discoveryTimeMs = System.currentTimeMillis() - discoveryStarted;
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery = stable.getDiscovery();

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "stable-region-scan-strict-correction-pseudo-arclength-network-assembly");
    output.addProperty("boundaryTopology", "hydrocarbon-water-two-to-three-phase");
    output.addProperty("continuationPointsPerDirection", continuationPointsPerDirection);
    output.addProperty("stableFlashEvaluationCount", stable.getStableRegionScan().getFlashEvaluations());
    output.addProperty("stableTransitionBracketCount", stable.getStableRegionScan().getBrackets().size());
    output.addProperty("stableScanFailureCount", stable.getStableRegionScan().getFailures().size());
    output.addProperty("initialCorrectedAnchorCount", discovery.getCorrectedAnchorCount());
    output.addProperty("initialBranchCount", discovery.getBranches().size());
    output.addProperty("initialIsolatedBranchCount", discovery.getIsolatedBranchCount());
    output.addProperty("withheldEndpointCandidateCount", discovery.getEndpointCandidateCount());
    output.addProperty("discoveryTimeMs", discoveryTimeMs);
    output.addProperty("stableScanTimeMs", stable.getStableScanTimeMs());
    output.addProperty("strictCorrectionTimeMs", stable.getCorrectionTimeMs());
    output.add("strictCorrectionDiagnostics", correctionDiagnostics(stable));
    long isolatedTopologyStarted = System.currentTimeMillis();
    IsolatedTopologyAudit isolatedTopology = isolatedRootDiagnostics(template, discovery);
    output.add("isolatedRootDiagnostics", isolatedTopology.json);
    output.addProperty("isolatedRootDiagnosticTimeMs", System.currentTimeMillis() - isolatedTopologyStarted);
    HydrocarbonWaterBoundaryAnchorDiscoverer.Result effectiveDiscovery = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        template).setClusteringControls(MAXIMUM_MATCHED_COMPOSITION_JUMP, 1.75, 0.5)
        .augmentWithIndependentAnchors(discovery, isolatedTopology.independentAnchors);
    output.addProperty("locallyRecoveredAnchorCount",
        effectiveDiscovery.getCorrectedAnchorCount() - discovery.getCorrectedAnchorCount());
    output.addProperty("correctedAnchorCount", effectiveDiscovery.getCorrectedAnchorCount());
    output.addProperty("branchCount", effectiveDiscovery.getBranches().size());
    output.addProperty("isolatedBranchCount", effectiveDiscovery.getIsolatedBranchCount());
    output.add("branches", branches(effectiveDiscovery));

    long networkStarted = System.currentTimeMillis();
    double[] sortedTemperatures = sortedUnique(temperaturesK);
    double[] sortedPressures = sortedUnique(pressuresBara);
    JsonObject network = network(template, effectiveDiscovery, continuationPointsPerDirection, sortedTemperatures[0],
        sortedTemperatures[sortedTemperatures.length - 1], sortedPressures[0],
        sortedPressures[sortedPressures.length - 1]);
    long networkTimeMs = System.currentTimeMillis() - networkStarted;
    JsonObject coverage = discoveryCoverage(temperaturesK, pressuresBara, stable, effectiveDiscovery);
    boolean branchNetworkEligible = network.get("engineeringEligible").getAsBoolean();
    boolean discoveryCoveragePassed = coverage.get("passed").getAsBoolean();
    network.addProperty("branchNetworkEligible", branchNetworkEligible);
    network.addProperty("discoveryCoveragePassed", discoveryCoveragePassed);
    network.addProperty("engineeringEligible", branchNetworkEligible && discoveryCoveragePassed);
    output.add("network", network);
    output.addProperty("networkTimeMs", networkTimeMs);
    output.add("discoveryCoverage", coverage);
    output.addProperty("engineeringEligible", branchNetworkEligible && discoveryCoveragePassed);

    long terminationStarted = System.currentTimeMillis();
    List<BoundaryEvidence> evidence = evidence(stable);
    JsonArray terminations = terminations(template, evidence);
    long terminationTimeMs = System.currentTimeMillis() - terminationStarted;
    output.add("terminations", terminations);
    int physicalTerminationCount = 0;
    for (int index = 0; index < terminations.size(); index++) {
      if (terminations.get(index).getAsJsonObject().get("physicalEndpoint").getAsBoolean()) {
        physicalTerminationCount++;
      }
    }
    int networkPhysicalEndpointCount = 0;
    JsonArray networkEndpoints = network.getAsJsonArray("endpoints");
    for (int index = 0; index < networkEndpoints.size(); index++) {
      SpecialPointType type = SpecialPointType
          .valueOf(networkEndpoints.get(index).getAsJsonObject().get("type").getAsString());
      if (type != SpecialPointType.DOMAIN_EXIT) {
        networkPhysicalEndpointCount++;
      }
    }
    physicalTerminationCount = Math.max(physicalTerminationCount, networkPhysicalEndpointCount);
    output.addProperty("physicalTerminationCount", physicalTerminationCount);
    output.addProperty("terminationTimeMs", terminationTimeMs);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    output.addProperty("summary",
        effectiveDiscovery.getCorrectedAnchorCount() + " strict globally stable anchors in "
            + effectiveDiscovery.getBranches().size() + " preserved branches; " + physicalTerminationCount
            + " independently refined physical terminations; " + network.get("attachedEndpointCount").getAsInt()
            + " endpoints attached to the global network. "
            + (branchNetworkEligible && discoveryCoveragePassed
                ? "Every preserved branch and physical endpoint passed the global topology gate."
                : branchNetworkEligible
                    ? "The discovered branch network is closed, but the probe grid is not dense enough to exclude "
                        + "missed isolated branches."
                    : "Network remains diagnostic until every branch and both endpoints pass the global topology "
                        + "gate."));
    return output;
  }

  /**
   * Runs only stable multiphase TP flashes and adjacent-region bracket refinement.
   *
   * <p>
   * This diagnostic intentionally omits zero-TPD correction, continuation, endpoint refinement, and network assembly.
   * It answers whether the selected model observes the stable regions needed to justify a full boundary audit, but it
   * can never produce an engineering curve.
   * </p>
   *
   * @param template fully configured non-reactive EOS system
   * @param temperaturesK temperature probes in kelvin
   * @param pressuresBara pressure probes in bara
   * @return stable-region inventory and adjacent 2P/GOW brackets
   */
  public static JsonObject runStableRegionTopologyScan(SystemInterface template, double[] temperaturesK,
      double[] pressuresBara) {
    if (template == null || temperaturesK == null || temperaturesK.length < 2 || pressuresBara == null
        || pressuresBara.length < 1) {
      throw new IllegalArgumentException("a thermodynamic template, two temperatures, and one pressure are required");
    }
    long started = System.currentTimeMillis();
    HydrocarbonWaterStableRegionTransitionScanner.Result scan = new HydrocarbonWaterStableRegionTransitionScanner(
        template).setNumericalControls(32, 0.25, 1.0e-9).scan(temperaturesK, pressuresBara);

    EnumMap<Region, Integer> regionCounts = new EnumMap<Region, Integer>(Region.class);
    JsonArray states = new JsonArray();
    for (HydrocarbonWaterStableRegionTransitionScanner.StableState state : scan.getStates()) {
      regionCounts.put(state.getRegion(), regionCounts.getOrDefault(state.getRegion(), 0) + 1);
      JsonObject row = new JsonObject();
      row.addProperty("temperatureK", state.getTemperatureK());
      row.addProperty("temperatureC", state.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", state.getPressureBara());
      row.addProperty("region", state.getRegion().getCode());
      row.addProperty("rawPhaseCount", state.getRawPhaseCount());
      row.addProperty("gasFraction", state.getPhaseFraction(CandidatePhase.GAS));
      row.addProperty("oilFraction", state.getPhaseFraction(CandidatePhase.OIL));
      row.addProperty("aqueousFraction", state.getPhaseFraction(CandidatePhase.AQUEOUS));
      states.add(row);
    }

    EnumMap<BoundaryFamily, Integer> bracketCounts = new EnumMap<BoundaryFamily, Integer>(BoundaryFamily.class);
    JsonArray brackets = new JsonArray();
    for (HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket : scan.getBrackets()) {
      bracketCounts.put(bracket.getFamily(), bracketCounts.getOrDefault(bracket.getFamily(), 0) + 1);
      JsonObject row = new JsonObject();
      row.addProperty("family", bracket.getFamily().name());
      row.addProperty("pressureBara", bracket.getPressureBara());
      row.addProperty("lowerTemperatureK", bracket.getLowerTemperatureState().getTemperatureK());
      row.addProperty("upperTemperatureK", bracket.getUpperTemperatureState().getTemperatureK());
      row.addProperty("temperatureWidthK", bracket.getTemperatureWidthK());
      row.addProperty("twoPhaseRegion", bracket.getTwoPhaseState().getRegion().getCode());
      row.addProperty("threePhaseRegion", bracket.getThreePhaseState().getRegion().getCode());
      row.addProperty("incipientPhaseFractionOnThreePhaseSide", bracket.getIncipientPhaseFractionOnThreePhaseSide());
      row.addProperty("retainedPhaseZeroSplitOnTwoPhaseSide", bracket.getRetainedPhaseZeroSplitOnTwoPhaseSide());
      row.addProperty("retainedPhaseZeroSplitOnThreePhaseSide", bracket.getRetainedPhaseZeroSplitOnThreePhaseSide());
      row.addProperty("retainedPhaseSplitJump", bracket.getRetainedPhaseSplitJump());
      row.addProperty("retainedPhaseZeroCompositionJump", bracket.getRetainedPhaseZeroCompositionJump());
      row.addProperty("retainedPhaseOneCompositionJump", bracket.getRetainedPhaseOneCompositionJump());
      row.addProperty("maximumRetainedPhaseCompositionJump", bracket.getMaximumRetainedPhaseCompositionJump());
      row.addProperty("bisectionCount", bracket.getBisections());
      brackets.add(row);
    }

    JsonArray failures = new JsonArray();
    for (HydrocarbonWaterStableRegionTransitionScanner.Failure failure : scan.getFailures()) {
      JsonObject row = new JsonObject();
      row.addProperty("temperatureK", failure.getTemperatureK());
      row.addProperty("temperatureC", failure.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", failure.getPressureBara());
      row.addProperty("reason", failure.getReason());
      failures.add(row);
    }

    JsonObject observedRegions = new JsonObject();
    for (Map.Entry<Region, Integer> entry : regionCounts.entrySet()) {
      observedRegions.addProperty(entry.getKey().getCode(), entry.getValue());
    }
    JsonObject observedBrackets = new JsonObject();
    for (Map.Entry<BoundaryFamily, Integer> entry : bracketCounts.entrySet()) {
      observedBrackets.addProperty(entry.getKey().name(), entry.getValue());
    }
    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "stable-multiphase-tp-region-scan-only");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("strictBoundaryCalculated", false);
    output.addProperty("requestedGridStateCount", temperaturesK.length * pressuresBara.length);
    output.addProperty("successfulGridStateCount", scan.getStates().size());
    output.addProperty("flashEvaluationCount", scan.getFlashEvaluations());
    output.addProperty("transitionBracketCount", scan.getBrackets().size());
    output.addProperty("failureCount", scan.getFailures().size());
    output.add("observedStableRegionCounts", observedRegions);
    output.add("transitionBracketCounts", observedBrackets);
    output.add("states", states);
    output.add("brackets", brackets);
    output.add("failures", failures);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /**
   * Runs a bounded fixed-pressure continuation diagnostic for one already observed regular boundary family.
   *
   * <p>
   * This diagnostic deliberately avoids the global network audit. It first rediscovers strict, globally stable anchors
   * at the supplied seed pressures and then asks the pressure-parameter tracer for a small number of new points. A
   * successful result proves only local continuation capability; it is not an engineering-quality envelope.
   * </p>
   *
   * @param template fully configured non-reactive EOS system
   * @param temperaturesK local stable-topology scan temperatures in kelvin
   * @param seedPressuresBara at least two pressure levels bracketing a regular segment
   * @param family explicit hydrocarbon-water boundary family
   * @param additionalPointCount requested new strict points beyond the two continuation seeds
   * @return auditable local continuation result
   */
  public static JsonObject runRegularPressureDiagnostic(SystemInterface template, double[] temperaturesK,
      double[] seedPressuresBara, BoundaryFamily family, int additionalPointCount) {
    return runRegularPressureDiagnostic(template, temperaturesK, seedPressuresBara, family, additionalPointCount, true);
  }

  /** Runs the regular-pressure diagnostic in an explicit increasing- or decreasing-pressure direction. */
  public static JsonObject runRegularPressureDiagnostic(SystemInterface template, double[] temperaturesK,
      double[] seedPressuresBara, BoundaryFamily family, int additionalPointCount, boolean increasingPressure) {
    if (template == null || temperaturesK == null || temperaturesK.length < 5 || seedPressuresBara == null
        || seedPressuresBara.length < 2 || family == null || additionalPointCount < 1 || additionalPointCount > 20) {
      throw new IllegalArgumentException(
          "a template, five temperatures, two seed pressures, a family, and 1-20 new points are required");
    }
    long started = System.currentTimeMillis();
    HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        template).setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).setStableScanControls(32, 0.25)
        .setMaximumCorrectionIntervalAttempts(STABLE_BRACKET_CORRECTION_INTERVAL_ATTEMPTS)
        .setClusteringControls(MAXIMUM_MATCHED_COMPOSITION_JUMP, 1.75, 0.5)
        .discoverFromStableRegionTransitions(temperaturesK, seedPressuresBara);

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "stable-region-seed-discovery-fixed-pressure-strict-continuation");
    output.addProperty("family", family.name());
    output.addProperty("pressureDirection", increasingPressure ? "INCREASING" : "DECREASING");
    output.addProperty("stableFlashEvaluationCount", stable.getStableRegionScan().getFlashEvaluations());
    output.addProperty("stableScanFailureCount", stable.getStableRegionScan().getFailures().size());
    output.addProperty("strictAnchorCount", stable.getDiscovery().getCorrectedAnchorCount());
    output.add("strictCorrectionDiagnostics", correctionDiagnostics(stable));
    output.addProperty("additionalPointCountRequested", additionalPointCount);

    Branch selected = selectRegularPressureSeedBranch(stable.getDiscovery(), family);
    if (selected == null) {
      output.addProperty("completedRequestedPoints", false);
      output.addProperty("newStrictPointCount", 0);
      output.addProperty("failure", "NO_TWO_POINT_COMPOSITION_CONTINUOUS_BRANCH_FOR_REQUESTED_FAMILY");
      output.add("points", new JsonArray());
      output.add("attempts", new JsonArray());
      output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
      return output;
    }

    List<AnchorPoint> ordered = new ArrayList<AnchorPoint>(selected.getPoints());
    ordered.sort(Comparator.comparingDouble(AnchorPoint::getPressureBara));
    AnchorPoint previous = increasingPressure ? ordered.get(ordered.size() - 2) : ordered.get(1);
    AnchorPoint current = increasingPressure ? ordered.get(ordered.size() - 1) : ordered.get(0);
    HydrocarbonWaterRegularPressureBoundaryTracer.Result trace = new HydrocarbonWaterRegularPressureBoundaryTracer(
        template, family).setStepControls(0.04, 0.0025, 0.08, 6)
        .setCorrectionControls(11, 8.0, 0.04, MAXIMUM_MATCHED_COMPOSITION_JUMP, 30.0)
        .trace(previous, current, additionalPointCount);

    output.addProperty("seedBranch", selected.getIdentifier());
    output.addProperty("completedRequestedPoints", trace.hasCompletedRequestedPoints());
    output.addProperty("newStrictPointCount", Math.max(0, trace.getPoints().size() - 2));
    output.addProperty("terminationReason",
        trace.getTerminationReason() == null ? null : trace.getTerminationReason().name());
    output.addProperty("failure", trace.getFailureMessage());
    output.addProperty("finalLogPressureStep", trace.getFinalLogPressureStep());
    output.add("points", regularPressurePoints(trace.getPoints()));
    output.add("continuationStates", anchorContinuationStates(trace.getPoints()));
    output.add("restartStates", lastAnchorContinuationStates(trace.getPoints()));
    output.add("attempts", regularPressureAttempts(trace.getAttempts()));
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /**
   * Corrects one explicitly proven stable 2P/GOW temperature bracket without repeating a broad topology scan.
   *
   * <p>
   * The two endpoint flashes are recomputed from the current thermodynamic model, refined to an adjacent stable-region
   * bracket, and then passed to the same specified-two-phase/zero-third-phase-TPD corrector and independent global
   * stability classifier used by the full audit. Supplying a bracket therefore reduces discovery cost but never
   * bypasses a physical quality gate.
   * </p>
   */
  public static JsonObject runStableBracketCorrectionDiagnostic(SystemInterface template, double lowerTemperatureK,
      double upperTemperatureK, double pressureBara, BoundaryFamily family) {
    if (template == null || family == null || !Double.isFinite(lowerTemperatureK) || !Double.isFinite(upperTemperatureK)
        || lowerTemperatureK < 50.0 || upperTemperatureK <= lowerTemperatureK || !Double.isFinite(pressureBara)
        || pressureBara <= 0.0) {
      throw new IllegalArgumentException(
          "a template, ordered physical temperature bracket, pressure, and family are required");
    }
    long started = System.currentTimeMillis();
    HydrocarbonWaterStableRegionTransitionScanner.Result scan = new HydrocarbonWaterStableRegionTransitionScanner(
        template).setNumericalControls(32, 1.0e-4, 1.0e-9)
        .scan(new double[] { lowerTemperatureK, upperTemperatureK }, new double[] { pressureBara });
    HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket selected = null;
    for (HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket : scan.getBrackets()) {
      if (bracket.getFamily() != family) {
        continue;
      }
      if (selected != null) {
        throw new IllegalStateException("explicit stable bracket resolved to multiple requested-family transitions");
      }
      selected = bracket;
    }

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "explicit-stable-bracket-strict-correction");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("family", family.name());
    output.addProperty("pressureBara", pressureBara);
    output.addProperty("requestedLowerTemperatureK", lowerTemperatureK);
    output.addProperty("requestedUpperTemperatureK", upperTemperatureK);
    output.addProperty("stableFlashEvaluationCount", scan.getFlashEvaluations());
    output.addProperty("stableScanFailureCount", scan.getFailures().size());
    if (selected == null) {
      output.addProperty("strictCorrectionConverged", false);
      output.addProperty("ordinaryRootCount", 0);
      output.addProperty("failure", "REQUESTED_FAMILY_STABLE_ADJACENCY_NOT_REPRODUCED");
      output.add("roots", new JsonArray());
      output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
      return output;
    }

    output.addProperty("refinedLowerTemperatureK", selected.getLowerTemperatureState().getTemperatureK());
    output.addProperty("refinedUpperTemperatureK", selected.getUpperTemperatureState().getTemperatureK());
    output.addProperty("refinedBracketWidthK", selected.getTemperatureWidthK());
    HydrocarbonWaterStableRegionBoundaryCorrector.Result correction = new HydrocarbonWaterStableRegionBoundaryCorrector(
        template).setNumericalControls(8, 80, 1.0e-7, 1.0e-8)
        .setMaximumIntervalAttempts(STABLE_BRACKET_CORRECTION_INTERVAL_ATTEMPTS).correct(selected);
    output.addProperty("strictCorrectionConverged", correction.isConverged());
    output.addProperty("ordinaryRootCount", correction.getAcceptedOrdinaryRootCount());
    output.addProperty("classificationCount", correction.getAllClassifications().size());
    output.addProperty("intervalAttemptCount", correction.getIntervalAttempts().size());
    output.addProperty("failure", correction.getFailureMessage());
    JsonArray roots = new JsonArray();
    for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : correction.getAllClassifications()) {
      TwoToThreePhaseBoundaryPointSolver.Result root = classification.getBoundaryRoot();
      JsonObject row = new JsonObject();
      row.addProperty("temperatureK", root.getTemperatureK());
      row.addProperty("temperatureC", root.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", root.getPressureBara());
      row.addProperty("beta", root.getBeta());
      row.addProperty("retainedFlashResidual", root.getRetainedFlashResidual());
      row.addProperty("tangentPlaneDistance", root.getTangentPlaneDistance());
      row.addProperty("stationarityResidual", root.getStationarityResidual());
      row.addProperty("classification", classification.getClassification().name());
      row.addProperty("ordinaryBoundaryPoint", classification.isOrdinaryBoundaryPoint());
      row.addProperty("classificationFailure", classification.getFailureMessage());
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = classification.getGlobalStabilityResult();
      row.addProperty("globalStabilityAccepted", stability != null && stability.isAccepted());
      row.addProperty("minimumNonTrivialTangentPlaneDistance",
          stability == null ? Double.NaN : stability.getMinimumNonTrivialTangentPlaneDistance());
      row.add("continuationState", continuationState(State.from(root)));
      roots.add(row);
    }
    output.add("roots", roots);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /**
   * Re-corrects a serialized strict boundary state with the currently configured thermodynamic model.
   *
   * <p>
   * The serialized state is used only as a numerical seed. All retained-phase equilibrium, third-phase stationarity,
   * zero TPD, endpoint, and global-stability equations are recomputed from {@code template}. This provides an auditable
   * EOS-homotopy entry path when a stable multiphase TP flash changes local minima discontinuously and therefore cannot
   * supply a continuous two-phase onset seed.
   * </p>
   */
  public static JsonObject runModelSeededFixedPressureDiagnostic(SystemInterface template, State seed,
      BoundaryFamily family, double minimumTemperatureK, double maximumTemperatureK) {
    if (template == null || seed == null || family == null || !Double.isFinite(minimumTemperatureK)
        || !Double.isFinite(maximumTemperatureK) || minimumTemperatureK < 50.0
        || maximumTemperatureK <= minimumTemperatureK || seed.getRetainedPhaseZero() != family.getRetainedPhaseZero()
        || seed.getRetainedPhaseOne() != family.getRetainedPhaseOne()
        || seed.getIncipientPhase() != family.getIncipientPhase()) {
      throw new IllegalArgumentException(
          "a template, topology-matching seed, and valid temperature interval are required");
    }
    long started = System.currentTimeMillis();
    TwoToThreePhaseBoundaryPointSolver.RootSet rootSet = new TwoToThreePhaseBoundaryPointSolver(template,
        family.getRetainedPhaseZero(), family.getRetainedPhaseOne(), family.getIncipientPhase())
        .setNumericalControls(48, 100, 1.0e-7, 1.0e-8).solveAll(seed.getPressureBara(), minimumTemperatureK,
            maximumTemperatureK, seed.getBeta(), seed.getPhaseZeroComposition(), seed.getPhaseOneComposition(),
            Arrays.asList(seed.getIncipientComposition()));
    HydrocarbonWaterBoundaryEndpointClassifier classifier = new HydrocarbonWaterBoundaryEndpointClassifier(template)
        .setTolerances(1.0e-6, 1.0e-5, 1.0e-8);
    JsonArray roots = new JsonArray();
    int ordinaryRootCount = 0;
    for (TwoToThreePhaseBoundaryPointSolver.Result root : rootSet.getRoots()) {
      HydrocarbonWaterBoundaryEndpointClassifier.Result classification = classifier.classify(root);
      JsonObject row = new JsonObject();
      row.addProperty("temperatureK", root.getTemperatureK());
      row.addProperty("temperatureC", root.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", root.getPressureBara());
      row.addProperty("beta", root.getBeta());
      row.addProperty("retainedFlashResidual", root.getRetainedFlashResidual());
      row.addProperty("tangentPlaneDistance", root.getTangentPlaneDistance());
      row.addProperty("stationarityResidual", root.getStationarityResidual());
      row.addProperty("classification", classification.getClassification().name());
      row.addProperty("ordinaryBoundaryPoint", classification.isOrdinaryBoundaryPoint());
      row.addProperty("classificationFailure", classification.getFailureMessage());
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = classification.getGlobalStabilityResult();
      row.addProperty("globalStabilityAccepted", stability != null && stability.isAccepted());
      row.addProperty("minimumNonTrivialTangentPlaneDistance",
          stability == null ? Double.NaN : stability.getMinimumNonTrivialTangentPlaneDistance());
      row.add("continuationState", continuationState(State.from(root)));
      roots.add(row);
      ordinaryRootCount += classification.isOrdinaryBoundaryPoint() ? 1 : 0;
    }

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "serialized-boundary-state-model-homotopy-fixed-pressure-recorrection");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("seedOnly", true);
    output.addProperty("family", family.name());
    output.addProperty("seedTemperatureK", seed.getTemperatureK());
    output.addProperty("seedPressureBara", seed.getPressureBara());
    output.addProperty("minimumTemperatureK", minimumTemperatureK);
    output.addProperty("maximumTemperatureK", maximumTemperatureK);
    output.addProperty("flashEvaluationCount", rootSet.getFlashEvaluations());
    output.addProperty("stabilityEvaluationCount", rootSet.getStabilityEvaluations());
    output.addProperty("usableStabilityStateCount", rootSet.getUsableStabilityStates());
    output.addProperty("minimumTangentPlaneDistance", rootSet.getMinimumTangentPlaneDistance());
    output.addProperty("maximumTangentPlaneDistance", rootSet.getMaximumTangentPlaneDistance());
    output.addProperty("strictRootCount", rootSet.getRoots().size());
    output.addProperty("ordinaryRootCount", ordinaryRootCount);
    output.addProperty("failure", rootSet.getFailureMessage());
    output.add("roots", roots);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /** Continues regular-pressure tracing from two serialized strict states without repeating stable-region discovery. */
  public static JsonObject runRegularPressureRestartDiagnostic(SystemInterface template, State previous, State current,
      BoundaryFamily family, int additionalPointCount) {
    if (template == null || previous == null || current == null || family == null || additionalPointCount < 1
        || additionalPointCount > 20 || !sameTopology(previous, current)
        || previous.getRetainedPhaseZero() != family.getRetainedPhaseZero()
        || previous.getRetainedPhaseOne() != family.getRetainedPhaseOne()
        || previous.getIncipientPhase() != family.getIncipientPhase()) {
      throw new IllegalArgumentException("two topology-matching strict states, a family, and 1-20 points are required");
    }
    long started = System.currentTimeMillis();
    AnchorPoint previousAnchor = AnchorPoint.from(family,
        TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(previous, 1.0e-9), template);
    AnchorPoint currentAnchor = AnchorPoint.from(family,
        TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(current, 1.0e-9), template);
    HydrocarbonWaterRegularPressureBoundaryTracer.Result trace = new HydrocarbonWaterRegularPressureBoundaryTracer(
        template, family).setStepControls(0.04, 0.0025, 0.08, 6)
        .setCorrectionControls(11, 8.0, 0.04, MAXIMUM_MATCHED_COMPOSITION_JUMP, 30.0)
        .trace(previousAnchor, currentAnchor, additionalPointCount);

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "strict-fixed-pressure-checkpoint-restart");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("family", family.name());
    output.addProperty("pressureDirection",
        current.getPressureBara() >= previous.getPressureBara() ? "INCREASING" : "DECREASING");
    output.addProperty("additionalPointCountRequested", additionalPointCount);
    output.addProperty("completedRequestedPoints", trace.hasCompletedRequestedPoints());
    output.addProperty("newStrictPointCount", Math.max(0, trace.getPoints().size() - 2));
    output.addProperty("terminationReason",
        trace.getTerminationReason() == null ? null : trace.getTerminationReason().name());
    output.addProperty("failure", trace.getFailureMessage());
    output.addProperty("finalLogPressureStep", trace.getFinalLogPressureStep());
    output.add("points", regularPressurePoints(trace.getPoints()));
    output.add("continuationStates", anchorContinuationStates(trace.getPoints()));
    output.add("restartStates", lastAnchorContinuationStates(trace.getPoints()));
    output.add("attempts", regularPressureAttempts(trace.getAttempts()));
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /**
   * Runs a bounded bidirectional hybrid continuation diagnostic for one discovered boundary branch.
   *
   * <p>
   * The supplied pressure levels are used only to establish a composition-continuous strict anchor sequence. Dense
   * regular sequences are extended with pressure as the local parameter; pseudo-arclength is reserved for a pressure
   * turning point. This diagnostic intentionally omits global branch discovery and network promotion, so it can be used
   * to test endpoint reachability without paying for unrelated topology scans.
   * </p>
   *
   * @param template fully configured non-reactive EOS system
   * @param temperaturesK local stable-topology scan temperatures in kelvin
   * @param seedPressuresBara at least three pressure levels on the target branch
   * @param family explicit hydrocarbon-water boundary family
   * @param additionalPointsPerDirection bounded extension request from each anchor end
   * @param minimumTemperatureK declared diagnostic domain minimum temperature
   * @param maximumTemperatureK declared diagnostic domain maximum temperature
   * @param minimumPressureBara declared diagnostic domain minimum pressure
   * @param maximumPressureBara declared diagnostic domain maximum pressure
   * @return strict branch evidence and continuation diagnostics; never an engineering-promoted envelope
   */
  public static JsonObject runHybridBranchDiagnostic(SystemInterface template, double[] temperaturesK,
      double[] seedPressuresBara, BoundaryFamily family, int additionalPointsPerDirection, double minimumTemperatureK,
      double maximumTemperatureK, double minimumPressureBara, double maximumPressureBara) {
    if (template == null || temperaturesK == null || temperaturesK.length < 5 || seedPressuresBara == null
        || seedPressuresBara.length < 3 || family == null || additionalPointsPerDirection < 1
        || additionalPointsPerDirection > 100 || !Double.isFinite(minimumTemperatureK)
        || !Double.isFinite(maximumTemperatureK) || maximumTemperatureK <= minimumTemperatureK
        || !Double.isFinite(minimumPressureBara) || minimumPressureBara <= 0.0 || !Double.isFinite(maximumPressureBara)
        || maximumPressureBara <= minimumPressureBara) {
      throw new IllegalArgumentException(
          "a template, five temperatures, three seed pressures, a family, 1-100 points, and a valid PT domain are required");
    }
    long started = System.currentTimeMillis();
    HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        template).setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).setStableScanControls(32, 0.25)
        .setMaximumCorrectionIntervalAttempts(STABLE_BRACKET_CORRECTION_INTERVAL_ATTEMPTS)
        .setClusteringControls(MAXIMUM_MATCHED_COMPOSITION_JUMP, 1.75, 0.5)
        .discoverFromStableRegionTransitions(temperaturesK, seedPressuresBara);

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "target-branch-dense-anchor-pressure-continuation-turning-point-pseudo-arclength");
    output.addProperty("family", family.name());
    output.addProperty("engineeringEligible", false);
    output.addProperty("diagnosticOnly", true);
    output.addProperty("stableFlashEvaluationCount", stable.getStableRegionScan().getFlashEvaluations());
    output.addProperty("stableScanFailureCount", stable.getStableRegionScan().getFailures().size());
    output.addProperty("strictAnchorCount", stable.getDiscovery().getCorrectedAnchorCount());
    output.add("strictCorrectionDiagnostics", correctionDiagnostics(stable));
    output.addProperty("additionalPointsPerDirectionRequested", additionalPointsPerDirection);

    Branch selected = selectRegularPressureSeedBranch(stable.getDiscovery(), family);
    if (selected == null || selected.getPoints().size() < 3) {
      output.addProperty("assembled", false);
      output.addProperty("failure", "NO_THREE_POINT_COMPOSITION_CONTINUOUS_BRANCH_FOR_REQUESTED_FAMILY");
      output.add("discoveredBranches", branches(stable.getDiscovery()));
      output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
      return output;
    }

    HydrocarbonWaterBoundaryBranchAssembler.Result assembly;
    try {
      assembly = new HydrocarbonWaterBoundaryBranchAssembler(template).setCorrectorControls(80, 1.0e-8, 2.0e-5)
          .setStepControls(1.0, 1.0e-4, 2.0, 16, MAXIMUM_MATCHED_COMPOSITION_JUMP)
          .setDomainBounds(minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara)
          .assembleBranches(Arrays.asList(selected), additionalPointsPerDirection);
    } catch (RuntimeException error) {
      output.addProperty("assembled", false);
      output.addProperty("seedBranch", selected.getIdentifier());
      output.addProperty("failure", diagnostic(error));
      output.add("anchors", anchorPoints(selected.getPoints()));
      output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
      return output;
    }

    HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch = assembly.getBranches().get(0);
    output.addProperty("assembled", true);
    output.addProperty("seedBranch", selected.getIdentifier());
    output.addProperty("anchorCount", branch.getAnchorCount());
    output.addProperty("regularPressurePrimarySelected",
        branch.getBackwardTrace() == null && branch.getForwardTrace() == null
            && branch.getRegularBackwardTrace() != null && branch.getRegularForwardTrace() != null);
    output.addProperty("strictEvidencePointCount", branch.getEvidencePoints().size());
    output.addProperty("newStrictPointCount", Math.max(0, branch.getEvidencePoints().size() - branch.getAnchorCount()));
    output.add("anchors", anchorPoints(selected.getPoints()));
    output.add("points", evidencePoints(branch));
    output.add("initialPseudoArcBackwardTrace", traceDiagnostic(branch.getBackwardTrace(), null));
    output.add("initialPseudoArcForwardTrace", traceDiagnostic(branch.getForwardTrace(), null));
    output.add("regularPressureBackwardTrace", regularTraceDiagnostic(branch.getRegularBackwardTrace()));
    output.add("regularPressureForwardTrace", regularTraceDiagnostic(branch.getRegularForwardTrace()));
    output.add("postRegularPseudoArcBackwardTrace",
        traceDiagnostic(branch.getPostRegularBackwardTrace(), branch.getBackwardTermination()));
    output.add("postRegularPseudoArcForwardTrace",
        traceDiagnostic(branch.getPostRegularForwardTrace(), branch.getForwardTermination()));
    output.addProperty("failure", (String) null);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /** Continues directly from two serialized strict pseudo-arclength states without repeating branch discovery. */
  public static JsonObject runPseudoArcRestartDiagnostic(SystemInterface template, State previous, State current,
      int additionalPointCount, double initialArcStep, double minimumArcStep, double finiteDifferenceStep,
      double minimumTemperatureK, double maximumTemperatureK, double minimumPressureBara, double maximumPressureBara) {
    return runPseudoArcRestartDiagnostic(template, previous, current, additionalPointCount, initialArcStep,
        minimumArcStep, Math.max(2.0, initialArcStep), finiteDifferenceStep, minimumTemperatureK, maximumTemperatureK,
        minimumPressureBara, maximumPressureBara);
  }

  /** Continues from two strict states with an explicit maximum pseudo-arclength step. */
  public static JsonObject runPseudoArcRestartDiagnostic(SystemInterface template, State previous, State current,
      int additionalPointCount, double initialArcStep, double minimumArcStep, double maximumArcStep,
      double finiteDifferenceStep, double minimumTemperatureK, double maximumTemperatureK, double minimumPressureBara,
      double maximumPressureBara) {
    if (template == null || previous == null || current == null || additionalPointCount < 1
        || additionalPointCount > 200 || !Double.isFinite(initialArcStep) || initialArcStep <= 0.0
        || !Double.isFinite(minimumArcStep) || minimumArcStep <= 0.0 || minimumArcStep > initialArcStep
        || !Double.isFinite(maximumArcStep) || maximumArcStep < initialArcStep || !Double.isFinite(finiteDifferenceStep)
        || finiteDifferenceStep <= 0.0 || !Double.isFinite(minimumTemperatureK) || !Double.isFinite(maximumTemperatureK)
        || maximumTemperatureK <= minimumTemperatureK || !Double.isFinite(minimumPressureBara)
        || minimumPressureBara <= 0.0 || !Double.isFinite(maximumPressureBara)
        || maximumPressureBara <= minimumPressureBara
        || previous.getRetainedPhaseZero() != current.getRetainedPhaseZero()
        || previous.getRetainedPhaseOne() != current.getRetainedPhaseOne()
        || previous.getIncipientPhase() != current.getIncipientPhase()) {
      throw new IllegalArgumentException(
          "two topology-matching states, 1-200 points, and a valid PT domain are required");
    }
    long started = System.currentTimeMillis();
    TwoToThreePhasePseudoArcLengthTracer.Result trace = new TwoToThreePhasePseudoArcLengthTracer(template,
        current.getRetainedPhaseZero(), current.getRetainedPhaseOne(), current.getIncipientPhase())
        .setCorrectorControls(80, 1.0e-8, finiteDifferenceStep)
        .setStepControls(initialArcStep, minimumArcStep, maximumArcStep, 16)
        .setMaximumCompositionJump(MAXIMUM_MATCHED_COMPOSITION_JUMP)
        .setDomainBounds(minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara)
        .trace(previous, current, additionalPointCount);

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "strict-pseudo-arclength-checkpoint-restart");
    output.addProperty("solverContractVersion", "independent-replay-v1");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("additionalPointCountRequested", additionalPointCount);
    output.addProperty("initialArcStep", initialArcStep);
    output.addProperty("minimumArcStep", minimumArcStep);
    output.addProperty("maximumArcStep", maximumArcStep);
    output.addProperty("finiteDifferenceStep", finiteDifferenceStep);
    output.addProperty("newStrictPointCount", trace.getAcceptedCorrections().size());
    output.addProperty("completedRequestedPoints", trace.hasCompletedRequestedPoints());
    HydrocarbonWaterBoundaryTerminationClassifier.Result termination = new HydrocarbonWaterBoundaryTerminationClassifier(
        template).classify(trace, TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(current, 1.0e-9));
    output.add("trace", traceDiagnostic(trace, termination));
    output.add("points", pseudoArcPoints(trace.getPoints()));
    output.add("acceptedCorrectionDiagnostics", pseudoArcCorrectionDiagnostics(trace));
    output.add("continuationStates", continuationStates(trace.getPoints()));
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /**
   * Probes both null-space tangent orientations at a stopped continuation point without promoting either candidate.
   *
   * <p>
   * This diagnostic distinguishes a correctable secant failure from a branch point. Every converged candidate retains
   * its forward-orientation cosine and independent global-stability evidence; a backward candidate is never accepted as
   * continuation progress.
   * </p>
   */
  public static JsonObject runLocalTangentProbeDiagnostic(SystemInterface template, State previous, State current,
      double[] arcSteps, double finiteDifferenceStep) {
    if (template == null || previous == null || current == null || arcSteps == null || arcSteps.length == 0
        || arcSteps.length > 16 || !Double.isFinite(finiteDifferenceStep) || finiteDifferenceStep <= 0.0
        || previous.getRetainedPhaseZero() != current.getRetainedPhaseZero()
        || previous.getRetainedPhaseOne() != current.getRetainedPhaseOne()
        || previous.getIncipientPhase() != current.getIncipientPhase()) {
      throw new IllegalArgumentException(
          "two topology-matching states, 1-16 arc steps, and a finite difference are required");
    }
    for (double arcStep : arcSteps) {
      if (!Double.isFinite(arcStep) || arcStep <= 0.0) {
        throw new IllegalArgumentException("local tangent probe arc steps must be finite and positive");
      }
    }
    long started = System.currentTimeMillis();
    TwoToThreePhaseArcLengthCorrector corrector = new TwoToThreePhaseArcLengthCorrector(template,
        current.getRetainedPhaseZero(), current.getRetainedPhaseOne(), current.getIncipientPhase())
        .setNumericalControls(80, 1.0e-8, finiteDifferenceStep);
    HydrocarbonWaterBoundaryGlobalStabilityGate stabilityGate = new HydrocarbonWaterBoundaryGlobalStabilityGate(
        template);
    JsonArray probes = new JsonArray();
    int convergedCount = 0;
    int forwardCount = 0;
    int stableForwardCount = 0;
    for (double arcStep : arcSteps) {
      for (int orientation : new int[] { 1, -1 }) {
        TwoToThreePhaseArcLengthCorrector.Result correction = corrector.correctFromLocalTangent(current, arcStep,
            orientation);
        JsonObject row = new JsonObject();
        row.addProperty("arcStep", arcStep);
        row.addProperty("orientation", orientation);
        row.addProperty("converged", correction.isConverged());
        row.addProperty("residualConverged", correction.isResidualConverged());
        row.addProperty("distinctPhases", correction.hasDistinctPhases());
        row.addProperty("maximumResidual", correction.getMaximumResidual());
        row.addProperty("thermodynamicMaximumResidual", correction.getThermodynamicMaximumResidual());
        row.addProperty("arcLengthResidual", correction.getArcLengthResidual());
        row.addProperty("iterations", correction.getIterations());
        row.addProperty("jacobianConditionNumber", correction.getJacobianConditionNumber());
        row.addProperty("failure", correction.getFailureMessage());
        if (correction.isConverged()) {
          convergedCount++;
          State candidate = correction.getState();
          double forwardCosine = TwoToThreePhasePseudoArcLengthTracer.forwardCosine(previous, current, candidate);
          double compositionJump = stateCompositionJump(current, candidate);
          boolean forward = Double.isFinite(forwardCosine) && forwardCosine > 0.0;
          row.addProperty("forwardCosine", forwardCosine);
          row.addProperty("compositionJump", compositionJump);
          row.addProperty("forward", forward);
          row.add("candidateState", continuationState(candidate));
          if (forward) {
            forwardCount++;
          }
          TwoToThreePhaseBoundaryPointSolver.Result root = TwoToThreePhaseBoundaryPointSolver.Result
              .fromContinuationState(candidate, correction.getThermodynamicMaximumResidual());
          HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = stabilityGate.evaluate(root);
          row.addProperty("globalStabilityAccepted", stability.isAccepted());
          row.addProperty("globalStabilityFailure", stability.getFailureMessage());
          row.addProperty("targetMatched", stability.isTargetMatched());
          row.addProperty("targetCompositionDistance", stability.getTargetCompositionDistance());
          row.addProperty("minimumNonTrivialTangentPlaneDistance",
              stability.getMinimumNonTrivialTangentPlaneDistance());
          if (forward && stability.isAccepted() && compositionJump <= MAXIMUM_MATCHED_COMPOSITION_JUMP) {
            stableForwardCount++;
          }
        }
        probes.add(row);
      }
    }
    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "two-sided-local-null-tangent-branch-probe");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("convergedCandidateCount", convergedCount);
    output.addProperty("forwardCandidateCount", forwardCount);
    output.addProperty("stableForwardCandidateCount", stableForwardCount);
    output.addProperty("forwardContinuationAvailable", stableForwardCount > 0);
    output.add("previousState", continuationState(previous));
    output.add("currentState", continuationState(current));
    output.add("probes", probes);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /** Audits the local stable-region topology around one serialized continuation state. */
  public static JsonObject runRestartStateTopologyDiagnostic(SystemInterface template, State state,
      BoundaryFamily family) {
    if (template == null || state == null || family == null
        || state.getRetainedPhaseZero() != family.getRetainedPhaseZero()
        || state.getRetainedPhaseOne() != family.getRetainedPhaseOne()
        || state.getIncipientPhase() != family.getIncipientPhase()) {
      throw new IllegalArgumentException("a template and topology-matching boundary state are required");
    }
    long started = System.currentTimeMillis();
    TwoToThreePhaseBoundaryPointSolver.Result root = TwoToThreePhaseBoundaryPointSolver.Result
        .fromContinuationState(state, 1.0e-9);
    AnchorPoint anchor = AnchorPoint.from(family, root, template);
    HydrocarbonWaterIsolatedBoundaryRootClassifier.Result diagnosis = new HydrocarbonWaterIsolatedBoundaryRootClassifier(
        template).setNeighborhoodControls(2, 0.01, 33, 0.02, 4.0)
        .setQualityControls(MAXIMUM_MATCHED_COMPOSITION_JUMP, 0.90).classify(anchor);
    HydrocarbonWaterBoundaryEndpointClassifier.Result endpoint = diagnosis.getEndpointClassification();

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "strict-restart-state-local-stable-topology-audit");
    output.addProperty("family", family.name());
    output.addProperty("temperatureK", state.getTemperatureK());
    output.addProperty("temperatureC", state.getTemperatureK() - 273.15);
    output.addProperty("pressureBara", state.getPressureBara());
    output.addProperty("beta", state.getBeta());
    output.addProperty("endpointClassification", endpoint.getClassification().name());
    output.addProperty("minimumPhaseDistance", endpoint.getMinimumPhaseDistance());
    output.addProperty("phaseZeroOneDistance", endpoint.getPhaseZeroOneDistance());
    output.addProperty("phaseZeroIncipientDistance", endpoint.getPhaseZeroIncipientDistance());
    output.addProperty("phaseOneIncipientDistance", endpoint.getPhaseOneIncipientDistance());
    output.addProperty("localTopology", diagnosis.getTopology().name());
    output.addProperty("independentBranchSupport", diagnosis.hasIndependentBranchSupport());
    output.addProperty("lowerPressureSupport", diagnosis.hasLowerPressureSupport());
    output.addProperty("higherPressureSupport", diagnosis.hasHigherPressureSupport());
    output.addProperty("successfulGridFraction", diagnosis.getSuccessfulGridFraction());
    output.addProperty("matchingStableBracketCount", diagnosis.getMatchingBracketCount());
    output.addProperty("matchingStrictRootCount", diagnosis.getMatchingRoots().size());
    output.addProperty("rejectedCompositionCount", diagnosis.getRejectedCompositionCount());
    output.addProperty("diagnostic", diagnosis.getDiagnostic());
    if (diagnosis.getStableRegionScan() != null) {
      output.addProperty("stableFlashEvaluationCount", diagnosis.getStableRegionScan().getFlashEvaluations());
      output.addProperty("stableScanFailureCount", diagnosis.getStableRegionScan().getFailures().size());
    }
    JsonArray roots = new JsonArray();
    for (TwoToThreePhaseBoundaryPointSolver.Result matching : diagnosis.getMatchingRoots()) {
      JsonObject row = new JsonObject();
      row.addProperty("temperatureK", matching.getTemperatureK());
      row.addProperty("temperatureC", matching.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", matching.getPressureBara());
      row.addProperty("beta", matching.getBeta());
      row.addProperty("retainedFlashResidual", matching.getRetainedFlashResidual());
      row.addProperty("tangentPlaneDistance", matching.getTangentPlaneDistance());
      row.addProperty("stationarityResidual", matching.getStationarityResidual());
      roots.add(row);
    }
    output.add("matchingStrictRoots", roots);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /** Discovers strict restartable anchor sequences at explicitly targeted PT levels without continuing them. */
  public static JsonObject runStrictBranchSeedDiagnostic(SystemInterface template, double[] temperaturesK,
      double[] pressuresBara, BoundaryFamily family) {
    if (template == null || temperaturesK == null || temperaturesK.length < 5 || pressuresBara == null
        || pressuresBara.length < 2 || family == null) {
      throw new IllegalArgumentException("a template, five temperatures, two pressures, and a family are required");
    }
    long started = System.currentTimeMillis();
    HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable = new HydrocarbonWaterBoundaryAnchorDiscoverer(
        template).setCorrectionControls(32, 80, 1.0e-5, 1.0e-8).setStableScanControls(32, 0.25)
        .setMaximumCorrectionIntervalAttempts(STABLE_BRACKET_CORRECTION_INTERVAL_ATTEMPTS)
        .setClusteringControls(MAXIMUM_MATCHED_COMPOSITION_JUMP, 1.75, 0.5)
        .discoverFromStableRegionTransitions(temperaturesK, pressuresBara);

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "targeted-stable-region-strict-branch-seed-discovery");
    output.addProperty("family", family.name());
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("stableFlashEvaluationCount", stable.getStableRegionScan().getFlashEvaluations());
    output.addProperty("stableScanFailureCount", stable.getStableRegionScan().getFailures().size());
    output.add("strictCorrectionDiagnostics", correctionDiagnostics(stable));
    JsonArray branchRows = new JsonArray();
    int familyAnchorCount = 0;
    for (Branch branch : stable.getDiscovery().getBranches()) {
      if (branch.getFamily() != family) {
        continue;
      }
      JsonObject row = new JsonObject();
      row.addProperty("identifier", branch.getIdentifier());
      row.addProperty("anchorCount", branch.getPoints().size());
      row.addProperty("continuationSeedEligible", branch.canSeedContinuation());
      row.add("anchors", anchorPoints(branch.getPoints()));
      JsonArray states = new JsonArray();
      for (AnchorPoint point : branch.getPoints()) {
        states.add(continuationState(point.toContinuationState()));
      }
      row.add("restartStates", states);
      branchRows.add(row);
      familyAnchorCount += branch.getPoints().size();
    }
    output.addProperty("branchCount", branchRows.size());
    output.addProperty("strictAnchorCount", familyAnchorCount);
    output.add("branches", branchRows);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /** Tracks the explicitly seeded competing stationary branch over serialized strict boundary states. */
  public static JsonObject runSecondaryStationaryBranchDiagnostic(SystemInterface template, List<State> orderedStates,
      State rejectedState, CandidatePhase stationaryPhase, double[] terminalCompositionSeed) {
    if (template == null || orderedStates == null || orderedStates.isEmpty() || rejectedState == null
        || stationaryPhase == null || terminalCompositionSeed == null) {
      throw new IllegalArgumentException(
          "a template, ordered states, rejected state, phase, and terminal stationary composition are required");
    }
    long started = System.currentTimeMillis();
    List<TwoToThreePhaseBoundaryPointSolver.Result> roots = new ArrayList<TwoToThreePhaseBoundaryPointSolver.Result>();
    State preceding = null;
    for (State state : orderedStates) {
      if (state == null || preceding != null && !sameTopology(preceding, state)) {
        throw new IllegalArgumentException("ordered secondary-branch states must share one topology");
      }
      roots.add(TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(state, 1.0e-9));
      preceding = state;
    }
    if (preceding != null && !sameTopology(preceding, rejectedState)) {
      throw new IllegalArgumentException("rejected secondary-branch state must match the ordered-state topology");
    }
    roots.add(TwoToThreePhaseBoundaryPointSolver.Result.fromContinuationState(rejectedState, 1.0e-9));
    HydrocarbonWaterSecondaryStationaryBranchTracker.Result tracked = new HydrocarbonWaterSecondaryStationaryBranchTracker(
        template).setNumericalControls(160, 1000, 1.0e-10, 2.0e-5, 0.15, MAXIMUM_BRANCH_SWITCH_COMPOSITION_JUMP, 1.0e-8)
        .trackBackward(roots, stationaryPhase, terminalCompositionSeed);

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "explicit-secondary-stationary-branch-backward-tracking");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("stationaryPhase", stationaryPhase.name());
    output.addProperty("inputBoundaryRootCount", roots.size());
    output.addProperty("trackedSampleCount", tracked.getSamples().size());
    output.addProperty("zeroTpdBracketCount", tracked.getZeroTpdBrackets().size());
    output.addProperty("branchLimitBracketCount", tracked.getBranchLimitBrackets().size());
    output.addProperty("targetBranchMergeBracketCount", tracked.getTargetBranchMergeBrackets().size());
    output.addProperty("reachesFirstBoundaryRoot", tracked.reachesFirstBoundaryRoot());
    output.addProperty("termination", tracked.getTerminationMessage());
    JsonArray samples = new JsonArray();
    for (HydrocarbonWaterSecondaryStationaryBranchTracker.Sample sample : tracked.getSamples()) {
      JsonObject row = new JsonObject();
      row.addProperty("temperatureK", sample.getBoundaryRoot().getTemperatureK());
      row.addProperty("temperatureC", sample.getBoundaryRoot().getTemperatureK() - 273.15);
      row.addProperty("pressureBara", sample.getBoundaryRoot().getPressureBara());
      row.addProperty("tracked", sample.isTracked());
      row.addProperty("collapsedStationaryBranch", sample.isCollapsedStationaryBranch());
      row.addProperty("tangentPlaneDistance", sample.getTangentPlaneDistance());
      row.addProperty("compositionJump", sample.getCompositionJump());
      row.addProperty("targetCompositionDistance", sample.getTargetCompositionDistance());
      row.addProperty("collapsedTargetBranch", sample.isCollapsedTargetBranch());
      row.addProperty("failure", sample.getFailureMessage());
      if (sample.getCandidate() != null) {
        row.add("candidate", stationaryCandidate(sample.getCandidate()));
      }
      samples.add(row);
    }
    output.add("samples", samples);
    JsonArray zeroBrackets = new JsonArray();
    for (HydrocarbonWaterSecondaryStationaryBranchTracker.ZeroTpdBracket bracket : tracked.getZeroTpdBrackets()) {
      JsonObject row = new JsonObject();
      row.add("first", secondarySampleEndpoint(bracket.getFirst()));
      row.add("second", secondarySampleEndpoint(bracket.getSecond()));
      zeroBrackets.add(row);
    }
    output.add("zeroTpdBrackets", zeroBrackets);
    JsonArray branchLimits = new JsonArray();
    for (HydrocarbonWaterSecondaryStationaryBranchTracker.BranchLimitBracket bracket : tracked
        .getBranchLimitBrackets()) {
      JsonObject row = new JsonObject();
      row.add("collapsed", secondarySampleEndpoint(bracket.getCollapsed()));
      row.add("nonTrivial", secondarySampleEndpoint(bracket.getNonTrivial()));
      branchLimits.add(row);
    }
    output.add("branchLimitBrackets", branchLimits);
    JsonArray targetMerges = new JsonArray();
    for (HydrocarbonWaterSecondaryStationaryBranchTracker.TargetBranchMergeBracket bracket : tracked
        .getTargetBranchMergeBrackets()) {
      JsonObject row = new JsonObject();
      row.add("merged", secondarySampleEndpoint(bracket.getMerged()));
      row.add("distinct", secondarySampleEndpoint(bracket.getDistinct()));
      targetMerges.add(row);
    }
    output.add("targetBranchMergeBrackets", targetMerges);
    JsonObject trimming = new JsonObject();
    trimming.addProperty("attempted", !tracked.getTargetBranchMergeBrackets().isEmpty());
    if (!tracked.getTargetBranchMergeBrackets().isEmpty()) {
      HydrocarbonWaterSecondaryStationaryBranchTracker.TargetBranchMergeBracket merge = tracked
          .getTargetBranchMergeBrackets().get(0);
      EndpointEvidence endpoint = EndpointEvidence.targetBranchMerge("tracked-target-branch-merge",
          boundaryFamily(orderedStates.get(0)), merge, End.END);
      List<EvidencePoint> evidence = new ArrayList<EvidencePoint>();
      for (TwoToThreePhaseBoundaryPointSolver.Result root : roots) {
        evidence.add(EvidencePoint.from(root));
      }
      TargetBranchMergeTrimResult trim = new HydrocarbonWaterBoundaryNetworkAssembler().trimTargetBranchMerge(evidence,
          endpoint);
      trimming.addProperty("accepted", trim.isAccepted());
      trimming.addProperty("failure", trim.getFailureMessage());
      trimming.addProperty("originalPointCount", evidence.size());
      trimming.addProperty("retainedPointCount", trim.getRetainedPoints().size());
      trimming.addProperty("trimmedPointCount", trim.getTrimmedPointCount());
      trimming.addProperty("excludedPointIndex", trim.getExcludedPointIndex());
      trimming.addProperty("transformedDistance", trim.getTransformedDistance());
      trimming.addProperty("endpointType", endpoint.getType().name());
      trimming.addProperty("endpointEvidenceSource", endpoint.getSource().name());
      trimming.addProperty("requestedBranchEnd", endpoint.getRequestedBranchEnd().name());
      trimming.addProperty("endpointTemperatureK", endpoint.getTemperatureK());
      trimming.addProperty("endpointTemperatureC", endpoint.getTemperatureK() - 273.15);
      trimming.addProperty("endpointPressureBara", endpoint.getPressureBara());
      JsonArray retainedPoints = new JsonArray();
      JsonArray retainedStates = new JsonArray();
      BoundaryFamily trimFamily = boundaryFamily(orderedStates.get(0));
      for (EvidencePoint point : trim.getRetainedPoints()) {
        JsonObject pointJson = new JsonObject();
        pointJson.addProperty("temperatureK", point.getTemperatureK());
        pointJson.addProperty("temperatureC", point.getTemperatureK() - 273.15);
        pointJson.addProperty("pressureBara", point.getPressureBara());
        retainedPoints.add(pointJson);
        State state = State.create(trimFamily.getRetainedPhaseZero(), trimFamily.getRetainedPhaseOne(),
            trimFamily.getIncipientPhase(), point.getTemperatureK(), point.getPressureBara(), point.getBeta(),
            point.getPhaseZeroComposition(), point.getPhaseOneComposition(), point.getIncipientComposition());
        retainedStates.add(continuationState(state));
      }
      trimming.add("retainedPoints", retainedPoints);
      trimming.add("retainedStates", retainedStates);
    }
    output.add("targetBranchMergeTrimming", trimming);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /** Assembles two independently traced, identity-compatible paths into one diagnostic closed loop. */
  public static JsonObject runClosedBoundaryLoopDiagnostic(List<State> firstPath, List<State> secondPath) {
    return runClosedBoundaryLoopDiagnostic(firstPath, secondPath, false);
  }

  /** Assembles either a target-branch-merge loop or an ordinary-to-ordinary seam loop. */
  public static JsonObject runClosedBoundaryLoopDiagnostic(List<State> firstPath, List<State> secondPath,
      boolean ordinaryClosure) {
    long started = System.currentTimeMillis();
    HydrocarbonWaterClosedBoundaryLoopAssembler assembler = new HydrocarbonWaterClosedBoundaryLoopAssembler();
    HydrocarbonWaterClosedBoundaryLoopAssembler.Result assembled = ordinaryClosure
        ? assembler.assembleOrdinaryClosure(firstPath, secondPath)
        : assembler.assemble(firstPath, secondPath);
    JsonObject output = new JsonObject();
    output.addProperty("status", assembled.isAccepted() ? "success" : "rejected");
    output.addProperty("method", "identity-compatible-two-path-closed-loop-assembly");
    output.addProperty("closureContract", ordinaryClosure ? "ORDINARY_STATE_SEAM" : "TARGET_BRANCH_MERGE");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("accepted", assembled.isAccepted());
    output.addProperty("failure", assembled.getFailureMessage());
    output.addProperty("firstPathStateCount", firstPath == null ? 0 : firstPath.size());
    output.addProperty("secondPathStateCount", secondPath == null ? 0 : secondPath.size());
    output.addProperty("closedLoopStateCount", assembled.getStates().size());
    output.addProperty("mergePointIndex", assembled.getMergePointIndex());
    output.addProperty("mergePointType",
        assembled.getMergePointType() == null ? null : assembled.getMergePointType().name());
    output.addProperty("ordinaryEndpointPtDistance", assembled.getOrdinaryEndpointPtDistance());
    output.addProperty("ordinaryEndpointCompositionDistance", assembled.getOrdinaryEndpointCompositionDistance());
    output.addProperty("mergeEndpointPtDistance", assembled.getMergeEndpointPtDistance());
    output.addProperty("mergeEndpointCompositionDistance", assembled.getMergeEndpointCompositionDistance());
    JsonArray states = new JsonArray();
    JsonArray points = new JsonArray();
    for (State state : assembled.getStates()) {
      states.add(continuationState(state));
      JsonObject point = new JsonObject();
      point.addProperty("temperatureK", state.getTemperatureK());
      point.addProperty("temperatureC", state.getTemperatureK() - 273.15);
      point.addProperty("pressureBara", state.getPressureBara());
      points.add(point);
    }
    output.add("closedLoopStates", states);
    output.add("closedLoopPoints", points);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /** Assembles a closed loop and optionally performs a full independent per-state thermodynamic verification. */
  public static JsonObject runClosedBoundaryLoopDiagnostic(SystemInterface template, List<State> firstPath,
      List<State> secondPath, boolean verifyEvidence) {
    return runClosedBoundaryLoopDiagnostic(template, firstPath, secondPath, verifyEvidence, false);
  }

  /** Assembles and optionally verifies a loop under an explicit closure contract. */
  public static JsonObject runClosedBoundaryLoopDiagnostic(SystemInterface template, List<State> firstPath,
      List<State> secondPath, boolean verifyEvidence, boolean ordinaryClosure) {
    JsonObject output = runClosedBoundaryLoopDiagnostic(firstPath, secondPath, ordinaryClosure);
    output.addProperty("evidenceVerificationRequested", verifyEvidence);
    if (!verifyEvidence || !output.get("accepted").getAsBoolean()) {
      return output;
    }
    long started = System.currentTimeMillis();
    HydrocarbonWaterClosedBoundaryLoopAssembler assembler = new HydrocarbonWaterClosedBoundaryLoopAssembler();
    HydrocarbonWaterClosedBoundaryLoopAssembler.Result assembled = ordinaryClosure
        ? assembler.assembleOrdinaryClosure(firstPath, secondPath)
        : assembler.assemble(firstPath, secondPath);
    HydrocarbonWaterClosedBoundaryLoopVerifier.Result verification = new HydrocarbonWaterClosedBoundaryLoopVerifier(
        template).verify(assembled.getStates());
    JsonObject evidenceVerification = new JsonObject();
    evidenceVerification.addProperty("uniquePointCount", verification.getPointResults().size());
    evidenceVerification.addProperty("acceptedUniquePointCount", verification.getAcceptedUniquePointCount());
    evidenceVerification.addProperty("allUniquePointsAccepted", verification.areAllUniquePointsAccepted());
    evidenceVerification.addProperty("qualityGateAccepted", verification.getQualityReport().isAccepted());
    evidenceVerification.addProperty("internalQualityEligible", verification.isInternalQualityEligible());
    evidenceVerification.addProperty("benchmarkPending", true);
    evidenceVerification.addProperty("engineeringEligible", false);
    evidenceVerification.addProperty("maximumRetainedFlashResidual", verification.getMaximumRetainedFlashResidual());
    evidenceVerification.addProperty("maximumTargetTangentPlaneDistance",
        verification.getMaximumTargetTangentPlaneDistance());
    evidenceVerification.addProperty("maximumTargetStationarityResidual",
        verification.getMaximumTargetStationarityResidual());
    evidenceVerification.addProperty("minimumNonTrivialTangentPlaneDistance",
        verification.getMinimumNonTrivialTangentPlaneDistance());
    JsonArray envelopeViolations = new JsonArray();
    for (String violation : verification.getQualityReport().getViolations()) {
      envelopeViolations.add(violation);
    }
    evidenceVerification.add("envelopeViolations", envelopeViolations);
    JsonArray branchReports = new JsonArray();
    for (TwoToThreePhaseBoundaryQualityGate.BranchReport report : verification.getQualityReport().getBranchReports()) {
      JsonObject row = new JsonObject();
      row.addProperty("identifier", report.getIdentifier());
      row.addProperty("pointCount", report.getPointCount());
      row.addProperty("maximumResidual", report.getMaximumResidual());
      row.addProperty("maximumCompositionJump", report.getMaximumCompositionJump());
      row.addProperty("minimumPhaseDistance", report.getMinimumPhaseDistance());
      row.addProperty("projectedIntersectionCount", report.getProjectedIntersectionCount());
      row.addProperty("identityIntersectionCount", report.getIdentityIntersectionCount());
      row.addProperty("minimumProjectedIntersectionStateDistance",
          report.getMinimumProjectedIntersectionStateDistance());
      row.addProperty("accepted", report.isAccepted());
      JsonArray violations = new JsonArray();
      for (String violation : report.getViolations()) {
        violations.add(violation);
      }
      row.add("violations", violations);
      branchReports.add(row);
    }
    evidenceVerification.add("branchReports", branchReports);
    JsonArray points = new JsonArray();
    for (HydrocarbonWaterClosedBoundaryLoopVerifier.PointResult point : verification.getPointResults()) {
      JsonObject row = new JsonObject();
      row.addProperty("index", point.getIndex());
      row.addProperty("temperatureK", point.getInputState().getTemperatureK());
      row.addProperty("temperatureC", point.getInputState().getTemperatureK() - 273.15);
      row.addProperty("pressureBara", point.getInputState().getPressureBara());
      row.addProperty("accepted", point.isAccepted());
      row.addProperty("failure", point.getGlobalStability().getFailureMessage());
      row.addProperty("targetMatched", point.getGlobalStability().isTargetMatched());
      row.addProperty("targetCompositionDistance", point.getGlobalStability().getTargetCompositionDistance());
      row.addProperty("minimumNonTrivialTangentPlaneDistance",
          point.getGlobalStability().getMinimumNonTrivialTangentPlaneDistance());
      row.addProperty("retainedFlashResidual", point.getEvidence().getRetainedFlashResidual());
      row.addProperty("targetTangentPlaneDistance", point.getEvidence().getTangentPlaneDistance());
      row.addProperty("targetStationarityResidual", point.getEvidence().getStationarityResidual());
      addCompositionDiagnostics(row, "phaseZero", point.getEvidence().getPhaseZeroComposition());
      addCompositionDiagnostics(row, "phaseOne", point.getEvidence().getPhaseOneComposition());
      addCompositionDiagnostics(row, "incipient", point.getEvidence().getIncipientComposition());
      points.add(row);
    }
    evidenceVerification.add("points", points);
    long verificationTimeMs = System.currentTimeMillis() - started;
    evidenceVerification.addProperty("computationTimeMs", verificationTimeMs);
    output.add("evidenceVerification", evidenceVerification);
    output.addProperty("internalQualityEligible", verification.isInternalQualityEligible());
    output.addProperty("benchmarkPending", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("computationTimeMs", output.get("computationTimeMs").getAsLong() + verificationTimeMs);
    return output;
  }

  /** Independently verifies an ordered domain-exit-to-spinodal open boundary branch. */
  public static JsonObject runOpenBoundaryVerificationDiagnostic(SystemInterface template, List<State> orderedStates,
      double minimumPressureBara, State firstRejectedHighState) {
    long started = System.currentTimeMillis();
    HydrocarbonWaterOpenBoundaryVerifier.Result verification = new HydrocarbonWaterOpenBoundaryVerifier(template)
        .verify(orderedStates, minimumPressureBara, firstRejectedHighState);
    JsonObject output = new JsonObject();
    output.addProperty("status", verification.isInternalQualityEligible() ? "success" : "rejected");
    output.addProperty("method", "independent-open-boundary-state-and-endpoint-verification");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("ordinaryStateCount", verification.getPointResults().size());
    output.addProperty("acceptedOrdinaryStateCount", verification.getAcceptedPointCount());
    output.addProperty("allOrdinaryStatesAccepted", verification.areAllOrdinaryPointsAccepted());
    output.addProperty("domainExitVerified", verification.isDomainExitVerified());
    output.addProperty("minimumPressureBara", verification.getMinimumPressureBara());
    output.addProperty("qualityGateAccepted", verification.getQualityReport().isAccepted());
    output.addProperty("internalQualityEligible", verification.isInternalQualityEligible());
    output.addProperty("benchmarkPending", true);
    output.addProperty("engineeringEligible", false);
    output.addProperty("maximumRetainedFlashResidual", verification.getMaximumRetainedFlashResidual());
    output.addProperty("maximumTargetTangentPlaneDistance", verification.getMaximumTargetTangentPlaneDistance());
    output.addProperty("maximumTargetStationarityResidual", verification.getMaximumTargetStationarityResidual());
    output.addProperty("minimumNonTrivialTangentPlaneDistance",
        verification.getMinimumNonTrivialTangentPlaneDistance());

    HydrocarbonWaterBoundaryTerminationClassifier.Result termination = verification.getHighTermination();
    JsonObject endpoint = new JsonObject();
    endpoint.addProperty("type", termination.getType().name());
    endpoint.addProperty("physicalEndpoint", termination.isPhysicalEndpoint());
    endpoint.addProperty("temperatureK", termination.getTemperatureK());
    endpoint.addProperty("temperatureC", termination.getTemperatureK() - 273.15);
    endpoint.addProperty("pressureBara", termination.getPressureBara());
    endpoint.addProperty("destabilizingPhase",
        termination.getDestabilizingPhase() == null ? null : termination.getDestabilizingPhase().name());
    endpoint.addProperty("qualityMeasure", termination.getQualityMeasure());
    endpoint.addProperty("diagnostic", termination.getDiagnostic());
    output.add("highPressureEndpoint", endpoint);

    JsonArray envelopeViolations = new JsonArray();
    for (String violation : verification.getQualityReport().getViolations()) {
      envelopeViolations.add(violation);
    }
    output.add("envelopeViolations", envelopeViolations);
    JsonArray branchReports = new JsonArray();
    for (TwoToThreePhaseBoundaryQualityGate.BranchReport report : verification.getQualityReport().getBranchReports()) {
      JsonObject row = new JsonObject();
      row.addProperty("identifier", report.getIdentifier());
      row.addProperty("pointCount", report.getPointCount());
      row.addProperty("maximumResidual", report.getMaximumResidual());
      row.addProperty("maximumCompositionJump", report.getMaximumCompositionJump());
      row.addProperty("minimumPhaseDistance", report.getMinimumPhaseDistance());
      row.addProperty("projectedIntersectionCount", report.getProjectedIntersectionCount());
      row.addProperty("identityIntersectionCount", report.getIdentityIntersectionCount());
      row.addProperty("minimumProjectedIntersectionStateDistance",
          report.getMinimumProjectedIntersectionStateDistance());
      row.addProperty("accepted", report.isAccepted());
      JsonArray violations = new JsonArray();
      for (String violation : report.getViolations()) {
        violations.add(violation);
      }
      row.add("violations", violations);
      branchReports.add(row);
    }
    output.add("branchReports", branchReports);

    JsonArray points = new JsonArray();
    for (HydrocarbonWaterOpenBoundaryVerifier.PointResult point : verification.getPointResults()) {
      JsonObject row = new JsonObject();
      row.addProperty("index", point.getIndex());
      row.addProperty("temperatureK", point.getInputState().getTemperatureK());
      row.addProperty("temperatureC", point.getInputState().getTemperatureK() - 273.15);
      row.addProperty("pressureBara", point.getInputState().getPressureBara());
      row.addProperty("accepted", point.isAccepted());
      row.addProperty("failure", point.getGlobalStability().getFailureMessage());
      row.addProperty("targetMatched", point.getGlobalStability().isTargetMatched());
      row.addProperty("targetCompositionDistance", point.getGlobalStability().getTargetCompositionDistance());
      row.addProperty("minimumNonTrivialTangentPlaneDistance",
          point.getGlobalStability().getMinimumNonTrivialTangentPlaneDistance());
      row.addProperty("retainedFlashResidual", point.getEvidence().getRetainedFlashResidual());
      row.addProperty("targetTangentPlaneDistance", point.getEvidence().getTangentPlaneDistance());
      row.addProperty("targetStationarityResidual", point.getEvidence().getStationarityResidual());
      addCompositionDiagnostics(row, "phaseZero", point.getEvidence().getPhaseZeroComposition());
      addCompositionDiagnostics(row, "phaseOne", point.getEvidence().getPhaseOneComposition());
      addCompositionDiagnostics(row, "incipient", point.getEvidence().getIncipientComposition());
      points.add(row);
    }
    output.add("points", points);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  /** Independently verifies the OW closed loop and GO open branch under one combined network quality gate. */
  public static JsonObject runFullNetworkVerificationDiagnostic(SystemInterface template, List<State> closedLoopStates,
      List<State> openBranchStates, double minimumPressureBara, State firstRejectedHighState) {
    long started = System.currentTimeMillis();
    HydrocarbonWaterFullNetworkVerifier.Result verification = new HydrocarbonWaterFullNetworkVerifier(template)
        .verify(closedLoopStates, openBranchStates, minimumPressureBara, firstRejectedHighState);
    JsonObject output = new JsonObject();
    output.addProperty("status", verification.isInternalQualityEligible() ? "success" : "rejected");
    output.addProperty("method", "independent-two-branch-first-scope-network-verification");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("requiredBoundaryFamilyCount", 2);
    output.addProperty("verifiedBoundaryFamilyCount", 2);
    output.addProperty("internalQualityEligible", verification.isInternalQualityEligible());
    output.addProperty("benchmarkPending", true);
    output.addProperty("engineeringEligible", false);

    HydrocarbonWaterClosedBoundaryLoopVerifier.Result closed = verification.getClosedLoop();
    JsonObject closedSummary = new JsonObject();
    closedSummary.addProperty("boundaryFamily", BoundaryFamily.OW_TO_GOW.name());
    closedSummary.addProperty("uniquePointCount", closed.getPointResults().size());
    closedSummary.addProperty("acceptedUniquePointCount", closed.getAcceptedUniquePointCount());
    closedSummary.addProperty("allUniquePointsAccepted", closed.areAllUniquePointsAccepted());
    closedSummary.addProperty("qualityGateAccepted", closed.getQualityReport().isAccepted());
    closedSummary.addProperty("internalQualityEligible", closed.isInternalQualityEligible());
    closedSummary.addProperty("maximumRetainedFlashResidual", closed.getMaximumRetainedFlashResidual());
    closedSummary.addProperty("maximumTargetTangentPlaneDistance", closed.getMaximumTargetTangentPlaneDistance());
    closedSummary.addProperty("maximumTargetStationarityResidual", closed.getMaximumTargetStationarityResidual());
    closedSummary.addProperty("minimumNonTrivialTangentPlaneDistance",
        closed.getMinimumNonTrivialTangentPlaneDistance());
    output.add("closedLoop", closedSummary);

    HydrocarbonWaterOpenBoundaryVerifier.Result open = verification.getOpenBranch();
    JsonObject openSummary = new JsonObject();
    openSummary.addProperty("boundaryFamily", BoundaryFamily.GO_TO_GOW.name());
    openSummary.addProperty("ordinaryPointCount", open.getPointResults().size());
    openSummary.addProperty("acceptedOrdinaryPointCount", open.getAcceptedPointCount());
    openSummary.addProperty("allOrdinaryPointsAccepted", open.areAllOrdinaryPointsAccepted());
    openSummary.addProperty("domainExitVerified", open.isDomainExitVerified());
    openSummary.addProperty("qualityGateAccepted", open.getQualityReport().isAccepted());
    openSummary.addProperty("internalQualityEligible", open.isInternalQualityEligible());
    openSummary.addProperty("maximumRetainedFlashResidual", open.getMaximumRetainedFlashResidual());
    openSummary.addProperty("maximumTargetTangentPlaneDistance", open.getMaximumTargetTangentPlaneDistance());
    openSummary.addProperty("maximumTargetStationarityResidual", open.getMaximumTargetStationarityResidual());
    openSummary.addProperty("minimumNonTrivialTangentPlaneDistance", open.getMinimumNonTrivialTangentPlaneDistance());
    openSummary.addProperty("highEndpointType", open.getHighTermination().getType().name());
    openSummary.addProperty("highEndpointPhysical", open.getHighTermination().isPhysicalEndpoint());
    openSummary.addProperty("highEndpointTemperatureK", open.getHighTermination().getTemperatureK());
    openSummary.addProperty("highEndpointPressureBara", open.getHighTermination().getPressureBara());
    output.add("openBranch", openSummary);

    JsonObject combined = new JsonObject();
    combined.addProperty("accepted", verification.getCombinedQualityReport().isAccepted());
    JsonArray violations = new JsonArray();
    for (String violation : verification.getCombinedQualityReport().getViolations()) {
      violations.add(violation);
    }
    combined.add("violations", violations);
    JsonArray branchReports = new JsonArray();
    for (TwoToThreePhaseBoundaryQualityGate.BranchReport report : verification.getCombinedQualityReport()
        .getBranchReports()) {
      JsonObject row = new JsonObject();
      row.addProperty("identifier", report.getIdentifier());
      row.addProperty("pointCount", report.getPointCount());
      row.addProperty("maximumResidual", report.getMaximumResidual());
      row.addProperty("maximumCompositionJump", report.getMaximumCompositionJump());
      row.addProperty("minimumPhaseDistance", report.getMinimumPhaseDistance());
      row.addProperty("projectedIntersectionCount", report.getProjectedIntersectionCount());
      row.addProperty("identityIntersectionCount", report.getIdentityIntersectionCount());
      row.addProperty("accepted", report.isAccepted());
      JsonArray reportViolations = new JsonArray();
      for (String violation : report.getViolations()) {
        reportViolations.add(violation);
      }
      row.add("violations", reportViolations);
      branchReports.add(row);
    }
    combined.add("branchReports", branchReports);
    output.add("combinedQualityGate", combined);
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  private static void addCompositionDiagnostics(JsonObject output, String prefix, double[] composition) {
    double sum = 0.0;
    double minimum = Double.POSITIVE_INFINITY;
    double maximum = Double.NEGATIVE_INFINITY;
    for (double value : composition) {
      sum += value;
      minimum = Math.min(minimum, value);
      maximum = Math.max(maximum, value);
    }
    output.addProperty(prefix + "CompositionSum", sum);
    output.addProperty(prefix + "CompositionMinimum", minimum);
    output.addProperty(prefix + "CompositionMaximum", maximum);
  }

  /** Applies strict seed gates before continuing from one locally supported multivalued branch to another. */
  public static JsonObject runBranchSwitchDiagnostic(SystemInterface template, State source, State targetPrevious,
      State targetCurrent, int additionalPointCount, double initialArcStep, double minimumArcStep,
      double finiteDifferenceStep, double minimumTemperatureK, double maximumTemperatureK, double minimumPressureBara,
      double maximumPressureBara) {
    if (template == null || source == null || targetPrevious == null || targetCurrent == null) {
      throw new IllegalArgumentException("a template, one source state, and two target-branch states are required");
    }
    long started = System.currentTimeMillis();
    TwoToThreePhaseBoundaryPointSolver.Result sourceRoot = TwoToThreePhaseBoundaryPointSolver.Result
        .fromContinuationState(source, 1.0e-9);
    TwoToThreePhaseBoundaryPointSolver.Result targetPreviousRoot = TwoToThreePhaseBoundaryPointSolver.Result
        .fromContinuationState(targetPrevious, 1.0e-9);
    TwoToThreePhaseBoundaryPointSolver.Result targetCurrentRoot = TwoToThreePhaseBoundaryPointSolver.Result
        .fromContinuationState(targetCurrent, 1.0e-9);
    HydrocarbonWaterBoundaryEndpointClassifier classifier = new HydrocarbonWaterBoundaryEndpointClassifier(template);
    HydrocarbonWaterBoundaryEndpointClassifier.Result sourceClassification = classifier.classify(sourceRoot);
    HydrocarbonWaterBoundaryEndpointClassifier.Result targetPreviousClassification = classifier
        .classify(targetPreviousRoot);
    HydrocarbonWaterBoundaryEndpointClassifier.Result targetCurrentClassification = classifier
        .classify(targetCurrentRoot);
    double switchCompositionJump = compositionJump(sourceRoot, targetPreviousRoot);
    double switchPtDistance = ptDistance(source, targetPrevious);
    double targetSeedCompositionJump = compositionJump(targetPreviousRoot, targetCurrentRoot);
    double targetSeedPtDistance = ptDistance(targetPrevious, targetCurrent);
    boolean topologyMatches = sameTopology(source, targetPrevious) && sameTopology(targetPrevious, targetCurrent);
    boolean accepted = topologyMatches && sourceClassification.isOrdinaryBoundaryPoint()
        && targetPreviousClassification.isOrdinaryBoundaryPoint()
        && targetCurrentClassification.isOrdinaryBoundaryPoint()
        && switchCompositionJump <= MAXIMUM_BRANCH_SWITCH_COMPOSITION_JUMP
        && switchPtDistance <= MAXIMUM_BRANCH_SWITCH_PT_DISTANCE
        && targetSeedCompositionJump <= MAXIMUM_BRANCH_SWITCH_COMPOSITION_JUMP
        && targetSeedPtDistance >= MINIMUM_BRANCH_SWITCH_TARGET_SEED_PT_DISTANCE
        && targetSeedPtDistance <= MAXIMUM_BRANCH_SWITCH_PT_DISTANCE;

    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "strict-multivalued-boundary-branch-switch");
    output.addProperty("diagnosticOnly", true);
    output.addProperty("engineeringEligible", false);
    JsonObject gate = new JsonObject();
    gate.addProperty("passed", accepted);
    gate.addProperty("topologyMatches", topologyMatches);
    gate.addProperty("sourceClassification", sourceClassification.getClassification().name());
    gate.addProperty("targetPreviousClassification", targetPreviousClassification.getClassification().name());
    gate.addProperty("targetCurrentClassification", targetCurrentClassification.getClassification().name());
    gate.addProperty("sourceGlobalStabilityAccepted", sourceClassification.getGlobalStabilityResult() != null
        && sourceClassification.getGlobalStabilityResult().isAccepted());
    gate.addProperty("targetPreviousGlobalStabilityAccepted",
        targetPreviousClassification.getGlobalStabilityResult() != null
            && targetPreviousClassification.getGlobalStabilityResult().isAccepted());
    gate.addProperty("targetCurrentGlobalStabilityAccepted",
        targetCurrentClassification.getGlobalStabilityResult() != null
            && targetCurrentClassification.getGlobalStabilityResult().isAccepted());
    gate.addProperty("switchCompositionJump", switchCompositionJump);
    gate.addProperty("targetSeedCompositionJump", targetSeedCompositionJump);
    gate.addProperty("maximumCompositionJumpLimit", MAXIMUM_BRANCH_SWITCH_COMPOSITION_JUMP);
    gate.addProperty("switchPtDistance", switchPtDistance);
    gate.addProperty("targetSeedPtDistance", targetSeedPtDistance);
    gate.addProperty("minimumTargetSeedPtDistance", MINIMUM_BRANCH_SWITCH_TARGET_SEED_PT_DISTANCE);
    gate.addProperty("ptDistanceLimit", MAXIMUM_BRANCH_SWITCH_PT_DISTANCE);
    output.add("seedQualityGate", gate);
    output.add("sourceState", continuationState(source));
    output.add("targetPreviousState", continuationState(targetPrevious));
    output.add("targetCurrentState", continuationState(targetCurrent));
    if (!accepted) {
      output.addProperty("switchAccepted", false);
      output.addProperty("failure", "branch-switch seed quality gate rejected the proposed connection");
      output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
      return output;
    }

    JsonObject continuation = runPseudoArcRestartDiagnostic(template, targetPrevious, targetCurrent,
        additionalPointCount, initialArcStep, minimumArcStep, finiteDifferenceStep, minimumTemperatureK,
        maximumTemperatureK, minimumPressureBara, maximumPressureBara);
    JsonObject continuationGate = branchSwitchContinuationQualityGate(targetPrevious, targetCurrent, continuation);
    boolean continuationAccepted = continuationGate.get("passed").getAsBoolean();
    output.addProperty("switchAccepted", continuationAccepted);
    output.addProperty("newStrictPointCount", continuation.get("newStrictPointCount").getAsInt());
    output.addProperty("usableNewStrictPointCount",
        continuationAccepted ? continuation.get("newStrictPointCount").getAsInt() : 0);
    output.addProperty("completedRequestedPoints", continuation.get("completedRequestedPoints").getAsBoolean());
    output.add("continuationQualityGate", continuationGate);
    output.add("continuation", continuation);
    if (!continuationAccepted) {
      output.addProperty("failure", "target-branch continuation did not make distinct forward progress");
    }
    output.addProperty("computationTimeMs", System.currentTimeMillis() - started);
    return output;
  }

  private static JsonObject branchSwitchContinuationQualityGate(State previous, State current,
      JsonObject continuation) {
    JsonObject gate = new JsonObject();
    double temperatureScale = Math.max(1.0, current.getTemperatureK());
    double seedTemperature = (current.getTemperatureK() - previous.getTemperatureK()) / temperatureScale;
    double seedPressure = Math.log(current.getPressureBara() / previous.getPressureBara());
    double seedNormSquared = seedTemperature * seedTemperature + seedPressure * seedPressure;
    double seedNorm = Math.sqrt(seedNormSquared);
    double maximumForwardAdvanceRatio = Double.NEGATIVE_INFINITY;
    double maximumDistanceRatio = 0.0;
    int distinctForwardPointCount = 0;
    JsonArray points = continuation.getAsJsonArray("points");
    double lastTemperatureK = current.getTemperatureK();
    double lastPressureBara = current.getPressureBara();
    for (int index = 2; points != null && index < points.size(); index++) {
      JsonObject point = points.get(index).getAsJsonObject();
      double temperatureK = point.get("temperatureK").getAsDouble();
      double pressureBara = point.get("pressureBara").getAsDouble();
      double deltaTemperature = (temperatureK - current.getTemperatureK()) / temperatureScale;
      double deltaPressure = Math.log(pressureBara / current.getPressureBara());
      double forwardRatio = seedNormSquared > 0.0
          ? (deltaTemperature * seedTemperature + deltaPressure * seedPressure) / seedNormSquared
          : Double.NEGATIVE_INFINITY;
      double distanceRatio = seedNorm > 0.0 ? Math.hypot(deltaTemperature, deltaPressure) / seedNorm : 0.0;
      maximumForwardAdvanceRatio = Math.max(maximumForwardAdvanceRatio, forwardRatio);
      maximumDistanceRatio = Math.max(maximumDistanceRatio, distanceRatio);
      double segmentTemperature = (temperatureK - lastTemperatureK) / temperatureScale;
      double segmentPressure = Math.log(pressureBara / lastPressureBara);
      if (seedNorm > 0.0 && Math.hypot(segmentTemperature, segmentPressure) / seedNorm >= 1.0e-3
          && forwardRatio > 0.0) {
        distinctForwardPointCount++;
      }
      lastTemperatureK = temperatureK;
      lastPressureBara = pressureBara;
    }
    boolean finite = Double.isFinite(seedNorm) && Double.isFinite(maximumForwardAdvanceRatio)
        && Double.isFinite(maximumDistanceRatio);
    boolean passed = finite && continuation.get("newStrictPointCount").getAsInt() > 0
        && maximumForwardAdvanceRatio >= MINIMUM_BRANCH_SWITCH_FORWARD_ADVANCE_RATIO && distinctForwardPointCount > 0;
    gate.addProperty("passed", passed);
    gate.addProperty("seedPtNorm", seedNorm);
    gate.addProperty("maximumForwardAdvanceRatio", maximumForwardAdvanceRatio);
    gate.addProperty("minimumForwardAdvanceRatio", MINIMUM_BRANCH_SWITCH_FORWARD_ADVANCE_RATIO);
    gate.addProperty("maximumDistanceRatio", maximumDistanceRatio);
    gate.addProperty("distinctForwardPointCount", distinctForwardPointCount);
    return gate;
  }

  private static boolean sameTopology(State first, State second) {
    return first.getRetainedPhaseZero() == second.getRetainedPhaseZero()
        && first.getRetainedPhaseOne() == second.getRetainedPhaseOne()
        && first.getIncipientPhase() == second.getIncipientPhase();
  }

  private static BoundaryFamily boundaryFamily(State state) {
    for (BoundaryFamily family : BoundaryFamily.values()) {
      if (family.getRetainedPhaseZero() == state.getRetainedPhaseZero()
          && family.getRetainedPhaseOne() == state.getRetainedPhaseOne()
          && family.getIncipientPhase() == state.getIncipientPhase()) {
        return family;
      }
    }
    throw new IllegalArgumentException("unsupported hydrocarbon-water continuation topology");
  }

  private static double ptDistance(State first, State second) {
    return Math.abs(second.getTemperatureK() - first.getTemperatureK()) / Math.max(1.0, first.getTemperatureK())
        + Math.abs(Math.log(second.getPressureBara() / first.getPressureBara()));
  }

  private static Branch selectRegularPressureSeedBranch(HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery,
      BoundaryFamily family) {
    Branch selected = null;
    for (Branch branch : discovery.getBranches()) {
      if (branch.getFamily() == family && branch.canSeedContinuation()
          && (selected == null || branch.getPoints().size() > selected.getPoints().size())) {
        selected = branch;
      }
    }
    return selected;
  }

  private static JsonArray regularPressurePoints(List<AnchorPoint> points) {
    JsonArray output = new JsonArray();
    for (int index = 0; index < points.size(); index++) {
      AnchorPoint point = points.get(index);
      JsonObject row = new JsonObject();
      row.addProperty("role", index < 2 ? "seed" : "continued");
      row.addProperty("temperatureK", point.getTemperatureK());
      row.addProperty("temperatureC", point.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", point.getPressureBara());
      row.addProperty("beta", point.getBeta());
      row.addProperty("retainedFlashResidual", point.getRetainedFlashResidual());
      row.addProperty("tangentPlaneDistance", point.getTangentPlaneDistance());
      row.addProperty("stationarityResidual", point.getStationarityResidual());
      row.addProperty("globalStabilityAccepted", point.getGlobalStabilityResult().isAccepted());
      output.add(row);
    }
    return output;
  }

  private static JsonArray anchorContinuationStates(List<AnchorPoint> points) {
    JsonArray output = new JsonArray();
    for (AnchorPoint point : points) {
      output.add(continuationState(point.toContinuationState()));
    }
    return output;
  }

  private static JsonArray lastAnchorContinuationStates(List<AnchorPoint> points) {
    JsonArray output = new JsonArray();
    int first = Math.max(0, points.size() - 2);
    for (int index = first; index < points.size(); index++) {
      output.add(continuationState(points.get(index).toContinuationState()));
    }
    return output;
  }

  private static JsonArray anchorPoints(List<AnchorPoint> points) {
    JsonArray output = new JsonArray();
    for (AnchorPoint point : points) {
      JsonObject row = new JsonObject();
      row.addProperty("temperatureK", point.getTemperatureK());
      row.addProperty("temperatureC", point.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", point.getPressureBara());
      row.addProperty("beta", point.getBeta());
      row.addProperty("retainedFlashResidual", point.getRetainedFlashResidual());
      row.addProperty("tangentPlaneDistance", point.getTangentPlaneDistance());
      row.addProperty("stationarityResidual", point.getStationarityResidual());
      row.addProperty("globalStabilityAccepted", point.getGlobalStabilityResult().isAccepted());
      output.add(row);
    }
    return output;
  }

  private static JsonArray evidencePoints(HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch) {
    JsonArray output = new JsonArray();
    for (EvidencePoint point : branch.getEvidencePoints()) {
      JsonObject row = new JsonObject();
      row.addProperty("temperatureK", point.getTemperatureK());
      row.addProperty("temperatureC", point.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", point.getPressureBara());
      row.addProperty("beta", point.getBeta());
      row.addProperty("retainedFlashResidual", point.getRetainedFlashResidual());
      row.addProperty("tangentPlaneDistance", point.getTangentPlaneDistance());
      row.addProperty("stationarityResidual", point.getStationarityResidual());
      row.addProperty("curve", branch.getFamily().name());
      output.add(row);
    }
    return output;
  }

  private static JsonArray pseudoArcPoints(List<State> points) {
    JsonArray output = new JsonArray();
    for (int index = 0; index < points.size(); index++) {
      State point = points.get(index);
      JsonObject row = new JsonObject();
      row.addProperty("role", index < 2 ? "seed" : "continued");
      row.addProperty("temperatureK", point.getTemperatureK());
      row.addProperty("temperatureC", point.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", point.getPressureBara());
      row.addProperty("beta", point.getBeta());
      output.add(row);
    }
    return output;
  }

  private static JsonArray pseudoArcCorrectionDiagnostics(TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    JsonArray output = new JsonArray();
    List<State> points = trace.getPoints();
    List<neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.Result> corrections = trace
        .getAcceptedCorrections();
    for (int correctionIndex = 0; correctionIndex < corrections.size(); correctionIndex++) {
      neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.Result correction = corrections
          .get(correctionIndex);
      JsonObject row = new JsonObject();
      row.addProperty("pointIndex", correctionIndex + 2);
      row.addProperty("arcStep", correction.getArcStep());
      row.addProperty("iterations", correction.getIterations());
      row.addProperty("jacobianConditionNumber", correction.getJacobianConditionNumber());
      row.addProperty("initialMaximumResidual", correction.getInitialMaximumResidual());
      row.addProperty("maximumResidual", correction.getMaximumResidual());
      row.addProperty("thermodynamicMaximumResidual", correction.getThermodynamicMaximumResidual());
      row.addProperty("arcLengthResidual", correction.getArcLengthResidual());
      if (correctionIndex + 2 < points.size()) {
        row.addProperty("forwardCosine", TwoToThreePhasePseudoArcLengthTracer.forwardCosine(points.get(correctionIndex),
            points.get(correctionIndex + 1), points.get(correctionIndex + 2)));
      }
      output.add(row);
    }
    return output;
  }

  private static JsonArray continuationStates(List<State> points) {
    JsonArray output = new JsonArray();
    for (State point : points) {
      output.add(continuationState(point));
    }
    return output;
  }

  private static JsonObject stationaryCandidate(Candidate candidate) {
    JsonObject output = new JsonObject();
    output.addProperty("seedPhase", candidate.getSeedPhase().name());
    output.addProperty("phase", candidate.getPhase().name());
    output.addProperty("tangentPlaneDistance", candidate.getTangentPlaneDistance());
    output.addProperty("stationarityResidual", candidate.getStationarityResidual());
    output.addProperty("stabilityFunction", candidate.getStabilityFunction());
    output.addProperty("iterations", candidate.getIterations());
    output.addProperty("converged", candidate.isConverged());
    output.addProperty("trivial", candidate.isTrivial());
    output.addProperty("failure", candidate.getFailureMessage());
    output.add("composition", numbers(candidate.getComposition()));
    return output;
  }

  private static JsonObject secondarySampleEndpoint(HydrocarbonWaterSecondaryStationaryBranchTracker.Sample sample) {
    JsonObject output = new JsonObject();
    output.addProperty("temperatureK", sample.getBoundaryRoot().getTemperatureK());
    output.addProperty("temperatureC", sample.getBoundaryRoot().getTemperatureK() - 273.15);
    output.addProperty("pressureBara", sample.getBoundaryRoot().getPressureBara());
    output.addProperty("tangentPlaneDistance", sample.getTangentPlaneDistance());
    output.addProperty("targetCompositionDistance", sample.getTargetCompositionDistance());
    output.addProperty("tracked", sample.isTracked());
    if (sample.getCandidate() != null) {
      output.add("candidate", stationaryCandidate(sample.getCandidate()));
    }
    return output;
  }

  private static JsonArray regularPressureAttempts(
      List<HydrocarbonWaterRegularPressureBoundaryTracer.Attempt> attempts) {
    JsonArray output = new JsonArray();
    for (HydrocarbonWaterRegularPressureBoundaryTracer.Attempt attempt : attempts) {
      JsonObject row = new JsonObject();
      row.addProperty("pressureBara", attempt.getPressureBara());
      row.addProperty("predictedTemperatureK", attempt.getPredictedTemperatureK());
      row.addProperty("logPressureStep", attempt.getLogPressureStep());
      row.addProperty("stableFlashEvaluations", attempt.getStableFlashEvaluations());
      row.addProperty("stableScanFailures", attempt.getStableScanFailures());
      row.addProperty("matchingBracketCount", attempt.getMatchingBracketCount());
      row.addProperty("strictOrdinaryRootCount", attempt.getStrictOrdinaryRootCount());
      row.addProperty("accepted", attempt.getAcceptedAnchor() != null);
      row.addProperty("failure", attempt.getFailureMessage());
      output.add(row);
    }
    return output;
  }

  private static IsolatedTopologyAudit isolatedRootDiagnostics(SystemInterface template,
      HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery) {
    JsonArray output = new JsonArray();
    List<AnchorPoint> independentAnchors = new ArrayList<AnchorPoint>();
    for (Branch branch : discovery.getBranches()) {
      if (!branch.isIsolated()) {
        continue;
      }
      AnchorPoint anchor = branch.getPoints().get(0);
      JsonObject row = new JsonObject();
      row.addProperty("identifier", branch.getIdentifier());
      row.addProperty("family", branch.getFamily().name());
      row.addProperty("temperatureK", anchor.getTemperatureK());
      row.addProperty("temperatureC", anchor.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", anchor.getPressureBara());
      try {
        HydrocarbonWaterIsolatedBoundaryRootClassifier.Result diagnosis = new HydrocarbonWaterIsolatedBoundaryRootClassifier(
            template).setNeighborhoodControls(2, 0.08, 17, 0.10, 12.0)
            .setQualityControls(MAXIMUM_MATCHED_COMPOSITION_JUMP, 0.90).classify(anchor);
        row.addProperty("topology", diagnosis.getTopology().name());
        row.addProperty("endpointClassification", diagnosis.getEndpointClassification().getClassification().name());
        row.addProperty("independentBranchSupport", diagnosis.hasIndependentBranchSupport());
        row.addProperty("lowerPressureSupport", diagnosis.hasLowerPressureSupport());
        row.addProperty("higherPressureSupport", diagnosis.hasHigherPressureSupport());
        row.addProperty("successfulGridFraction", diagnosis.getSuccessfulGridFraction());
        row.addProperty("matchingStableBracketCount", diagnosis.getMatchingBracketCount());
        row.addProperty("strictCorrectionCount", diagnosis.getBracketCorrections().size());
        row.addProperty("matchingStrictRootCount", diagnosis.getMatchingRoots().size());
        row.addProperty("rejectedCompositionCount", diagnosis.getRejectedCompositionCount());
        row.addProperty("diagnostic", diagnosis.getDiagnostic());
        if (diagnosis.getStableRegionScan() != null) {
          row.addProperty("stableFlashEvaluationCount", diagnosis.getStableRegionScan().getFlashEvaluations());
          row.addProperty("stableScanFailureCount", diagnosis.getStableRegionScan().getFailures().size());
          row.addProperty("allStableBracketCount", diagnosis.getStableRegionScan().getBrackets().size());
        }
        JsonArray roots = new JsonArray();
        for (TwoToThreePhaseBoundaryPointSolver.Result root : diagnosis.getMatchingRoots()) {
          if (diagnosis.hasIndependentBranchSupport()) {
            independentAnchors.add(AnchorPoint.from(branch.getFamily(), root, template));
          }
          JsonObject rootJson = new JsonObject();
          rootJson.addProperty("temperatureK", root.getTemperatureK());
          rootJson.addProperty("temperatureC", root.getTemperatureK() - 273.15);
          rootJson.addProperty("pressureBara", root.getPressureBara());
          rootJson.addProperty("beta", root.getBeta());
          rootJson.addProperty("retainedFlashResidual", root.getRetainedFlashResidual());
          rootJson.addProperty("tangentPlaneDistance", root.getTangentPlaneDistance());
          rootJson.addProperty("stationarityResidual", root.getStationarityResidual());
          roots.add(rootJson);
        }
        row.add("matchingStrictRoots", roots);
        JsonArray corrections = new JsonArray();
        for (HydrocarbonWaterIsolatedBoundaryRootClassifier.BracketCorrection correction : diagnosis
            .getBracketCorrections()) {
          JsonObject correctionJson = new JsonObject();
          correctionJson.addProperty("pressureBara", correction.getBracket().getPressureBara());
          correctionJson.addProperty("bracketMidpointTemperatureK", correction.getBracket().getMidpointTemperatureK());
          correctionJson.addProperty("bracketWidthK", correction.getBracket().getTemperatureWidthK());
          correctionJson.addProperty("strictCorrectionConverged", correction.getCorrection().isConverged());
          correctionJson.addProperty("ordinaryRootCount", correction.getCorrection().getAcceptedOrdinaryRootCount());
          correctionJson.addProperty("compositionContinuousRootCount",
              correction.getAcceptedCompositionContinuousRootCount());
          correctionJson.addProperty("failure", correction.getCorrection().getFailureMessage());
          corrections.add(correctionJson);
        }
        row.add("bracketCorrections", corrections);
      } catch (RuntimeException error) {
        row.addProperty("topology", "INCONCLUSIVE");
        row.addProperty("independentBranchSupport", false);
        row.addProperty("diagnostic", diagnostic(error));
        row.add("matchingStrictRoots", new JsonArray());
        row.add("bracketCorrections", new JsonArray());
      }
      output.add(row);
    }
    return new IsolatedTopologyAudit(output, independentAnchors);
  }

  private static JsonArray correctionDiagnostics(
      HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable) {
    JsonArray output = new JsonArray();
    for (int index = 0; index < stable.getCorrections().size(); index++) {
      HydrocarbonWaterStableRegionBoundaryCorrector.Result correction = stable.getCorrections().get(index);
      JsonObject row = new JsonObject();
      row.addProperty("family", correction.getBracket().getFamily().name());
      row.addProperty("pressureBara", correction.getBracket().getPressureBara());
      row.addProperty("bracketMidpointTemperatureK", correction.getBracket().getMidpointTemperatureK());
      row.addProperty("bracketWidthK", correction.getBracket().getTemperatureWidthK());
      row.addProperty("elapsedTimeMs", stable.getCorrectionTimesMs().get(index));
      row.addProperty("intervalAttemptCount", correction.getIntervalAttempts().size());
      row.addProperty("acceptedOrdinaryRootCount", correction.getAcceptedOrdinaryRootCount());
      row.addProperty("classificationCount", correction.getAllClassifications().size());
      row.addProperty("classificationEvaluationCount", correction.getClassificationEvaluationCount());
      row.addProperty("classificationReuseCount", correction.getClassificationReuseCount());
      row.addProperty("failure", correction.getFailureMessage());
      JsonArray intervalAttempts = new JsonArray();
      for (HydrocarbonWaterStableRegionBoundaryCorrector.IntervalAttempt attempt : correction.getIntervalAttempts()) {
        JsonObject attemptJson = new JsonObject();
        attemptJson.addProperty("seedSource", attempt.getSeedSource());
        attemptJson.addProperty("minimumTemperatureK", attempt.getMinimumTemperatureK());
        attemptJson.addProperty("maximumTemperatureK", attempt.getMaximumTemperatureK());
        attemptJson.addProperty("acceptedOrdinaryRootCount", attempt.getAcceptedOrdinaryRootCount());
        TwoToThreePhaseBoundaryPointSolver.RootSet roots = attempt.getRootSet();
        attemptJson.addProperty("strictRootCount", roots == null ? 0 : roots.getRoots().size());
        attemptJson.addProperty("usableStabilityStateCount", roots == null ? 0 : roots.getUsableStabilityStates());
        attemptJson.addProperty("minimumTangentPlaneDistance",
            roots == null ? Double.NaN : roots.getMinimumTangentPlaneDistance());
        attemptJson.addProperty("maximumTangentPlaneDistance",
            roots == null ? Double.NaN : roots.getMaximumTangentPlaneDistance());
        attemptJson.addProperty("failure", roots == null ? "NO_ROOT_SET" : roots.getFailureMessage());
        intervalAttempts.add(attemptJson);
      }
      row.add("intervalAttempts", intervalAttempts);
      output.add(row);
    }
    return output;
  }

  private static JsonObject discoveryCoverage(double[] temperaturesK, double[] pressuresBara,
      HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable,
      HydrocarbonWaterBoundaryAnchorDiscoverer.Result effectiveDiscovery) {
    double[] temperatures = sortedUnique(temperaturesK);
    double[] pressures = sortedUnique(pressuresBara);
    List<String> violations = new ArrayList<String>();
    if (temperatures.length < 8) {
      violations.add("INSUFFICIENT_TEMPERATURE_PROBES:" + temperatures.length + "/8");
    }
    if (pressures.length < 8) {
      violations.add("INSUFFICIENT_PRESSURE_PROBES:" + pressures.length + "/8");
    }
    double temperatureGapFraction = maximumLinearGapFraction(temperatures);
    double pressureGapFraction = maximumLogGapFraction(pressures);
    if (temperatureGapFraction > 0.25) {
      violations.add("TEMPERATURE_PROBE_GAP_FRACTION:" + temperatureGapFraction);
    }
    if (pressureGapFraction > 0.25) {
      violations.add("PRESSURE_PROBE_GAP_FRACTION:" + pressureGapFraction);
    }
    int requestedGridStates = temperatures.length * pressures.length;
    double successfulFraction = requestedGridStates == 0 ? 0.0
        : (double) stable.getStableRegionScan().getStates().size() / requestedGridStates;
    if (successfulFraction < 0.98) {
      violations.add("STABLE_GRID_SUCCESS_FRACTION:" + successfulFraction);
    }
    for (Branch branch : effectiveDiscovery.getBranches()) {
      if (branch.canSeedContinuation() && branch.getPoints().size() < 3) {
        violations
            .add("INSUFFICIENT_BRANCH_ANCHORS:" + branch.getIdentifier() + ":" + branch.getPoints().size() + "/3");
      }
    }
    if (stable.getCorrections().size() != stable.getStableRegionScan().getBrackets().size()) {
      violations.add("TRANSITION_CORRECTION_INVENTORY_MISMATCH");
    }

    Set<Region> observedRegions = EnumSet.noneOf(Region.class);
    for (neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterStableRegionTransitionScanner.StableState state : stable
        .getStableRegionScan().getStates()) {
      observedRegions.add(state.getRegion());
    }
    Map<BoundaryFamily, Integer> familyCounts = new EnumMap<BoundaryFamily, Integer>(BoundaryFamily.class);
    for (BoundaryFamily family : BoundaryFamily.values()) {
      familyCounts.put(family, 0);
    }
    for (neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterStableRegionTransitionScanner.TransitionBracket bracket : stable
        .getStableRegionScan().getBrackets()) {
      familyCounts.put(bracket.getFamily(), familyCounts.get(bracket.getFamily()) + 1);
    }

    JsonObject output = new JsonObject();
    output.addProperty("passed", violations.isEmpty());
    output.addProperty("temperatureProbeCount", temperatures.length);
    output.addProperty("pressureProbeCount", pressures.length);
    output.addProperty("temperatureMinimumK", temperatures[0]);
    output.addProperty("temperatureMaximumK", temperatures[temperatures.length - 1]);
    output.addProperty("pressureMinimumBara", pressures[0]);
    output.addProperty("pressureMaximumBara", pressures[pressures.length - 1]);
    output.addProperty("maximumTemperatureGapFraction", temperatureGapFraction);
    output.addProperty("maximumLogPressureGapFraction", pressureGapFraction);
    output.addProperty("requestedGridStateCount", requestedGridStates);
    output.addProperty("successfulGridStateCount", stable.getStableRegionScan().getStates().size());
    output.addProperty("successfulGridStateFraction", successfulFraction);
    output.add("violations", strings(violations));
    JsonArray regions = new JsonArray();
    for (Region region : observedRegions) {
      regions.add(region.getCode());
    }
    output.add("observedStableRegions", regions);
    JsonObject bracketCounts = new JsonObject();
    for (BoundaryFamily family : BoundaryFamily.values()) {
      bracketCounts.addProperty(family.name(), familyCounts.get(family));
    }
    output.add("transitionBracketCounts", bracketCounts);
    return output;
  }

  private static double[] sortedUnique(double[] values) {
    double[] sorted = values.clone();
    Arrays.sort(sorted);
    List<Double> unique = new ArrayList<Double>();
    for (double value : sorted) {
      if (unique.isEmpty() || Math.abs(value - unique.get(unique.size() - 1)) > 1.0e-12 * Math.max(1.0, value)) {
        unique.add(value);
      }
    }
    double[] result = new double[unique.size()];
    for (int index = 0; index < result.length; index++) {
      result[index] = unique.get(index);
    }
    return result;
  }

  private static double maximumLinearGapFraction(double[] values) {
    double span = values[values.length - 1] - values[0];
    if (!(span > 0.0)) {
      return Double.POSITIVE_INFINITY;
    }
    double maximum = 0.0;
    for (int index = 1; index < values.length; index++) {
      maximum = Math.max(maximum, (values[index] - values[index - 1]) / span);
    }
    return maximum;
  }

  private static double maximumLogGapFraction(double[] values) {
    double span = Math.log(values[values.length - 1] / values[0]);
    if (!(span > 0.0)) {
      return Double.POSITIVE_INFINITY;
    }
    double maximum = 0.0;
    for (int index = 1; index < values.length; index++) {
      maximum = Math.max(maximum, Math.log(values[index] / values[index - 1]) / span);
    }
    return maximum;
  }

  private static JsonObject network(SystemInterface template, HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery,
      int continuationPointsPerDirection, double minimumTemperatureK, double maximumTemperatureK,
      double minimumPressureBara, double maximumPressureBara) {
    JsonObject output = new JsonObject();
    HydrocarbonWaterBoundaryBranchAssembler.Result assembly;
    try {
      assembly = new HydrocarbonWaterBoundaryBranchAssembler(template).setCorrectorControls(80, 1.0e-8, 2.0e-5)
          .setStepControls(1.0, 1.0e-4, 2.0, 16, MAXIMUM_MATCHED_COMPOSITION_JUMP)
          .setDomainBounds(minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara)
          .assemble(discovery, continuationPointsPerDirection);
    } catch (RuntimeException error) {
      output.addProperty("attempted", true);
      output.addProperty("engineeringEligible", false);
      output.addProperty("branchCount", discovery.getBranches().size());
      output.addProperty("endpointCount", 0);
      output.addProperty("attachedEndpointCount", 0);
      JsonArray violations = new JsonArray();
      violations.add("BRANCH_ASSEMBLY_FAILED:" + diagnostic(error));
      output.add("violations", violations);
      output.add("branches", new JsonArray());
      output.add("endpoints", new JsonArray());
      return output;
    }

    List<EndpointEvidence> endpoints = independentlySolvedEndpoints(template, assembly, minimumTemperatureK,
        maximumTemperatureK, minimumPressureBara, maximumPressureBara);
    HydrocarbonWaterBoundaryNetworkAssembler.Result network = new HydrocarbonWaterBoundaryNetworkAssembler()
        .assemble(assembly, endpoints);
    output.addProperty("attempted", true);
    output.addProperty("engineeringEligible", network.isEngineeringEligible());
    output.addProperty("branchCount", network.getBranches().size());
    output.addProperty("endpointCount", network.getEndpoints().size());
    int attachedEndpointCount = 0;
    for (EndpointEvidence endpoint : network.getEndpoints()) {
      attachedEndpointCount += endpoint.getAttachmentCount() > 0 ? 1 : 0;
    }
    output.addProperty("attachedEndpointCount", attachedEndpointCount);
    output.addProperty("originalEndpointCandidateCount", network.getOriginalEndpointCandidateCount());
    output.add("violations", strings(network.getViolations()));
    output.add("branches", networkBranches(network.getBranches()));
    output.add("endpoints", networkEndpoints(network.getEndpoints()));
    JsonObject quality = new JsonObject();
    quality.addProperty("calculated", network.getQualityReport() != null);
    quality.addProperty("accepted", network.getQualityReport() != null && network.getQualityReport().isAccepted());
    quality.add("violations",
        network.getQualityReport() == null ? new JsonArray() : strings(network.getQualityReport().getViolations()));
    JsonArray branchReports = new JsonArray();
    if (network.getQualityReport() != null) {
      for (neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.BranchReport report : network
          .getQualityReport().getBranchReports()) {
        JsonObject branchReport = new JsonObject();
        branchReport.addProperty("identifier", report.getIdentifier());
        branchReport.addProperty("accepted", report.isAccepted());
        branchReport.addProperty("pointCount", report.getPointCount());
        branchReport.addProperty("maximumResidual", report.getMaximumResidual());
        branchReport.addProperty("maximumCompositionJump", report.getMaximumCompositionJump());
        branchReport.addProperty("minimumPhaseDistance", report.getMinimumPhaseDistance());
        branchReport.add("violations", strings(report.getViolations()));
        branchReports.add(branchReport);
      }
    }
    quality.add("branchReports", branchReports);
    output.add("qualityGate", quality);
    return output;
  }

  private static List<EndpointEvidence> independentlySolvedEndpoints(SystemInterface template,
      HydrocarbonWaterBoundaryBranchAssembler.Result assembly, double minimumTemperatureK, double maximumTemperatureK,
      double minimumPressureBara, double maximumPressureBara) {
    List<EndpointEvidence> endpoints = new ArrayList<EndpointEvidence>();
    for (HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch : assembly.getBranches()) {
      if (branch.isIsolated()) {
        continue;
      }
      if (branch.getPostRegularBackwardTrace() != null) {
        addTraceEndpoint(template, branch, branch.getPostRegularBackwardTrace(), branch.getBackwardTermination(),
            "start", endpoints, minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara);
      } else if (!addRegularTraceEndpoint(branch, branch.getRegularBackwardTrace(), "start", endpoints,
          minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara)) {
        addTraceEndpoint(template, branch, branch.getBackwardTrace(), branch.getBackwardTermination(), "start",
            endpoints, minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara);
      }
      if (branch.getPostRegularForwardTrace() != null) {
        addTraceEndpoint(template, branch, branch.getPostRegularForwardTrace(), branch.getForwardTermination(), "end",
            endpoints, minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara);
      } else if (!addRegularTraceEndpoint(branch, branch.getRegularForwardTrace(), "end", endpoints,
          minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara)) {
        addTraceEndpoint(template, branch, branch.getForwardTrace(), branch.getForwardTermination(), "end", endpoints,
            minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara);
      }
    }
    return endpoints;
  }

  private static boolean addRegularTraceEndpoint(HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch,
      HydrocarbonWaterRegularPressureBoundaryTracer.Result trace, String end, List<EndpointEvidence> endpoints,
      double minimumTemperatureK, double maximumTemperatureK, double minimumPressureBara, double maximumPressureBara) {
    if (trace == null || trace.getPoints().size() <= 2) {
      return false;
    }
    try {
      endpoints.add(EndpointEvidence.domainExit(branch.getIdentifier() + "-" + end + "-domain", branch.getFamily(),
          trace, minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara, 1.0e-5));
    } catch (RuntimeException error) {
      // A locally extended regular segment without a proven endpoint remains an open network end.
    }
    return true;
  }

  private static void addTraceEndpoint(SystemInterface template,
      HydrocarbonWaterBoundaryBranchAssembler.AssembledBranch branch,
      neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhasePseudoArcLengthTracer.Result trace,
      HydrocarbonWaterBoundaryTerminationClassifier.Result termination, String end, List<EndpointEvidence> endpoints,
      double minimumTemperatureK, double maximumTemperatureK, double minimumPressureBara, double maximumPressureBara) {
    if (trace == null || trace.getPoints().isEmpty() || termination != null && termination.isPhysicalEndpoint()) {
      return;
    }
    try {
      endpoints.add(EndpointEvidence.domainExit(branch.getIdentifier() + "-" + end + "-domain", branch.getFamily(),
          trace, minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara, 1.0e-5));
      return;
    } catch (RuntimeException error) {
      // Continue with independent critical-endpoint equations.
    }
    if (trace.hasCompletedRequestedPoints()) {
      return;
    }
    neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.Result terminalCorrection = trace
        .getTerminalCorrection();
    if (terminalCorrection == null
        || !terminalCorrection.isResidualConverged() && terminalCorrection.hasDistinctPhases()) {
      return;
    }
    neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State last = trace
        .getPoints().get(trace.getPoints().size() - 1);
    for (CandidatePhase criticalPhase : new CandidatePhase[] { branch.getFamily().getRetainedPhaseZero(),
        branch.getFamily().getRetainedPhaseOne() }) {
      HydrocarbonWaterCriticalEndpointSolver.Result critical = new HydrocarbonWaterCriticalEndpointSolver(template,
          branch.getFamily().getRetainedPhaseZero(), branch.getFamily().getRetainedPhaseOne(),
          branch.getFamily().getIncipientPhase(), criticalPhase)
          .setNumericalControls(24, 1.0e-6, 1.0e-9, 1.0e-3, 2.0e-4, 1.0e-7, 1.0e-4).solve(last);
      if (critical.isPhysicalEndpoint()) {
        endpoints.add(EndpointEvidence.criticalEndpoint(branch.getIdentifier() + "-" + end + "-cep", branch.getFamily(),
            critical));
        return;
      }
    }
  }

  private static JsonArray networkBranches(List<NetworkBranch> branches) {
    JsonArray output = new JsonArray();
    for (NetworkBranch branch : branches) {
      JsonObject row = new JsonObject();
      row.addProperty("identifier", branch.getIdentifier());
      row.addProperty("family", branch.getFamily().name());
      row.addProperty("boundaryCode", branch.getFamily().getCode());
      row.addProperty("isolated", branch.isIsolated());
      row.addProperty("generatedSeedPointCount", branch.getGeneratedSeedPointCount());
      row.addProperty("internalBridgePointCount", branch.getInternalBridgePointCount());
      row.addProperty("trimmedPointCount", branch.getTrimmedPointCount());
      row.addProperty("correctedPointCount", branch.getPoints().size());
      row.add("start", attachment(branch.getStartAttachment()));
      row.add("end", attachment(branch.getEndAttachment()));
      row.add("backwardTrace", traceDiagnostic(branch.getBackwardTrace(),
          branch.getRegularBackwardTrace() == null ? branch.getBackwardTermination() : null));
      row.add("forwardTrace", traceDiagnostic(branch.getForwardTrace(),
          branch.getRegularForwardTrace() == null ? branch.getForwardTermination() : null));
      row.add("regularPressureBackwardTrace", regularTraceDiagnostic(branch.getRegularBackwardTrace()));
      row.add("regularPressureForwardTrace", regularTraceDiagnostic(branch.getRegularForwardTrace()));
      row.add("postRegularPseudoArcBackwardTrace",
          traceDiagnostic(branch.getPostRegularBackwardTrace(), branch.getBackwardTermination()));
      row.add("postRegularPseudoArcForwardTrace",
          traceDiagnostic(branch.getPostRegularForwardTrace(), branch.getForwardTermination()));
      JsonArray points = new JsonArray();
      for (EvidencePoint point : branch.getPoints()) {
        JsonObject pointJson = new JsonObject();
        pointJson.addProperty("temperatureK", point.getTemperatureK());
        pointJson.addProperty("temperatureC", point.getTemperatureK() - 273.15);
        pointJson.addProperty("pressureBara", point.getPressureBara());
        pointJson.addProperty("beta", point.getBeta());
        pointJson.addProperty("retainedFlashResidual", point.getRetainedFlashResidual());
        pointJson.addProperty("tangentPlaneDistance", point.getTangentPlaneDistance());
        pointJson.addProperty("stationarityResidual", point.getStationarityResidual());
        pointJson.addProperty("curve", branch.getFamily().name());
        points.add(pointJson);
      }
      row.add("points", points);
      output.add(row);
    }
    return output;
  }

  private static JsonObject traceDiagnostic(
      neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhasePseudoArcLengthTracer.Result trace,
      HydrocarbonWaterBoundaryTerminationClassifier.Result termination) {
    JsonObject output = new JsonObject();
    output.addProperty("attempted", trace != null);
    if (trace == null) {
      return output;
    }
    output.addProperty("acceptedCorrectionCount", trace.getAcceptedCorrections().size());
    output.addProperty("attemptedCorrectionCount", trace.getAttemptedCorrections());
    output.addProperty("rejectedCorrectionCount", trace.getRejectedCorrections());
    output.addProperty("localTangentRescueAttemptCount", trace.getLocalTangentRescueAttempts());
    output.addProperty("acceptedLocalTangentRescueCount", trace.getAcceptedLocalTangentRescues());
    output.addProperty("completedRequestedPoints", trace.hasCompletedRequestedPoints());
    output.addProperty("finalArcStep", trace.getFinalArcStep());
    output.addProperty("terminationReason",
        trace.getTerminationReason() == null ? null : trace.getTerminationReason().name());
    output.addProperty("failure", trace.getFailureMessage());
    output.addProperty("suggestedRestartArcStep", trace.getFinalArcStep());
    JsonArray restartStates = new JsonArray();
    int restartStart = Math.max(0, trace.getPoints().size() - 2);
    for (int index = restartStart; index < trace.getPoints().size(); index++) {
      restartStates.add(continuationState(trace.getPoints().get(index)));
    }
    output.add("restartStates", restartStates);
    output.addProperty("physicalTermination", termination != null && termination.isPhysicalEndpoint());
    output.addProperty("terminationType", termination == null ? null : termination.getType().name());
    output.addProperty("terminationDiagnostic", termination == null ? null : termination.getDiagnostic());
    if (trace.getTerminalCorrection() != null) {
      neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.Result correction = trace
          .getTerminalCorrection();
      JsonObject correctionJson = new JsonObject();
      correctionJson.addProperty("converged", correction.isConverged());
      correctionJson.addProperty("residualConverged", correction.isResidualConverged());
      correctionJson.addProperty("distinctPhases", correction.hasDistinctPhases());
      correctionJson.addProperty("initialMaximumResidual", correction.getInitialMaximumResidual());
      correctionJson.addProperty("maximumResidual", correction.getMaximumResidual());
      correctionJson.addProperty("thermodynamicMaximumResidual", correction.getThermodynamicMaximumResidual());
      correctionJson.addProperty("arcLengthResidual", correction.getArcLengthResidual());
      if (correction.getEquationResiduals() != null) {
        correctionJson.add("equationResiduals", equationResidualsJson(correction.getEquationResiduals()));
      }
      if (correction.getIndependentReplayEquationResiduals() != null) {
        correctionJson.add("independentReplayEquationResiduals",
            equationResidualsJson(correction.getIndependentReplayEquationResiduals()));
      }
      if (correction.getReturnedStateReplayEquationResiduals() != null) {
        correctionJson.add("returnedStateReplayEquationResiduals",
            equationResidualsJson(correction.getReturnedStateReplayEquationResiduals()));
      }
      if (correction.getSolveFingerprint() != null) {
        correctionJson.add("solveFingerprint", stateFingerprintJson(correction.getSolveFingerprint()));
      }
      if (correction.getReturnedStateReplayFingerprint() != null) {
        correctionJson.add("returnedStateReplayFingerprint",
            stateFingerprintJson(correction.getReturnedStateReplayFingerprint()));
      }
      if (correction.getIndependentReplayFingerprint() != null) {
        correctionJson.add("independentReplayFingerprint",
            stateFingerprintJson(correction.getIndependentReplayFingerprint()));
      }
      correctionJson.addProperty("independentlyReproducible", correction.isIndependentlyReproducible());
      correctionJson.addProperty("arcStep", correction.getArcStep());
      correctionJson.addProperty("iterations", correction.getIterations());
      correctionJson.addProperty("jacobianConditionNumber", correction.getJacobianConditionNumber());
      correctionJson.addProperty("failure", correction.getFailureMessage());
      if (trace.getPoints().size() >= 2 && correction.getState().getIncipientComposition().length > 0) {
        correctionJson.addProperty("forwardCosine",
            TwoToThreePhasePseudoArcLengthTracer.forwardCosine(trace.getPoints().get(trace.getPoints().size() - 2),
                trace.getPoints().get(trace.getPoints().size() - 1), correction.getState()));
      }
      output.add("terminalCorrection", correctionJson);
    }
    if (!trace.getGlobalStabilityEvidence().isEmpty()) {
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result terminalEvidence = trace.getGlobalStabilityEvidence()
          .get(trace.getGlobalStabilityEvidence().size() - 1);
      JsonObject stability = new JsonObject();
      stability.addProperty("accepted", terminalEvidence.isAccepted());
      stability.addProperty("failure", terminalEvidence.getFailureMessage());
      stability.addProperty("temperatureK", terminalEvidence.getBoundaryRoot().getTemperatureK());
      stability.addProperty("pressureBara", terminalEvidence.getBoundaryRoot().getPressureBara());
      stability.addProperty("minimumNonTrivialTangentPlaneDistance",
          terminalEvidence.getMinimumNonTrivialTangentPlaneDistance());
      stability.addProperty("targetCompositionDistance", terminalEvidence.getTargetCompositionDistance());
      stability.addProperty("targetMatched", terminalEvidence.isTargetMatched());
      stability.addProperty("replayedBoundaryResidual", terminalEvidence.getReplayedBoundaryResidual());
      if (terminalEvidence.getReplayedEquationResiduals() != null) {
        stability.add("replayedEquationResiduals",
            equationResidualsJson(terminalEvidence.getReplayedEquationResiduals()));
      }
      if (terminalEvidence.getRetainedFlash() != null) {
        stability.addProperty("retainedFlashMaximumResidual", terminalEvidence.getRetainedFlash().getMaximumResidual());
      }
      Candidate target = terminalEvidence.getTargetCandidate();
      if (target != null) {
        stability.addProperty("targetTangentPlaneDistance", target.getTangentPlaneDistance());
        stability.addProperty("targetStationarityResidual", target.getStationarityResidual());
        stability.addProperty("targetIterations", target.getIterations());
      }
      SeedEvaluation targetSeed = terminalEvidence.getTargetSeedEvaluation();
      if (targetSeed != null) {
        stability.addProperty("targetSeedTangentPlaneDistance", targetSeed.getTangentPlaneDistance());
        stability.addProperty("targetSeedMaximumBoundaryResidual", targetSeed.getMaximumBoundaryResidual());
        stability.addProperty("targetSeedMaximumStationarityDeparture", targetSeed.getMaximumStationarityDeparture());
        stability.addProperty("maximumReferenceLogFugacitySpread", targetSeed.getMaximumReferenceLogFugacitySpread());
      }
      SeedEvaluation correctorSlotSeed = terminalEvidence.getTargetCorrectorSlotSeedEvaluation();
      if (correctorSlotSeed != null) {
        stability.addProperty("targetCorrectorSlotSeedTangentPlaneDistance",
            correctorSlotSeed.getTangentPlaneDistance());
        stability.addProperty("targetCorrectorSlotSeedMaximumBoundaryResidual",
            correctorSlotSeed.getMaximumBoundaryResidual());
        stability.addProperty("targetCorrectorSlotSeedMaximumStationarityDeparture",
            correctorSlotSeed.getMaximumStationarityDeparture());
      }
      stability.addProperty("trialCount", terminalEvidence.getTrials().size());
      stability.addProperty("seededCandidateCount", terminalEvidence.getSeededCandidates().size());
      stability.add("rejectedBoundaryState", continuationState(State.from(terminalEvidence.getBoundaryRoot())));
      Candidate mostUnstable = terminalEvidence.getMostUnstableCandidate();
      if (mostUnstable != null) {
        JsonObject candidate = new JsonObject();
        candidate.addProperty("seedPhase", mostUnstable.getSeedPhase().name());
        candidate.addProperty("phase", mostUnstable.getPhase().name());
        candidate.addProperty("tangentPlaneDistance", mostUnstable.getTangentPlaneDistance());
        candidate.addProperty("stationarityResidual", mostUnstable.getStationarityResidual());
        candidate.addProperty("stabilityFunction", mostUnstable.getStabilityFunction());
        candidate.addProperty("iterations", mostUnstable.getIterations());
        candidate.addProperty("converged", mostUnstable.isConverged());
        candidate.addProperty("trivial", mostUnstable.isTrivial());
        candidate.add("composition", numbers(mostUnstable.getComposition()));
        stability.add("mostUnstableCandidate", candidate);
      }
      JsonArray retainedPhases = new JsonArray();
      for (HydrocarbonWaterBoundaryGlobalStabilityGate.RetainedPhaseLocalStability retained : terminalEvidence
          .getRetainedPhaseStability()) {
        JsonObject retainedJson = new JsonObject();
        retainedJson.addProperty("phase", retained.getPhase().name());
        retainedJson.addProperty("accepted", retained.isAccepted());
        retainedJson.addProperty("failure", retained.getFailureMessage());
        IncipientPhaseStationarityJacobianAnalyzer.Result curvature = retained.getCurvature();
        if (curvature != null) {
          retainedJson.addProperty("finiteEvaluation", curvature.isFiniteEvaluation());
          retainedJson.addProperty("hasRealBifurcationMode", curvature.hasRealBifurcationMode());
          retainedJson.addProperty("bifurcationEigenvalue", curvature.getBifurcationEigenvalue());
          retainedJson.addProperty("symmetricMinimumEigenvalue", curvature.getMinimumEigenvalue());
          retainedJson.addProperty("minimumSingularValue", curvature.getMinimumSingularValue());
          retainedJson.addProperty("maximumResidual", curvature.getMaximumResidual());
          retainedJson.addProperty("maximumAntisymmetry", curvature.getMaximumAntisymmetry());
          retainedJson.addProperty("finiteDifferenceStep", curvature.getFiniteDifferenceStep());
        }
        retainedPhases.add(retainedJson);
      }
      stability.add("retainedPhases", retainedPhases);
      output.add("terminalGlobalStability", stability);
    }
    return output;
  }

  private static JsonObject equationResidualsJson(TwoToThreePhaseArcLengthCorrector.EquationResiduals residuals) {
    JsonObject output = new JsonObject();
    output.addProperty("retainedFugacityMaximumResidual", residuals.getRetainedFugacityMaximumResidual());
    output.addProperty("retainedMaterialBalanceResidual", residuals.getRetainedMaterialBalanceResidual());
    output.addProperty("incipientFugacityMaximumResidual", residuals.getIncipientFugacityMaximumResidual());
    output.addProperty("incipientNormalizationResidual", residuals.getIncipientNormalizationResidual());
    output.addProperty("thermodynamicMaximumResidual", residuals.getThermodynamicMaximumResidual());
    output.addProperty("retainedFugacityWorstComponent", residuals.getRetainedFugacityWorstComponent());
    output.addProperty("incipientFugacityWorstComponent", residuals.getIncipientFugacityWorstComponent());
    return output;
  }

  private static JsonObject stateFingerprintJson(TwoToThreePhaseArcLengthCorrector.StateFingerprint fingerprint) {
    JsonObject output = new JsonObject();
    output.addProperty("componentCount", fingerprint.getComponentCount());
    JsonArray names = new JsonArray();
    String[] componentNames = fingerprint.getComponentNames();
    for (int index = 0; index < componentNames.length; index++) {
      names.add(componentNames[index]);
    }
    output.add("componentNames", names);
    output.add("overallComposition", numbers(fingerprint.getOverallComposition()));
    output.addProperty("overallCompositionSum", fingerprint.getOverallCompositionSum());
    output.addProperty("beta", fingerprint.getBeta());
    output.addProperty("temperatureK", fingerprint.getTemperatureK());
    output.addProperty("pressureBara", fingerprint.getPressureBara());
    output.add("equilibriumRatios", numbers(fingerprint.getEquilibriumRatios()));
    output.addProperty("equilibriumRatioDefinition", fingerprint.getEquilibriumRatioDefinition());
    output.add("unnormalizedPhaseZeroComposition", numbers(fingerprint.getUnnormalizedPhaseZeroComposition()));
    output.add("unnormalizedPhaseOneComposition", numbers(fingerprint.getUnnormalizedPhaseOneComposition()));
    output.add("incipientComposition", numbers(fingerprint.getIncipientComposition()));
    output.add("materialBalanceByComponent", numbers(fingerprint.getMaterialBalanceByComponent()));
    output.addProperty("rachfordRiceResidual", fingerprint.getRachfordRiceResidual());
    JsonArray identities = new JsonArray();
    String[] phaseModelIdentities = fingerprint.getPhaseModelIdentities();
    for (int index = 0; index < phaseModelIdentities.length; index++) {
      identities.add(phaseModelIdentities[index]);
    }
    output.add("phaseModelIdentities", identities);
    return output;
  }

  private static JsonObject continuationState(State state) {
    JsonObject output = new JsonObject();
    output.addProperty("retainedPhaseZero", state.getRetainedPhaseZero().name());
    output.addProperty("retainedPhaseOne", state.getRetainedPhaseOne().name());
    output.addProperty("incipientPhase", state.getIncipientPhase().name());
    output.addProperty("temperatureK", state.getTemperatureK());
    output.addProperty("pressureBara", state.getPressureBara());
    output.addProperty("beta", state.getBeta());
    if (state.getSolverVariables() != null) {
      // Carried so a restart from this file re-enters the solver on the exact variables, instead of reconstructing
      // ln(K) from compositions that have already been truncated at the EOS composition floor.
      output.add("solverVariables", numbers(state.getSolverVariables()));
    }
    output.add("phaseZeroComposition", numbers(state.getPhaseZeroComposition()));
    output.add("phaseOneComposition", numbers(state.getPhaseOneComposition()));
    output.add("incipientComposition", numbers(state.getIncipientComposition()));
    return output;
  }

  private static JsonArray numbers(double[] values) {
    JsonArray output = new JsonArray();
    for (double value : values) {
      output.add(value);
    }
    return output;
  }

  private static JsonObject regularTraceDiagnostic(HydrocarbonWaterRegularPressureBoundaryTracer.Result trace) {
    JsonObject output = new JsonObject();
    output.addProperty("attempted", trace != null);
    if (trace == null) {
      return output;
    }
    output.addProperty("acceptedPointCount", Math.max(0, trace.getPoints().size() - 2));
    output.addProperty("attemptCount", trace.getAttempts().size());
    output.addProperty("completedRequestedPoints", trace.hasCompletedRequestedPoints());
    output.addProperty("terminationReason",
        trace.getTerminationReason() == null ? null : trace.getTerminationReason().name());
    output.addProperty("finalLogPressureStep", trace.getFinalLogPressureStep());
    output.addProperty("failure", trace.getFailureMessage());
    output.add("attempts", regularPressureAttempts(trace.getAttempts()));
    return output;
  }

  private static JsonObject attachment(EndpointAttachment attachment) {
    JsonObject output = new JsonObject();
    output.addProperty("attached", attachment != null);
    if (attachment != null) {
      output.addProperty("endpointIdentifier", attachment.getEndpointIdentifier());
      output.addProperty("type", attachment.getType().name());
      output.addProperty("temperatureK", attachment.getTemperatureK());
      output.addProperty("temperatureC", attachment.getTemperatureK() - 273.15);
      output.addProperty("pressureBara", attachment.getPressureBara());
      output.addProperty("evidenceSource", attachment.getSource().name());
      output.addProperty("transformedDistance", attachment.getTransformedDistance());
    }
    return output;
  }

  private static JsonArray networkEndpoints(List<EndpointEvidence> endpoints) {
    JsonArray output = new JsonArray();
    for (EndpointEvidence endpoint : endpoints) {
      JsonObject row = new JsonObject();
      row.addProperty("identifier", endpoint.getIdentifier());
      row.addProperty("type", endpoint.getType().name());
      row.addProperty("temperatureK", endpoint.getTemperatureK());
      row.addProperty("temperatureC", endpoint.getTemperatureK() - 273.15);
      row.addProperty("pressureBara", endpoint.getPressureBara());
      row.addProperty("evidenceSource", endpoint.getSource().name());
      row.addProperty("requiredDegree", endpoint.getRequiredDegree());
      row.addProperty("attachmentCount", endpoint.getAttachmentCount());
      row.addProperty("requestedBranchEnd",
          endpoint.getRequestedBranchEnd() == null ? null : endpoint.getRequestedBranchEnd().name());
      JsonArray families = new JsonArray();
      for (BoundaryFamily family : endpoint.getCompatibleFamilies()) {
        families.add(family.name());
      }
      row.add("compatibleFamilies", families);
      output.add(row);
    }
    return output;
  }

  private static JsonArray strings(List<String> values) {
    JsonArray output = new JsonArray();
    for (String value : values) {
      output.add(value);
    }
    return output;
  }

  private static String diagnostic(RuntimeException error) {
    return error.getMessage() == null ? error.getClass().getSimpleName() : error.getMessage();
  }

  private static JsonArray branches(HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery) {
    JsonArray output = new JsonArray();
    for (Branch branch : discovery.getBranches()) {
      JsonObject branchJson = new JsonObject();
      branchJson.addProperty("identifier", branch.getIdentifier());
      branchJson.addProperty("family", branch.getFamily().name());
      branchJson.addProperty("boundaryCode", branch.getFamily().getCode());
      branchJson.addProperty("continuationSeedEligible", branch.canSeedContinuation());
      JsonArray points = new JsonArray();
      for (AnchorPoint anchor : branch.getPoints()) {
        JsonObject point = new JsonObject();
        point.addProperty("temperatureK", anchor.getTemperatureK());
        point.addProperty("temperatureC", anchor.getTemperatureK() - 273.15);
        point.addProperty("pressureBara", anchor.getPressureBara());
        point.addProperty("beta", anchor.getBeta());
        point.addProperty("retainedFlashResidual", anchor.getRetainedFlashResidual());
        point.addProperty("tangentPlaneDistance", anchor.getTangentPlaneDistance());
        point.addProperty("stationarityResidual", anchor.getStationarityResidual());
        points.add(point);
      }
      branchJson.add("anchors", points);
      output.add(branchJson);
    }
    return output;
  }

  private static List<BoundaryEvidence> evidence(
      HydrocarbonWaterBoundaryAnchorDiscoverer.StableDiscoveryResult stable) {
    List<BoundaryEvidence> output = new ArrayList<BoundaryEvidence>();
    for (HydrocarbonWaterStableRegionBoundaryCorrector.Result correction : stable.getCorrections()) {
      BoundaryFamily family = correction.getBracket().getFamily();
      for (HydrocarbonWaterBoundaryEndpointClassifier.Result classification : correction.getAllClassifications()) {
        if (classification.getGlobalStabilityResult() != null) {
          output.add(new BoundaryEvidence(family, classification.getBoundaryRoot(),
              classification.getGlobalStabilityResult()));
        }
      }
    }
    output.sort(Comparator.comparing((BoundaryEvidence item) -> item.family.name())
        .thenComparingDouble(item -> item.root.getPressureBara())
        .thenComparingDouble(item -> item.root.getTemperatureK()));
    return output;
  }

  private static JsonArray terminations(SystemInterface template, List<BoundaryEvidence> evidence) {
    JsonArray output = new JsonArray();
    HydrocarbonWaterBoundaryTerminationClassifier classifier = new HydrocarbonWaterBoundaryTerminationClassifier(
        template).setNumericalControls(40, 80, 1.0e-8, 2.0e-5, 2.0e-4, 1.0e-4, 1.0e-6);
    for (BoundaryEvidence rejected : evidence) {
      if (rejected.stability.isAccepted() || mostUnstableRetainedPhase(rejected.stability) == null) {
        continue;
      }
      BoundaryEvidence accepted = nearestAccepted(evidence, rejected);
      if (accepted == null) {
        output.add(unresolvedTermination(rejected, "no composition-continuous accepted boundary root was found"));
        continue;
      }
      HydrocarbonWaterBoundaryTerminationClassifier.Result termination;
      try {
        termination = classifier.classify(accepted.root, rejected.stability);
      } catch (RuntimeException error) {
        output.add(unresolvedTermination(rejected, error.getMessage()));
        continue;
      }
      JsonObject row = termination(accepted, rejected, termination);
      if (!containsEquivalentTermination(output, row)) {
        output.add(row);
      }
    }
    return output;
  }

  private static BoundaryEvidence nearestAccepted(List<BoundaryEvidence> evidence, BoundaryEvidence rejected) {
    BoundaryEvidence selected = null;
    double selectedScore = Double.POSITIVE_INFINITY;
    for (BoundaryEvidence candidate : evidence) {
      if (candidate.family != rejected.family || !candidate.stability.isAccepted()) {
        continue;
      }
      double compositionJump = compositionJump(candidate.root, rejected.root);
      if (compositionJump > MAXIMUM_MATCHED_COMPOSITION_JUMP) {
        continue;
      }
      double pressureDistance = Math.abs(Math.log(candidate.root.getPressureBara() / rejected.root.getPressureBara()));
      double score = compositionJump + 0.1 * pressureDistance;
      if (score < selectedScore) {
        selected = candidate;
        selectedScore = score;
      }
    }
    return selected;
  }

  private static CandidatePhase mostUnstableRetainedPhase(
      HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability) {
    CandidatePhase selected = null;
    double minimumEigenvalue = 0.0;
    for (HydrocarbonWaterBoundaryGlobalStabilityGate.RetainedPhaseLocalStability retained : stability
        .getRetainedPhaseStability()) {
      if (!retained.isAccepted() && retained.getMinimumEigenvalue() < minimumEigenvalue) {
        minimumEigenvalue = retained.getMinimumEigenvalue();
        selected = retained.getPhase();
      }
    }
    return selected;
  }

  private static JsonObject termination(BoundaryEvidence accepted, BoundaryEvidence rejected,
      HydrocarbonWaterBoundaryTerminationClassifier.Result termination) {
    JsonObject row = new JsonObject();
    row.addProperty("family", rejected.family.name());
    row.addProperty("boundaryCode", rejected.family.getCode());
    row.addProperty("type", termination.getType().name());
    row.addProperty("physicalEndpoint", termination.isPhysicalEndpoint());
    row.addProperty("temperatureK", termination.getTemperatureK());
    row.addProperty("temperatureC", termination.getTemperatureK() - 273.15);
    row.addProperty("pressureBara", termination.getPressureBara());
    row.addProperty("destabilizingPhase",
        termination.getDestabilizingPhase() == null ? null : termination.getDestabilizingPhase().name());
    row.addProperty("qualityMeasure", termination.getQualityMeasure());
    row.addProperty("acceptedBracketPressureBara", accepted.root.getPressureBara());
    row.addProperty("rejectedBracketPressureBara", rejected.root.getPressureBara());
    row.addProperty("diagnostic", termination.getDiagnostic());
    if (termination.getSpinodalResult() != null) {
      row.addProperty("pressureBracketWidthBara", termination.getSpinodalResult().getPressureBracketWidthBara());
      row.addProperty("refinementIterations", termination.getSpinodalResult().getIterations());
    }
    return row;
  }

  private static JsonObject unresolvedTermination(BoundaryEvidence rejected, String diagnostic) {
    JsonObject row = new JsonObject();
    row.addProperty("family", rejected.family.name());
    row.addProperty("boundaryCode", rejected.family.getCode());
    row.addProperty("type", "UNRESOLVED_GLOBAL_STABILITY_LIMIT");
    row.addProperty("physicalEndpoint", false);
    row.addProperty("temperatureK", rejected.root.getTemperatureK());
    row.addProperty("temperatureC", rejected.root.getTemperatureK() - 273.15);
    row.addProperty("pressureBara", rejected.root.getPressureBara());
    CandidatePhase unstable = mostUnstableRetainedPhase(rejected.stability);
    row.addProperty("destabilizingPhase", unstable == null ? null : unstable.name());
    row.addProperty("diagnostic", diagnostic == null ? "unresolved global-stability limit" : diagnostic);
    return row;
  }

  private static boolean containsEquivalentTermination(JsonArray rows, JsonObject candidate) {
    for (int index = 0; index < rows.size(); index++) {
      JsonObject existing = rows.get(index).getAsJsonObject();
      if (existing.get("family").getAsString().equals(candidate.get("family").getAsString())
          && existing.get("type").getAsString().equals(candidate.get("type").getAsString()) && Math
              .abs(existing.get("pressureBara").getAsDouble() - candidate.get("pressureBara").getAsDouble()) < 1.0e-3) {
        return true;
      }
    }
    return false;
  }

  private static double compositionJump(TwoToThreePhaseBoundaryPointSolver.Result first,
      TwoToThreePhaseBoundaryPointSolver.Result second) {
    return Math.max(distance(first.getPhaseZeroComposition(), second.getPhaseZeroComposition()),
        Math.max(distance(first.getPhaseOneComposition(), second.getPhaseOneComposition()),
            distance(first.getIncipientComposition(), second.getIncipientComposition())));
  }

  private static double stateCompositionJump(State first, State second) {
    return Math.max(distance(first.getPhaseZeroComposition(), second.getPhaseZeroComposition()),
        Math.max(distance(first.getPhaseOneComposition(), second.getPhaseOneComposition()),
            distance(first.getIncipientComposition(), second.getIncipientComposition())));
  }

  private static double distance(double[] first, double[] second) {
    if (first.length != second.length) {
      return Double.POSITIVE_INFINITY;
    }
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }

  private static final class IsolatedTopologyAudit {
    private final JsonArray json;
    private final List<AnchorPoint> independentAnchors;

    private IsolatedTopologyAudit(JsonArray json, List<AnchorPoint> independentAnchors) {
      this.json = json;
      this.independentAnchors = new ArrayList<AnchorPoint>(independentAnchors);
    }
  }

  private static final class BoundaryEvidence {
    private final BoundaryFamily family;
    private final TwoToThreePhaseBoundaryPointSolver.Result root;
    private final HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability;

    private BoundaryEvidence(BoundaryFamily family, TwoToThreePhaseBoundaryPointSolver.Result root,
        HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability) {
      this.family = family;
      this.root = root;
      this.stability = stability;
    }
  }
}
