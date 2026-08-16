package neqsim.mcp.runners;

import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import neqsim.thermo.mixingrule.HVMixingRulesInterface;
import neqsim.process.util.monitor.FluidResponse;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemElectrolyteCPAstatoil;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkCPAstatoil;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermo.util.readwrite.JsonFluidReadWrite;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.flashops.TPflash;
import neqsim.thermodynamicoperations.flashops.reactiveflash.ReactiveMultiphaseTPflash;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterContinuationChainAssembler;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterContinuationChainAssembler.Segment;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseBoundaryPointSolver;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoHydrocarbonPhaseEnvelopeSolver;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseArcLengthCorrector.State;

/**
 * Runner for the field-fluid validation pages that require electrolyte chemistry or a characterized heavy-oil
 * sub-fluid.
 *
 * <p>
 * The input composition defines one mole of the final field fluid. When a PVTsim fluid definition is supplied, its
 * normalized composition is treated as the heavy-oil sub-fluid and is scaled by the requested {@code heavyOil}
 * fraction. Explicit water, carbon dioxide, nitrogen, and ammonia fractions are then blended into the same
 * thermodynamic system.
 * </p>
 */
public final class FieldFluidRunner {
  private static final Gson GSON = new GsonBuilder().setPrettyPrinting().serializeSpecialFloatingPointValues().create();
  private static final double TRACE_MOLES = 1.0e-20;
  private static final double MIN_CONTINUOUS_PRESSURE_RATIO = 0.25;
  private static final double MAX_ENGINEERING_ENVELOPE_PRESSURE_BARA = 10000.0;
  private static final double DEFAULT_PHASE_MAP_MIN_TEMPERATURE_C = -60.0;
  private static final double DEFAULT_PHASE_MAP_MAX_TEMPERATURE_C = 500.0;
  private static final double DEFAULT_PHASE_MAP_MIN_PRESSURE_BARA = 0.05;
  private static final double DEFAULT_PHASE_MAP_MAX_PRESSURE_BARA = 500.0;
  private static final int DEFAULT_PHASE_MAP_TEMPERATURE_POINTS = 41;
  private static final int DEFAULT_PHASE_MAP_PRESSURE_POINTS = 25;
  private static final int PHASE_MAP_BISECTION_ITERATIONS = 8;
  private static final String[] REACTIVE_SPECIES = new String[] { "NH4+", "HCO3-", "CO3--", "H3O+", "OH-" };
  /**
   * Highest pressure the envelope continuation may be extended to when a branch is still climbing at the requested
   * ceiling. The phase-map pressure range is a plotting window and must not silently truncate the envelope, but a cubic
   * EOS should not be trusted for engineering phase behaviour far beyond this.
   */
  private static final double ENVELOPE_PRESSURE_CEILING_LIMIT_BARA = 2000.0;
  /** Pressures sampled along the hydrocarbon dew branch of a water-containing feed. */
  private static final int DEW_BRANCH_PRESSURE_POINTS = 28;
  private static final double DEW_BRANCH_MINIMUM_TEMPERATURE_K = 173.15;
  private static final double DEW_BRANCH_MAXIMUM_TEMPERATURE_K = 823.15;
  /** Species carried by the aqueous phase, excluded from the hydrocarbon phase-envelope basis. */
  private static final java.util.Set<String> AQUEOUS_BASIS_SPECIES = new java.util.HashSet<String>(
      java.util.Arrays.asList("water", "ammonia"));

  private FieldFluidRunner() {
  }

  /**
   * Runs one state point and a definition-compatible phase boundary. Characterized non-reactive heavy fluids use an
   * Michelsen predictor-corrector continuation with a stable-phase topology fallback; simple feeds use a dew-pressure
   * scan over the requested temperature range.
   *
   * @param json field-fluid request JSON
   * @return JSON result containing state properties, speciation, and calculated boundary points
   */
  public static String run(String json) {
    long started = System.currentTimeMillis();
    try {
      JsonObject input = JsonParser.parseString(json).getAsJsonObject();
      double temperatureC = requiredNumber(input, "temperatureC");
      double pressureBara = requiredPositiveNumber(input, "pressureBara");
      double minimumC = optionalNumber(input, "temperatureMinC", temperatureC);
      double maximumC = optionalNumber(input, "temperatureMaxC", temperatureC);
      int pointCount = Math.max(11, Math.min(optionalInteger(input, "pointCount", 31), 81));
      boolean reactive = input.has("reactive") && input.get("reactive").getAsBoolean();
      String requestedModel = input.has("eos") ? input.get("eos").getAsString()
          : input.has("model") ? input.get("model").getAsString() : reactive ? "Electrolyte-CPA" : "SRK";
      if (maximumC < minimumC) {
        throw new IllegalArgumentException("temperatureMaxC must be greater than temperatureMinC");
      }

      JsonObject composition = input.getAsJsonObject("components");
      if (composition == null || composition.size() == 0) {
        throw new IllegalArgumentException("components must be a non-empty object");
      }
      JsonObject fluidDefinition = input.has("fluidDefinition") && input.get("fluidDefinition").isJsonObject()
          ? input.getAsJsonObject("fluidDefinition")
          : null;
      String importedEos = fluidDefinition == null ? null
          : fluidDefinition.has("eos") ? fluidDefinition.get("eos").getAsString().toUpperCase(Locale.ROOT) : "SRK";
      String sourceEos = fluidDefinition == null ? null
          : fluidDefinition.has("sourceEos") ? fluidDefinition.get("sourceEos").getAsString().toUpperCase(Locale.ROOT)
              : importedEos;
      boolean modelTransferredFromSource = fluidDefinition != null && !reactive && !sourceEos.equals(importedEos);

      SystemInterface state = buildSystem(composition, fluidDefinition, requestedModel, reactive, temperatureC + 273.15,
          pressureBara);
      ReactiveFlashAudit reactiveAudit = runTpFlash(state, reactive);
      state.initProperties();

      JsonArray boundary = new JsonArray();
      JsonArray phaseMap = new JsonArray();
      List<String> boundaryWarnings = new ArrayList<String>();
      int rejectedBoundaryPoints = 0;
      String boundaryMethod = "not-requested";
      String boundaryTopology = "not-calculated";
      boolean fullEnvelopeCalculated = false;
      boolean boundaryDefinitionCompatible = fluidDefinition == null;
      double aqueousBasisFraction = 0.0;
      JsonObject hydrocarbonDewBranchQuality = null;
      PhaseMapAudit phaseMapAudit = null;
      TwoHydrocarbonPhaseEnvelopeSolver.Result twoHydrocarbonEnvelope = null;
      if (maximumC > minimumC && fluidDefinition != null && !reactive) {
        boundaryMethod = "two-hydrocarbon-phase-envelope";
        boundaryTopology = "hydrocarbon-one-to-two-phase";
        double mapMinimumC = Math.min(minimumC,
            optionalNumber(input, "phaseMapTemperatureMinC", DEFAULT_PHASE_MAP_MIN_TEMPERATURE_C));
        double mapMaximumC = Math.max(maximumC,
            optionalNumber(input, "phaseMapTemperatureMaxC", DEFAULT_PHASE_MAP_MAX_TEMPERATURE_C));
        double mapMinimumPressure = optionalNumber(input, "phaseMapPressureMinBara",
            DEFAULT_PHASE_MAP_MIN_PRESSURE_BARA);
        double mapMaximumPressure = optionalNumber(input, "phaseMapPressureMaxBara",
            DEFAULT_PHASE_MAP_MAX_PRESSURE_BARA);
        int mapTemperaturePoints = Math.max(21,
            Math.min(optionalInteger(input, "phaseMapTemperaturePoints", DEFAULT_PHASE_MAP_TEMPERATURE_POINTS), 81));
        int mapPressurePoints = Math.max(13,
            Math.min(optionalInteger(input, "phaseMapPressurePoints", DEFAULT_PHASE_MAP_PRESSURE_POINTS), 41));
        if (!(mapMinimumPressure > 0.0) || mapMaximumPressure <= mapMinimumPressure) {
          throw new IllegalArgumentException("phase-envelope pressure range must be positive and increasing");
        }
        aqueousBasisFraction = aqueousFraction(composition);
        JsonObject envelopeComposition = aqueousBasisFraction > 0.0 ? hydrocarbonBasisComposition(composition)
            : composition;
        SystemInterface envelopeSystem = buildBoundarySystem(envelopeComposition, fluidDefinition, requestedModel,
            temperatureC + 273.15, pressureBara);
        boolean huronVidalEnvelope = fluidDefinition.has("polarModel")
            && "HV".equalsIgnoreCase(fluidDefinition.get("polarModel").getAsString());
        boolean allowExpensiveTopologyFallback = !huronVidalEnvelope || (input.has("allowExpensiveTopologyFallback")
            && input.get("allowExpensiveTopologyFallback").getAsBoolean());
        twoHydrocarbonEnvelope = new TwoHydrocarbonPhaseEnvelopeSolver(envelopeSystem)
            .setTemperatureRange(mapMinimumC + 273.15, mapMaximumC + 273.15)
            .setPressureRange(mapMinimumPressure, mapMaximumPressure)
            .setTopologySeedGrid(mapTemperaturePoints, mapPressurePoints).setMaximumSteps(10.0, 10.0)
            .setMaximumContinuationIterations(800).setTopologyFallbackEnabled(allowExpensiveTopologyFallback)
            .setPressureCeilingExtensionLimit(ENVELOPE_PRESSURE_CEILING_LIMIT_BARA).setStabilityAnalysisEnabled(false)
            .solve();
        boundary = twoHydrocarbonEnvelopeRows(twoHydrocarbonEnvelope);
        rejectedBoundaryPoints = twoHydrocarbonEnvelope.getFailedFlashCount();
        boundaryMethod = twoHydrocarbonEnvelope.getMethod();
        fullEnvelopeCalculated = twoHydrocarbonEnvelope.isEnvelopeClosed();
        boundaryDefinitionCompatible = twoHydrocarbonEnvelope.hasFinitePhysicalValues() && boundary.size() >= 8
            && (!twoHydrocarbonEnvelope.isIterationLimitReached() || fullEnvelopeCalculated);
        if (twoHydrocarbonEnvelope.getFailureMessage() != null) {
          boundaryWarnings.add(twoHydrocarbonEnvelope.getFailureMessage());
        }
        if (twoHydrocarbonEnvelope.isIterationLimitReached() && !fullEnvelopeCalculated) {
          boundaryWarnings
              .add("Michelsen continuation reached its hard iteration limit before closing both physical branches");
        }
        if (!fullEnvelopeCalculated && twoHydrocarbonEnvelope.getTerminationSummary().contains("PRESSURE_CEILING")) {
          boundaryWarnings.add("The traced boundary ran into the maximum envelope pressure of " + mapMaximumPressure
              + " bara while still climbing, so it is a truncated piece of the envelope rather than a closed one");
        }
        if (huronVidalEnvelope && !allowExpensiveTopologyFallback && boundary.size() < 8) {
          boundaryWarnings.add(
              "Native HV continuation did not produce a complete boundary; the expensive TP-topology fallback is disabled for interactive calculation");
        }
      } else if (maximumC > minimumC && fluidDefinition != null) {
        boundaryMethod = "reactive-heavy-envelope-unavailable";
        boundaryWarnings.add(
            "A reaction-coupled full PT envelope is not available for imported heavy-oil pseudo-components; the reactive state point remains valid");
      } else if (maximumC > minimumC) {
        boundaryMethod = "temperature-scan-dew-pressure";
        boundaryTopology = "vapor-liquid";
        for (int index = 0; index < pointCount; index++) {
          double fraction = (double) index / (pointCount - 1);
          double pointTemperatureC = minimumC + fraction * (maximumC - minimumC);
          try {
            SystemInterface point = buildBoundarySystem(composition, fluidDefinition, requestedModel,
                pointTemperatureC + 273.15, 1.0);
            ThermodynamicOperations operations = new ThermodynamicOperations(point);
            operations.dewPointPressureFlash();
            double pointPressure = point.getPressure();
            if (Double.isFinite(pointPressure) && pointPressure > 0.0
                && pointPressure <= MAX_ENGINEERING_ENVELOPE_PRESSURE_BARA) {
              if (boundary.size() > 0) {
                double previousPressure = boundary.get(boundary.size() - 1).getAsJsonObject().get("pressure_bara")
                    .getAsDouble();
                if (pointPressure < previousPressure * MIN_CONTINUOUS_PRESSURE_RATIO) {
                  rejectedBoundaryPoints++;
                  boundaryWarnings
                      .add("Boundary point at " + pointTemperatureC + " C was rejected: pressure changed from "
                          + previousPressure + " bara to " + pointPressure + " bara outside the continuity gate");
                  continue;
                }
              }
              JsonObject row = new JsonObject();
              row.addProperty("temperature_C", pointTemperatureC);
              row.addProperty("temperature_K", pointTemperatureC + 273.15);
              row.addProperty("pressure_bara", pointPressure);
              row.addProperty("curve",
                  normalizedModel(requestedModel, reactive).toLowerCase(Locale.ROOT).replace('-', '_') + "-feed-dew");
              boundary.add(row);
            }
          } catch (Exception error) {
            boundaryWarnings.add("Boundary point at " + pointTemperatureC + " C failed: " + error.getMessage());
          }
        }
      }

      JsonObject hydrocarbonWaterBoundaryAudit = null;
      if (input.has("hydrocarbonWaterBoundaryAudit")) {
        if (!input.get("hydrocarbonWaterBoundaryAudit").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterBoundaryAudit must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "hydrocarbon-water boundary audit requires a non-reactive characterized fluid definition");
        }
        JsonObject auditRequest = input.getAsJsonObject("hydrocarbonWaterBoundaryAudit");
        double[] auditTemperaturesK = requiredNumberArray(auditRequest, "temperaturesC", 273.15);
        double[] auditPressuresBara = requiredNumberArray(auditRequest, "pressuresBara", 0.0);
        int continuationPointsPerDirection = Math.max(1,
            Math.min(optionalInteger(auditRequest, "continuationPointsPerDirection", 80), 400));
        SystemInterface auditSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            auditTemperaturesK[0], auditPressuresBara[0]);
        hydrocarbonWaterBoundaryAudit = HydrocarbonWaterBoundaryAuditRunner.run(auditSystem, auditTemperaturesK,
            auditPressuresBara, continuationPointsPerDirection);
      }

      JsonObject hydrocarbonWaterStableRegionTopologyScan = null;
      if (input.has("hydrocarbonWaterStableRegionTopologyScan")) {
        if (!input.get("hydrocarbonWaterStableRegionTopologyScan").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterStableRegionTopologyScan must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "stable-region topology scan requires a non-reactive characterized fluid definition");
        }
        JsonObject scanRequest = input.getAsJsonObject("hydrocarbonWaterStableRegionTopologyScan");
        double[] scanTemperaturesK = requiredNumberArray(scanRequest, "temperaturesC", 273.15);
        double[] scanPressuresBara = requiredNumberArray(scanRequest, "pressuresBara", 0.0);
        SystemInterface scanSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            scanTemperaturesK[0], scanPressuresBara[0]);
        hydrocarbonWaterStableRegionTopologyScan = HydrocarbonWaterBoundaryAuditRunner
            .runStableRegionTopologyScan(scanSystem, scanTemperaturesK, scanPressuresBara);
      }

      JsonObject hydrocarbonWaterModelSeededBoundaryDiagnostic = null;
      if (input.has("hydrocarbonWaterModelSeededBoundaryDiagnostic")) {
        if (!input.get("hydrocarbonWaterModelSeededBoundaryDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterModelSeededBoundaryDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "model-seeded boundary diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject diagnosticRequest = input.getAsJsonObject("hydrocarbonWaterModelSeededBoundaryDiagnostic");
        State seedState = continuationState(diagnosticRequest, "seedState");
        HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily family = HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily
            .valueOf(diagnosticRequest.get("family").getAsString());
        double minimumTemperatureK = optionalNumber(diagnosticRequest, "minimumTemperatureC",
            seedState.getTemperatureK() - 273.15 - 50.0) + 273.15;
        double maximumTemperatureK = optionalNumber(diagnosticRequest, "maximumTemperatureC",
            seedState.getTemperatureK() - 273.15 + 50.0) + 273.15;
        SystemInterface diagnosticSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            seedState.getTemperatureK(), seedState.getPressureBara());
        hydrocarbonWaterModelSeededBoundaryDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runModelSeededFixedPressureDiagnostic(diagnosticSystem, seedState, family, minimumTemperatureK,
                maximumTemperatureK);
      }

      JsonObject hydrocarbonWaterStableBracketCorrectionDiagnostic = null;
      if (input.has("hydrocarbonWaterStableBracketCorrectionDiagnostic")) {
        if (!input.get("hydrocarbonWaterStableBracketCorrectionDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterStableBracketCorrectionDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "stable-bracket correction requires a non-reactive characterized fluid definition");
        }
        JsonObject diagnosticRequest = input.getAsJsonObject("hydrocarbonWaterStableBracketCorrectionDiagnostic");
        double lowerTemperatureK = requiredNumber(diagnosticRequest, "lowerTemperatureC") + 273.15;
        double upperTemperatureK = requiredNumber(diagnosticRequest, "upperTemperatureC") + 273.15;
        double bracketPressureBara = requiredPositiveNumber(diagnosticRequest, "pressureBara");
        HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily family = HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily
            .valueOf(diagnosticRequest.get("family").getAsString());
        SystemInterface diagnosticSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            lowerTemperatureK, bracketPressureBara);
        hydrocarbonWaterStableBracketCorrectionDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runStableBracketCorrectionDiagnostic(diagnosticSystem, lowerTemperatureK, upperTemperatureK,
                bracketPressureBara, family);
      }

      JsonObject hydrocarbonWaterRegularPressureDiagnostic = null;
      if (input.has("hydrocarbonWaterRegularPressureDiagnostic")) {
        if (!input.get("hydrocarbonWaterRegularPressureDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterRegularPressureDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "regular-pressure boundary diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject diagnosticRequest = input.getAsJsonObject("hydrocarbonWaterRegularPressureDiagnostic");
        double[] diagnosticTemperaturesK = requiredNumberArray(diagnosticRequest, "temperaturesC", 273.15);
        double[] diagnosticPressuresBara = requiredNumberArray(diagnosticRequest, "seedPressuresBara", 0.0);
        HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily family = HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily
            .valueOf(diagnosticRequest.get("family").getAsString());
        int additionalPointCount = Math.max(1,
            Math.min(optionalInteger(diagnosticRequest, "additionalPointCount", 1), 20));
        String pressureDirection = diagnosticRequest.has("pressureDirection")
            ? diagnosticRequest.get("pressureDirection").getAsString().toUpperCase(Locale.ROOT)
            : "INCREASING";
        if (!"INCREASING".equals(pressureDirection) && !"DECREASING".equals(pressureDirection)) {
          throw new IllegalArgumentException("pressureDirection must be INCREASING or DECREASING");
        }
        SystemInterface diagnosticSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            diagnosticTemperaturesK[0], diagnosticPressuresBara[0]);
        hydrocarbonWaterRegularPressureDiagnostic = HydrocarbonWaterBoundaryAuditRunner.runRegularPressureDiagnostic(
            diagnosticSystem, diagnosticTemperaturesK, diagnosticPressuresBara, family, additionalPointCount,
            "INCREASING".equals(pressureDirection));
      }

      JsonObject hydrocarbonWaterRegularPressureRestartDiagnostic = null;
      if (input.has("hydrocarbonWaterRegularPressureRestartDiagnostic")) {
        if (!input.get("hydrocarbonWaterRegularPressureRestartDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterRegularPressureRestartDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "regular-pressure restart diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject restartRequest = input.getAsJsonObject("hydrocarbonWaterRegularPressureRestartDiagnostic");
        State previousState = continuationState(restartRequest, "previousState");
        State currentState = continuationState(restartRequest, "currentState");
        HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily family = HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily
            .valueOf(restartRequest.get("family").getAsString());
        int additionalPointCount = Math.max(1,
            Math.min(optionalInteger(restartRequest, "additionalPointCount", 5), 20));
        SystemInterface restartSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            currentState.getTemperatureK(), currentState.getPressureBara());
        hydrocarbonWaterRegularPressureRestartDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runRegularPressureRestartDiagnostic(restartSystem, previousState, currentState, family,
                additionalPointCount);
      }

      JsonObject hydrocarbonWaterHybridBranchDiagnostic = null;
      if (input.has("hydrocarbonWaterHybridBranchDiagnostic")) {
        if (!input.get("hydrocarbonWaterHybridBranchDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterHybridBranchDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "hybrid branch diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject diagnosticRequest = input.getAsJsonObject("hydrocarbonWaterHybridBranchDiagnostic");
        double[] diagnosticTemperaturesK = requiredNumberArray(diagnosticRequest, "temperaturesC", 273.15);
        double[] diagnosticPressuresBara = requiredNumberArray(diagnosticRequest, "seedPressuresBara", 0.0);
        HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily family = HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily
            .valueOf(diagnosticRequest.get("family").getAsString());
        int additionalPointsPerDirection = Math.max(1,
            Math.min(optionalInteger(diagnosticRequest, "additionalPointsPerDirection", 2), 100));
        double minimumTemperatureK = optionalNumber(diagnosticRequest, "minimumTemperatureC", -100.0) + 273.15;
        double maximumTemperatureK = optionalNumber(diagnosticRequest, "maximumTemperatureC", 700.0) + 273.15;
        double minimumPressureBara = optionalNumber(diagnosticRequest, "minimumPressureBara", 0.01);
        double maximumPressureBara = optionalNumber(diagnosticRequest, "maximumPressureBara", 1000.0);
        SystemInterface diagnosticSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            diagnosticTemperaturesK[0], diagnosticPressuresBara[0]);
        hydrocarbonWaterHybridBranchDiagnostic = HydrocarbonWaterBoundaryAuditRunner.runHybridBranchDiagnostic(
            diagnosticSystem, diagnosticTemperaturesK, diagnosticPressuresBara, family, additionalPointsPerDirection,
            minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara);
      }

      JsonObject hydrocarbonWaterPseudoArcRestartDiagnostic = null;
      if (input.has("hydrocarbonWaterPseudoArcRestartDiagnostic")) {
        if (!input.get("hydrocarbonWaterPseudoArcRestartDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterPseudoArcRestartDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "pseudo-arclength restart diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject restartRequest = input.getAsJsonObject("hydrocarbonWaterPseudoArcRestartDiagnostic");
        State previousState = continuationState(restartRequest, "previousState");
        State currentState = continuationState(restartRequest, "currentState");
        int additionalPointCount = Math.max(1,
            Math.min(optionalInteger(restartRequest, "additionalPointCount", 5), 200));
        double initialArcStep = optionalNumber(restartRequest, "initialArcStep", 0.25);
        double minimumArcStep = optionalNumber(restartRequest, "minimumArcStep", 1.0e-4);
        double maximumArcStep = optionalNumber(restartRequest, "maximumArcStep", Math.max(2.0, initialArcStep));
        double finiteDifferenceStep = optionalNumber(restartRequest, "finiteDifferenceStep", 2.0e-5);
        double minimumTemperatureK = optionalNumber(restartRequest, "minimumTemperatureC", -100.0) + 273.15;
        double maximumTemperatureK = optionalNumber(restartRequest, "maximumTemperatureC", 700.0) + 273.15;
        double minimumPressureBara = optionalNumber(restartRequest, "minimumPressureBara", 0.01);
        double maximumPressureBara = optionalNumber(restartRequest, "maximumPressureBara", 1000.0);
        SystemInterface restartSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            currentState.getTemperatureK(), currentState.getPressureBara());
        hydrocarbonWaterPseudoArcRestartDiagnostic = HydrocarbonWaterBoundaryAuditRunner.runPseudoArcRestartDiagnostic(
            restartSystem, previousState, currentState, additionalPointCount, initialArcStep, minimumArcStep,
            maximumArcStep, finiteDifferenceStep, minimumTemperatureK, maximumTemperatureK, minimumPressureBara,
            maximumPressureBara);
      }

      JsonObject hydrocarbonWaterLocalTangentProbeDiagnostic = null;
      if (input.has("hydrocarbonWaterLocalTangentProbeDiagnostic")) {
        if (!input.get("hydrocarbonWaterLocalTangentProbeDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterLocalTangentProbeDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "local tangent probe requires a non-reactive characterized fluid definition");
        }
        JsonObject probeRequest = input.getAsJsonObject("hydrocarbonWaterLocalTangentProbeDiagnostic");
        State previousState = continuationState(probeRequest, "previousState");
        State currentState = continuationState(probeRequest, "currentState");
        double[] arcSteps = requiredNumberArray(probeRequest, "arcSteps", 0.0);
        double finiteDifferenceStep = optionalNumber(probeRequest, "finiteDifferenceStep", 2.0e-5);
        SystemInterface probeSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            currentState.getTemperatureK(), currentState.getPressureBara());
        hydrocarbonWaterLocalTangentProbeDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runLocalTangentProbeDiagnostic(probeSystem, previousState, currentState, arcSteps, finiteDifferenceStep);
      }

      JsonObject hydrocarbonWaterRestartStateTopologyDiagnostic = null;
      if (input.has("hydrocarbonWaterRestartStateTopologyDiagnostic")) {
        if (!input.get("hydrocarbonWaterRestartStateTopologyDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterRestartStateTopologyDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "restart-state topology diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject topologyRequest = input.getAsJsonObject("hydrocarbonWaterRestartStateTopologyDiagnostic");
        State topologyState = continuationState(topologyRequest, "state");
        HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily family = HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily
            .valueOf(topologyRequest.get("family").getAsString());
        SystemInterface topologySystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            topologyState.getTemperatureK(), topologyState.getPressureBara());
        hydrocarbonWaterRestartStateTopologyDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runRestartStateTopologyDiagnostic(topologySystem, topologyState, family);
      }

      JsonObject hydrocarbonWaterStrictBranchSeedDiagnostic = null;
      if (input.has("hydrocarbonWaterStrictBranchSeedDiagnostic")) {
        if (!input.get("hydrocarbonWaterStrictBranchSeedDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterStrictBranchSeedDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "strict branch-seed diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject seedRequest = input.getAsJsonObject("hydrocarbonWaterStrictBranchSeedDiagnostic");
        double[] seedTemperaturesK = requiredNumberArray(seedRequest, "temperaturesC", 273.15);
        double[] seedPressuresBara = requiredNumberArray(seedRequest, "pressuresBara", 0.0);
        HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily family = HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily
            .valueOf(seedRequest.get("family").getAsString());
        SystemInterface seedSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            seedTemperaturesK[0], seedPressuresBara[0]);
        hydrocarbonWaterStrictBranchSeedDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runStrictBranchSeedDiagnostic(seedSystem, seedTemperaturesK, seedPressuresBara, family);
      }

      JsonObject hydrocarbonWaterBranchSwitchDiagnostic = null;
      if (input.has("hydrocarbonWaterBranchSwitchDiagnostic")) {
        if (!input.get("hydrocarbonWaterBranchSwitchDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterBranchSwitchDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "branch-switch diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject switchRequest = input.getAsJsonObject("hydrocarbonWaterBranchSwitchDiagnostic");
        State sourceState = continuationState(switchRequest, "sourceState");
        State targetPreviousState = continuationState(switchRequest, "targetPreviousState");
        State targetCurrentState = continuationState(switchRequest, "targetCurrentState");
        int additionalPointCount = Math.max(1,
            Math.min(optionalInteger(switchRequest, "additionalPointCount", 5), 200));
        double initialArcStep = optionalNumber(switchRequest, "initialArcStep", 0.25);
        double minimumArcStep = optionalNumber(switchRequest, "minimumArcStep", 1.0e-4);
        double finiteDifferenceStep = optionalNumber(switchRequest, "finiteDifferenceStep", 2.0e-5);
        double minimumTemperatureK = optionalNumber(switchRequest, "minimumTemperatureC", -100.0) + 273.15;
        double maximumTemperatureK = optionalNumber(switchRequest, "maximumTemperatureC", 700.0) + 273.15;
        double minimumPressureBara = optionalNumber(switchRequest, "minimumPressureBara", 0.01);
        double maximumPressureBara = optionalNumber(switchRequest, "maximumPressureBara", 1000.0);
        SystemInterface switchSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            targetCurrentState.getTemperatureK(), targetCurrentState.getPressureBara());
        hydrocarbonWaterBranchSwitchDiagnostic = HydrocarbonWaterBoundaryAuditRunner.runBranchSwitchDiagnostic(
            switchSystem, sourceState, targetPreviousState, targetCurrentState, additionalPointCount, initialArcStep,
            minimumArcStep, finiteDifferenceStep, minimumTemperatureK, maximumTemperatureK, minimumPressureBara,
            maximumPressureBara);
      }

      JsonObject hydrocarbonWaterSecondaryStationaryBranchDiagnostic = null;
      if (input.has("hydrocarbonWaterSecondaryStationaryBranchDiagnostic")) {
        if (!input.get("hydrocarbonWaterSecondaryStationaryBranchDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterSecondaryStationaryBranchDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "secondary-stationary branch diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject secondaryRequest = input.getAsJsonObject("hydrocarbonWaterSecondaryStationaryBranchDiagnostic");
        List<State> orderedStates = new ArrayList<State>();
        JsonObject checkpointChainDiagnostic = null;
        if (secondaryRequest.has("orderedStateSegments")) {
          if (!secondaryRequest.get("orderedStateSegments").isJsonArray()) {
            throw new IllegalArgumentException("orderedStateSegments must be an array of named state segments");
          }
          JsonArray serializedSegments = secondaryRequest.getAsJsonArray("orderedStateSegments");
          if (serializedSegments.size() < 2 || serializedSegments.size() > 40) {
            throw new IllegalArgumentException("orderedStateSegments must contain 2-40 segments");
          }
          List<Segment> segments = new ArrayList<Segment>();
          for (int segmentIndex = 0; segmentIndex < serializedSegments.size(); segmentIndex++) {
            if (!serializedSegments.get(segmentIndex).isJsonObject()) {
              throw new IllegalArgumentException("orderedStateSegments entries must be objects");
            }
            JsonObject serializedSegment = serializedSegments.get(segmentIndex).getAsJsonObject();
            String identifier = serializedSegment.has("identifier") ? serializedSegment.get("identifier").getAsString()
                : "checkpoint-segment-" + segmentIndex;
            if (!serializedSegment.has("states") || !serializedSegment.get("states").isJsonArray()) {
              throw new IllegalArgumentException("orderedStateSegments entries require a states array");
            }
            JsonArray serializedStates = serializedSegment.getAsJsonArray("states");
            List<State> states = continuationStates(serializedStates, "orderedStateSegments[" + segmentIndex + "]");
            segments.add(new Segment(identifier, states));
          }
          HydrocarbonWaterContinuationChainAssembler.Result checkpointChain = new HydrocarbonWaterContinuationChainAssembler()
              .assemble(segments);
          checkpointChainDiagnostic = checkpointChainDiagnostic(checkpointChain);
          if (!checkpointChain.isSingleConnectedChain()) {
            throw new IllegalArgumentException(
                "orderedStateSegments do not form one identity-connected checkpoint chain: "
                    + checkpointChain.getViolations());
          }
          orderedStates.addAll(checkpointChain.getChains().get(0).getStates());
        } else {
          if (!secondaryRequest.has("orderedStates") || !secondaryRequest.get("orderedStates").isJsonArray()) {
            throw new IllegalArgumentException("orderedStates must be an array of continuation states");
          }
          orderedStates.addAll(continuationStates(secondaryRequest.getAsJsonArray("orderedStates"), "orderedStates"));
        }
        State rejectedState = continuationState(secondaryRequest, "rejectedState");
        CandidatePhase stationaryPhase = CandidatePhase.valueOf(secondaryRequest.get("stationaryPhase").getAsString());
        double[] terminalComposition = requiredNumberArray(secondaryRequest, "terminalComposition", 0.0);
        SystemInterface secondarySystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            rejectedState.getTemperatureK(), rejectedState.getPressureBara());
        hydrocarbonWaterSecondaryStationaryBranchDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runSecondaryStationaryBranchDiagnostic(secondarySystem, orderedStates, rejectedState, stationaryPhase,
                terminalComposition);
        if (checkpointChainDiagnostic != null) {
          hydrocarbonWaterSecondaryStationaryBranchDiagnostic.add("checkpointChainAssembly", checkpointChainDiagnostic);
        }
      }

      JsonObject hydrocarbonWaterClosedBoundaryLoopDiagnostic = null;
      if (input.has("hydrocarbonWaterClosedBoundaryLoopDiagnostic")) {
        if (!input.get("hydrocarbonWaterClosedBoundaryLoopDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterClosedBoundaryLoopDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "closed hydrocarbon-water boundary-loop diagnostic requires a non-reactive characterized fluid definition");
        }
        JsonObject loopRequest = input.getAsJsonObject("hydrocarbonWaterClosedBoundaryLoopDiagnostic");
        if (!loopRequest.has("firstPath") || !loopRequest.get("firstPath").isJsonArray()
            || !loopRequest.has("secondPath") || !loopRequest.get("secondPath").isJsonArray()) {
          throw new IllegalArgumentException(
              "closed boundary-loop diagnostic requires firstPath and secondPath arrays");
        }
        List<State> firstPath = continuationStates(loopRequest.getAsJsonArray("firstPath"), "firstPath");
        List<State> secondPath = continuationStates(loopRequest.getAsJsonArray("secondPath"), "secondPath");
        boolean verifyEvidence = loopRequest.has("verifyEvidence") && loopRequest.get("verifyEvidence").getAsBoolean();
        String closureContract = loopRequest.has("closureContract")
            ? loopRequest.get("closureContract").getAsString().toUpperCase(Locale.ROOT)
            : "TARGET_BRANCH_MERGE";
        if (!"TARGET_BRANCH_MERGE".equals(closureContract) && !"ORDINARY_STATE_SEAM".equals(closureContract)) {
          throw new IllegalArgumentException(
              "closed boundary-loop closureContract must be TARGET_BRANCH_MERGE or ORDINARY_STATE_SEAM");
        }
        State referenceState = firstPath.get(0);
        SystemInterface loopSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            referenceState.getTemperatureK(), referenceState.getPressureBara());
        hydrocarbonWaterClosedBoundaryLoopDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runClosedBoundaryLoopDiagnostic(loopSystem, firstPath, secondPath, verifyEvidence,
                "ORDINARY_STATE_SEAM".equals(closureContract));
      }

      JsonObject hydrocarbonWaterOpenBoundaryVerificationDiagnostic = null;
      if (input.has("hydrocarbonWaterOpenBoundaryVerificationDiagnostic")) {
        if (!input.get("hydrocarbonWaterOpenBoundaryVerificationDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterOpenBoundaryVerificationDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "open hydrocarbon-water boundary verification requires a non-reactive characterized fluid definition");
        }
        JsonObject verificationRequest = input.getAsJsonObject("hydrocarbonWaterOpenBoundaryVerificationDiagnostic");
        List<State> orderedStates = new ArrayList<State>();
        JsonObject checkpointChainDiagnostic = null;
        if (verificationRequest.has("orderedStateSegments")
            && verificationRequest.get("orderedStateSegments").isJsonArray()) {
          JsonArray serializedSegments = verificationRequest.getAsJsonArray("orderedStateSegments");
          List<Segment> segments = new ArrayList<Segment>();
          for (int index = 0; index < serializedSegments.size(); index++) {
            JsonObject serialized = serializedSegments.get(index).getAsJsonObject();
            String identifier = serialized.has("identifier") ? serialized.get("identifier").getAsString()
                : "open-checkpoint-segment-" + index;
            segments.add(new Segment(identifier,
                continuationStates(serialized.getAsJsonArray("states"), identifier + ".states")));
          }
          HydrocarbonWaterContinuationChainAssembler.Result checkpointChain = new HydrocarbonWaterContinuationChainAssembler()
              .assemble(segments);
          checkpointChainDiagnostic = checkpointChainDiagnostic(checkpointChain);
          if (!checkpointChain.isSingleConnectedChain()) {
            throw new IllegalArgumentException(
                "orderedStateSegments do not form one identity-connected checkpoint chain: "
                    + checkpointChain.getViolations());
          }
          orderedStates.addAll(checkpointChain.getChains().get(0).getStates());
        } else if (verificationRequest.has("orderedStates") && verificationRequest.get("orderedStates").isJsonArray()) {
          orderedStates
              .addAll(continuationStates(verificationRequest.getAsJsonArray("orderedStates"), "orderedStates"));
        } else {
          throw new IllegalArgumentException(
              "open boundary verification requires orderedStates or orderedStateSegments");
        }
        State rejectedHighState = continuationState(verificationRequest, "firstRejectedHighState");
        double minimumPressureBara = requiredPositiveNumber(verificationRequest, "minimumPressureBara");
        State referenceState = orderedStates.get(0);
        SystemInterface verificationSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            referenceState.getTemperatureK(), referenceState.getPressureBara());
        hydrocarbonWaterOpenBoundaryVerificationDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runOpenBoundaryVerificationDiagnostic(verificationSystem, orderedStates, minimumPressureBara,
                rejectedHighState);
        if (checkpointChainDiagnostic != null) {
          hydrocarbonWaterOpenBoundaryVerificationDiagnostic.add("checkpointChainAssembly", checkpointChainDiagnostic);
        }
      }

      JsonObject hydrocarbonWaterFullNetworkVerificationDiagnostic = null;
      if (input.has("hydrocarbonWaterFullNetworkVerificationDiagnostic")) {
        if (!input.get("hydrocarbonWaterFullNetworkVerificationDiagnostic").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterFullNetworkVerificationDiagnostic must be an object");
        }
        if (fluidDefinition == null || reactive) {
          throw new IllegalArgumentException(
              "full hydrocarbon-water network verification requires a non-reactive characterized fluid definition");
        }
        JsonObject verificationRequest = input.getAsJsonObject("hydrocarbonWaterFullNetworkVerificationDiagnostic");
        List<State> closedLoopStates = continuationStates(verificationRequest.getAsJsonArray("closedLoopStates"),
            "closedLoopStates");
        JsonArray serializedSegments = verificationRequest.getAsJsonArray("openStateSegments");
        if (serializedSegments == null || serializedSegments.size() < 2) {
          throw new IllegalArgumentException("full network verification requires at least two openStateSegments");
        }
        List<Segment> segments = new ArrayList<Segment>();
        for (int index = 0; index < serializedSegments.size(); index++) {
          JsonObject serialized = serializedSegments.get(index).getAsJsonObject();
          String identifier = serialized.has("identifier") ? serialized.get("identifier").getAsString()
              : "full-network-open-segment-" + index;
          segments.add(
              new Segment(identifier, continuationStates(serialized.getAsJsonArray("states"), identifier + ".states")));
        }
        HydrocarbonWaterContinuationChainAssembler.Result openChain = new HydrocarbonWaterContinuationChainAssembler()
            .assemble(segments);
        if (!openChain.isSingleConnectedChain()) {
          throw new IllegalArgumentException(
              "openStateSegments do not form one identity-connected chain: " + openChain.getViolations());
        }
        List<State> openStates = openChain.getChains().get(0).getStates();
        State rejectedHighState = continuationState(verificationRequest, "firstRejectedHighState");
        double minimumPressureBara = requiredPositiveNumber(verificationRequest, "minimumPressureBara");
        State referenceState = closedLoopStates.get(0);
        SystemInterface verificationSystem = buildBoundarySystem(composition, fluidDefinition, requestedModel,
            referenceState.getTemperatureK(), referenceState.getPressureBara());
        hydrocarbonWaterFullNetworkVerificationDiagnostic = HydrocarbonWaterBoundaryAuditRunner
            .runFullNetworkVerificationDiagnostic(verificationSystem, closedLoopStates, openStates, minimumPressureBara,
                rejectedHighState);
        hydrocarbonWaterFullNetworkVerificationDiagnostic.add("openCheckpointChainAssembly",
            checkpointChainDiagnostic(openChain));
      }

      JsonObject hydrocarbonWaterBoundaryRegression = null;
      if (input.has("hydrocarbonWaterBoundaryRegression")) {
        if (!input.get("hydrocarbonWaterBoundaryRegression").isJsonObject()) {
          throw new IllegalArgumentException("hydrocarbonWaterBoundaryRegression must be an object");
        }
        hydrocarbonWaterBoundaryRegression = HydrocarbonWaterBoundaryRegressionRunner
            .run(input.getAsJsonObject("hydrocarbonWaterBoundaryRegression"));
      }

      JsonObject result = new JsonObject();
      result.addProperty("status", "success");
      String model = reactive
          ? (fluidDefinition == null ? "Electrolyte-CPA" : "Electrolyte-CPA + imported heavy-oil pseudo-components")
          : (fluidDefinition == null ? normalizedModel(requestedModel, false) : importedEos);
      result.addProperty("model", model);
      boolean huronVidalModel = fluidDefinition != null && fluidDefinition.has("polarModel")
          && "HV".equalsIgnoreCase(fluidDefinition.get("polarModel").getAsString());
      String activeMixingRule = huronVidalModel ? "Huron-Vidal (NRTL)" : "Classic";
      result.addProperty("mixingRule", activeMixingRule);
      if (fluidDefinition != null) {
        result.addProperty("sourceModel", sourceEos);
        result.addProperty("candidateModel", model);
      }
      result.addProperty("method",
          reactive ? "Electrolyte-CPA TP flash with aqueous reaction set and charge-balance audit"
              : (fluidDefinition == null ? model + " non-reactive multiphase TP flash"
                  : "Characterized " + model
                      + " multiphase TP flash with topology-preserving 2-HC phase-envelope tracing"));
      result.add("fluid",
          JsonParser.parseString(GSON.toJson(new FluidResponse("field-fluid", state))).getAsJsonObject());
      JsonObject modelIdentity = HydrocarbonWaterModelIdentity.describe(state, model, activeMixingRule,
          fluidDefinition);
      result.add("modelIdentity", modelIdentity);

      JsonObject flash = new JsonObject();
      flash.addProperty("model", model);
      flash.addProperty("flashType", "TP");
      flash.addProperty("numberOfPhases", state.getNumberOfPhases());
      flash.addProperty("reactiveEquilibriumConverged", reactiveAudit.converged);
      flash.addProperty("reactiveIterations", reactiveAudit.iterations);
      flash.addProperty("independentReactions", reactiveAudit.reactions);
      JsonArray phases = new JsonArray();
      for (int index = 0; index < state.getNumberOfPhases(); index++) {
        phases.add(state.getPhase(index).getPhaseTypeName());
      }
      flash.add("phases", phases);
      result.add("flash", flash);
      result.add("envelope", boundary);
      result.add("phaseMap", phaseMap);
      result.add("criticalPoint", criticalPoint(twoHydrocarbonEnvelope));
      result.add("speciation", buildSpeciation(state, reactive));
      result.addProperty("boundaryTopology", boundaryTopology);
      if (hydrocarbonWaterBoundaryAudit != null) {
        result.add("hydrocarbonWaterBoundaryAudit", hydrocarbonWaterBoundaryAudit);
      }
      if (hydrocarbonWaterStableRegionTopologyScan != null) {
        result.add("hydrocarbonWaterStableRegionTopologyScan", hydrocarbonWaterStableRegionTopologyScan);
      }
      if (hydrocarbonWaterRegularPressureDiagnostic != null) {
        result.add("hydrocarbonWaterRegularPressureDiagnostic", hydrocarbonWaterRegularPressureDiagnostic);
      }
      if (hydrocarbonWaterModelSeededBoundaryDiagnostic != null) {
        result.add("hydrocarbonWaterModelSeededBoundaryDiagnostic", hydrocarbonWaterModelSeededBoundaryDiagnostic);
      }
      if (hydrocarbonWaterStableBracketCorrectionDiagnostic != null) {
        result.add("hydrocarbonWaterStableBracketCorrectionDiagnostic",
            hydrocarbonWaterStableBracketCorrectionDiagnostic);
      }
      if (hydrocarbonWaterRegularPressureRestartDiagnostic != null) {
        result.add("hydrocarbonWaterRegularPressureRestartDiagnostic",
            hydrocarbonWaterRegularPressureRestartDiagnostic);
      }
      if (hydrocarbonWaterHybridBranchDiagnostic != null) {
        result.add("hydrocarbonWaterHybridBranchDiagnostic", hydrocarbonWaterHybridBranchDiagnostic);
      }
      if (hydrocarbonWaterPseudoArcRestartDiagnostic != null) {
        result.add("hydrocarbonWaterPseudoArcRestartDiagnostic", hydrocarbonWaterPseudoArcRestartDiagnostic);
      }
      if (hydrocarbonWaterLocalTangentProbeDiagnostic != null) {
        result.add("hydrocarbonWaterLocalTangentProbeDiagnostic", hydrocarbonWaterLocalTangentProbeDiagnostic);
      }
      if (hydrocarbonWaterRestartStateTopologyDiagnostic != null) {
        result.add("hydrocarbonWaterRestartStateTopologyDiagnostic", hydrocarbonWaterRestartStateTopologyDiagnostic);
      }
      if (hydrocarbonWaterStrictBranchSeedDiagnostic != null) {
        result.add("hydrocarbonWaterStrictBranchSeedDiagnostic", hydrocarbonWaterStrictBranchSeedDiagnostic);
      }
      if (hydrocarbonWaterBranchSwitchDiagnostic != null) {
        result.add("hydrocarbonWaterBranchSwitchDiagnostic", hydrocarbonWaterBranchSwitchDiagnostic);
      }
      if (hydrocarbonWaterSecondaryStationaryBranchDiagnostic != null) {
        result.add("hydrocarbonWaterSecondaryStationaryBranchDiagnostic",
            hydrocarbonWaterSecondaryStationaryBranchDiagnostic);
      }
      if (hydrocarbonWaterClosedBoundaryLoopDiagnostic != null) {
        result.add("hydrocarbonWaterClosedBoundaryLoopDiagnostic", hydrocarbonWaterClosedBoundaryLoopDiagnostic);
      }
      if (hydrocarbonWaterOpenBoundaryVerificationDiagnostic != null) {
        result.add("hydrocarbonWaterOpenBoundaryVerificationDiagnostic",
            hydrocarbonWaterOpenBoundaryVerificationDiagnostic);
      }
      if (hydrocarbonWaterFullNetworkVerificationDiagnostic != null) {
        result.add("hydrocarbonWaterFullNetworkVerificationDiagnostic",
            hydrocarbonWaterFullNetworkVerificationDiagnostic);
      }
      if (hydrocarbonWaterBoundaryRegression != null) {
        result.add("hydrocarbonWaterBoundaryRegression", hydrocarbonWaterBoundaryRegression);
      }

      JsonObject quality = new JsonObject();
      boolean boundaryCoveragePassed = boundary.size() >= 8;
      boolean boundaryContinuityPassed = boundaryCoveragePassed
          && (twoHydrocarbonEnvelope != null || phaseMapAudit != null || fullEnvelopeCalculated
              ? boundaryHasFinitePositiveValues(boundary)
              : boundaryIsContinuous(boundary));
      double maximumBoundaryPressure = maximumBoundaryPressure(boundary);
      double minimumComparablePressure = Math.max(1.0, pressureBara * 0.1);
      boolean boundaryPressureScalePassed = boundaryCoveragePassed
          && maximumBoundaryPressure >= minimumComparablePressure;
      // A stability-contour fallback is a marching-squares trace of a grid criterion, not a continuation-traced
      // envelope. Scoring one against a reference envelope produces a large temperature difference that looks
      // exactly like an inequivalent model, so it must never be advertised as comparable. This enforces the
      // frozen project rule that grid contours stay diagnostic and never become the engineering boundary.
      boolean continuationTracedBoundary = !boundaryMethod.contains("marching-squares")
          && !boundaryMethod.contains("stability-contour");
      boolean boundaryComparisonEligible = boundaryCoveragePassed && boundaryContinuityPassed
          && boundaryPressureScalePassed && boundaryDefinitionCompatible && continuationTracedBoundary;
      quality.addProperty("continuationTracedBoundary", continuationTracedBoundary);
      quality.addProperty("verdict", boundaryComparisonEligible ? "passed" : "warning");
      quality.addProperty("statePointCalculated", true);
      quality.addProperty("reactiveEquilibriumConverged", reactiveAudit.converged);
      quality.addProperty("boundaryPointCount", boundary.size());
      quality.addProperty("requestedBoundaryPointCount", pointCount);
      quality.addProperty("rejectedBoundaryPointCount", rejectedBoundaryPoints);
      quality.addProperty("boundaryMethod", boundaryMethod);
      quality.addProperty("boundaryTopology", boundaryTopology);
      quality.addProperty("fullEnvelopeCalculated", fullEnvelopeCalculated);
      quality.addProperty("phaseMapCalculated", phaseMapAudit != null);
      quality.addProperty("boundaryDefinitionCompatible", boundaryDefinitionCompatible);
      quality.addProperty("boundaryContinuityPassed", boundaryContinuityPassed);
      quality.addProperty("boundaryPressureScalePassed", boundaryPressureScalePassed);
      quality.addProperty("boundaryMaximumPressureBara", maximumBoundaryPressure);
      quality.addProperty("boundaryMinimumComparablePressureBara", minimumComparablePressure);
      quality.addProperty("boundaryComparisonEligible", boundaryComparisonEligible);
      quality.addProperty("hydrocarbonWaterBoundaryAuditCalculated", hydrocarbonWaterBoundaryAudit != null);
      quality.addProperty("hydrocarbonWaterStableRegionTopologyScanCalculated",
          hydrocarbonWaterStableRegionTopologyScan != null);
      quality.addProperty("hydrocarbonWaterStableBracketCorrectionDiagnosticCalculated",
          hydrocarbonWaterStableBracketCorrectionDiagnostic != null);
      quality.addProperty("hydrocarbonWaterRegularPressureRestartDiagnosticCalculated",
          hydrocarbonWaterRegularPressureRestartDiagnostic != null);
      quality.addProperty("hydrocarbonWaterHybridBranchDiagnosticCalculated",
          hydrocarbonWaterHybridBranchDiagnostic != null);
      quality.addProperty("hydrocarbonWaterPseudoArcRestartDiagnosticCalculated",
          hydrocarbonWaterPseudoArcRestartDiagnostic != null);
      quality.addProperty("hydrocarbonWaterLocalTangentProbeDiagnosticCalculated",
          hydrocarbonWaterLocalTangentProbeDiagnostic != null);
      quality.addProperty("hydrocarbonWaterRestartStateTopologyDiagnosticCalculated",
          hydrocarbonWaterRestartStateTopologyDiagnostic != null);
      quality.addProperty("hydrocarbonWaterStrictBranchSeedDiagnosticCalculated",
          hydrocarbonWaterStrictBranchSeedDiagnostic != null);
      quality.addProperty("hydrocarbonWaterBranchSwitchDiagnosticCalculated",
          hydrocarbonWaterBranchSwitchDiagnostic != null);
      quality.addProperty("hydrocarbonWaterSecondaryStationaryBranchDiagnosticCalculated",
          hydrocarbonWaterSecondaryStationaryBranchDiagnostic != null);
      quality.addProperty("hydrocarbonWaterClosedBoundaryLoopDiagnosticCalculated",
          hydrocarbonWaterClosedBoundaryLoopDiagnostic != null);
      quality.addProperty("hydrocarbonWaterOpenBoundaryVerificationDiagnosticCalculated",
          hydrocarbonWaterOpenBoundaryVerificationDiagnostic != null);
      quality.addProperty("hydrocarbonWaterFullNetworkVerificationDiagnosticCalculated",
          hydrocarbonWaterFullNetworkVerificationDiagnostic != null);
      quality.addProperty("hydrocarbonWaterBoundaryRegressionCalculated", hydrocarbonWaterBoundaryRegression != null);
      if (hydrocarbonWaterBoundaryRegression != null) {
        quality.addProperty("hydrocarbonWaterBoundaryRegressionComparisonEligible",
            hydrocarbonWaterBoundaryRegression.get("comparisonEligible").getAsBoolean());
        quality.addProperty("hydrocarbonWaterBoundaryRegressionAccepted",
            hydrocarbonWaterBoundaryRegression.get("accepted").getAsBoolean());
        quality.addProperty("hydrocarbonWaterBoundaryRegressionThresholdsFrozen",
            hydrocarbonWaterBoundaryRegression.get("thresholdsFrozen").getAsBoolean());
        quality.addProperty("hydrocarbonWaterBoundaryRegressionBenchmarkPending",
            hydrocarbonWaterBoundaryRegression.get("benchmarkPending").getAsBoolean());
      }
      if (hydrocarbonWaterBoundaryAudit != null) {
        quality.addProperty("hydrocarbonWaterPhysicalTerminationCount",
            hydrocarbonWaterBoundaryAudit.get("physicalTerminationCount").getAsInt());
        quality.addProperty("hydrocarbonWaterAuditEngineeringEligible",
            hydrocarbonWaterBoundaryAudit.get("engineeringEligible").getAsBoolean());
      }
      quality.addProperty("modelTransferredFromSource", modelTransferredFromSource);
      if (phaseMapAudit != null) {
        quality.addProperty("phaseMapTemperaturePointCount", phaseMapAudit.temperaturePointCount);
        quality.addProperty("phaseMapPressurePointCount", phaseMapAudit.pressurePointCount);
        quality.addProperty("phaseMapCoarsePointCount", phaseMapAudit.samples.size());
        quality.addProperty("phaseMapFlashCount", phaseMapAudit.flashCount);
        quality.addProperty("phaseMapFailedFlashCount", phaseMapAudit.failedFlashCount);
        quality.addProperty("phaseMapHydrocarbonComponentCount", phaseMapAudit.hydrocarbonComponentCount);
        quality.addProperty("phaseMapSelectedHydrocarbonCellCount", phaseMapAudit.selectedHydrocarbonCellCount);
        quality.addProperty("phaseMapExcludedHydrocarbonCellCount", phaseMapAudit.excludedHydrocarbonCellCount);
        quality.addProperty("phaseMapTemperatureMinimumC", phaseMapAudit.temperatureMinimumC);
        quality.addProperty("phaseMapTemperatureMaximumC", phaseMapAudit.temperatureMaximumC);
        quality.addProperty("phaseMapPressureMinimumBara", phaseMapAudit.pressureMinimumBara);
        quality.addProperty("phaseMapPressureMaximumBara", phaseMapAudit.pressureMaximumBara);
      }
      if (hydrocarbonDewBranchQuality != null) {
        for (Map.Entry<String, JsonElement> entry : hydrocarbonDewBranchQuality.entrySet()) {
          quality.add(entry.getKey(), entry.getValue());
        }
        quality.addProperty("boundaryComparableRange", "hydrocarbon dew branch only");
      }
      if (twoHydrocarbonEnvelope != null) {
        quality.addProperty("envelopeSolver", twoHydrocarbonEnvelope.getMethod());
        quality.addProperty("envelopeSegmentCount", twoHydrocarbonEnvelope.getSegments().size());
        quality.addProperty("envelopeClosed", twoHydrocarbonEnvelope.isEnvelopeClosed());
        quality.addProperty("envelopeCriticalPointCount", twoHydrocarbonEnvelope.getCriticalPointCount());
        quality.addProperty("envelopeTermination", twoHydrocarbonEnvelope.getTerminationSummary());
        quality.addProperty("envelopeBasis", "hydrocarbon-only");
        quality.addProperty("envelopeBasisExcludedAqueousFraction", aqueousBasisFraction);
        quality.addProperty("continuationIterationLimitReached", twoHydrocarbonEnvelope.isIterationLimitReached());
        quality.addProperty("stabilityFlashCount", twoHydrocarbonEnvelope.getFlashCount());
        quality.addProperty("stabilityFailedFlashCount", twoHydrocarbonEnvelope.getFailedFlashCount());
        quality.addProperty("threePhaseBoundaryPointCount", twoHydrocarbonEnvelope.getThreePhaseStabilityPointCount());
        quality.addProperty("topologyFallbackEnabled", !huronVidalModel || (input.has("allowExpensiveTopologyFallback")
            && input.get("allowExpensiveTopologyFallback").getAsBoolean()));
      }
      quality.addProperty("engineeringReviewRequired", true);
      quality.addProperty("summary", boundaryComparisonEligible ? boundary.size()
          + " calculated definition-compatible phase-boundary points passed coverage, pressure-scale, and continuity checks; "
          + rejectedBoundaryPoints + " invalid points were excluded"
          : "State point calculated; a definition-compatible phase envelope did not pass coverage and pressure-scale gates, so it is not eligible for curve comparison");
      result.add("qualityGate", quality);

      JsonArray limitations = new JsonArray();
      if (reactive) {
        limitations.add(
            "Electrolyte reactions are solved only in the aqueous phase; imported heavy-oil pseudo-components remain non-reactive.");
        limitations.add(
            "NH3/CO2 speciation must be reviewed against laboratory pH or ion-analysis data before design freeze.");
        limitations.add(fluidDefinition == null
            ? "The plotted phase boundary uses the selected Electrolyte-CPA feed species without reaction coupling; the state point and aqueous speciation use the converged reactive flash."
            : "The plotted phase boundary uses the imported heavy-oil EOS without reaction coupling; the state point and aqueous speciation use Electrolyte-CPA with imported pseudo-components.");
      }
      if (fluidDefinition != null) {
        if (modelTransferredFromSource) {
          limitations.add("The equation of state was changed from " + sourceEos + " to " + importedEos
              + "; pseudo-component properties, volume shifts, and common Kij were transferred without retuning, so an independent benchmark is required.");
        }
        limitations.add(
            "Imported pseudo-component critical properties, volume shifts, and common Kij are retained; new water/NH3 interactions use the selected NeqSim mixing-rule database defaults.");
        if (twoHydrocarbonEnvelope != null) {
          if (aqueousBasisFraction > 0.0) {
            limitations.add("The 2-HC phase envelope is traced on the normalized hydrocarbon sub-composition; "
                + "aqueous species remain in the reported state-point flash but are excluded from the envelope basis.");
          }
          limitations.add(
              "The characterized-fluid 2-HC boundary uses Michelsen predictor-corrector continuation when it initializes; otherwise it traces the stable gas-and-hydrocarbon-liquid coexistence contour without connecting unrelated TP states.");
          limitations.add(
              "All disconnected stable 2-HC boundary segments are preserved. Aqueous-phase appearance is audited separately and does not redefine the PVTsim 2-HC boundary.");
          limitations.add(
              "The stability-contour fallback returns only bracketed transitions inside its automatic domain and never extrapolates or joins missing branches.");
        }
      }
      int warningLimit = Math.min(boundaryWarnings.size(), 3);
      for (int index = 0; index < warningLimit; index++) {
        limitations.add(boundaryWarnings.get(index));
      }
      result.add("limitations", limitations);
      result.addProperty("computationTimeMs", System.currentTimeMillis() - started);

      JsonObject data = new JsonObject();
      data.add("flash", flash.deepCopy());
      data.add("fluid", result.getAsJsonObject("fluid").deepCopy());
      data.add("envelope", boundary.deepCopy());
      data.add("phaseMap", phaseMap.deepCopy());
      data.add("criticalPoint", criticalPoint(twoHydrocarbonEnvelope));
      data.add("modelIdentity", modelIdentity.deepCopy());
      if (hydrocarbonWaterBoundaryAudit != null) {
        data.add("hydrocarbonWaterBoundaryAudit", hydrocarbonWaterBoundaryAudit.deepCopy());
      }
      if (hydrocarbonWaterStableRegionTopologyScan != null) {
        data.add("hydrocarbonWaterStableRegionTopologyScan", hydrocarbonWaterStableRegionTopologyScan.deepCopy());
      }
      if (hydrocarbonWaterRegularPressureDiagnostic != null) {
        data.add("hydrocarbonWaterRegularPressureDiagnostic", hydrocarbonWaterRegularPressureDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterModelSeededBoundaryDiagnostic != null) {
        data.add("hydrocarbonWaterModelSeededBoundaryDiagnostic",
            hydrocarbonWaterModelSeededBoundaryDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterStableBracketCorrectionDiagnostic != null) {
        data.add("hydrocarbonWaterStableBracketCorrectionDiagnostic",
            hydrocarbonWaterStableBracketCorrectionDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterRegularPressureRestartDiagnostic != null) {
        data.add("hydrocarbonWaterRegularPressureRestartDiagnostic",
            hydrocarbonWaterRegularPressureRestartDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterHybridBranchDiagnostic != null) {
        data.add("hydrocarbonWaterHybridBranchDiagnostic", hydrocarbonWaterHybridBranchDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterPseudoArcRestartDiagnostic != null) {
        data.add("hydrocarbonWaterPseudoArcRestartDiagnostic", hydrocarbonWaterPseudoArcRestartDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterLocalTangentProbeDiagnostic != null) {
        data.add("hydrocarbonWaterLocalTangentProbeDiagnostic", hydrocarbonWaterLocalTangentProbeDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterRestartStateTopologyDiagnostic != null) {
        data.add("hydrocarbonWaterRestartStateTopologyDiagnostic",
            hydrocarbonWaterRestartStateTopologyDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterStrictBranchSeedDiagnostic != null) {
        data.add("hydrocarbonWaterStrictBranchSeedDiagnostic", hydrocarbonWaterStrictBranchSeedDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterBranchSwitchDiagnostic != null) {
        data.add("hydrocarbonWaterBranchSwitchDiagnostic", hydrocarbonWaterBranchSwitchDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterSecondaryStationaryBranchDiagnostic != null) {
        data.add("hydrocarbonWaterSecondaryStationaryBranchDiagnostic",
            hydrocarbonWaterSecondaryStationaryBranchDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterClosedBoundaryLoopDiagnostic != null) {
        data.add("hydrocarbonWaterClosedBoundaryLoopDiagnostic",
            hydrocarbonWaterClosedBoundaryLoopDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterOpenBoundaryVerificationDiagnostic != null) {
        data.add("hydrocarbonWaterOpenBoundaryVerificationDiagnostic",
            hydrocarbonWaterOpenBoundaryVerificationDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterFullNetworkVerificationDiagnostic != null) {
        data.add("hydrocarbonWaterFullNetworkVerificationDiagnostic",
            hydrocarbonWaterFullNetworkVerificationDiagnostic.deepCopy());
      }
      if (hydrocarbonWaterBoundaryRegression != null) {
        data.add("hydrocarbonWaterBoundaryRegression", hydrocarbonWaterBoundaryRegression.deepCopy());
      }
      result.add("data", data);
      return GSON.toJson(result);
    } catch (Exception error) {
      JsonObject result = new JsonObject();
      result.addProperty("status", "error");
      result.addProperty("message", "Field-fluid calculation failed: " + error.getMessage());
      return GSON.toJson(result);
    }
  }

  private static JsonArray twoHydrocarbonEnvelopeRows(TwoHydrocarbonPhaseEnvelopeSolver.Result envelope) {
    JsonArray output = new JsonArray();
    int segmentIndex = 0;
    for (TwoHydrocarbonPhaseEnvelopeSolver.Segment segment : envelope.getSegments()) {
      segmentIndex++;
      String branch = segment.getPhysicalType().name().toLowerCase(Locale.ROOT);
      String methodPrefix = envelope.getMethod().startsWith("michelsen") ? "michelsen" : "stability";
      String curve = methodPrefix + "-2HC-" + branch + "-" + segmentIndex;
      double[] temperatures = segment.getTemperaturesK();
      double[] pressures = segment.getPressuresBara();
      for (int pointIndex = 0; pointIndex < temperatures.length; pointIndex++) {
        JsonObject row = new JsonObject();
        row.addProperty("temperature_C", temperatures[pointIndex] - 273.15);
        row.addProperty("temperature_K", temperatures[pointIndex]);
        row.addProperty("pressure_bara", pressures[pointIndex]);
        row.addProperty("curve", curve);
        row.addProperty("branch", branch);
        row.addProperty("segmentIndex", segmentIndex);
        row.addProperty("pointIndex", pointIndex);
        row.addProperty("boundaryKind", "stable-two-hydrocarbon-phase-boundary");
        output.add(row);
      }
    }
    return output;
  }

  private static JsonObject criticalPoint(TwoHydrocarbonPhaseEnvelopeSolver.Result envelope) {
    if (envelope == null || envelope.getCriticalPointCount() < 1) {
      return null;
    }
    double[] critical = envelope.getCriticalPoint();
    if (critical.length < 2 || !Double.isFinite(critical[0]) || !Double.isFinite(critical[1]) || critical[0] < 50.0
        || critical[1] <= 0.0) {
      return null;
    }
    JsonObject output = new JsonObject();
    output.addProperty("temperature_C", critical[0] - 273.15);
    output.addProperty("temperature_K", critical[0]);
    output.addProperty("pressure_bara", critical[1]);
    output.addProperty("method", envelope.getMethod().startsWith("michelsen") ? "thermodynamic-critical-point"
        : "minimum-phase-composition-distance-estimate");
    return output;
  }

  private static PhaseMapAudit calculateHydrocarbonPhaseMap(JsonObject composition, JsonObject fluidDefinition,
      String requestedModel, double minimumTemperatureC, double maximumTemperatureC, double minimumPressureBara,
      double maximumPressureBara, int temperaturePointCount, int pressurePointCount, double targetTemperatureC,
      double targetPressureBara) {
    SystemInterface template = buildBoundarySystem(composition, fluidDefinition, requestedModel,
        minimumTemperatureC + 273.15, minimumPressureBara);
    JsonArray samples = new JsonArray();
    JsonArray boundary = new JsonArray();
    List<String> warnings = new ArrayList<String>();
    int flashCount = 0;
    int failedFlashCount = 0;
    double logarithmicPressureSpan = Math.log(maximumPressureBara / minimumPressureBara);

    PhaseMapState[][] states = new PhaseMapState[temperaturePointCount][pressurePointCount];
    for (int temperatureIndex = 0; temperatureIndex < temperaturePointCount; temperatureIndex++) {
      double temperatureFraction = (double) temperatureIndex / (temperaturePointCount - 1);
      double temperatureC = minimumTemperatureC + temperatureFraction * (maximumTemperatureC - minimumTemperatureC);
      for (int pressureIndex = 0; pressureIndex < pressurePointCount; pressureIndex++) {
        double pressureFraction = (double) pressureIndex / (pressurePointCount - 1);
        double pressureBara = minimumPressureBara * Math.exp(pressureFraction * logarithmicPressureSpan);
        flashCount++;
        try {
          PhaseMapState phaseState = evaluatePhaseMapState(template, temperatureC, pressureBara);
          states[temperatureIndex][pressureIndex] = phaseState;
        } catch (Exception error) {
          failedFlashCount++;
          states[temperatureIndex][pressureIndex] = PhaseMapState.failed(temperatureC, pressureBara);
          warnings.add(
              "Phase-map flash failed at " + temperatureC + " C and " + pressureBara + " bara: " + error.getMessage());
        }
      }
    }

    HydrocarbonComponentSelection selection = selectHydrocarbonComponent(states, targetTemperatureC, targetPressureBara,
        minimumTemperatureC, maximumTemperatureC, minimumPressureBara, maximumPressureBara);
    for (int temperatureIndex = 0; temperatureIndex < temperaturePointCount; temperatureIndex++) {
      for (int pressureIndex = 0; pressureIndex < pressurePointCount; pressureIndex++) {
        samples.add(phaseMapSample(states[temperatureIndex][pressureIndex],
            selection.selected[temperatureIndex][pressureIndex]));
      }
    }

    for (int temperatureIndex = 0; temperatureIndex < temperaturePointCount; temperatureIndex++) {
      int firstSelectedPressureIndex = -1;
      int lastSelectedPressureIndex = -1;
      for (int pressureIndex = 0; pressureIndex < pressurePointCount; pressureIndex++) {
        if (selection.selected[temperatureIndex][pressureIndex]) {
          if (firstSelectedPressureIndex < 0) {
            firstSelectedPressureIndex = pressureIndex;
          }
          lastSelectedPressureIndex = pressureIndex;
        }
      }
      if (firstSelectedPressureIndex < 0) {
        continue;
      }
      if (firstSelectedPressureIndex > 0) {
        BoundaryRefinement refinement = refineHydrocarbonBoundary(template,
            states[temperatureIndex][firstSelectedPressureIndex - 1],
            states[temperatureIndex][firstSelectedPressureIndex], warnings);
        flashCount += refinement.flashCount;
        failedFlashCount += refinement.failedFlashCount;
        if (refinement.row != null) {
          refinement.row.addProperty("curve", "adaptive-2HC-entry");
          boundary.add(refinement.row);
        }
      }
      if (lastSelectedPressureIndex < pressurePointCount - 1) {
        BoundaryRefinement refinement = refineHydrocarbonBoundary(template,
            states[temperatureIndex][lastSelectedPressureIndex],
            states[temperatureIndex][lastSelectedPressureIndex + 1], warnings);
        flashCount += refinement.flashCount;
        failedFlashCount += refinement.failedFlashCount;
        if (refinement.row != null) {
          refinement.row.addProperty("curve", "adaptive-2HC-exit");
          boundary.add(refinement.row);
        }
      }
    }
    return new PhaseMapAudit(boundary, samples, warnings, flashCount, failedFlashCount, temperaturePointCount,
        pressurePointCount, minimumTemperatureC, maximumTemperatureC, minimumPressureBara, maximumPressureBara,
        selection.componentCount, selection.selectedCellCount, selection.excludedCellCount);
  }

  private static HydrocarbonComponentSelection selectHydrocarbonComponent(PhaseMapState[][] states,
      double targetTemperatureC, double targetPressureBara, double minimumTemperatureC, double maximumTemperatureC,
      double minimumPressureBara, double maximumPressureBara) {
    int temperaturePointCount = states.length;
    int pressurePointCount = states[0].length;
    boolean[][] visited = new boolean[temperaturePointCount][pressurePointCount];
    List<List<int[]>> components = new ArrayList<List<int[]>>();
    int[][] neighborOffsets = new int[][] { { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } };
    for (int temperatureIndex = 0; temperatureIndex < temperaturePointCount; temperatureIndex++) {
      for (int pressureIndex = 0; pressureIndex < pressurePointCount; pressureIndex++) {
        if (visited[temperatureIndex][pressureIndex] || !states[temperatureIndex][pressureIndex].hydrocarbonTwoPhase) {
          continue;
        }
        List<int[]> component = new ArrayList<int[]>();
        ArrayDeque<int[]> queue = new ArrayDeque<int[]>();
        queue.add(new int[] { temperatureIndex, pressureIndex });
        visited[temperatureIndex][pressureIndex] = true;
        while (!queue.isEmpty()) {
          int[] current = queue.removeFirst();
          component.add(current);
          for (int[] offset : neighborOffsets) {
            int nextTemperatureIndex = current[0] + offset[0];
            int nextPressureIndex = current[1] + offset[1];
            if (nextTemperatureIndex < 0 || nextTemperatureIndex >= temperaturePointCount || nextPressureIndex < 0
                || nextPressureIndex >= pressurePointCount || visited[nextTemperatureIndex][nextPressureIndex]
                || !states[nextTemperatureIndex][nextPressureIndex].hydrocarbonTwoPhase) {
              continue;
            }
            visited[nextTemperatureIndex][nextPressureIndex] = true;
            queue.addLast(new int[] { nextTemperatureIndex, nextPressureIndex });
          }
        }
        components.add(component);
      }
    }

    boolean[][] selected = new boolean[temperaturePointCount][pressurePointCount];
    if (components.isEmpty()) {
      return new HydrocarbonComponentSelection(selected, 0, 0, 0);
    }
    double temperatureScale = Math.max(maximumTemperatureC - minimumTemperatureC, 1.0);
    double logarithmicPressureScale = Math.max(Math.log(maximumPressureBara / minimumPressureBara), 1.0);
    double targetLogPressure = Math.log(Math.max(targetPressureBara, minimumPressureBara));
    List<int[]> selectedComponent = null;
    double bestDistance = Double.POSITIVE_INFINITY;
    for (List<int[]> component : components) {
      double componentDistance = Double.POSITIVE_INFINITY;
      for (int[] cell : component) {
        PhaseMapState state = states[cell[0]][cell[1]];
        double temperatureDistance = (state.temperatureC - targetTemperatureC) / temperatureScale;
        double pressureDistance = (Math.log(state.pressureBara) - targetLogPressure) / logarithmicPressureScale;
        componentDistance = Math.min(componentDistance,
            temperatureDistance * temperatureDistance + pressureDistance * pressureDistance);
      }
      if (componentDistance < bestDistance) {
        bestDistance = componentDistance;
        selectedComponent = component;
      }
    }
    for (int[] cell : selectedComponent) {
      selected[cell[0]][cell[1]] = true;
    }
    int totalHydrocarbonCellCount = 0;
    for (List<int[]> component : components) {
      totalHydrocarbonCellCount += component.size();
    }
    return new HydrocarbonComponentSelection(selected, components.size(), selectedComponent.size(),
        totalHydrocarbonCellCount - selectedComponent.size());
  }

  private static BoundaryRefinement refineHydrocarbonBoundary(SystemInterface template, PhaseMapState lower,
      PhaseMapState upper, List<String> warnings) {
    if (!lower.valid || !upper.valid || lower.hydrocarbonTwoPhase == upper.hydrocarbonTwoPhase) {
      return new BoundaryRefinement(null, 0, 0);
    }
    int flashCount = 0;
    for (int iteration = 0; iteration < PHASE_MAP_BISECTION_ITERATIONS; iteration++) {
      double midpointPressure = Math.sqrt(lower.pressureBara * upper.pressureBara);
      flashCount++;
      try {
        PhaseMapState midpoint = evaluatePhaseMapState(template, lower.temperatureC, midpointPressure);
        if (midpoint.hydrocarbonTwoPhase == lower.hydrocarbonTwoPhase) {
          lower = midpoint;
        } else {
          upper = midpoint;
        }
      } catch (Exception error) {
        warnings.add("Phase-map boundary refinement failed at " + lower.temperatureC + " C and " + midpointPressure
            + " bara: " + error.getMessage());
        return new BoundaryRefinement(null, flashCount, 1);
      }
    }
    JsonObject row = new JsonObject();
    row.addProperty("temperature_C", lower.temperatureC);
    row.addProperty("temperature_K", lower.temperatureC + 273.15);
    row.addProperty("pressure_bara", Math.sqrt(lower.pressureBara * upper.pressureBara));
    row.addProperty("boundaryKind", "target-connected-gas-oil-coexistence-transition");
    row.addProperty("relativePressureBracket", upper.pressureBara / lower.pressureBara - 1.0);
    row.add("phasesBelow", toJsonArray(lower.phases));
    row.add("phasesAbove", toJsonArray(upper.phases));
    return new BoundaryRefinement(row, flashCount, 0);
  }

  private static PhaseMapState evaluatePhaseMapState(SystemInterface template, double temperatureC,
      double pressureBara) {
    SystemInterface point = template.clone();
    point.setTemperature(temperatureC + 273.15);
    point.setPressure(pressureBara);
    runTpFlash(point, false);
    List<String> phases = new ArrayList<String>();
    boolean gas = false;
    boolean oil = false;
    boolean aqueous = false;
    for (int phaseIndex = 0; phaseIndex < point.getNumberOfPhases(); phaseIndex++) {
      String phaseName = point.getPhase(phaseIndex).getPhaseTypeName();
      String normalized = phaseName.toLowerCase(Locale.ROOT);
      if (!phases.contains(phaseName)) {
        phases.add(phaseName);
      }
      gas = gas || normalized.contains("gas") || normalized.contains("vapour") || normalized.contains("vapor");
      aqueous = aqueous || normalized.contains("aqueous") || normalized.contains("water");
      oil = oil || normalized.contains("oil")
          || (normalized.contains("liquid") && !normalized.contains("aqueous") && !normalized.contains("water"));
    }
    return new PhaseMapState(temperatureC, pressureBara, true, gas, oil, aqueous, phases);
  }

  private static JsonObject phaseMapSample(PhaseMapState state, boolean selectedHydrocarbonComponent) {
    JsonObject sample = new JsonObject();
    sample.addProperty("temperature_C", state.temperatureC);
    sample.addProperty("temperature_K", state.temperatureC + 273.15);
    sample.addProperty("pressure_bara", state.pressureBara);
    sample.addProperty("region", state.region());
    sample.addProperty("hydrocarbonTwoPhase", state.hydrocarbonTwoPhase);
    sample.addProperty("selectedHydrocarbonComponent", selectedHydrocarbonComponent);
    sample.add("phases", toJsonArray(state.phases));
    return sample;
  }

  private static JsonArray toJsonArray(List<String> values) {
    JsonArray output = new JsonArray();
    for (String value : values) {
      output.add(value);
    }
    return output;
  }

  private static final class PhaseMapState {
    private final double temperatureC;
    private final double pressureBara;
    private final boolean valid;
    private final boolean gas;
    private final boolean oil;
    private final boolean aqueous;
    private final boolean hydrocarbonTwoPhase;
    private final List<String> phases;

    private PhaseMapState(double temperatureC, double pressureBara, boolean valid, boolean gas, boolean oil,
        boolean aqueous, List<String> phases) {
      this.temperatureC = temperatureC;
      this.pressureBara = pressureBara;
      this.valid = valid;
      this.gas = gas;
      this.oil = oil;
      this.aqueous = aqueous;
      this.hydrocarbonTwoPhase = gas && oil;
      this.phases = phases;
    }

    private static PhaseMapState failed(double temperatureC, double pressureBara) {
      return new PhaseMapState(temperatureC, pressureBara, false, false, false, false, new ArrayList<String>());
    }

    private String region() {
      if (gas && oil && aqueous) {
        return "gas-oil-aqueous";
      }
      if (gas && oil) {
        return "gas-oil";
      }
      if (gas && aqueous) {
        return "gas-aqueous";
      }
      if (oil && aqueous) {
        return "oil-aqueous";
      }
      if (gas) {
        return "gas";
      }
      if (oil) {
        return "oil";
      }
      if (aqueous) {
        return "aqueous";
      }
      return "unclassified";
    }
  }

  private static final class PhaseMapAudit {
    private final JsonArray boundary;
    private final JsonArray samples;
    private final List<String> warnings;
    private final int flashCount;
    private final int failedFlashCount;
    private final int temperaturePointCount;
    private final int pressurePointCount;
    private final double temperatureMinimumC;
    private final double temperatureMaximumC;
    private final double pressureMinimumBara;
    private final double pressureMaximumBara;
    private final int hydrocarbonComponentCount;
    private final int selectedHydrocarbonCellCount;
    private final int excludedHydrocarbonCellCount;

    private PhaseMapAudit(JsonArray boundary, JsonArray samples, List<String> warnings, int flashCount,
        int failedFlashCount, int temperaturePointCount, int pressurePointCount, double temperatureMinimumC,
        double temperatureMaximumC, double pressureMinimumBara, double pressureMaximumBara,
        int hydrocarbonComponentCount, int selectedHydrocarbonCellCount, int excludedHydrocarbonCellCount) {
      this.boundary = boundary;
      this.samples = samples;
      this.warnings = warnings;
      this.flashCount = flashCount;
      this.failedFlashCount = failedFlashCount;
      this.temperaturePointCount = temperaturePointCount;
      this.pressurePointCount = pressurePointCount;
      this.temperatureMinimumC = temperatureMinimumC;
      this.temperatureMaximumC = temperatureMaximumC;
      this.pressureMinimumBara = pressureMinimumBara;
      this.pressureMaximumBara = pressureMaximumBara;
      this.hydrocarbonComponentCount = hydrocarbonComponentCount;
      this.selectedHydrocarbonCellCount = selectedHydrocarbonCellCount;
      this.excludedHydrocarbonCellCount = excludedHydrocarbonCellCount;
    }
  }

  private static final class HydrocarbonComponentSelection {
    private final boolean[][] selected;
    private final int componentCount;
    private final int selectedCellCount;
    private final int excludedCellCount;

    private HydrocarbonComponentSelection(boolean[][] selected, int componentCount, int selectedCellCount,
        int excludedCellCount) {
      this.selected = selected;
      this.componentCount = componentCount;
      this.selectedCellCount = selectedCellCount;
      this.excludedCellCount = excludedCellCount;
    }
  }

  private static final class BoundaryRefinement {
    private final JsonObject row;
    private final int flashCount;
    private final int failedFlashCount;

    private BoundaryRefinement(JsonObject row, int flashCount, int failedFlashCount) {
      this.row = row;
      this.flashCount = flashCount;
      this.failedFlashCount = failedFlashCount;
    }
  }

  private static boolean boundaryIsContinuous(JsonArray boundary) {
    for (int index = 1; index < boundary.size(); index++) {
      double previous = boundary.get(index - 1).getAsJsonObject().get("pressure_bara").getAsDouble();
      double current = boundary.get(index).getAsJsonObject().get("pressure_bara").getAsDouble();
      if (!Double.isFinite(previous) || !Double.isFinite(current) || previous <= 0.0 || current <= 0.0
          || current < previous * MIN_CONTINUOUS_PRESSURE_RATIO) {
        return false;
      }
    }
    return true;
  }

  private static boolean boundaryHasFinitePositiveValues(JsonArray boundary) {
    for (JsonElement element : boundary) {
      JsonObject point = element.getAsJsonObject();
      double temperature = point.get("temperature_K").getAsDouble();
      double pressure = point.get("pressure_bara").getAsDouble();
      if (!Double.isFinite(temperature) || !Double.isFinite(pressure) || temperature < 50.0 || pressure <= 0.0
          || pressure > MAX_ENGINEERING_ENVELOPE_PRESSURE_BARA) {
        return false;
      }
    }
    return boundary.size() > 0;
  }

  private static double maximumBoundaryPressure(JsonArray boundary) {
    double maximum = 0.0;
    for (JsonElement element : boundary) {
      maximum = Math.max(maximum, element.getAsJsonObject().get("pressure_bara").getAsDouble());
    }
    return maximum;
  }

  private static double[] requiredNumberArray(JsonObject input, String name, double offset) {
    if (input == null || !input.has(name) || !input.get(name).isJsonArray()) {
      throw new IllegalArgumentException(name + " must be a numeric array");
    }
    JsonArray values = input.getAsJsonArray(name);
    if (values.size() < 2) {
      throw new IllegalArgumentException(name + " must contain at least two values");
    }
    double[] result = new double[values.size()];
    for (int index = 0; index < values.size(); index++) {
      if (!values.get(index).isJsonPrimitive() || !values.get(index).getAsJsonPrimitive().isNumber()) {
        throw new IllegalArgumentException(name + " must contain only numbers");
      }
      result[index] = values.get(index).getAsDouble() + offset;
      if (!Double.isFinite(result[index])) {
        throw new IllegalArgumentException(name + " contains a non-finite value");
      }
    }
    return result;
  }

  private static State continuationState(JsonObject input, String name) {
    if (input == null || !input.has(name) || !input.get(name).isJsonObject()) {
      throw new IllegalArgumentException(name + " must be a continuation-state object");
    }
    return continuationStateValue(input.getAsJsonObject(name), name);
  }

  private static List<State> continuationStates(JsonArray serializedStates, String name) {
    if (serializedStates == null || serializedStates.size() == 0 || serializedStates.size() > 400) {
      throw new IllegalArgumentException(name + " must contain 1-400 continuation states");
    }
    List<State> states = new ArrayList<State>();
    for (int index = 0; index < serializedStates.size(); index++) {
      if (!serializedStates.get(index).isJsonObject()) {
        throw new IllegalArgumentException(name + " must contain only continuation-state objects");
      }
      states.add(continuationStateValue(serializedStates.get(index).getAsJsonObject(), name + "[" + index + "]"));
    }
    return states;
  }

  private static JsonObject checkpointChainDiagnostic(HydrocarbonWaterContinuationChainAssembler.Result result) {
    JsonObject output = new JsonObject();
    output.addProperty("singleConnectedChain", result.isSingleConnectedChain());
    output.addProperty("chainCount", result.getChains().size());
    output.addProperty("joinCount", result.getJoins().size());
    JsonArray violations = new JsonArray();
    for (String violation : result.getViolations()) {
      violations.add(violation);
    }
    output.add("violations", violations);
    JsonArray chains = new JsonArray();
    for (HydrocarbonWaterContinuationChainAssembler.Chain chain : result.getChains()) {
      JsonObject row = new JsonObject();
      row.addProperty("stateCount", chain.getStates().size());
      JsonArray sources = new JsonArray();
      for (String source : chain.getSourceIdentifiers()) {
        sources.add(source);
      }
      row.add("sourceIdentifiers", sources);
      chains.add(row);
    }
    output.add("chains", chains);
    JsonArray joins = new JsonArray();
    for (HydrocarbonWaterContinuationChainAssembler.Join join : result.getJoins()) {
      JsonObject row = new JsonObject();
      row.addProperty("identifier", join.getIdentifier());
      row.addProperty("firstSegmentIdentifier", join.getFirstSegmentIdentifier());
      row.addProperty("secondSegmentIdentifier", join.getSecondSegmentIdentifier());
      row.addProperty("overlapStateCount", join.getOverlapStateCount());
      row.addProperty("maximumPtDistance", join.getMaximumPtDistance());
      row.addProperty("maximumCompositionDistance", join.getMaximumCompositionDistance());
      row.addProperty("reversedFirst", join.isReversedFirst());
      row.addProperty("reversedSecond", join.isReversedSecond());
      joins.add(row);
    }
    output.add("joins", joins);
    return output;
  }

  private static State continuationStateValue(JsonObject state, String name) {
    if (state == null) {
      throw new IllegalArgumentException(name + " must be a continuation-state object");
    }
    double[] phaseZero = requiredNumberArray(state, "phaseZeroComposition", 0.0);
    double[] phaseOne = requiredNumberArray(state, "phaseOneComposition", 0.0);
    double[] incipient = requiredNumberArray(state, "incipientComposition", 0.0);
    if (phaseZero.length != phaseOne.length || phaseZero.length != incipient.length) {
      throw new IllegalArgumentException(name + " phase-composition lengths must match");
    }
    double[] solverVariables = null;
    if (state.has("solverVariables") && state.get("solverVariables").isJsonArray()) {
      solverVariables = requiredNumberArray(state, "solverVariables", 0.0);
      if (solverVariables.length != 2 * phaseZero.length + 3) {
        throw new IllegalArgumentException(name + " solverVariables must hold 2 * componentCount + 3 entries");
      }
    }
    return State.create(CandidatePhase.valueOf(state.get("retainedPhaseZero").getAsString()),
        CandidatePhase.valueOf(state.get("retainedPhaseOne").getAsString()),
        CandidatePhase.valueOf(state.get("incipientPhase").getAsString()), requiredNumber(state, "temperatureK"),
        requiredPositiveNumber(state, "pressureBara"), requiredNumber(state, "beta"), phaseZero, phaseOne, incipient,
        solverVariables);
  }

  private static SystemInterface buildSystem(JsonObject composition, JsonObject fluidDefinition, String requestedModel,
      boolean reactive, double temperatureK, double pressureBara) {
    JsonObject definition = fluidDefinition == null ? null : fluidDefinition.deepCopy();
    SystemInterface system;
    if (definition == null) {
      String model = normalizedModel(requestedModel, reactive);
      if ("PR".equals(model)) {
        system = new SystemPrEos(temperatureK, pressureBara);
      } else if ("SRK".equals(model)) {
        system = new SystemSrkEos(temperatureK, pressureBara);
      } else if ("CPA".equals(model)) {
        system = new SystemSrkCPAstatoil(temperatureK, pressureBara);
      } else {
        system = new SystemElectrolyteCPAstatoil(temperatureK, pressureBara);
      }
    } else {
      if (reactive) {
        definition.addProperty("eos", "Electrolyte-CPA");
      }
      system = JsonFluidReadWrite.readString(definition.toString());
      system.setTemperature(temperatureK);
      system.setPressure(pressureBara);
    }

    double heavyFraction = fraction(composition, "heavyOil");
    if (definition == null && heavyFraction > 1.0e-12) {
      throw new IllegalArgumentException("heavyOil requires a complete PVTsim fluid definition");
    }

    Map<String, Double> explicit = new LinkedHashMap<String, Double>();
    addExplicit(explicit, "water", fraction(composition, "water"));
    addExplicit(explicit, "CO2", fraction(composition, "CO2"));
    addExplicit(explicit, "nitrogen", fraction(composition, "nitrogen"));
    addExplicit(explicit, "ammonia", fraction(composition, "ammonia"));

    for (Map.Entry<String, Double> entry : explicit.entrySet()) {
      if (entry.getValue() > 0.0 && !system.hasComponent(entry.getKey())) {
        system.addComponent(entry.getKey(), TRACE_MOLES);
      }
    }
    if (reactive) {
      for (String species : REACTIVE_SPECIES) {
        if (!system.hasComponent(species)) {
          system.addComponent(species, TRACE_MOLES);
        }
      }
    }

    if (reactive) {
      system.chemicalReactionInit();
      system.createDatabase(true);
      system.setMixingRule(10);
      applyImportedKij(system, definition);
    } else if (definition != null) {
      // "polarModel" decides two unrelated things at once: which mixing rule gets built here, and - over in
      // JsonFluidReadWrite - whether pseudo components keep their database names. An explicit "mixingRule" states
      // the first without disturbing the second, which matters because the kij lookups match on component names.
      // Absent the field, behaviour is exactly what it was.
      boolean huronVidal = definition.has("polarModel")
          && "HV".equalsIgnoreCase(definition.get("polarModel").getAsString());
      if (definition.has("mixingRule")) {
        String requested = definition.get("mixingRule").getAsString();
        if ("classic".equalsIgnoreCase(requested)) {
          huronVidal = false;
        } else if ("HV".equalsIgnoreCase(requested)) {
          huronVidal = true;
        } else {
          throw new IllegalArgumentException("unsupported mixingRule '" + requested + "'; expected classic or HV");
        }
      }
      if (huronVidal) {
        system.setMixingRule("HV", "NRTL");
      } else {
        system.setMixingRule(2);
      }
      applyImportedKij(system, definition);
      if (huronVidal) {
        // Huron-Vidal parameters mean nothing to the classic rule; feeding them there would only mislead.
        applyImportedHuronVidalParameters(system, definition);
      }
    } else {
      system.createDatabase(true);
      String model = normalizedModel(requestedModel, false);
      system.setMixingRule("CPA".equals(model) || "Electrolyte-CPA".equals(model) ? 10 : 2);
    }

    if (definition != null && !reactive) {
      JsonFluidReadWrite.applyComponentAlphaParameters(system, definition);
    }

    double[] target = new double[system.getNumberOfComponents()];
    if (definition != null) {
      JsonArray imported = definition.getAsJsonArray("components");
      for (JsonElement element : imported) {
        JsonObject component = element.getAsJsonObject();
        String name = mappedName(component.get("name").getAsString(),
            component.has("isPseudo") && component.get("isPseudo").getAsBoolean());
        int index = componentIndex(system, name);
        if (index >= 0) {
          target[index] += heavyFraction * component.get("moleFraction").getAsDouble();
        }
      }
    }
    for (Map.Entry<String, Double> entry : explicit.entrySet()) {
      int index = componentIndex(system, entry.getKey());
      if (index >= 0) {
        target[index] += entry.getValue();
      }
    }
    if (reactive) {
      for (String species : REACTIVE_SPECIES) {
        int index = componentIndex(system, species);
        if (index >= 0) {
          target[index] = Math.max(target[index], TRACE_MOLES);
        }
      }
    }
    system.setMolarComposition(target);
    system.setMultiPhaseCheck(true);
    return system;
  }

  private static SystemInterface buildBoundarySystem(JsonObject composition, JsonObject fluidDefinition,
      String requestedModel, double temperatureK, double pressureBara) {
    return buildSystem(composition, fluidDefinition, requestedModel, false, temperatureK, pressureBara);
  }

  private static String normalizedModel(String requestedModel, boolean reactive) {
    if (reactive) {
      return "Electrolyte-CPA";
    }
    String normalized = String.valueOf(requestedModel).replaceAll("[^A-Za-z0-9]", "").toUpperCase(Locale.ROOT);
    if (normalized.equals("PR") || normalized.contains("PENGROBINSON") || normalized.contains("PR78")) {
      return "PR";
    }
    if (normalized.equals("SRK") || normalized.contains("SOAVEREDLICHKWONG")) {
      return "SRK";
    }
    if (normalized.contains("ELECTROLYTECPA")) {
      return "Electrolyte-CPA";
    }
    if (normalized.contains("CPA")) {
      return "CPA";
    }
    throw new IllegalArgumentException("Unsupported field-fluid model: " + requestedModel
        + "; supported models are SRK, PR, CPA, and Electrolyte-CPA");
  }

  private static ReactiveFlashAudit runTpFlash(SystemInterface system, boolean reactive) {
    if (reactive) {
      system.setMaxNumberOfPhases(3);
      ReactiveMultiphaseTPflash flash = new ReactiveMultiphaseTPflash(system);
      flash.setMaxNumberOfPhases(3);
      flash.run();
      if (!flash.isConverged()) {
        throw new IllegalStateException("Reactive multiphase TP flash did not converge");
      }
      return new ReactiveFlashAudit(true, flash.getTotalIterations(), flash.getNumberOfReactions());
    }
    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();
    if (operations.getOperation() instanceof TPflash) {
      TPflash flash = (TPflash) operations.getOperation();
      if (!flash.isLastMultiphaseSolveAccepted()) {
        throw new IllegalStateException("Non-reactive multiphase TP flash failed with "
            + flash.getLastMultiphaseSolveStatus() + ": " + flash.getLastMultiphaseSolveMessage());
      }
    }
    return new ReactiveFlashAudit(true, 0, 0);
  }

  private static final class ReactiveFlashAudit {
    private final boolean converged;
    private final int iterations;
    private final int reactions;

    private ReactiveFlashAudit(boolean converged, int iterations, int reactions) {
      this.converged = converged;
      this.iterations = iterations;
      this.reactions = reactions;
    }
  }

  private static JsonObject buildSpeciation(SystemInterface system, boolean reactive) {
    JsonObject output = new JsonObject();
    output.addProperty("enabled", reactive);
    if (!reactive) {
      return output;
    }
    PhaseInterface aqueous = null;
    for (int index = 0; index < system.getNumberOfPhases(); index++) {
      if ("aqueous".equalsIgnoreCase(system.getPhase(index).getPhaseTypeName())) {
        aqueous = system.getPhase(index);
        break;
      }
    }
    if (aqueous == null) {
      output.addProperty("aqueousPhasePresent", false);
      output.addProperty("chargeBalance", 0.0);
      return output;
    }
    output.addProperty("aqueousPhasePresent", true);
    double charge = 0.0;
    JsonObject moleFractions = new JsonObject();
    for (int index = 0; index < aqueous.getNumberOfComponents(); index++) {
      double value = aqueous.getComponent(index).getx();
      double ionicCharge = aqueous.getComponent(index).getIonicCharge();
      charge += value * ionicCharge;
      if (value > 1.0e-16) {
        moleFractions.addProperty(aqueous.getComponent(index).getComponentName(), value);
      }
    }
    output.addProperty("chargeBalance", charge);
    try {
      output.addProperty("pH", aqueous.getpH());
    } catch (Exception ignored) {
      output.add("pH", null);
    }
    output.add("aqueousMoleFractions", moleFractions);
    output.addProperty("chargeBalancePassed", Math.abs(charge) <= 1.0e-6);
    return output;
  }

  private static void applyImportedKij(SystemInterface system, JsonObject definition) {
    if (definition == null || !definition.has("binaryInteractionCoefficients")) {
      return;
    }
    JsonArray coefficients = definition.getAsJsonArray("binaryInteractionCoefficients");
    clearUnspecifiedImportedKij(system, definition, coefficients);
    for (JsonElement element : coefficients) {
      JsonObject coefficient = element.getAsJsonObject();
      int first = componentIndex(system, mappedName(coefficient.get("i").getAsString(), false));
      int second = componentIndex(system, mappedName(coefficient.get("j").getAsString(), false));
      if (first < 0 || second < 0) {
        first = componentIndex(system, coefficient.get("i").getAsString());
        second = componentIndex(system, coefficient.get("j").getAsString());
      }
      if (first < 0 || second < 0) {
        continue;
      }
      double kij = coefficient.get("kij").getAsDouble();
      for (int phase = 0; phase < system.getMaxNumberOfPhases(); phase++) {
        if (system.getPhase(phase) instanceof PhaseEosInterface) {
          PhaseEosInterface eosPhase = (PhaseEosInterface) system.getPhase(phase);
          eosPhase.getEosMixingRule().setBinaryInteractionParameter(first, second, kij);
          eosPhase.getEosMixingRule().setBinaryInteractionParameter(second, first, kij);
        }
      }
    }
  }

  /**
   * Zero every interaction parameter between imported components that the source fluid does not specify.
   *
   * <p>
   * A PVTsim fluid writes a Kij only where it is non-zero; every pair it omits is zero by definition. NeqSim fills
   * unspecified pairs from its own database and correlations, so without this the imported fluid quietly runs on a
   * different parameter set than the one it was tuned as. Pairs involving components the caller added on top of the
   * import - water, ammonia - are left at the NeqSim defaults, which is what the reported limitations describe.
   * </p>
   *
   * @param system the built system
   * @param definition the imported fluid definition
   * @param coefficients the declared interaction coefficients
   */
  private static void clearUnspecifiedImportedKij(SystemInterface system, JsonObject definition,
      JsonArray coefficients) {
    java.util.Set<String> declared = new java.util.HashSet<String>();
    for (JsonElement element : coefficients) {
      JsonObject coefficient = element.getAsJsonObject();
      declared.add(pairKey(coefficient.get("i").getAsString(), coefficient.get("j").getAsString()));
    }

    Map<Integer, String> importedIndices = new LinkedHashMap<Integer, String>();
    for (JsonElement element : definition.getAsJsonArray("components")) {
      JsonObject component = element.getAsJsonObject();
      String sourceName = component.get("name").getAsString();
      int index = componentIndex(system,
          mappedName(sourceName, component.has("isPseudo") && component.get("isPseudo").getAsBoolean()));
      if (index < 0) {
        index = componentIndex(system, sourceName);
      }
      if (index >= 0) {
        importedIndices.put(Integer.valueOf(index), sourceName);
      }
    }

    for (Map.Entry<Integer, String> first : importedIndices.entrySet()) {
      for (Map.Entry<Integer, String> second : importedIndices.entrySet()) {
        if (first.getKey().intValue() >= second.getKey().intValue()
            || declared.contains(pairKey(first.getValue(), second.getValue()))) {
          continue;
        }
        for (int phase = 0; phase < system.getMaxNumberOfPhases(); phase++) {
          if (system.getPhase(phase) instanceof PhaseEosInterface) {
            PhaseEosInterface eosPhase = (PhaseEosInterface) system.getPhase(phase);
            eosPhase.getEosMixingRule().setBinaryInteractionParameter(first.getKey().intValue(),
                second.getKey().intValue(), 0.0);
            eosPhase.getEosMixingRule().setBinaryInteractionParameter(second.getKey().intValue(),
                first.getKey().intValue(), 0.0);
          }
        }
      }
    }
  }

  private static String pairKey(String first, String second) {
    return first.compareTo(second) <= 0 ? first + "|" + second : second + "|" + first;
  }

  private static void applyImportedHuronVidalParameters(SystemInterface system, JsonObject definition) {
    if (definition == null || !definition.has("huronVidalInteractionParameters")
        || !definition.get("huronVidalInteractionParameters").isJsonArray()) {
      return;
    }
    for (JsonElement element : definition.getAsJsonArray("huronVidalInteractionParameters")) {
      JsonObject parameter = element.getAsJsonObject();
      int first = componentIndex(system, mappedName(parameter.get("i").getAsString(), false));
      int second = componentIndex(system, mappedName(parameter.get("j").getAsString(), false));
      if (first < 0 || second < 0) {
        first = componentIndex(system, parameter.get("i").getAsString());
        second = componentIndex(system, parameter.get("j").getAsString());
      }
      if (first < 0 || second < 0) {
        continue;
      }
      double forward = parameter.get("forward").getAsDouble();
      double reverse = parameter.get("reverse").getAsDouble();
      if (!Double.isFinite(forward) || !Double.isFinite(reverse)) {
        continue;
      }
      for (int phaseIndex = 0; phaseIndex < system.getMaxNumberOfPhases(); phaseIndex++) {
        if (!(system.getPhase(phaseIndex) instanceof PhaseEosInterface)) {
          continue;
        }
        PhaseEosInterface phase = (PhaseEosInterface) system.getPhase(phaseIndex);
        if (phase.getEosMixingRule() instanceof HVMixingRulesInterface) {
          HVMixingRulesInterface rule = (HVMixingRulesInterface) phase.getEosMixingRule();
          rule.setHVDijParameter(first, second, forward);
          rule.setHVDijParameter(second, first, reverse);
        }
      }
    }
  }

  /**
   * Traces the hydrocarbon dew branch of a water-containing feed, pressure by pressure.
   *
   * <p>
   * The water stays in the mixture: it is part of the fluid whose boundary is being asked about, and removing it traces
   * a different mixture. At each pressure the incipient-phase solve is given the full temperature window, because a
   * window centred on the previous point stops the solver bracketing a root at all. Only the dew branch comes out of
   * this - the low-temperature branch of the same curve is a second root the solve does not return, and the caller is
   * told so rather than shown a filled-in guess.
   * </p>
   *
   * @param composition the requested overall composition, water included
   * @param fluidDefinition the imported characterized fluid
   * @param requestedModel the equation of state
   * @param minimumPressureBara lowest pressure to trace
   * @param maximumPressureBara highest pressure to trace
   * @param stateTemperatureK the state-point temperature, used to seed each build
   * @return boundary rows, possibly empty
   */
  private static JsonArray traceHydrocarbonDewBranchWithWater(JsonObject composition, JsonObject fluidDefinition,
      String requestedModel, double minimumPressureBara, double maximumPressureBara, double stateTemperatureK) {
    JsonArray rows = new JsonArray();
    double floor = Math.max(minimumPressureBara, 1.0e-3);
    double ceiling = Math.max(floor * 1.001, maximumPressureBara);
    for (int index = 0; index < DEW_BRANCH_PRESSURE_POINTS; index++) {
      double fraction = (double) index / (DEW_BRANCH_PRESSURE_POINTS - 1);
      double pressure = floor * Math.pow(ceiling / floor, fraction);
      try {
        SystemInterface system = buildBoundarySystem(composition, fluidDefinition, requestedModel, stateTemperatureK,
            pressure);
        IncipientPhaseBoundaryPointSolver.Result result = new IncipientPhaseBoundaryPointSolver(system,
            IncipientPhaseStabilityAnalyzer.CandidatePhase.GAS, IncipientPhaseStabilityAnalyzer.CandidatePhase.OIL)
            .solve(pressure, DEW_BRANCH_MINIMUM_TEMPERATURE_K, DEW_BRANCH_MAXIMUM_TEMPERATURE_K);
        if (!result.isConverged()) {
          continue;
        }
        double temperatureK = result.getTemperatureK();
        if (!Double.isFinite(temperatureK) || temperatureK <= 0.0) {
          continue;
        }
        JsonObject row = new JsonObject();
        row.addProperty("temperature_C", temperatureK - 273.15);
        row.addProperty("temperature_K", temperatureK);
        row.addProperty("pressure_bara", pressure);
        row.addProperty("curve", "hydrocarbon-dew-with-water");
        rows.add(row);
      } catch (Exception ignored) {
        // A pressure with no reachable root is simply absent; it is not a failure of the whole trace.
      }
    }
    return rows;
  }

  /**
   * Build the composition the hydrocarbon phase envelope is traced on.
   *
   * <p>
   * Water and ammonia are aqueous-phase species: including them turns the two-phase boundary into an aqueous dew line
   * instead of the hydrocarbon envelope. The remaining fractions are renormalized so the traced fluid is the
   * hydrocarbon sub-composition at its own basis.
   * </p>
   *
   * @param composition the requested overall composition
   * @return the hydrocarbon sub-composition, renormalized to sum to the same total
   */
  private static JsonObject hydrocarbonBasisComposition(JsonObject composition) {
    JsonObject hydrocarbon = new JsonObject();
    double retained = 0.0;
    for (Map.Entry<String, JsonElement> entry : composition.entrySet()) {
      if (AQUEOUS_BASIS_SPECIES.contains(entry.getKey())) {
        continue;
      }
      double value = entry.getValue().getAsDouble();
      if (value > 0.0) {
        hydrocarbon.addProperty(entry.getKey(), value);
        retained += value;
      }
    }
    if (retained <= 0.0) {
      throw new IllegalArgumentException(
          "a hydrocarbon phase envelope needs a non-aqueous fraction; the requested composition is aqueous only");
    }
    for (Map.Entry<String, JsonElement> entry : new ArrayList<Map.Entry<String, JsonElement>>(hydrocarbon.entrySet())) {
      hydrocarbon.addProperty(entry.getKey(), entry.getValue().getAsDouble() / retained);
    }
    return hydrocarbon;
  }

  /**
   * Sum of the aqueous-phase species excluded from the hydrocarbon envelope basis.
   *
   * @param composition the requested overall composition
   * @return the excluded mole fraction, normalized by the composition total
   */
  private static double aqueousFraction(JsonObject composition) {
    double aqueous = 0.0;
    double total = 0.0;
    for (Map.Entry<String, JsonElement> entry : composition.entrySet()) {
      double value = entry.getValue().getAsDouble();
      if (value <= 0.0) {
        continue;
      }
      total += value;
      if (AQUEOUS_BASIS_SPECIES.contains(entry.getKey())) {
        aqueous += value;
      }
    }
    return total > 0.0 ? aqueous / total : 0.0;
  }

  private static void addExplicit(Map<String, Double> explicit, String name, double value) {
    if (value > 0.0) {
      explicit.put(name, value);
    }
  }

  private static double fraction(JsonObject composition, String name) {
    return composition.has(name) ? composition.get(name).getAsDouble() : 0.0;
  }

  private static int componentIndex(SystemInterface system, String name) {
    for (int index = 0; index < system.getNumberOfComponents(); index++) {
      String actualName = system.getComponent(index).getComponentName();
      if (name.equalsIgnoreCase(actualName) || (name + "_PC").equalsIgnoreCase(actualName)
          || (actualName.endsWith("_PC") && name.equalsIgnoreCase(actualName.substring(0, actualName.length() - 3)))) {
        return index;
      }
    }
    return -1;
  }

  private static String mappedName(String name, boolean pseudo) {
    if (pseudo) {
      return name;
    }
    if ("C1".equalsIgnoreCase(name))
      return "methane";
    if ("C2".equalsIgnoreCase(name))
      return "ethane";
    if ("C3".equalsIgnoreCase(name))
      return "propane";
    if ("iC4".equalsIgnoreCase(name))
      return "i-butane";
    if ("nC4".equalsIgnoreCase(name) || "C4".equalsIgnoreCase(name))
      return "n-butane";
    if ("iC5".equalsIgnoreCase(name))
      return "i-pentane";
    if ("nC5".equalsIgnoreCase(name) || "C5".equalsIgnoreCase(name))
      return "n-pentane";
    if ("C6".equalsIgnoreCase(name))
      return "n-hexane";
    if ("N2".equalsIgnoreCase(name))
      return "nitrogen";
    if ("H2O".equalsIgnoreCase(name))
      return "water";
    if ("NH3".equalsIgnoreCase(name))
      return "ammonia";
    return name;
  }

  private static double requiredNumber(JsonObject input, String name) {
    if (!input.has(name))
      throw new IllegalArgumentException(name + " is required");
    double value = input.get(name).getAsDouble();
    if (!Double.isFinite(value))
      throw new IllegalArgumentException(name + " must be finite");
    return value;
  }

  private static double requiredPositiveNumber(JsonObject input, String name) {
    double value = requiredNumber(input, name);
    if (value <= 0.0)
      throw new IllegalArgumentException(name + " must be positive");
    return value;
  }

  private static double optionalNumber(JsonObject input, String name, double fallback) {
    return input.has(name) ? requiredNumber(input, name) : fallback;
  }

  private static int optionalInteger(JsonObject input, String name, int fallback) {
    return input.has(name) ? input.get(name).getAsInt() : fallback;
  }
}
