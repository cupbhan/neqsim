package neqsim.mcp.runners;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.AcceptanceCriteria;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.CriticalPointMetrics;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.Curve;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.DirectionalMetrics;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.Metrics;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.ModelFingerprint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.Point;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryRegression.ReferenceSource;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.BoundaryDefinition;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;

/** JSON adapter for auditable paper/PVTsim same-topology hydrocarbon-water boundary regression. */
public final class HydrocarbonWaterBoundaryRegressionRunner {
  private HydrocarbonWaterBoundaryRegressionRunner() {
  }

  /**
   * Runs one model-, composition-, and topology-locked boundary comparison.
   *
   * @param request regression request containing {@code reference}, {@code candidate}, and optional frozen
   * {@code criteria}
   * @return machine-readable contract, metric, and threshold verdict
   */
  public static JsonObject run(JsonObject request) {
    if (request == null || !request.has("reference") || !request.get("reference").isJsonObject()
        || !request.has("candidate") || !request.get("candidate").isJsonObject()) {
      throw new IllegalArgumentException("boundary regression requires reference and candidate objects");
    }
    Curve reference = curve(request.getAsJsonObject("reference"), "reference");
    Curve candidate = curve(request.getAsJsonObject("candidate"), "candidate");
    AcceptanceCriteria criteria = request.has("criteria") && request.get("criteria").isJsonObject()
        ? criteria(request.getAsJsonObject("criteria"))
        : null;
    HydrocarbonWaterBoundaryRegression.Result result = new HydrocarbonWaterBoundaryRegression().compare(reference,
        candidate, criteria);
    return serialize(reference, candidate, result);
  }

  private static Curve curve(JsonObject json, String label) {
    String identifier = requiredString(json, "identifier", label);
    ReferenceSource source = ReferenceSource.valueOf(requiredString(json, "source", label).toUpperCase(Locale.ROOT));
    JsonObject model = requiredObject(json, "modelFingerprint", label);
    ModelFingerprint fingerprint = new ModelFingerprint(requiredString(model, "equationOfState", label + " model"),
        requiredString(model, "volumeTranslation", label + " model"),
        requiredString(model, "mixingRule", label + " model"),
        requiredString(model, "compositionFingerprint", label + " model"),
        requiredString(model, "parameterFingerprint", label + " model"));
    BoundaryDefinition topology = topology(requiredString(json, "boundaryTopology", label));
    boolean topologyResolved = json.has("topologyResolved") && json.get("topologyResolved").getAsBoolean();
    String topologyEvidence = requiredString(json, "topologyEvidence", label);
    if (!json.has("points") || !json.get("points").isJsonArray()) {
      throw new IllegalArgumentException(label + " points must be an array");
    }
    List<Point> points = points(json.getAsJsonArray("points"), label + " points");
    Point criticalPoint = json.has("criticalPoint") && json.get("criticalPoint").isJsonObject()
        ? point(json.getAsJsonObject("criticalPoint"), label + " critical point")
        : null;
    return new Curve(identifier, source, fingerprint, topology, topologyResolved, topologyEvidence, points,
        criticalPoint);
  }

  private static AcceptanceCriteria criteria(JsonObject json) {
    return new AcceptanceCriteria(requiredNumber(json, "minimumReferenceCoverageFraction", "criteria"),
        requiredNumber(json, "minimumCandidateCoverageFraction", "criteria"),
        requiredNumber(json, "maximumRelativeTemperatureDifference", "criteria"),
        requiredNumber(json, "maximumRmsRelativeTemperatureDifference", "criteria"),
        requiredNumber(json, "maximumP95RelativeTemperatureDifference", "criteria"),
        requiredNumber(json, "maximumSymmetricRelativeTemperatureDifference", "criteria"),
        json.has("criticalPointRequired") && json.get("criticalPointRequired").getAsBoolean(),
        requiredNumber(json, "maximumCriticalTemperatureDifference", "criteria"),
        requiredNumber(json, "maximumCriticalPressureDifference", "criteria"));
  }

  private static BoundaryDefinition topology(String value) {
    String normalized = value.trim().toUpperCase(Locale.ROOT).replace('-', '_').replace('>', '_');
    if (normalized.equals("OW_TO_GOW") || normalized.equals("OIL_AQUEOUS_TO_GAS_OIL_AQUEOUS")) {
      return new BoundaryDefinition(Region.OIL_AQUEOUS, Region.GAS_OIL_AQUEOUS);
    }
    if (normalized.equals("GW_TO_GOW") || normalized.equals("GAS_AQUEOUS_TO_GAS_OIL_AQUEOUS")) {
      return new BoundaryDefinition(Region.GAS_AQUEOUS, Region.GAS_OIL_AQUEOUS);
    }
    if (normalized.equals("GO_TO_GOW") || normalized.equals("GAS_OIL_TO_GAS_OIL_AQUEOUS")) {
      return new BoundaryDefinition(Region.GAS_OIL, Region.GAS_OIL_AQUEOUS);
    }
    throw new IllegalArgumentException("unsupported two-to-three-phase boundary topology " + value);
  }

  private static List<Point> points(JsonArray array, String label) {
    List<Point> points = new ArrayList<Point>();
    for (int index = 0; index < array.size(); index++) {
      if (!array.get(index).isJsonObject()) {
        throw new IllegalArgumentException(label + "[" + index + "] must be an object");
      }
      points.add(point(array.get(index).getAsJsonObject(), label + "[" + index + "]"));
    }
    return points;
  }

  private static Point point(JsonObject json, String label) {
    double temperatureK;
    if (json.has("temperatureK")) {
      temperatureK = requiredNumber(json, "temperatureK", label);
    } else if (json.has("temperatureC")) {
      temperatureK = requiredNumber(json, "temperatureC", label) + 273.15;
    } else if (json.has("temperature_C")) {
      temperatureK = requiredNumber(json, "temperature_C", label) + 273.15;
    } else {
      throw new IllegalArgumentException(label + " requires temperatureK or temperatureC");
    }
    double pressureBara = json.has("pressureBara") ? requiredNumber(json, "pressureBara", label)
        : requiredNumber(json, "pressure_bara", label);
    return new Point(temperatureK, pressureBara);
  }

  private static JsonObject serialize(Curve reference, Curve candidate,
      HydrocarbonWaterBoundaryRegression.Result result) {
    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "bidirectional-same-pressure-temperature-regression");
    output.addProperty("referenceIdentifier", reference.getIdentifier());
    output.addProperty("referenceSource", reference.getSource().name());
    output.addProperty("candidateIdentifier", candidate.getIdentifier());
    output.addProperty("candidateSource", candidate.getSource().name());
    output.addProperty("boundaryTopology", reference.getBoundaryDefinition().getCode());
    output.addProperty("referenceTopologyEvidence", reference.getTopologyEvidence());
    output.addProperty("candidateTopologyEvidence", candidate.getTopologyEvidence());
    output.addProperty("comparisonEligible", result.isComparisonEligible());
    output.addProperty("thresholdsFrozen", result.isThresholdsFrozen());
    output.addProperty("benchmarkPending", result.isBenchmarkPending());
    output.addProperty("accepted", result.isAccepted());
    output.addProperty("engineeringEligible", false);
    output.addProperty("releaseRule",
        "benchmark acceptance must be conjoined with the independent internal topology and thermodynamic quality gate");
    JsonArray violations = new JsonArray();
    for (String violation : result.getViolations()) {
      violations.add(violation);
    }
    output.add("violations", violations);
    if (result.getMetrics() != null) {
      output.add("metrics", metrics(result.getMetrics()));
    }
    return output;
  }

  private static JsonObject metrics(Metrics metrics) {
    JsonObject output = new JsonObject();
    output.add("referenceToCandidate", directional(metrics.getReferenceToCandidate()));
    output.add("candidateToReference", directional(metrics.getCandidateToReference()));
    output.addProperty("symmetricMaximumRelativeTemperatureDifference",
        metrics.getSymmetricMaximumRelativeTemperatureDifference());
    CriticalPointMetrics critical = metrics.getCriticalPointMetrics();
    JsonObject criticalJson = new JsonObject();
    criticalJson.addProperty("available", critical.isAvailable());
    if (critical.isAvailable()) {
      criticalJson.addProperty("absoluteTemperatureDifferenceK", critical.getAbsoluteTemperatureDifferenceK());
      criticalJson.addProperty("absolutePressureDifferenceBara", critical.getAbsolutePressureDifferenceBara());
      criticalJson.addProperty("relativeTemperatureDifference", critical.getRelativeTemperatureDifference());
      criticalJson.addProperty("relativePressureDifference", critical.getRelativePressureDifference());
    }
    output.add("criticalPoint", criticalJson);
    return output;
  }

  private static JsonObject directional(DirectionalMetrics metrics) {
    JsonObject output = new JsonObject();
    output.addProperty("sourcePointCount", metrics.getSourcePointCount());
    output.addProperty("matchedPointCount", metrics.getMatchedPointCount());
    output.addProperty("pressureCoverageFraction", metrics.getCoverageFraction());
    output.addProperty("maximumAbsoluteTemperatureDifferenceK", metrics.getMaximumAbsoluteTemperatureDifferenceK());
    output.addProperty("rmsAbsoluteTemperatureDifferenceK", metrics.getRmsAbsoluteTemperatureDifferenceK());
    output.addProperty("p95AbsoluteTemperatureDifferenceK", metrics.getP95AbsoluteTemperatureDifferenceK());
    output.addProperty("maximumRelativeTemperatureDifference", metrics.getMaximumRelativeTemperatureDifference());
    output.addProperty("rmsRelativeTemperatureDifference", metrics.getRmsRelativeTemperatureDifference());
    output.addProperty("p95RelativeTemperatureDifference", metrics.getP95RelativeTemperatureDifference());
    return output;
  }

  private static JsonObject requiredObject(JsonObject json, String name, String label) {
    if (!json.has(name) || !json.get(name).isJsonObject()) {
      throw new IllegalArgumentException(label + " " + name + " must be an object");
    }
    return json.getAsJsonObject(name);
  }

  private static String requiredString(JsonObject json, String name, String label) {
    if (!json.has(name) || json.get(name).isJsonNull()) {
      throw new IllegalArgumentException(label + " " + name + " is required");
    }
    String value = json.get(name).getAsString().trim();
    if (value.isEmpty()) {
      throw new IllegalArgumentException(label + " " + name + " is required");
    }
    return value;
  }

  private static double requiredNumber(JsonObject json, String name, String label) {
    JsonElement element = json.get(name);
    if (element == null || element.isJsonNull()) {
      throw new IllegalArgumentException(label + " " + name + " is required");
    }
    double value = element.getAsDouble();
    if (!Double.isFinite(value)) {
      throw new IllegalArgumentException(label + " " + name + " must be finite");
    }
    return value;
  }
}
