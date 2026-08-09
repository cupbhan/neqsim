package neqsim.mcp.runners;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.BoundaryDefinition;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterTopologyBenchmark;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterTopologyBenchmark.Candidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterTopologyBenchmark.Reference;

/** Reads, verifies, and applies a topology-only hydrocarbon-water paper benchmark manifest. */
public final class HydrocarbonWaterTopologyBenchmarkRunner {
  private static final long MAXIMUM_MANIFEST_BYTES = 5L * 1024L * 1024L;

  private HydrocarbonWaterTopologyBenchmarkRunner() {
  }

  /**
   * Runs a topology-only regression after independently hashing all local source artifacts.
   *
   * @param request request containing an absolute {@code referenceManifestPath} and a candidate topology object
   * @return auditable source-contract and topology verdict
   */
  public static JsonObject run(JsonObject request) {
    if (request == null || !request.has("candidate") || !request.get("candidate").isJsonObject()) {
      throw new IllegalArgumentException("paper topology benchmark requires a candidate object");
    }
    Path manifestPath = Paths.get(requiredString(request, "referenceManifestPath", "request")).toAbsolutePath()
        .normalize();
    JsonObject manifest = readManifest(manifestPath);
    ReferenceAudit audit = auditReference(manifestPath, manifest);
    Reference reference = reference(manifest, audit.isVerified());
    Candidate candidate = candidate(request.getAsJsonObject("candidate"));
    HydrocarbonWaterTopologyBenchmark.Result result = new HydrocarbonWaterTopologyBenchmark().compare(reference,
        candidate);
    return serialize(manifestPath, manifest, audit, reference, candidate, result);
  }

  private static JsonObject readManifest(Path path) {
    try {
      if (!Files.isRegularFile(path) || !Files.isReadable(path)) {
        throw new IllegalArgumentException("reference manifest is not a readable regular file");
      }
      long size = Files.size(path);
      if (size <= 0L || size > MAXIMUM_MANIFEST_BYTES) {
        throw new IllegalArgumentException("reference manifest size is outside the allowed range");
      }
      JsonElement parsed = JsonParser.parseString(new String(Files.readAllBytes(path), StandardCharsets.UTF_8));
      if (!parsed.isJsonObject()) {
        throw new IllegalArgumentException("reference manifest must contain one JSON object");
      }
      return parsed.getAsJsonObject();
    } catch (IOException error) {
      throw new IllegalArgumentException("unable to read reference manifest: " + error.getMessage(), error);
    }
  }

  private static ReferenceAudit auditReference(Path manifestPath, JsonObject manifest) {
    List<String> violations = new ArrayList<String>();
    if (!"hydrocarbon-water-paper-topology-benchmark-v1".equals(optionalString(manifest, "schemaVersion"))) {
      violations.add("MANIFEST_SCHEMA_UNSUPPORTED");
    }
    JsonArray artifactArray = manifest.has("artifacts") && manifest.get("artifacts").isJsonArray()
        ? manifest.getAsJsonArray("artifacts")
        : new JsonArray();
    JsonArray artifactAudit = new JsonArray();
    if (artifactArray.size() == 0) {
      violations.add("SOURCE_ARTIFACTS_MISSING");
    }
    Path base = manifestPath.getParent();
    for (JsonElement element : artifactArray) {
      JsonObject output = new JsonObject();
      if (!element.isJsonObject()) {
        violations.add("SOURCE_ARTIFACT_DECLARATION_INVALID");
        output.addProperty("verified", false);
        artifactAudit.add(output);
        continue;
      }
      JsonObject artifact = element.getAsJsonObject();
      String relativePath = optionalString(artifact, "path");
      String expectedHash = optionalString(artifact, "sha256").toLowerCase(Locale.ROOT);
      Path resolved = base.resolve(relativePath).normalize();
      boolean contained = resolved.startsWith(base);
      boolean readable = contained && Files.isRegularFile(resolved) && Files.isReadable(resolved);
      String actualHash = readable ? sha256(resolved) : "";
      long expectedBytes = optionalLong(artifact, "bytes", -1L);
      long actualBytes = readable ? size(resolved) : -1L;
      boolean verified = artifact.has("verified") && artifact.get("verified").getAsBoolean() && readable
          && expectedHash.matches("[0-9a-f]{64}") && expectedHash.equals(actualHash) && expectedBytes == actualBytes;
      output.addProperty("role", optionalString(artifact, "role"));
      output.addProperty("path", relativePath);
      output.addProperty("containedByManifestDirectory", contained);
      output.addProperty("expectedSha256", expectedHash);
      output.addProperty("actualSha256", actualHash);
      output.addProperty("expectedBytes", expectedBytes);
      output.addProperty("actualBytes", actualBytes);
      output.addProperty("verified", verified);
      artifactAudit.add(output);
      if (!verified) {
        violations.add("SOURCE_ARTIFACT_INTEGRITY_FAILED");
      }
    }

    JsonObject fluid = requiredObject(manifest, "referenceFluid", "manifest");
    JsonArray components = fluid.has("components") && fluid.get("components").isJsonArray()
        ? fluid.getAsJsonArray("components")
        : new JsonArray();
    double compositionSum = 0.0;
    for (JsonElement element : components) {
      if (!element.isJsonObject() || !element.getAsJsonObject().has("moleFraction")) {
        violations.add("REFERENCE_COMPOSITION_INVALID");
        continue;
      }
      double fraction = element.getAsJsonObject().get("moleFraction").getAsDouble();
      if (!Double.isFinite(fraction) || fraction < 0.0) {
        violations.add("REFERENCE_COMPOSITION_INVALID");
      } else {
        compositionSum += fraction;
      }
    }
    if (components.size() == 0 || Math.abs(compositionSum - 1.0) > 1.0e-10) {
      violations.add("REFERENCE_COMPOSITION_NOT_NORMALIZED");
    }

    auditDigitizedCoordinates(manifest, violations);
    return new ReferenceAudit(violations.isEmpty(), manifestSha256(manifestPath), compositionSum, artifactAudit,
        violations);
  }

  private static void auditDigitizedCoordinates(JsonObject manifest, List<String> violations) {
    JsonObject digitization = requiredObject(manifest, "figureDigitization", "manifest");
    JsonObject axes = requiredObject(digitization, "axisCalibration", "figureDigitization");
    JsonObject temperature = requiredObject(axes, "temperature", "axisCalibration");
    JsonObject pressure = requiredObject(axes, "pressure", "axisCalibration");
    double minimumTemperature = requiredNumber(temperature, "physicalMinimumK", "temperature axis");
    double maximumTemperature = requiredNumber(temperature, "physicalMaximumK", "temperature axis");
    double minimumPressure = requiredNumber(pressure, "physicalMinimumBar", "pressure axis");
    double maximumPressure = requiredNumber(pressure, "physicalMaximumBar", "pressure axis");
    JsonObject points = requiredObject(digitization, "points", "figureDigitization");
    JsonArray boundaries = requiredObject(manifest, "topologyReference", "manifest")
        .getAsJsonArray("requiredTwoToThreePhaseBoundaries");
    for (JsonElement boundaryElement : boundaries) {
      String code = requiredString(boundaryElement.getAsJsonObject(), "code", "boundary");
      if (!points.has(code) || !points.get(code).isJsonArray() || points.getAsJsonArray(code).size() == 0) {
        violations.add("DIGITIZED_BOUNDARY_POINTS_MISSING");
        continue;
      }
      for (JsonElement pointElement : points.getAsJsonArray(code)) {
        if (!pointElement.isJsonObject()) {
          violations.add("DIGITIZED_BOUNDARY_POINT_INVALID");
          continue;
        }
        JsonObject point = pointElement.getAsJsonObject();
        double pointTemperature = requiredNumber(point, "temperatureK", code + " point");
        double pointPressure = requiredNumber(point, "pressureBar", code + " point");
        if (pointTemperature < minimumTemperature || pointTemperature > maximumTemperature
            || pointPressure < minimumPressure || pointPressure > maximumPressure) {
          violations.add("DIGITIZED_BOUNDARY_POINT_OUTSIDE_AXIS");
        }
      }
    }
  }

  private static Reference reference(JsonObject manifest, boolean sourceContractVerified) {
    JsonObject policy = requiredObject(manifest, "comparisonPolicy", "manifest");
    boolean topologyOnly = "TOPOLOGY_ONLY".equals(optionalString(policy, "scope").toUpperCase(Locale.ROOT));
    JsonObject topology = requiredObject(manifest, "topologyReference", "manifest");
    JsonArray requiredRegionArray = topology.has("requiredRegionsForTwoToThreePhaseBenchmark")
        ? requiredArray(topology, "requiredRegionsForTwoToThreePhaseBenchmark", "topologyReference")
        : requiredArray(topology, "observedRegions", "topologyReference");
    Set<Region> regions = regions(requiredRegionArray);
    Set<String> boundaries = new LinkedHashSet<String>();
    for (JsonElement element : topology.getAsJsonArray("requiredTwoToThreePhaseBoundaries")) {
      boundaries.add(boundaryCode(requiredString(element.getAsJsonObject(), "code", "reference boundary")));
    }
    return new Reference(requiredString(manifest, "benchmarkId", "manifest"), sourceContractVerified, topologyOnly,
        regions, boundaries);
  }

  private static Candidate candidate(JsonObject json) {
    Set<Region> regions = regions(requiredArray(json, "observedRegions", "candidate"));
    Set<String> boundaries = new LinkedHashSet<String>();
    for (JsonElement element : requiredArray(json, "observedBoundaryTopologies", "candidate")) {
      boundaries.add(boundaryCode(element.getAsString()));
    }
    return new Candidate(requiredString(json, "identifier", "candidate"),
        json.has("topologyEvidenceResolved") && json.get("topologyEvidenceResolved").getAsBoolean(),
        json.has("internalQualityEligible") && json.get("internalQualityEligible").getAsBoolean(), regions, boundaries);
  }

  private static JsonObject serialize(Path manifestPath, JsonObject manifest, ReferenceAudit audit, Reference reference,
      Candidate candidate, HydrocarbonWaterTopologyBenchmark.Result result) {
    JsonObject output = new JsonObject();
    output.addProperty("status", "success");
    output.addProperty("method", "verified-paper-source-topology-adjacency-regression");
    output.addProperty("referenceManifestPath", manifestPath.toString());
    output.addProperty("referenceManifestSha256", audit.getManifestSha256());
    output.addProperty("referenceIdentifier", reference.getIdentifier());
    output.addProperty("referenceDoi", optionalString(manifest, "doi"));
    output.addProperty("candidateIdentifier", candidate.getIdentifier());
    output.addProperty("sourceContractVerified", audit.isVerified());
    output.addProperty("referenceCompositionSum", audit.getCompositionSum());
    output.add("artifactAudit", audit.getArtifactAudit());
    output.add("sourceContractViolations", strings(audit.getViolations()));
    output.addProperty("topologyComparisonEligible", result.isComparisonEligible());
    output.addProperty("topologyAccepted", result.isTopologyAccepted());
    output.addProperty("candidateInternalQualityEligible", candidate.isInternalQualityEligible());
    output.addProperty("numericComparisonEligible", result.isNumericComparisonEligible());
    output.addProperty("numericThresholdsFrozen", false);
    output.addProperty("numericBenchmarkPending", result.isNumericBenchmarkPending());
    output.addProperty("engineeringEligible", result.isEngineeringEligible());
    output.addProperty("releaseRule",
        "paper topology acceptance is supporting evidence only; exact-fluid exact-model numeric regression and the independent internal quality gate are still required");
    output.add("requiredRegions", regionCodes(reference.getRequiredRegions()));
    output.add("observedRegions", regionCodes(candidate.getObservedRegions()));
    output.add("missingRegions", regionCodes(result.getMissingRegions()));
    output.add("requiredBoundaryTopologies", strings(reference.getRequiredBoundaryCodes()));
    output.add("observedBoundaryTopologies", strings(candidate.getObservedBoundaryCodes()));
    output.add("missingBoundaryTopologies", strings(result.getMissingBoundaryCodes()));
    output.add("violations", strings(result.getViolations()));
    return output;
  }

  private static Set<Region> regions(JsonArray array) {
    Set<Region> output = EnumSet.noneOf(Region.class);
    for (JsonElement element : array) {
      output.add(region(element.getAsString()));
    }
    return output;
  }

  private static Region region(String value) {
    String token = value.trim().toUpperCase(Locale.ROOT).replaceAll("[^A-Z]", "");
    if (token.equals("G") || token.equals("V") || token.equals("GAS")) {
      return Region.GAS;
    }
    if (token.equals("O") || token.equals("L") || token.equals("OIL")) {
      return Region.OIL;
    }
    if (token.equals("W") || token.equals("AQ") || token.equals("AQUEOUS")) {
      return Region.AQUEOUS;
    }
    if (token.equals("GO") || token.equals("LV") || token.equals("GASOIL")) {
      return Region.GAS_OIL;
    }
    if (token.equals("GW") || token.equals("VW") || token.equals("GASAQUEOUS")) {
      return Region.GAS_AQUEOUS;
    }
    if (token.equals("OW") || token.equals("LW") || token.equals("OILAQUEOUS")) {
      return Region.OIL_AQUEOUS;
    }
    if (token.equals("GOW") || token.equals("LVW") || token.equals("GASOILAQUEOUS")) {
      return Region.GAS_OIL_AQUEOUS;
    }
    throw new IllegalArgumentException("unsupported phase-region token " + value);
  }

  private static String boundaryCode(String value) {
    String token = value.trim().toUpperCase(Locale.ROOT).replaceAll("[^A-Z]", "");
    BoundaryDefinition boundary;
    if (token.equals("OWTOGOW") || token.equals("OWGOW")) {
      boundary = new BoundaryDefinition(Region.OIL_AQUEOUS, Region.GAS_OIL_AQUEOUS);
    } else if (token.equals("GWTOGOW") || token.equals("GWGOW")) {
      boundary = new BoundaryDefinition(Region.GAS_AQUEOUS, Region.GAS_OIL_AQUEOUS);
    } else if (token.equals("GOTOGOW") || token.equals("GOGOW")) {
      boundary = new BoundaryDefinition(Region.GAS_OIL, Region.GAS_OIL_AQUEOUS);
    } else {
      throw new IllegalArgumentException("unsupported two-to-three-phase boundary token " + value);
    }
    return boundary.getCode();
  }

  private static String manifestSha256(Path path) {
    return sha256(path);
  }

  private static String sha256(Path path) {
    try {
      MessageDigest digest = MessageDigest.getInstance("SHA-256");
      byte[] bytes = Files.readAllBytes(path);
      StringBuilder output = new StringBuilder();
      for (byte value : digest.digest(bytes)) {
        output.append(String.format(Locale.ROOT, "%02x", value & 0xff));
      }
      return output.toString();
    } catch (IOException | NoSuchAlgorithmException error) {
      throw new IllegalArgumentException("unable to hash " + path + ": " + error.getMessage(), error);
    }
  }

  private static long size(Path path) {
    try {
      return Files.size(path);
    } catch (IOException error) {
      return -1L;
    }
  }

  private static JsonObject requiredObject(JsonObject json, String name, String label) {
    if (!json.has(name) || !json.get(name).isJsonObject()) {
      throw new IllegalArgumentException(label + " " + name + " must be an object");
    }
    return json.getAsJsonObject(name);
  }

  private static JsonArray requiredArray(JsonObject json, String name, String label) {
    if (!json.has(name) || !json.get(name).isJsonArray()) {
      throw new IllegalArgumentException(label + " " + name + " must be an array");
    }
    return json.getAsJsonArray(name);
  }

  private static String requiredString(JsonObject json, String name, String label) {
    String value = optionalString(json, name);
    if (value.isEmpty()) {
      throw new IllegalArgumentException(label + " " + name + " is required");
    }
    return value;
  }

  private static String optionalString(JsonObject json, String name) {
    return json.has(name) && !json.get(name).isJsonNull() ? json.get(name).getAsString().trim() : "";
  }

  private static long optionalLong(JsonObject json, String name, long fallback) {
    return json.has(name) && !json.get(name).isJsonNull() ? json.get(name).getAsLong() : fallback;
  }

  private static double requiredNumber(JsonObject json, String name, String label) {
    if (!json.has(name) || json.get(name).isJsonNull()) {
      throw new IllegalArgumentException(label + " " + name + " is required");
    }
    double value = json.get(name).getAsDouble();
    if (!Double.isFinite(value)) {
      throw new IllegalArgumentException(label + " " + name + " must be finite");
    }
    return value;
  }

  private static JsonArray strings(Iterable<String> values) {
    JsonArray output = new JsonArray();
    for (String value : values) {
      output.add(value);
    }
    return output;
  }

  private static JsonArray regionCodes(Iterable<Region> values) {
    JsonArray output = new JsonArray();
    for (Region value : values) {
      output.add(value.getCode());
    }
    return output;
  }

  private static final class ReferenceAudit {
    private final boolean verified;
    private final String manifestSha256;
    private final double compositionSum;
    private final JsonArray artifactAudit;
    private final List<String> violations;

    private ReferenceAudit(boolean verified, String manifestSha256, double compositionSum, JsonArray artifactAudit,
        List<String> violations) {
      this.verified = verified;
      this.manifestSha256 = manifestSha256;
      this.compositionSum = compositionSum;
      this.artifactAudit = artifactAudit;
      this.violations = violations;
    }

    private boolean isVerified() {
      return verified;
    }

    private String getManifestSha256() {
      return manifestSha256;
    }

    private double getCompositionSum() {
      return compositionSum;
    }

    private JsonArray getArtifactAudit() {
      return artifactAudit.deepCopy();
    }

    private List<String> getViolations() {
      return violations;
    }
  }
}
