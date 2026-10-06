package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;

/**
 * Validated enthalpy and viscosity evaluation on the equilibrium compositions of an explicit flash.
 *
 * @author NeqSim development team
 * @version 1.0
 */
public final class ExplicitPropertyCalculator {
  static final String PROFILE = "native-properties-v1";
  static final String ENTHALPY = "native-enthalpy-sdk-convention-v1";

  /**
   * Read a required finite numeric parameter.
   *
   * @param value JSON object containing the parameter
   * @param key parameter or property name
   * @return evaluated value
   */
  public static double n(JsonNode value, String key) {
    return ExplicitFlash.number(value, key);
  }

  /**
   * Check finite values against combined absolute and relative tolerances.
   *
   * @param a BWR polynomial coefficients
   * @param b BWR exponential coefficients
   * @param absolute absolute tolerance
   * @param relative relative tolerance
   * @return evaluated value
   */
  public static boolean close(double a, double b, double absolute, double relative) {
    return Double.isFinite(a) && Double.isFinite(b) && Math.abs(a - b) <= absolute + relative * Math.abs(a);
  }

  /**
   * Record an unavailable property and its reason.
   *
   * @param info property availability metadata
   * @param key parameter or property name
   * @param code unavailability reason code
   * @return evaluated value
   */
  public static ObjectNode missing(ObjectNode info, String key, String code) {
    return info.putObject(key).put("status", "unavailable").put("code", code);
  }

  /**
   * Attach validated enthalpy, viscosity and availability metadata to an equilibrium phase.
   *
   * @param input explicit model and calculation settings
   * @param point temperature in degrees Celsius and pressure in MPa
   * @param phase phase type and mole-fraction composition
   */
  public static void properties(JsonNode input, JsonNode point, ObjectNode phase) {
    ObjectNode info = phase.putObject("properties");
    phase.putNull("enthalpyJMol");
    phase.putNull("viscosityPaS");
    String blocked = null;
    for (JsonNode c : input.path("components"))
      if (c.has("saltIonNumber"))
        blocked = "UNSUPPORTED_SALT";
    double tc = n(point, "temperatureC"), p = n(point, "pressureMPa");
    if (blocked == null && !(tc >= 25 && tc <= 150 && p >= 1 && p <= 20))
      blocked = "OUTSIDE_PROPERTY_RANGE";
    if (blocked == null && !input.path("nativeProperties").path("enabled").asBoolean())
      blocked = "PARAMETERS_MISSING";
    if (blocked != null) {
      missing(info, "enthalpy", blocked);
      missing(info, "viscosity", blocked);
      return;
    }
    try {
      if (!PROFILE.equals(input.path("nativeProperties").path("profile").asText())
          || input.path("propertyAuditVersion").asInt() != 4)
        throw new IllegalArgumentException("Unknown native property profile");
      ObjectNode values = ExplicitProperties.phase(input, point, phase);
      phase.set("propertyDiagnostics", values);
      double t = tc + 273.15, ideal = 0;
      JsonNode xs = phase.path("composition");
      for (int i = 0; i < xs.size(); i++) {
        JsonNode cp = input.path("components").get(i).path("cpJMolK");
        for (int j = 0; j < 4; j++)
          ideal += xs.get(i).asDouble() * cp.get(j).asDouble() * (Math.pow(t, j + 1) - Math.pow(273.15, j + 1))
              / (j + 1);
      }
      double residual = n(values, "residualFiniteDifferenceJMol");
      double h = ideal + residual * 8.3147295 / ExplicitFlash.R + n(values, "penelouxEnthalpyCorrectionJMol");
      if (Double.isFinite(h) && close(residual, n(values, "residualHalfStepJMol"), .01, 1e-6)
          && close(n(phase, "densityKgM3"), n(values, "densityKgM3"), .01, 1e-4)
          && close(n(phase, "z"), n(values, "z"), 1e-6, 1e-4)
          && close(ideal, n(values, "enthalpyBoundCpJMol") - n(values, "residualAnalyticFrozenKijJMol"), 1e-6, 1e-8)) {
        phase.put("enthalpyJMol", h);
        info.putObject("enthalpy").put("status", "available").put("method", ENTHALPY)
            .put("referenceTemperatureK", 273.15).put("basis", "explicit-Cp-and-EOS-residual");
      } else
        missing(info, "enthalpy", "ENTHALPY_CHECK_FAILED");
      boolean water = phase.path("type").asText().equals("aqueous");
      JsonNode viscosity = values.path(water ? "waterBaseline" : "nativeCsp");
      String expected = input.path("transport").path(water ? "waterCandidate" : "nativeCandidate").asText();
      double mu = viscosity.path("viscosityPaS").asDouble(Double.NaN);
      if (expected.equals(viscosity.path("method").asText()) && Double.isFinite(mu) && mu > 0) {
        phase.put("viscosityPaS", mu);
        ObjectNode v = info.putObject("viscosity").put("status", "available").put("method", expected)
            .put("basis", water ? "pure-water-approximation" : "hydrocarbon-csp").put("mixtureProperty", !water);
        if (water)
          v.put("waterMoleFraction", n(viscosity, "waterMoleFraction"));
      } else
        missing(info, "viscosity", water ? "WATER_MIXTURE_UNSUPPORTED" : "VISCOSITY_UNSUPPORTED");
    } catch (Exception error) {
      phase.putNull("enthalpyJMol");
      phase.putNull("viscosityPaS");
      missing(info, "enthalpy", "PROPERTY_CALCULATION_FAILED");
      missing(info, "viscosity", "PROPERTY_CALCULATION_FAILED");
      phase.put("propertyError", error.toString());
    }
  }
}
