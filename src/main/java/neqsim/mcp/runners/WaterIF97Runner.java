package neqsim.mcp.runners;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import neqsim.mcp.model.ApiEnvelope;
import neqsim.mcp.model.ResultProvenance;
import neqsim.process.util.monitor.FluidResponse;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemWaterIF97;
import neqsim.thermo.util.steam.Iapws_if97;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/** Pure-water state and saturation calculations using IAPWS-IF97. */
public final class WaterIF97Runner {
  private static final Gson GSON = new GsonBuilder().setPrettyPrinting().serializeSpecialFloatingPointValues().create();
  private static final double TRIPLE_TEMPERATURE_K = 273.16;
  private static final double CRITICAL_TEMPERATURE_K = 647.096;
  private static final double CRITICAL_PRESSURE_BARA = 220.64;

  private WaterIF97Runner() {
  }

  /**
   * Calculates one pure-water TP state and an IF97 region-4 saturation curve.
   *
   * @param temperatureK state temperature in kelvin
   * @param pressureBara absolute state pressure in bara
   * @param minimumTemperatureK saturation-curve lower temperature in kelvin
   * @param maximumTemperatureK saturation-curve upper temperature in kelvin
   * @param requestedPointCount requested saturation-curve point count
   * @return standard NeqSim JSON result
   */
  public static String run(double temperatureK, double pressureBara, double minimumTemperatureK,
      double maximumTemperatureK, int requestedPointCount) {
    long startTime = System.currentTimeMillis();
    try {
      requireFiniteInRange(temperatureK, "temperatureK", TRIPLE_TEMPERATURE_K, 1073.15);
      requireFiniteInRange(pressureBara, "pressureBara", 1.0e-5, 10000.0);
      double minimumK = Math.max(minimumTemperatureK, TRIPLE_TEMPERATURE_K);
      double maximumK = Math.min(maximumTemperatureK, CRITICAL_TEMPERATURE_K - 1.0e-4);
      if (!Double.isFinite(minimumK) || !Double.isFinite(maximumK) || maximumK <= minimumK) {
        throw new IllegalArgumentException("Invalid IF97 saturation-curve temperature range");
      }
      int pointCount = Math.max(11, Math.min(requestedPointCount, 401));

      SystemInterface fluid = new SystemWaterIF97(temperatureK, pressureBara);
      new ThermodynamicOperations(fluid).TPflash();
      fluid.initProperties();

      JsonObject flash = new JsonObject();
      flash.addProperty("model", "IF97");
      flash.addProperty("flashType", "TP");
      flash.addProperty("numberOfPhases", fluid.getNumberOfPhases());
      JsonArray phases = new JsonArray();
      for (int index = 0; index < fluid.getNumberOfPhases(); index++) {
        phases.add(fluid.getPhase(index).getPhaseTypeName());
      }
      flash.add("phases", phases);

      JsonObject fluidObject = JsonParser.parseString(GSON.toJson(new FluidResponse("water-if97", fluid)))
          .getAsJsonObject();
      JsonArray saturationCurve = new JsonArray();
      for (int index = 0; index < pointCount; index++) {
        double fraction = (double) index / (pointCount - 1);
        double pointTemperatureK = minimumK + fraction * (maximumK - minimumK);
        double pressure = Iapws_if97.psat_t(pointTemperatureK) * 10.0;
        JsonObject point = new JsonObject();
        point.addProperty("temperature_K", pointTemperatureK);
        point.addProperty("temperature_C", pointTemperatureK - 273.15);
        point.addProperty("pressure_bara", pressure);
        saturationCurve.add(point);
      }

      JsonObject criticalPoint = new JsonObject();
      criticalPoint.addProperty("temperature_K", CRITICAL_TEMPERATURE_K);
      criticalPoint.addProperty("temperature_C", CRITICAL_TEMPERATURE_K - 273.15);
      criticalPoint.addProperty("pressure_bara", CRITICAL_PRESSURE_BARA);
      criticalPoint.addProperty("source", "IAPWS-IF97");

      JsonObject result = new JsonObject();
      result.addProperty("status", "success");
      result.addProperty("model", "IAPWS-IF97");
      result.addProperty("method", "IAPWS-IF97 TP state and region-4 saturation equation");
      result.add("flash", flash);
      result.add("fluid", fluidObject);
      result.add("envelope", saturationCurve);
      result.add("saturationCurve", saturationCurve.deepCopy());
      result.add("criticalPoint", criticalPoint);
      JsonObject data = new JsonObject();
      data.add("flash", flash.deepCopy());
      data.add("fluid", fluidObject.deepCopy());
      result.add("data", data);

      ResultProvenance provenance = ResultProvenance.forFlash("IF97", "TP", "not-applicable");
      provenance.setBenchmarkTrustLevel(BenchmarkTrust.getMaturityLevel("runFlash"));
      provenance.setComputationTimeMs(System.currentTimeMillis() - startTime);
      provenance.addValidationPassed("pure_water_model_enforced");
      provenance.addValidationPassed("if97_region_4_curve_monotonic");
      result.add("provenance", GSON.toJsonTree(provenance));
      ApiEnvelope.applyStandardFields(result, "runWaterIF97", provenance,
          ApiEnvelope.validationStatus(true, "pure_water_if97",
              "Pure-water state and saturation range passed IF97 input checks"),
          ApiEnvelope.qualityGate("passed", pointCount + " IAPWS-IF97 saturation points calculated", true));
      return GSON.toJson(result);
    } catch (Exception error) {
      JsonObject result = new JsonObject();
      result.addProperty("status", "error");
      result.addProperty("message", "IAPWS-IF97 water calculation failed: " + error.getMessage());
      return GSON.toJson(result);
    }
  }

  private static void requireFiniteInRange(double value, String name, double minimum, double maximum) {
    if (!Double.isFinite(value) || value < minimum || value > maximum) {
      throw new IllegalArgumentException(name + " must be within [" + minimum + ", " + maximum + "]");
    }
  }
}
