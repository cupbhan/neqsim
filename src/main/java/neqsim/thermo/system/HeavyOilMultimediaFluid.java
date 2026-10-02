package neqsim.thermo.system;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import neqsim.thermo.mixingrule.EosMixingRuleType;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.util.readwrite.JsonFluidReadWrite;

/**
 * Factory for the versioned development heavy-oil fluid containing H2O, CO2, N2, and NH3.
 *
 * <p>
 * The fluid package is the software-delivery baseline validated against PVTsim Nova using nonreactive SRK, Peneloux
 * volume translation, and the tabulated classical mixing rule. It is a development fluid derived from a PVTsim
 * demonstration heavy-end distribution and is not a field-qualified reservoir-fluid description.
 * </p>
 *
 * <p>
 * The released envelope covers 20-250 degrees Celsius and 1-100 bara. The delivery gate accepts at least 98 percent
 * coarse-grid topology agreement, mean physical-boundary deviation up to 1.5 bara, and maximum physical-boundary
 * deviation up to 5 bara. The largest observed deviation is 4.25 bara on the gas-disappearance boundary near 80 degrees
 * Celsius. NH3 is nonreactive: this class makes no pH, ion-speciation, or reaction-heat claim.
 * </p>
 */
public final class HeavyOilMultimediaFluid {
  /** Versioned PVTsim/NeqSim compatibility contract identifier. */
  public static final String CONTRACT_ID = "HeavyOil-DEV-H2O65-CO27-N23-NH3-SRK-v1";

  /** Minimum temperature of the released software envelope in degrees Celsius. */
  public static final double MIN_TEMPERATURE_C = 20.0;

  /** Maximum temperature of the released software envelope in degrees Celsius. */
  public static final double MAX_TEMPERATURE_C = 250.0;

  /** Minimum pressure of the released software envelope in bara. */
  public static final double MIN_PRESSURE_BARA = 1.0;

  /** Maximum pressure of the released software envelope in bara. */
  public static final double MAX_PRESSURE_BARA = 100.0;

  /** Released minimum coarse-grid topology agreement in percent. */
  public static final double MIN_TOPOLOGY_AGREEMENT_PCT = 98.0;

  /** Released maximum mean physical-boundary deviation in bara. */
  public static final double MAX_MEAN_BOUNDARY_DEVIATION_BARA = 1.5;

  /** Released maximum physical-boundary deviation in bara. */
  public static final double MAX_BOUNDARY_DEVIATION_BARA = 5.0;

  private static final String RESOURCE_ROOT = "/neqsim/thermo/fluid/";

  private static final String[] HEAVY_COMPONENTS = {"C7", "C8", "C9", "C10-C12", "C13-C15", "C16-C18", "C19-C22",
      "C23-C26", "C27-C31", "C32-C38", "C39-C49", "C50-C80"};

  private static final double[] CO2_WATER_KNOTS_K = toKelvin(20.0, 203.0, 220.0, 222.0, 250.0);
  private static final double[] CO2_WATER_VALUES = {0.10, 0.10, 0.143, 0.10, 0.10};
  private static final double[] CO2_OTHER_KNOTS_K = toKelvin(20.0, 180.0, 200.0, 203.0, 220.0, 222.0, 250.0);
  private static final double[] CO2_OTHER_VALUES = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  private static final double[] N2_WATER_KNOTS_K = toKelvin(20.0, 50.0, 80.0, 150.0, 220.0, 250.0);
  private static final double[] N2_WATER_VALUES = {-1.50, -1.50, -1.60, -1.90, -1.90, -2.20};
  private static final double[] N2_HEAVY_KNOTS_K = toKelvin(20.0, 50.0, 80.0, 150.0, 220.0, 250.0);
  private static final double[] N2_HEAVY_VALUES = {-0.10, -0.225, -0.10, -0.12, -0.40, -0.40};

  /** Supported, independently validated NH3 dose anchors. */
  public enum Nh3Dose {
    /** 0.1 mol% NH3 dose. */
    NH3_0P1(0.1, "heavy-oil-multimedia-nh3-0p1-v1.json", 0.0, 0.0),
    /** 1 mol% NH3 dose. */
    NH3_1(1.0, "heavy-oil-multimedia-nh3-1-v1.json", -0.11, -0.18),
    /** 5 mol% NH3 dose. */
    NH3_5(5.0, "heavy-oil-multimedia-nh3-5-v1.json", -0.125, -0.16);

    private final double molPercent;
    private final String resourceName;
    private final double waterKij;
    private final double heavyKij;

    Nh3Dose(double molPercent, String resourceName, double waterKij, double heavyKij) {
      this.molPercent = molPercent;
      this.resourceName = resourceName;
      this.waterKij = waterKij;
      this.heavyKij = heavyKij;
    }

    /**
     * Get the NH3 amount in mol%.
     *
     * @return NH3 dose in mol%
     */
    public double getMolPercent() {
      return molPercent;
    }

    private static Nh3Dose fromMolPercent(double molPercent) {
      for (Nh3Dose dose : values()) {
        if (Math.abs(dose.molPercent - molPercent) < 1.0e-12) {
          return dose;
        }
      }
      throw new IllegalArgumentException(
          "Unsupported NH3 dose " + molPercent + " mol%. Supported doses are 0.1, 1.0, and 5.0 mol%.");
    }
  }

  private HeavyOilMultimediaFluid() {
  }

  /**
   * Create a released heavy-oil multimedia fluid.
   *
   * @param dose validated NH3 dose anchor
   * @param temperatureK temperature in Kelvin
   * @param pressureBara pressure in bara
   * @return configured nonreactive SRK fluid
   */
  public static SystemInterface create(Nh3Dose dose, double temperatureK, double pressureBara) {
    if (dose == null) {
      throw new IllegalArgumentException("NH3 dose cannot be null.");
    }
    String json = readResource(dose.resourceName);
    JsonObject definition = JsonParser.parseString(json).getAsJsonObject();
    SystemInterface fluid = JsonFluidReadWrite.readString(json);
    fluid.setTemperature(temperatureK);
    fluid.setPressure(pressureBara);
    configureCompatibilityRule(fluid, definition, dose);
    fluid.useVolumeCorrection(true);
    fluid.setMultiPhaseCheck(true);
    fluid.init(0);
    return fluid;
  }

  /**
   * Create a released heavy-oil multimedia fluid using a numeric NH3 dose.
   *
   * @param nh3MolPercent NH3 amount; must be 0.1, 1.0, or 5.0 mol%
   * @param temperatureK temperature in Kelvin
   * @param pressureBara pressure in bara
   * @return configured nonreactive SRK fluid
   */
  public static SystemInterface create(double nh3MolPercent, double temperatureK, double pressureBara) {
    return create(Nh3Dose.fromMolPercent(nh3MolPercent), temperatureK, pressureBara);
  }

  /**
   * Check whether conditions are inside the released software-validation envelope.
   *
   * @param temperatureC temperature in degrees Celsius
   * @param pressureBara pressure in bara
   * @return true when both temperature and pressure are in the released envelope
   */
  public static boolean isInsideReleasedEnvelope(double temperatureC, double pressureBara) {
    return temperatureC >= MIN_TEMPERATURE_C && temperatureC <= MAX_TEMPERATURE_C && pressureBara >= MIN_PRESSURE_BARA
        && pressureBara <= MAX_PRESSURE_BARA;
  }

  /**
   * State whether this development fluid is field-qualified.
   *
   * @return always false for contract version 1
   */
  public static boolean isFieldQualified() {
    return false;
  }

  /**
   * State whether NH3/CO2/H2O reaction chemistry is enabled.
   *
   * @return always false for the nonreactive SRK delivery
   */
  public static boolean isReactive() {
    return false;
  }

  private static void configureCompatibilityRule(SystemInterface fluid, JsonObject definition, Nh3Dose dose) {
    fluid.setMixingRule(EosMixingRuleType.CLASSIC_T_TABULATED);
    Map<String, Integer> indices = componentIndices(fluid);
    applyBaseBinaryInteractions(fluid, definition, indices);
    int co2 = requiredIndex(indices, "CO2");
    int water = requiredIndex(indices, "H2O");
    int n2 = requiredIndex(indices, "N2");
    int nh3 = requiredIndex(indices, "NH3");
    for (int phaseIndex = 0; phaseIndex < fluid.getMaxNumberOfPhases(); phaseIndex++) {
      EosMixingRulesInterface rule = ((PhaseEosInterface) fluid.getPhase(phaseIndex)).getEosMixingRule();
      rule.setBinaryInteractionParameterTemperatureTable(co2, water, CO2_WATER_KNOTS_K, CO2_WATER_VALUES);
      for (Map.Entry<String, Integer> entry : indices.entrySet()) {
        if (!"CO2".equals(entry.getKey()) && !"H2O".equals(entry.getKey())) {
          rule.setBinaryInteractionParameterTemperatureTable(co2, entry.getValue(), CO2_OTHER_KNOTS_K,
              CO2_OTHER_VALUES);
        }
      }
      rule.setBinaryInteractionParameterTemperatureTable(n2, water, N2_WATER_KNOTS_K, N2_WATER_VALUES);
      for (String name : HEAVY_COMPONENTS) {
        rule.setBinaryInteractionParameterTemperatureTable(n2, requiredIndex(indices, name), N2_HEAVY_KNOTS_K,
            N2_HEAVY_VALUES);
        rule.setBinaryInteractionParameter(nh3, requiredIndex(indices, name), dose.heavyKij);
      }
      rule.setBinaryInteractionParameter(nh3, water, dose.waterKij);
      rule.setBinaryInteractionParameter(nh3, co2, 0.0);
      rule.setBinaryInteractionParameter(nh3, n2, 0.0);
    }
  }

  private static void applyBaseBinaryInteractions(SystemInterface fluid, JsonObject definition,
      Map<String, Integer> indices) {
    JsonArray pairs = definition.getAsJsonArray("binaryInteractionCoefficients");
    if (pairs == null) {
      return;
    }
    for (int phaseIndex = 0; phaseIndex < fluid.getMaxNumberOfPhases(); phaseIndex++) {
      EosMixingRulesInterface rule = ((PhaseEosInterface) fluid.getPhase(phaseIndex)).getEosMixingRule();
      for (JsonElement element : pairs) {
        JsonObject pair = element.getAsJsonObject();
        int first = requiredIndex(indices, canonicalName(pair.get("i").getAsString()));
        int second = requiredIndex(indices, canonicalName(pair.get("j").getAsString()));
        rule.setBinaryInteractionParameter(first, second, pair.get("kij").getAsDouble());
      }
    }
  }

  private static Map<String, Integer> componentIndices(SystemInterface fluid) {
    Map<String, Integer> indices = new LinkedHashMap<String, Integer>();
    for (int index = 0; index < fluid.getPhase(0).getNumberOfComponents(); index++) {
      String runtimeName = fluid.getPhase(0).getComponent(index).getComponentName();
      indices.put(canonicalName(runtimeName), index);
    }
    return indices;
  }

  private static int requiredIndex(Map<String, Integer> indices, String componentName) {
    Integer index = indices.get(componentName);
    if (index == null) {
      throw new IllegalStateException("Missing component in heavy-oil multimedia fluid: " + componentName);
    }
    return index.intValue();
  }

  private static String canonicalName(String name) {
    String normalized = name.endsWith("_PC") ? name.substring(0, name.length() - 3) : name;
    String lower = normalized.toLowerCase(Locale.ROOT);
    if ("water".equals(lower) || "h2o".equals(lower)) {
      return "H2O";
    }
    if ("carbon dioxide".equals(lower) || "co2".equals(lower)) {
      return "CO2";
    }
    if ("nitrogen".equals(lower) || "n2".equals(lower)) {
      return "N2";
    }
    if ("ammonia".equals(lower) || "nh3".equals(lower)) {
      return "NH3";
    }
    return normalized;
  }

  private static String readResource(String resourceName) {
    InputStream input = HeavyOilMultimediaFluid.class.getResourceAsStream(RESOURCE_ROOT + resourceName);
    if (input == null) {
      throw new IllegalStateException("Missing bundled heavy-oil multimedia fluid resource: " + resourceName);
    }
    try (InputStream stream = input; ByteArrayOutputStream output = new ByteArrayOutputStream()) {
      byte[] buffer = new byte[8192];
      int count;
      while ((count = stream.read(buffer)) >= 0) {
        output.write(buffer, 0, count);
      }
      return new String(output.toByteArray(), StandardCharsets.UTF_8);
    } catch (IOException exception) {
      throw new IllegalStateException("Failed to read bundled heavy-oil multimedia fluid resource: " + resourceName,
          exception);
    }
  }

  private static double[] toKelvin(double... temperaturesC) {
    double[] temperaturesK = new double[temperaturesC.length];
    for (int index = 0; index < temperaturesC.length; index++) {
      temperaturesK[index] = temperaturesC[index] + 273.15;
    }
    return temperaturesK;
  }
}
