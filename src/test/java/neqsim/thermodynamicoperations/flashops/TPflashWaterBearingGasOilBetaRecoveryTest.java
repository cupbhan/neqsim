package neqsim.thermodynamicoperations.flashops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.util.readwrite.JsonFluidReadWrite;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

/**
 * Checks conservation near water-bearing heavy-oil bubble points after speculative third-phase trials.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
class TPflashWaterBearingGasOilBetaRecoveryTest {
  /**
   * Requires a valid gas/oil split even when all water is dissolved in those two phases.
   *
   * @param waterFraction overall water mole fraction
   * @param temperatureC temperature in degrees Celsius
   * @param pressureBara absolute pressure in bar
   * @throws IOException if the public development fluid resource cannot be read
   */
  @ParameterizedTest
  @CsvSource({"0.1,240,9.65", "0.1,240,9.70", "0.1,250,10.90", "0.3,250,21.95", "0.3,250,22.00"})
  void coldFlashConservesFeedWithoutAnAqueousPhase(double waterFraction, double temperatureC, double pressureBara)
      throws IOException {
    SystemInterface fluid = createFluid(waterFraction, temperatureC, pressureBara);
    TPflash operation = new TPflash(fluid);
    operation.run();
    fluid.init(1);
    assertEquilibrium(fluid);
    assertTrue(operation.isLastMultiphaseSolveAccepted(), operation.getLastMultiphaseSolveMessage());
    assertEquals(2, fluid.getNumberOfPhases());
    assertTrue(fluid.hasPhaseType(PhaseType.GAS));
    assertTrue(fluid.hasPhaseType(PhaseType.OIL));
    assertFalse(fluid.hasPhaseType(PhaseType.AQUEOUS));
  }

  /**
   * Checks the public API and pressure round trips across the nearby phase boundary.
   *
   * @param waterFraction overall water mole fraction
   * @param temperatureC temperature in degrees Celsius
   * @param pressureBara absolute pressure in bar
   * @throws IOException if the public development fluid resource cannot be read
   */
  @ParameterizedTest
  @CsvSource({"0.1,240,9.65", "0.1,240,9.70", "0.1,250,10.90", "0.3,250,21.95", "0.3,250,22.00"})
  void repeatedAndNeighbouringFlashesRemainBalanced(double waterFraction, double temperatureC, double pressureBara)
      throws IOException {
    SystemInterface fluid = createFluid(waterFraction, temperatureC, pressureBara);
    ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
    operations.TPflash();
    fluid.init(1);
    assertEquilibrium(fluid);
    double originalGasBeta = fluid.getBeta(fluid.getPhaseNumberOfPhase("gas"));
    for (double offset : new double[] {0.0, -0.02, 0.0, 0.05, 0.0}) {
      fluid.setPressure(pressureBara + offset);
      operations.TPflash();
      fluid.init(1);
      assertEquilibrium(fluid);
      if (offset == 0.0) {
        assertEquals(2, fluid.getNumberOfPhases());
        assertEquals(originalGasBeta, fluid.getBeta(fluid.getPhaseNumberOfPhase("gas")), 1.0e-8);
      }
    }
  }

  /**
   * Derives an explicit zero-kij water/heavy-oil contract without modifying the frozen NH3 preset.
   *
   * @param waterFraction overall water mole fraction
   * @param temperatureC temperature in degrees Celsius
   * @param pressureBara absolute pressure in bar
   * @return a fresh SRK fluid with multiphase stability checks enabled
   * @throws IOException if the public development fluid resource cannot be read
   */
  private SystemInterface createFluid(double waterFraction, double temperatureC, double pressureBara)
      throws IOException {
    JsonObject source;
    try (InputStream stream = getClass()
        .getResourceAsStream("/neqsim/thermo/fluid/heavy-oil-multimedia-nh3-5-v1.json")) {
      assertNotNull(stream);
      try (InputStreamReader reader = new InputStreamReader(stream, StandardCharsets.UTF_8)) {
        source = JsonParser.parseReader(reader).getAsJsonObject();
      }
    }
    JsonArray components = new JsonArray();
    double heavyTotal = 0.0;
    for (JsonElement element : source.getAsJsonArray("components")) {
      JsonObject component = element.getAsJsonObject();
      if (component.get("isPseudo").getAsBoolean()) {
        heavyTotal += component.get("moleFraction").getAsDouble();
      }
    }
    for (JsonElement element : source.getAsJsonArray("components")) {
      JsonObject component = element.getAsJsonObject().deepCopy();
      if ("H2O".equals(component.get("name").getAsString())) {
        component.addProperty("moleFraction", waterFraction);
        components.add(component);
      } else if (component.get("isPseudo").getAsBoolean()) {
        component.addProperty("moleFraction",
            component.get("moleFraction").getAsDouble() / heavyTotal * (1.0 - waterFraction));
        components.add(component);
      }
    }
    JsonArray interactions = new JsonArray();
    for (int first = 0; first < components.size(); first++) {
      for (int second = first + 1; second < components.size(); second++) {
        JsonObject pair = new JsonObject();
        pair.addProperty("i", components.get(first).getAsJsonObject().get("name").getAsString());
        pair.addProperty("j", components.get(second).getAsJsonObject().get("name").getAsString());
        pair.addProperty("kij", 0.0);
        interactions.add(pair);
      }
    }
    JsonObject contract = new JsonObject();
    contract.addProperty("eos", "SRK");
    contract.add("components", components);
    contract.add("binaryInteractionCoefficients", interactions);
    SystemInterface fluid = JsonFluidReadWrite.readString(contract.toString());
    fluid.setTemperature(temperatureC + 273.15);
    fluid.setPressure(pressureBara);
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }

  /**
   * Checks phase amounts, normalized compositions, every component balance and fugacity equality.
   *
   * @param fluid flashed state to validate
   */
  private void assertEquilibrium(SystemInterface fluid) {
    double betaTotal = 0.0;
    for (int phase = 0; phase < fluid.getNumberOfPhases(); phase++) {
      double beta = fluid.getBeta(phase);
      assertTrue(Double.isFinite(beta) && beta >= 0.0 && beta <= 1.0);
      betaTotal += beta;
      double compositionTotal = 0.0;
      for (int component = 0; component < fluid.getNumberOfComponents(); component++) {
        double fraction = fluid.getPhase(phase).getComponent(component).getx();
        assertTrue(Double.isFinite(fraction) && fraction >= 0.0 && fraction <= 1.0);
        compositionTotal += fraction;
      }
      assertEquals(1.0, compositionTotal, 1.0e-10, "phase composition sum");
    }
    assertEquals(1.0, betaTotal, 1.0e-10, "phase fraction sum");
    for (int component = 0; component < fluid.getNumberOfComponents(); component++) {
      double reconstructed = 0.0;
      double reference = logFugacity(fluid, 0, component);
      for (int phase = 0; phase < fluid.getNumberOfPhases(); phase++) {
        reconstructed += fluid.getBeta(phase) * fluid.getPhase(phase).getComponent(component).getx();
        assertEquals(reference, logFugacity(fluid, phase, component), 1.0e-8, "log fugacity equality");
      }
      assertEquals(fluid.getPhase(0).getComponent(component).getz(), reconstructed, 1.0e-10,
          "component balance: " + fluid.getPhase(0).getComponent(component).getComponentName());
    }
  }

  /**
   * Computes log fugacity with the common pressure term omitted.
   *
   * @param fluid fluid state
   * @param phase phase index
   * @param component component index
   * @return logarithm of mole fraction times fugacity coefficient
   */
  private double logFugacity(SystemInterface fluid, int phase, int component) {
    return Math.log(Math.max(fluid.getPhase(phase).getComponent(component).getx(), Double.MIN_NORMAL))
        + Math.log(fluid.getPhase(phase).getComponent(component).getFugacityCoefficient());
  }
}
