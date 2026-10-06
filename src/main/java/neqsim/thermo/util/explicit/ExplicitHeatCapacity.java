package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import neqsim.thermo.system.SystemInterface;

/**
 * Fixed-composition isobaric heat capacity for the explicit native enthalpy convention.
 *
 * @author NeqSim development team
 * @version 1.0
 */
final class ExplicitHeatCapacity {
  /** Prevent construction of the stateless evaluator. */
  private ExplicitHeatCapacity() {
  }

  /**
   * Differentiate enthalpy with independent step sizes while retaining pressure, composition and EOS root.
   *
   * @param input explicit component and property parameters
   * @param point requested temperature and pressure
   * @param phase previously validated equilibrium phase
   * @param info property availability metadata
   */
  static void attach(JsonNode input, JsonNode point, ObjectNode phase, ObjectNode info) {
    try {
      double t = ExplicitFlash.number(point, "temperatureC") + 273.15;
      double p = ExplicitFlash.number(point, "pressureMPa") * 10;
      double h = ExplicitFlash.number(phase, "enthalpyJMol");
      double coarse = derivative(input, phase, t, p, h, 0.05);
      double fine = derivative(input, phase, t, p, h, 0.025);
      if (!(fine > 0) || !ExplicitPropertyCalculator.close(fine, coarse, 0.03, 3e-4)) {
        ExplicitPropertyCalculator.missing(info, "cp", "CP_DERIVATIVE_NOT_CONVERGED");
        return;
      }
      phase.put("cpJMolK", fine);
      info.putObject("cp").put("status", "available").put("method", "native-enthalpy-temperature-derivative-v1")
          .put("basis", "fixed-phase-composition-and-pressure").put("temperatureStepK", 0.025)
          .put("stepDifferenceJMolK", Math.abs(fine - coarse));
    } catch (Exception error) {
      ExplicitPropertyCalculator.missing(info, "cp", "CP_CALCULATION_FAILED");
    }
  }

  /**
   * Evaluate a second-order derivative without crossing the qualified thermal range.
   *
   * @param input explicit model
   * @param phase fixed phase and composition
   * @param t temperature in kelvin
   * @param p pressure in bara
   * @param h previously validated enthalpy in J/mol
   * @param step perturbation in kelvin
   * @return isobaric heat capacity in J/(mol K)
   */
  private static double derivative(JsonNode input, JsonNode phase, double t, double p, double h, double step) {
    if (t - step < 298.15) {
      return (-3 * h + 4 * enthalpy(input, phase, t + step, p) - enthalpy(input, phase, t + 2 * step, p)) / (2 * step);
    }
    if (t + step > 423.15) {
      return (3 * h - 4 * enthalpy(input, phase, t - step, p) + enthalpy(input, phase, t - 2 * step, p)) / (2 * step);
    }
    return (enthalpy(input, phase, t + step, p) - enthalpy(input, phase, t - step, p)) / (2 * step);
  }

  /**
   * Replay the validated enthalpy convention on an isolated fixed-composition state.
   *
   * @param input explicit model
   * @param phase fixed equilibrium phase
   * @param t temperature in kelvin
   * @param p pressure in bara
   * @return molar enthalpy in J/mol
   */
  private static double enthalpy(JsonNode input, JsonNode phase, double t, double p) {
    SystemInterface fluid = ExplicitProperties.state(input, phase.path("composition"), phase.path("type").asText(), t,
        p);
    double ideal = 0, correction = 0;
    for (int i = 0; i < input.path("components").size(); i++) {
      JsonNode c = input.path("components").get(i);
      double x = phase.path("composition").get(i).asDouble();
      for (int j = 0; j < 4; j++) {
        ideal += x * c.path("cpJMolK").get(j).asDouble() * (Math.pow(t, j + 1) - Math.pow(273.15, j + 1)) / (j + 1);
      }
      correction += x * (-ExplicitFlash.number(c, "penelouxLmol") + 288.15 * ExplicitFlash.number(c, "penelouxTLmolK"))
          * 1e-3 * p * 1e5;
    }
    double z = fluid.getPhase(0).getZ();
    double baseZ = ExplicitFlash.number(phase, "untranslatedZ");
    if (!Double.isFinite(z) || Math.abs(z - baseZ) > 0.05) {
      throw new IllegalStateException("Fixed-phase root changed during heat-capacity evaluation");
    }
    double derivative = (ExplicitProperties.evaluatePhi(fluid, input, t + 0.01)
        - ExplicitProperties.evaluatePhi(fluid, input, t - 0.01)) / 0.02;
    double half = (ExplicitProperties.evaluatePhi(fluid, input, t + 0.005)
        - ExplicitProperties.evaluatePhi(fluid, input, t - 0.005)) / 0.01;
    double residual = -ExplicitFlash.R * t * t * derivative;
    if (!ExplicitPropertyCalculator.close(residual, -ExplicitFlash.R * t * t * half, 0.01, 1e-6)) {
      throw new IllegalStateException("Enthalpy derivative failed its step convergence check");
    }
    return ideal + residual * 8.3147295 / ExplicitFlash.R + correction;
  }
}
