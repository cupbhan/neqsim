package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ObjectNode;

/**
 * IF97 region-1 density and industrial IAPWS 2008 viscosity, excluding critical enhancement and mixture corrections.
 *
 * @author NeqSim development team
 * @version 1.0
 */
public final class PureWaterViscosity {
  static final String METHOD = "iapws2008-if97-region1-pure-water-v1";
  // Region-1 Gibbs terms with I > 0; I=0 terms have zero pressure derivative.
  static final double[][] REGION1 = {{1, -9, .28319080123804e-3}, {1, -7, -.60706301565874e-3},
      {1, -1, -.18990068218419e-1}, {1, 0, -.32529748770505e-1}, {1, 1, -.21841717175414e-1},
      {1, 3, -.52838357969930e-4}, {2, -3, -.47184321073267e-3}, {2, 0, -.30001780793026e-3},
      {2, 1, .47661393906987e-4}, {2, 3, -.44141845330846e-5}, {2, 17, -.72694996297594e-15},
      {3, -4, -.31679644845054e-4}, {3, 0, -.28270797985312e-5}, {3, 6, -.85205128120103e-9},
      {4, -5, -.22425281908000e-5}, {4, -2, -.65171222895601e-6}, {4, 10, -.14341729937924e-12},
      {5, -8, -.40516996860117e-6}, {8, -11, -.12734301741641e-8}, {8, -6, -.17424871230634e-9},
      {21, -29, -.68762131295531e-18}, {23, -31, .14478307828521e-19}, {29, -38, .26335781662795e-22},
      {30, -39, -.11947622640071e-22}, {31, -40, .18228094581404e-23}, {32, -41, -.93537087292458e-25}};
  static final double[][] VISCOSITY = {{0, 0, .520094}, {1, 0, .0850895}, {2, 0, -1.08374}, {3, 0, -.289555},
      {0, 1, .222531}, {1, 1, .999115}, {2, 1, 1.88797}, {3, 1, 1.26613}, {5, 1, .120573}, {0, 2, -.281378},
      {1, 2, -.906851}, {2, 2, -.772479}, {3, 2, -.489837}, {4, 2, -.257040}, {0, 3, .161913}, {1, 3, .257399},
      {0, 4, -.0325372}, {3, 4, .0698452}, {4, 5, .00872102}, {3, 6, -.00435673}, {5, 6, -.000593264}};

  /**
   * Evaluate IF97 saturation pressure in MPa.
   *
   * @param t absolute temperature in kelvin
   * @return evaluated value
   */
  public static double saturationPressureMPa(double t) {
    if (!(t >= 273.15 && t <= 647.096))
      throw new IllegalArgumentException("Water saturation temperature outside IF97 range");
    double theta = t - .23855557567849 / (t - 650.17534844798);
    double a = theta * theta + 1167.0521452767 * theta - 724213.16703206;
    double b = -17.073846940092 * theta * theta + 12020.824702470 * theta - 3232555.0322333;
    double c = 14.915108613530 * theta * theta - 4823.2657361591 * theta + 405113.40542057;
    return Math.pow(2 * c / (-b + Math.sqrt(b * b - 4 * a * c)), 4);
  }

  /**
   * Evaluate stable IF97 region-1 liquid density in kg/m3.
   *
   * @param t absolute temperature in kelvin
   * @param pMPa pressure in MPa
   * @return evaluated value
   */
  public static double densityKgM3(double t, double pMPa) {
    if (!(t >= 273.16 && t <= 623.15 && pMPa <= 100 && pMPa >= saturationPressureMPa(t)))
      throw new IllegalArgumentException("Pure-water baseline outside stable IF97 region 1");
    double pi = pMPa / 16.53, tau = 1386 / t, gammaP = 0;
    for (double[] term : REGION1)
      gammaP -= term[2] * term[0] * Math.pow(7.1 - pi, term[0] - 1) * Math.pow(tau - 1.222, term[1]);
    return 16.53e6 / (461.526 * t * gammaP);
  }

  /**
   * Evaluate industrial IAPWS water viscosity in Pa s without critical enhancement.
   *
   * @param t absolute temperature in kelvin
   * @param rhoKgM3 water mass density in kg/m3
   * @return evaluated value
   */
  public static double viscosityPaS(double t, double rhoKgM3) {
    if (!(Double.isFinite(t) && t > 0 && Double.isFinite(rhoKgM3) && rhoKgM3 >= 0))
      throw new IllegalArgumentException("Invalid water temperature or density");
    double tr = t / 647.096, dr = rhoKgM3 / 322;
    double mu0 = 100 * Math.sqrt(tr) / (1.67752 + 2.20462 / tr + .6366564 / (tr * tr) - .241605 / (tr * tr * tr));
    double sum = 0;
    for (double[] term : VISCOSITY)
      sum += term[2] * Math.pow(1 / tr - 1, term[0]) * Math.pow(dr - 1, term[1]);
    return 1e-6 * mu0 * Math.exp(dr * sum);
  }

  /**
   * Evaluate the requested model and return its diagnostics.
   *
   * @param input explicit model and calculation settings
   * @param state temperature in degrees Celsius and pressure in MPa
   * @param phase phase type and mole-fraction composition
   * @return calculated values and diagnostics
   */
  public static ObjectNode calculate(JsonNode input, JsonNode state, JsonNode phase) {
    if (!phase.path("type").asText().equals("aqueous"))
      throw new IllegalArgumentException("Pure-water baseline requires an aqueous phase");
    JsonNode components = input.path("components"), xs = phase.path("composition");
    if (xs.size() != components.size())
      throw new IllegalArgumentException("Water composition size mismatch");
    double water = 0, total = 0;
    for (int i = 0; i < components.size(); i++) {
      JsonNode component = components.get(i), amount = xs.get(i);
      if (!amount.isNumber() || !Double.isFinite(amount.asDouble()) || amount.asDouble() < 0)
        throw new IllegalArgumentException("Invalid water phase composition");
      double x = amount.asDouble();
      total += x;
      String name = component.path("name").asText();
      if (name.equals("H2O"))
        water += x;
      else if (x > 1e-12 && (component.has("saltIonNumber") || !(component.path("componentType").asText().equals("HC")
          || name.equals("N2") || name.equals("CO2") || name.equals("H2S"))))
        throw new IllegalArgumentException("Salt or inhibitor mixture viscosity is not the pure-water baseline");
    }
    if (Math.abs(total - 1) > 1e-8 || water < .99)
      throw new IllegalArgumentException("Pure-water baseline requires at least 99 mol% water");
    double t = ExplicitFlash.number(state, "temperatureC") + 273.15;
    double p = ExplicitFlash.number(state, "pressureMPa");
    double rho = densityKgM3(t, p), mu = viscosityPaS(t, rho);
    if (!(Double.isFinite(mu) && mu > 0))
      throw new IllegalStateException("Invalid pure-water viscosity");
    return ExplicitFlash.JSON.createObjectNode().put("method", METHOD)
        .put("basis", "pure water at phase temperature and pressure").put("mixtureProperty", false)
        .put("waterMoleFraction", water).put("densityKgM3", rho).put("viscosityPaS", mu)
        .put("criticalEnhancementApplied", false);
  }
}
