package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.util.ArrayList;
import java.util.List;

/**
 * Hydrocarbon corresponding-states viscosity with explicit mixture factors and versioned methane reference profiles.
 *
 * @author NeqSim development team
 * @version 1.0
 */
public final class CspViscosity {
  /**
   * Copy the published dilute-reference coefficients.
   *
   * @return independent coefficient array
   */
  public static double[] diluteCoefficients() {
    return GV.clone();
  }

  /**
   * Copy the high-temperature dense-reference coefficients.
   *
   * @return independent coefficient array
   */
  public static double[] highDensityCoefficients() {
    return J.clone();
  }

  /**
   * Copy the low-temperature dense-reference coefficients.
   *
   * @return independent coefficient array
   */
  public static double[] lowDensityCoefficients() {
    return K.clone();
  }

  /** CSP profile with a 100 K floor for the initial reference density. */
  public static final String DENSITY_FLOOR_METHOD = "native-csp-mccarty1974-densityfloor100-v2";
  /** Opt-in SDK parity hypothesis with a 95 K reference blend center. */
  public static final String TRANSITION_METHOD = "native-csp-mccarty1974-sdktransition95-v3";
  /** Methane BWR gas constant in L atm/(mol K). */
  public static final double R_ATM = 0.08205616;
  /** Methane molar mass in g/mol. */
  public static final double MW = 16.042;
  private static final double GAMMA = 0.0096;
  private static final double TC = 190.555, PC_ATM = 45.387, RHOC = 10.15;
  static final double[] N = {-1.8439486666e-2, 1.0510162064, -1.6057820303e1, 8.4844027562e2, -4.2738409106e4,
      7.6565285254e-4, -4.8360724197e-1, 8.5195473835e1, -1.6607434721e4, -3.7521074532e-5, 2.8616309259e-2,
      -2.8685285973, 1.1906973942e-4, -8.5315715699e-3, 3.8365063841, 2.4986828379e-5, 5.7974531455e-6,
      -7.1648329297e-3, 1.2577853784e-4, 2.2240102466e4, -1.4800512328e6, 5.0498054887e1, 1.6428375992e6,
      2.1325387196e-1, 3.7791273422e1, -1.1857016815e-5, -3.1630780767e1, -4.1006782941e-6, 1.4870043284e-3,
      3.1512261532e-9, -2.1670774745e-6, 2.4000551079e-5};
  static final double[] GV = {-2.090975e5, 2.647269e5, -1.472818e5, 4.716740e4, -9.491872e3, 1.219979e3, -9.627993e1,
      4.274152, -8.141531e-2};
  static final double[] J = {-10.35060586, 17.571599671, -3019.3918656, 188.73011594, 0.042903609488, 145.29023444,
      6127.6818706};
  static final double[] K = {-9.74602, 18.0834, -4126.66, 44.6055, 0.976544, 81.8134, 15649.9};

  /**
   * Evaluate the methane BWR polynomial coefficients.
   *
   * @param t absolute temperature in kelvin
   * @return evaluated value
   */
  public static double[] polynomial(double t) {
    return new double[] {R_ATM * t, N[0] * t + N[1] * Math.sqrt(t) + N[2] + N[3] / t + N[4] / (t * t),
        N[5] * t + N[6] + N[7] / t + N[8] / (t * t), N[9] * t + N[10] + N[11] / t, N[12], N[13] / t + N[14] / (t * t),
        N[15] / t, N[16] / t + N[17] / (t * t), N[18] / (t * t)};
  }

  /**
   * Evaluate the methane BWR exponential coefficients.
   *
   * @param t absolute temperature in kelvin
   * @return evaluated value
   */
  public static double[] exponential(double t) {
    return new double[] {N[19] / (t * t) + N[20] / Math.pow(t, 3), N[21] / (t * t) + N[22] / Math.pow(t, 4),
        N[23] / (t * t) + N[24] / Math.pow(t, 3), N[25] / (t * t) + N[26] / Math.pow(t, 4),
        N[27] / (t * t) + N[28] / Math.pow(t, 3), N[29] / (t * t) + N[30] / Math.pow(t, 3) + N[31] / Math.pow(t, 4)};
  }

  /**
   * Evaluate methane reference pressure in atmospheres.
   *
   * @param rho methane molar density in mol/L
   * @param a BWR polynomial coefficients
   * @param b BWR exponential coefficients
   * @return evaluated value
   */
  public static double pressure(double rho, double[] a, double[] b) {
    double p = 0, q = 0;
    for (int i = a.length - 1; i >= 0; i--)
      p = (p + a[i]) * rho;
    for (int i = b.length - 1; i >= 0; i--)
      q = q * rho * rho + b[i];
    return p + q * rho * rho * rho * Math.exp(-GAMMA * rho * rho);
  }

  /**
   * Integrate the exponential BWR density term using a stable series.
   *
   * @param rho methane molar density in mol/L
   * @param m nonnegative integral order
   * @return evaluated value
   */
  public static double expIntegral(double rho, int m) {
    double u = GAMMA * rho * rho, term = 1, series = 1;
    for (int k = 1; k < 200; k++) {
      term *= u / (m + 1 + k);
      series += term;
      if (Math.abs(term) < 1e-15 * Math.abs(series))
        break;
    }
    return Math.pow(rho, 2 * m + 2) * Math.exp(-u) * series / (2 * (m + 1));
  }

  /**
   * Return the mole-fraction-weighted log fugacity coefficient.
   *
   * @param rho methane molar density in mol/L
   * @param p pressure in bar
   * @param a BWR polynomial coefficients
   * @param b BWR exponential coefficients
   * @return evaluated value
   */
  public static double logPhi(double rho, double p, double[] a, double[] b) {
    double integral = 0;
    for (int i = 1; i < a.length; i++)
      integral += a[i] * Math.pow(rho, i) / i;
    for (int i = 0; i < b.length; i++)
      integral += b[i] * expIntegral(rho, i);
    double z = p / (rho * a[0]);
    return integral / a[0] + z - 1 - Math.log(z);
  }

  /**
   * Solve methane reference molar density in mol/L using mechanically stable roots.
   *
   * @param t absolute temperature in kelvin
   * @param pAtm pressure in atmospheres
   * @return evaluated value
   */
  public static double density(double t, double pAtm) {
    if (!(t >= 45 && t <= 1000 && pAtm > 0 && pAtm <= 10000))
      throw new IllegalArgumentException("Methane reference state outside diagnostic solver limits");
    double[] a = polynomial(t), b = exponential(t);
    List<Double> roots = new ArrayList<>();
    double left = 1e-12, fLeft = pressure(left, a, b) - pAtm;
    while (left < 45) {
      double right = Math.min(45, left < .1 ? left * 1.5 : left + .025);
      double fRight = pressure(right, a, b) - pAtm;
      if (fLeft <= 0 && fRight >= 0) {
        double lo = left, hi = right;
        for (int i = 0; i < 80; i++) {
          double mid = (lo + hi) / 2;
          if (pressure(mid, a, b) < pAtm)
            lo = mid;
          else
            hi = mid;
        }
        roots.add((lo + hi) / 2);
      }
      left = right;
      fLeft = fRight;
    }
    if (roots.isEmpty())
      throw new IllegalStateException("No stable methane BWR root");
    double rho = roots.get(0), best = logPhi(rho, pAtm, a, b);
    if (t < 90.69) {
      // CSP extrapolates a dense fluid below methane freezing. MBWR has a
      // spurious intermediate-density Gibbs minimum there; continue the dense
      // outer root, not an equilibrium solid/vapour prediction.
      rho = roots.get(roots.size() - 1);
    } else
      for (double candidate : roots) {
        double g = logPhi(candidate, pAtm, a, b);
        if (g < best) {
          best = g;
          rho = candidate;
        }
      }
    if (Math.abs(pressure(rho, a, b) - pAtm) > 1e-6 * Math.max(1, pAtm))
      throw new IllegalStateException("Methane pressure residual too large");
    return rho;
  }

  /**
   * Evaluate the dense methane reference contribution.
   *
   * @param t absolute temperature in kelvin
   * @param massDensity methane mass density in g/cm3
   * @param q dense-reference correlation coefficients
   * @return evaluated value
   */
  public static double denseContribution(double t, double massDensity, double[] q) {
    double theta = (massDensity - RHOC * MW / 1000) / (RHOC * MW / 1000);
    double e = Math.pow(massDensity, .1) * (q[1] + q[2] / Math.pow(t, 1.5))
        + theta * Math.sqrt(massDensity) * (q[4] + q[5] / t + q[6] / (t * t));
    return Math.exp(q[0] + q[3] / t) * Math.expm1(e);
  }

  /**
   * Evaluate methane reference viscosity in Pa s.
   *
   * @param t absolute temperature in kelvin
   * @param pAtm pressure in atmospheres
   * @return evaluated value
   */
  public static double referenceViscosity(double t, double pAtm) {
    return referenceViscosity(t, pAtm, 90.69);
  }

  /**
   * Evaluate methane reference viscosity in Pa s.
   *
   * @param t absolute temperature in kelvin
   * @param pAtm pressure in atmospheres
   * @param transitionK reference blend center in kelvin
   * @return evaluated value
   */
  public static double referenceViscosity(double t, double pAtm, double transitionK) {
    double rho = density(t, pAtm) * MW / 1000, dilute = 0;
    for (int i = 0; i < GV.length; i++)
      dilute += GV[i] * Math.pow(t, (i - 3) / 3.0);
    double linear = (1.696985927 - .133372346 * Math.pow(1.4 - Math.log(t / 168), 2)) * rho;
    double weight = (Math.tanh((t - transitionK) / 5) + 1) / 2;
    return (dilute + linear + weight * denseContribution(t, rho, J) + (1 - weight) * denseContribution(t, rho, K))
        * 1e-7;
  }

  /**
   * Evaluate the requested model and return its diagnostics.
   *
   * @param input explicit model and calculation settings
   * @param point temperature in degrees Celsius and pressure in MPa
   * @param phase phase type and mole-fraction composition
   * @return calculated values and diagnostics
   */
  public static ObjectNode calculate(JsonNode input, JsonNode point, JsonNode phase) {
    if (phase.path("type").asText().equals("aqueous"))
      throw new IllegalArgumentException("Aqueous viscosity law is not the hydrocarbon CSP branch");
    JsonNode components = input.path("components"), xs = phase.path("composition");
    if (!xs.isArray() || xs.size() != components.size() || xs.size() == 0) {
      throw new IllegalArgumentException("CSP composition size mismatch");
    }
    double compositionTotal = 0;
    for (JsonNode x : xs) {
      if (!x.isNumber() || !Double.isFinite(x.asDouble()) || x.asDouble() < 0) {
        throw new IllegalArgumentException("Invalid CSP phase composition");
      }
      compositionTotal += x.asDouble();
    }
    if (Math.abs(compositionTotal - 1) > 1e-8) {
      throw new IllegalArgumentException("CSP phase composition must sum to one");
    }
    double t = ExplicitFlash.number(point, "temperatureC") + 273.15;
    double pAtm = ExplicitFlash.number(point, "pressureMPa") * 1e6 / 101325;
    double mn = 0, mwNumerator = 0, sumV = 0, sumVT = 0;
    for (int i = 0; i < components.size(); i++) {
      JsonNode ci = components.get(i);
      double xi = xs.get(i).asDouble();
      double mni = ExplicitFlash.number(ci, "mwKgMol") * 1000;
      mn += xi * mni;
      mwNumerator += xi * mni * ExplicitFlash.number(ci, "mwWeightKgMol") * 1000;
      for (int j = 0; j < components.size(); j++) {
        JsonNode cj = components.get(j);
        double ti = ExplicitFlash.number(ci, "tcK"), tj = ExplicitFlash.number(cj, "tcK");
        double pi = ExplicitFlash.number(ci, "pcBar") / 1.01325, pj = ExplicitFlash.number(cj, "pcBar") / 1.01325;
        double v = xi * xs.get(j).asDouble() * Math.pow(Math.cbrt(ti / pi) + Math.cbrt(tj / pj), 3);
        sumV += v;
        sumVT += v * Math.sqrt(ti * tj);
      }
    }
    double mw = mwNumerator / mn, tcMix = sumVT / sumV, pcMix = 8 * sumVT / (sumV * sumV);
    JsonNode normal = input.path("transport").path("normalFactors"),
        heavy = input.path("transport").path("heavyFactors");
    if (normal.size() != 2 || heavy.size() != 2)
      throw new IllegalArgumentException("Explicit CSP factors required");
    for (JsonNode values : new JsonNode[] {normal, heavy})
      for (JsonNode f : values)
        if (!f.isNumber() || !Double.isFinite(f.asDouble()) || f.asDouble() <= 0)
          throw new IllegalArgumentException("CSP factors must be finite and positive");
    double power = 2.303 * normal.get(1).asDouble();
    double mMix = mn + 1.304e-4 * normal.get(0).asDouble() * (Math.pow(mw, power) - Math.pow(mn, power));
    String method = input.path("transport").path("nativeCandidate").asText();
    if (!(method.isEmpty() || method.equals("native-csp-mccarty1974-v1") || method.equals(DENSITY_FLOOR_METHOD)
        || method.equals(TRANSITION_METHOD))) {
      throw new IllegalArgumentException("Unknown CSP profile: " + method);
    }
    boolean sdkTransition = method.equals(TRANSITION_METHOD);
    boolean densityFloor = sdkTransition || method.equals(DENSITY_FLOOR_METHOD);
    double initialTemperature = t * TC / tcMix;
    // Independent heavy-factor perturbations isolate the SDK blend weight.
    // Below 100 K its implied alpha-density remains the BWR density at 100 K;
    // the final methane viscosity temperature is not clamped. Keep v1 intact.
    double densityTemperature = densityFloor ? Math.max(100, initialTemperature) : initialTemperature;
    double rho = density(densityTemperature, pAtm * PC_ATM / pcMix), reduced = rho / RHOC;
    double alpha = 1 + 7.378e-3 * Math.pow(reduced, 1.847) * Math.pow(mMix, .5173);
    double alpha0 = 1 + 7.378e-3 * Math.pow(reduced, 1.847) * Math.pow(MW, .5173);
    double t0 = t * TC / tcMix * alpha0 / alpha, p0 = pAtm * PC_ATM / pcMix * alpha0 / alpha;
    double v3 = .2252 * t / mn + .9738, v4 = .5354 * v3 - .1170;
    double mHeavy = mn
        * Math.pow(Math.max(1.5, mw / mn) / (v3 * heavy.get(0).asDouble()), v4 * heavy.get(1).asDouble());
    double sign = t > 564.49 ? -1 : 1;
    double muHeavy = Math.pow(10, -.07995 - sign * .01101 * mHeavy - 371.8 / t + sign * 6.215 * mHeavy / t)
        * Math.exp(.00384 * (Math.pow(pAtm, .822574) - 1) / .822574) * 1e-3;
    double mu = muHeavy, csp = Double.NaN;
    if (t0 >= 50) {
      // 95 K is a separately identified SDK parity hypothesis, not methane's
      // physical freezing point or a numerical value specified by the manual.
      csp = referenceViscosity(t0, p0, sdkTransition ? 95 : 90.69) * Math.pow(tcMix / TC, -1.0 / 6)
          * Math.pow(pcMix / PC_ATM, 2.0 / 3) * Math.sqrt(mMix / MW) * alpha / alpha0;
      double w = t0 > 75 ? 1 : (Math.tanh((t0 - 65) / 5) + 1) / 2;
      mu = w * csp + (1 - w) * muHeavy;
    }
    if (!Double.isFinite(mu) || mu <= 0)
      throw new IllegalStateException("Invalid native CSP candidate");
    ObjectNode out = ExplicitFlash.JSON.createObjectNode();
    out.put("method",
        sdkTransition ? TRANSITION_METHOD : densityFloor ? DENSITY_FLOOR_METHOD : "native-csp-mccarty1974-v1")
        .put("viscosityPaS", mu).put("numberAverageGmol", mn).put("weightAverageGmol", mw).put("cspMolarMassGmol", mMix)
        .put("heavyMolarMassGmol", mHeavy).put("mixtureCriticalTemperatureK", tcMix)
        .put("mixtureCriticalPressureAtm", pcMix).put("methaneTemperatureK", t0).put("methanePressureAtm", p0)
        .put("methaneInitialDensityMolL", rho).put("heavyViscosityPaS", muHeavy)
        .put("denseReferenceExtrapolation", t0 < 90.69).put("branch", t0 < 50 ? "heavy" : t0 > 75 ? "csp" : "blend");
    if (Double.isFinite(csp))
      out.put("cspViscosityPaS", csp);
    if (densityFloor)
      out.put("initialReferenceTemperatureK", initialTemperature).put("initialDensityTemperatureK", densityTemperature)
          .put("initialDensityTemperatureFloorK", 100);
    if (sdkTransition)
      out.put("referenceTransitionCenterK", 95).put("referenceTransitionBasis",
          "SDK pure-methane control hypothesis; not a physical freezing point");
    return out;
  }
}
