package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import neqsim.physicalproperties.system.PhysicalProperties;
import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.component.attractiveeosterm.AttractiveTermPr1978;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.mixingrule.HVMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;

/**
 * Fixed-composition enthalpy derivatives and transport diagnostics for explicit fluid models.
 *
 * @author NeqSim development team
 * @version 1.0
 */
public final class ExplicitProperties {
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
   * Create an isolated fixed-composition system with explicitly bound heat capacities.
   *
   * @param input explicit model and calculation settings
   * @param composition component mole fractions in model order
   * @param type gas, oil or aqueous phase
   * @param t absolute temperature in kelvin
   * @param p pressure in bar
   * @return isolated thermodynamic system
   */
  public static SystemInterface state(JsonNode input, JsonNode composition, String type, double t, double p) {
    String eos = input.path("eos").asText();
    boolean srk = eos.startsWith("Srk");
    JsonNode components = input.path("components");
    SystemInterface fluid = srk ? new SystemSrkEos(t, p) : new SystemPrEos(t, p);
    for (int i = 0; i < components.size(); i++) {
      JsonNode a = components.get(i);
      if (a.has("saltIonNumber"))
        throw new IllegalArgumentException("Salt thermal/transport basis is unverified");
      fluid.addTBPfraction("probe_" + i, composition.get(i).asDouble(), n(a, "mwKgMol"), 1.0);
    }
    for (int ph = 0; ph < fluid.getMaxNumberOfPhases(); ph++) {
      if (fluid.getPhase(ph) == null)
        continue;
      for (int i = 0; i < components.size(); i++) {
        JsonNode a = components.get(i);
        double x = composition.get(i).asDouble();
        String template = a.path("templateName").asText("methane");
        ComponentEos c = srk ? new ExplicitFlash.NativeSrk(template, x, i, n(a, "omegaB"))
            : new ExplicitFlash.NativePr(template, x, i, n(a, "omegaB"));
        c.setComponentName(a.path("engineName").asText(a.path("name").asText()));
        if (a.has("componentType"))
          c.setComponentType(a.path("componentType").asText());
        c.setTC(n(a, "tcK"));
        c.setPC(n(a, "pcBar"));
        c.setMolarMass(n(a, "mwKgMol"));
        c.setAcentricFactor(n(a, "acentricFactor"));
        c.setCriticalVolume(n(a, "vcCm3Mol"));
        c.setNormalBoilingPoint(n(a, "tbK"));
        c.setOmegaA(n(a, "omegaA"));
        c.setVolumeCorrection(0);
        c.setVolumeCorrectionT(0);
        if (eos.contains("78"))
          c.setAttractiveParameter(new AttractiveTermPr1978(c));
        else
          c.setAttractiveTerm(srk ? 0 : 1);
        if (a.path("alpha").asText().equals("MathiasCopeman")) {
          double[] mc = {n(a, "mc1"), n(a, "mc2"), n(a, "mc3")};
          if (srk) {
            c.setMatiascopemanParams(mc);
            c.setAttractiveTerm(4);
          } else {
            for (int j = 0; j < 3; j++)
              c.setMatiascopemanParamsPR(j, mc[j]);
            c.setAttractiveTerm(13);
          }
        }
        JsonNode cp = a.path("cpJMolK");
        if (cp.size() != 4)
          throw new IllegalArgumentException("Four explicit Cp coefficients required");
        for (JsonNode coefficient : cp) {
          if (!coefficient.isNumber() || !Double.isFinite(coefficient.asDouble())) {
            throw new IllegalArgumentException("Cp coefficients must be finite numbers");
          }
        }
        c.setCpA(cp.get(0).asDouble());
        c.setCpB(cp.get(1).asDouble());
        c.setCpC(cp.get(2).asDouble());
        // Pinned rc.6 Component.getCpD() adds 1e-10. Bind the observable coefficient.
        c.setCpD(cp.get(3).asDouble() - 1e-10);
        c.setCpE(0);
        c.setUseIdealGasEnthalpyOfFormation(false);
        fluid.getPhase(ph).getcomponentArray()[i] = c;
      }
    }
    boolean useHv = input.path("mixingRule").asText().equals("HV");
    fluid.setMixingRule(useHv ? 4 : 2);
    fluid.useVolumeCorrection(false);
    fluid.init(0);
    for (int ph = 0; ph < fluid.getMaxNumberOfPhases(); ph++) {
      if (fluid.getPhase(ph) == null)
        continue;
      EosMixingRulesInterface rule = ((PhaseEosInterface) fluid.getPhase(ph)).getEosMixingRule();
      if (useHv)
        for (int i = 0; i < components.size(); i++)
          for (int j = 0; j < components.size(); j++) {
            HVMixingRulesInterface hv = (HVMixingRulesInterface) rule;
            hv.setClassicOrHV(i, j, "Classic");
            hv.setHVDijParameter(i, j, 0);
            hv.setHVDijTParameter(i, j, 0);
            hv.setHValphaParameter(i, j, 0);
          }
      for (JsonNode pair : input.path("interactions")) {
        int i = pair.path("i").asInt(), j = pair.path("j").asInt();
        rule.setBinaryInteractionParameter(i, j, n(pair, "kij") + n(pair, "kijT") * (t - 288.15));
        if (useHv && pair.path("mixingRule").asText().equals("HuronVidal")) {
          HVMixingRulesInterface hv = (HVMixingRulesInterface) rule;
          hv.setClassicOrHV(i, j, "HV");
          hv.setClassicOrHV(j, i, "HV");
          hv.setHValphaParameter(i, j, n(pair, "alpha"));
          hv.setHVDijParameter(i, j, n(pair, "gij"));
          hv.setHVDijParameter(j, i, n(pair, "gji"));
          hv.setHVDijTParameter(i, j, n(pair, "gijT"));
          hv.setHVDijTParameter(j, i, n(pair, "gjiT"));
        }
      }
    }
    fluid.setNumberOfPhases(1);
    fluid.setBeta(0, 1.0);
    fluid.setPhaseType(0, type.equals("gas") ? PhaseType.GAS : type.equals("oil") ? PhaseType.OIL : PhaseType.AQUEOUS);
    fluid.init(3);
    return fluid;
  }

  /**
   * Evaluate perturbed log fugacity on a clone while rebinding temperature-dependent interactions.
   *
   * @param base unmodified base system
   * @param input explicit model and calculation settings
   * @param t absolute temperature in kelvin
   * @return evaluated value
   */
  public static double evaluatePhi(SystemInterface base, JsonNode input, double t) {
    SystemInterface fluid = base.clone();
    fluid.setTemperature(t);
    for (int ph = 0; ph < fluid.getMaxNumberOfPhases(); ph++) {
      if (fluid.getPhase(ph) == null)
        continue;
      EosMixingRulesInterface rule = ((PhaseEosInterface) fluid.getPhase(ph)).getEosMixingRule();
      for (JsonNode pair : input.path("interactions"))
        rule.setBinaryInteractionParameter(pair.path("i").asInt(), pair.path("j").asInt(),
            n(pair, "kij") + n(pair, "kijT") * (t - 288.15));
    }
    fluid.init(3);
    return logPhi(fluid);
  }

  /**
   * Return the mole-fraction-weighted log fugacity coefficient.
   *
   * @param fluid thermodynamic system
   * @return evaluated value
   */
  public static double logPhi(SystemInterface fluid) {
    PhaseInterface phase = fluid.getPhase(0);
    double total = 0;
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      ComponentInterface c = phase.getComponent(i);
      total += c.getx() * c.getLogFugacityCoefficient();
    }
    return total;
  }

  /**
   * Evaluate fixed-composition density, enthalpy derivatives and transport candidates.
   *
   * @param input explicit model and calculation settings
   * @param point temperature in degrees Celsius and pressure in MPa
   * @param reference fixed phase type and composition
   * @return calculated values and diagnostics
   */
  public static ObjectNode phase(JsonNode input, JsonNode point, JsonNode reference) {
    double t = n(point, "temperatureC") + 273.15, p = n(point, "pressureMPa") * 10;
    String type = reference.path("type").asText();
    JsonNode x = reference.path("composition");
    SystemInterface fluid = state(input, x, type, t, p);
    PhaseInterface phase = fluid.getPhase(0);
    ObjectNode out = ExplicitFlash.JSON.createObjectNode();
    out.put("type", type);
    out.set("composition", x.deepCopy());
    double translation = 0, hCorrection = 0;
    for (int i = 0; i < x.size(); i++) {
      JsonNode c = input.path("components").get(i);
      translation += x.get(i).asDouble() * (n(c, "penelouxLmol") + n(c, "penelouxTLmolK") * (t - 288.15)) * 1e-3;
      hCorrection += x.get(i).asDouble() * (-n(c, "penelouxLmol") + 288.15 * n(c, "penelouxTLmolK")) * 1e-3 * p * 1e5;
    }
    double volume = phase.getZ() * ExplicitFlash.R * t / (p * 1e5) - translation;
    if (!(volume > 0))
      throw new IllegalStateException("Nonpositive translated volume");
    out.put("densityKgM3", phase.getMolarMass() / volume);
    out.put("z", volume * p * 1e5 / (ExplicitFlash.R * t));
    out.put("enthalpyBoundCpJMol", phase.getEnthalpy("J/mol"));
    ArrayNode audit = out.putArray("parameterAudit");
    for (int i = 0; i < x.size(); i++) {
      ComponentInterface c = phase.getComponent(i);
      JsonNode expected = input.path("components").get(i);
      double error = Math.abs(c.getx() - x.get(i).asDouble());
      if (error > 1e-10 || c.getMolarMass() != n(expected, "mwKgMol"))
        throw new IllegalStateException("Fixed phase composition or molar mass changed");
      double[] cp = {c.getCpA(), c.getCpB(), c.getCpC(), c.getCpD(), c.getCpE()};
      for (int j = 0; j < 5; j++)
        if (Math.abs(cp[j] - (j == 4 ? 0 : expected.path("cpJMolK").get(j).asDouble())) > 1e-20)
          throw new IllegalStateException("Bound Cp coefficient changed: " + i + ":" + j + " actual=" + cp[j]
              + " expected=" + (j == 4 ? 0 : expected.path("cpJMolK").get(j).asDouble()));
      audit.addObject().put("name", expected.path("name").asText()).put("idealEnthalpyJMol", c.getHID(t));
    }
    out.put("residualAnalyticFrozenKijJMol", phase.getHresTP() / phase.getNumberOfMolesInPhase());
    out.put("penelouxEnthalpyCorrectionJMol", hCorrection);
    for (double step : new double[] {0.01, 0.005}) {
      double derivative = (evaluatePhi(fluid, input, t + step) - evaluatePhi(fluid, input, t - step)) / (2 * step);
      out.put(step == 0.01 ? "residualFiniteDifferenceJMol" : "residualHalfStepJMol",
          -ExplicitFlash.R * t * t * derivative);
    }
    for (String method : new String[] {"PFCT", "PFCT-Heavy-Oil"}) {
      try {
        PhysicalProperties properties = phase.getPhysicalProperties();
        properties.setViscosityModel(method);
        double mu = properties.getViscosityModel().calcViscosity();
        if (!Double.isFinite(mu) || mu <= 0)
          throw new IllegalStateException("Invalid viscosity");
        out.put(method.equals("PFCT") ? "pfctViscosityPaS" : "pfctHeavyOilViscosityPaS", mu);
        out.put(method.equals("PFCT") ? "pfctClass" : "pfctHeavyOilClass",
            properties.getViscosityModel().getClass().getName());
      } catch (Exception error) {
        out.put(method + "Error", error.toString());
      }
    }
    if (input.path("propertyAuditVersion").asInt() >= 2) {
      try {
        out.set("nativeCsp", CspViscosity.calculate(input, point, reference));
      } catch (Exception error) {
        out.put("nativeCspError", error.toString());
      }
    }
    if (input.path("propertyAuditVersion").asInt() >= 3 && type.equals("aqueous")) {
      try {
        out.set("waterBaseline", PureWaterViscosity.calculate(input, point, reference));
      } catch (Exception error) {
        out.put("waterBaselineError", error.toString());
      }
    }
    return out;
  }

}
