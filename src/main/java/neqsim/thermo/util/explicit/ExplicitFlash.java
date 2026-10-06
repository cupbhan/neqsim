package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.component.ComponentPR;
import neqsim.thermo.component.ComponentSrk;
import neqsim.thermo.component.attractiveeosterm.AttractiveTermPr1978;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.mixingrule.HVMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Explicit-component PR/SRK equilibrium with audited mixing rules, volume translation and NaCl mole-basis conversion.
 *
 * @author NeqSim development team
 * @version 1.0
 */
public final class ExplicitFlash {
  static final ObjectMapper JSON = new ObjectMapper();
  static final double R = neqsim.thermo.ThermodynamicConstantsInterface.R;

  /**
   * Cubic EOS component retaining an explicit OmegaB.
   *
   * @author NeqSim development team
   * @version 1.0
   */
  static final class NativeSrk extends ComponentSrk {
    private static final long serialVersionUID = 1000L;
    double omegaB;

    /**
     * Construct an explicitly parameterized component.
     *
     * @param template database allocation template
     * @param z mole amount
     * @param i component index
     * @param b dimensionless OmegaB
     */
    NativeSrk(String template, double z, int i, double b) {
      super(template, z, z, i);
      omegaB = b;
    }

    /**
     * Evaluate the explicitly specified EOS covolume coefficient.
     *
     * @return evaluated value
     */
    @Override
    public double calcb() {
      return omegaB * R * getTC() / getPC();
    }
  }

  /**
   * Cubic EOS component retaining an explicit OmegaB.
   *
   * @author NeqSim development team
   * @version 1.0
   */
  static final class NativePr extends ComponentPR {
    private static final long serialVersionUID = 1000L;
    double omegaB;

    /**
     * Construct an explicitly parameterized component.
     *
     * @param template database allocation template
     * @param z mole amount
     * @param i component index
     * @param b dimensionless OmegaB
     */
    NativePr(String template, double z, int i, double b) {
      super(template, z, z, i);
      omegaB = b;
    }

    /**
     * Evaluate the explicitly specified EOS covolume coefficient.
     *
     * @return evaluated value
     */
    @Override
    public double calcb() {
      return omegaB * R * getTC() / getPC();
    }
  }

  /**
   * Read a required finite numeric parameter.
   *
   * @param n JSON object containing the parameter
   * @param key parameter or property name
   * @return evaluated value
   */
  public static double number(JsonNode n, String key) {
    if (!n.path(key).isNumber() || !Double.isFinite(n.path(key).asDouble()))
      throw new IllegalArgumentException("Missing/nonfinite parameter " + key);
    return n.path(key).asDouble();
  }

  /**
   * Evaluate the requested model and return its diagnostics.
   *
   * @param input explicit model and calculation settings
   * @param point temperature in degrees Celsius and pressure in MPa
   * @return calculated values and diagnostics
   */
  public static ObjectNode calculate(JsonNode input, JsonNode point) {
    for (JsonNode component : input.path("components"))
      if (component.has("saltIonNumber"))
        return calculateBrine(input, point);
    return calculateEos(input, point, false);
  }

  /**
   * Flash explicit EOS parameters or evaluate a fixed aqueous state.
   *
   * @param input explicit model and calculation settings
   * @param point temperature in degrees Celsius and pressure in MPa
   * @param fixedAqueous whether to evaluate the specified aqueous composition without a flash
   * @return calculated values and diagnostics
   */
  public static ObjectNode calculateEos(JsonNode input, JsonNode point, boolean fixedAqueous) {
    double t = number(point, "temperatureC") + 273.15;
    double p = number(point, "pressureMPa") * 10;
    String eos = input.path("eos").asText();
    boolean srk = eos.startsWith("Srk");
    JsonNode components = input.path("components");
    SystemInterface fluid = srk ? new SystemSrkEos(t, p) : new SystemPrEos(t, p);
    // Unique placeholders allocate the native EOS phase arrays. Every thermodynamic
    // parameter used below is then supplied explicitly, including alpha and OmegaA/B.
    for (int i = 0; i < components.size(); i++)
      fluid.addTBPfraction("native_" + i, number(components.get(i), "z"), number(components.get(i), "mwKgMol"), 1.0); // Allocation
                                                                                                                      // only;
                                                                                                                      // replace
                                                                                                                      // component
                                                                                                                      // objects
                                                                                                                      // below.
    for (int ph = 0; ph < fluid.getMaxNumberOfPhases(); ph++) {
      if (fluid.getPhase(ph) == null)
        continue;
      for (int i = 0; i < components.size(); i++) {
        JsonNode a = components.get(i);
        double z = number(a, "z");
        String template = a.path("templateName").asText("methane");
        ComponentEos c = srk ? new NativeSrk(template, z, i, number(a, "omegaB"))
            : new NativePr(template, z, i, number(a, "omegaB"));
        c.setComponentName(a.path("engineName").asText(a.path("name").asText()));
        if (a.has("componentType"))
          c.setComponentType(a.path("componentType").asText());
        c.setTC(number(a, "tcK"));
        c.setPC(number(a, "pcBar"));
        c.setMolarMass(number(a, "mwKgMol"));
        c.setAcentricFactor(number(a, "acentricFactor"));
        if (a.has("vcCm3Mol"))
          c.setCriticalVolume(number(a, "vcCm3Mol"));
        c.setNormalBoilingPoint(number(a, "tbK"));
        c.setOmegaA(number(a, "omegaA"));
        c.setVolumeCorrection(0);
        c.setVolumeCorrectionT(0);
        if (eos.contains("78"))
          c.setAttractiveParameter(new AttractiveTermPr1978(c));
        else
          c.setAttractiveTerm(srk ? 0 : 1);
        if (a.path("alpha").asText().equals("MathiasCopeman")) {
          double[] coefficients = {number(a, "mc1"), number(a, "mc2"), number(a, "mc3")};
          if (srk) {
            c.setMatiascopemanParams(coefficients);
            c.setAttractiveTerm(4);
          } else {
            for (int j = 0; j < 3; j++)
              c.setMatiascopemanParamsPR(j, coefficients[j]);
            c.setAttractiveTerm(13);
          }
        }
        fluid.getPhase(ph).getcomponentArray()[i] = c;
      }
    }
    boolean useHv = input.path("mixingRule").asText().equals("HV");
    fluid.setMixingRule(useHv ? 4 : 2);
    fluid.useVolumeCorrection(false);
    fluid.setMultiPhaseCheck(true);
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
        double kij = number(pair, "kij") + number(pair, "kijT") * (t - 288.15);
        rule.setBinaryInteractionParameter(i, j, kij);
        if (useHv && pair.path("mixingRule").asText().equals("HuronVidal")) {
          HVMixingRulesInterface hv = (HVMixingRulesInterface) rule;
          hv.setClassicOrHV(i, j, "HV");
          hv.setClassicOrHV(j, i, "HV");
          hv.setHValphaParameter(i, j, number(pair, "alpha"));
          hv.setHVDijParameter(i, j, number(pair, "gij"));
          hv.setHVDijParameter(j, i, number(pair, "gji"));
          hv.setHVDijTParameter(i, j, number(pair, "gijT"));
          hv.setHVDijTParameter(j, i, number(pair, "gjiT"));
        }
      }
    }
    if (fixedAqueous) {
      fluid.setNumberOfPhases(1);
      fluid.setBeta(0, 1.0);
      fluid.setPhaseType(0, neqsim.thermo.phase.PhaseType.AQUEOUS);
    } else
      new ThermodynamicOperations(fluid).TPflash();
    fluid.init(1);
    ObjectNode result = (ObjectNode) point.deepCopy();
    result.put("status", "ok");
    ArrayNode phases = result.putArray("phases");
    double balanceError = 0;
    double fugacityError = 0;
    for (int i = 0; i < components.size(); i++) {
      double actual = 0;
      for (int ph = 0; ph < fluid.getNumberOfPhases(); ph++)
        actual += fluid.getBeta(ph) * fluid.getPhase(ph).getComponent(i).getx();
      balanceError = Math.max(balanceError, Math.abs(actual - number(components.get(i), "z")));
      // Constrained nonvolatile salt has no gas/oil fugacity-equality equation.
      if (!fluid.getPhase(0).getComponent(i).isIsIon() && number(components.get(i), "z") > 1e-12)
        for (int ph = 1; ph < fluid.getNumberOfPhases(); ph++) {
          ComponentInterface first = fluid.getPhase(0).getComponent(i);
          ComponentInterface other = fluid.getPhase(ph).getComponent(i);
          double residual = Math.abs(Math
              .log((other.getx() * other.getFugacityCoefficient()) / (first.getx() * first.getFugacityCoefficient())));
          if (!Double.isFinite(residual))
            throw new IllegalStateException("Invalid fugacity residual");
          fugacityError = Math.max(fugacityError, residual);
        }
    }
    result.put("maxMaterialBalanceError", balanceError);
    result.put("maxLogFugacityResidual", fugacityError);
    if (!Double.isFinite(balanceError) || balanceError > 1e-7)
      throw new IllegalStateException("Material balance failed: " + balanceError);
    if (fugacityError > 1e-6)
      throw new IllegalStateException("Fugacity equilibrium failed: " + fugacityError);
    for (int ph = 0; ph < fluid.getNumberOfPhases(); ph++) {
      PhaseInterface phase = fluid.getPhase(ph);
      double volumeM3Mol = phase.getMolarVolume() * 1e-5;
      double correctionM3Mol = 0;
      for (int i = 0; i < components.size(); i++) {
        JsonNode c = components.get(i);
        correctionM3Mol += phase.getComponent(i).getx()
            * (number(c, "penelouxLmol") + number(c, "penelouxTLmolK") * (t - 288.15)) * 1e-3;
      }
      double translatedVolume = volumeM3Mol - correctionM3Mol;
      if (!(translatedVolume > 0))
        throw new IllegalStateException("Nonpositive translated volume");
      ObjectNode row = phases.addObject();
      row.put("type", phase.getType().toString().toLowerCase(java.util.Locale.ROOT));
      row.put("moleFraction", fluid.getBeta(ph));
      row.put("densityKgM3", phase.getMolarMass() / translatedVolume);
      row.put("molarMassKgMol", phase.getMolarMass());
      row.put("z", translatedVolume * p * 1e5 / (R * t));
      row.put("untranslatedZ", phase.getZ());
      row.put("mixtureA", phase.getA() * p / Math.pow(phase.getNumberOfMolesInPhase() * R * t, 2));
      row.put("mixtureB", phase.getB() * p / (phase.getNumberOfMolesInPhase() * R * t));
      EosMixingRulesInterface actualRule = ((PhaseEosInterface) phase).getEosMixingRule();
      for (JsonNode pair : input.path("interactions")) {
        int i = pair.path("i").asInt(), j = pair.path("j").asInt();
        double expectedKij = number(pair, "kij") + number(pair, "kijT") * (t - 288.15);
        if (Math.abs(actualRule.getBinaryInteractionParameter(i, j) - expectedKij) > 1e-12)
          throw new IllegalStateException("Post-flash kij differs at phase " + ph + " pair " + i + "," + j);
        if (useHv && pair.path("mixingRule").asText().equals("HuronVidal")) {
          HVMixingRulesInterface hv = (HVMixingRulesInterface) actualRule;
          if (Math.abs(hv.getHValphaParameter(i, j) - number(pair, "alpha")) > 1e-12
              || Math.abs(hv.getHVDijParameter(i, j) - number(pair, "gij")) > 1e-12
              || Math.abs(hv.getHVDijParameter(j, i) - number(pair, "gji")) > 1e-12
              || Math.abs(hv.getHVDijTParameter(i, j) - number(pair, "gijT")) > 1e-12
              || Math.abs(hv.getHVDijTParameter(j, i) - number(pair, "gjiT")) > 1e-12)
            throw new IllegalStateException("Post-flash HV parameters differ");
        }
      }
      ArrayNode fractions = row.putArray("composition");
      for (int i = 0; i < components.size(); i++)
        fractions.add(phase.getComponent(i).getx());
    }
    // Audit post-flash parameters so initialization/clone regressions cannot silently replace inputs.
    ArrayNode audit = result.putArray("parameterAudit");
    for (int i = 0; i < components.size(); i++) {
      ComponentEos c = (ComponentEos) fluid.getPhase(0).getComponent(i);
      JsonNode expected = components.get(i);
      if (expected.path("componentType").asText().equals("ion") && !c.isIsIon())
        throw new IllegalStateException("Nonvolatile salt constraint was lost");
      double oa = c.calca() * c.getPC() / Math.pow(R * c.getTC(), 2);
      double ob = c.calcb() * c.getPC() / (R * c.getTC());
      if (Math.abs(oa - number(expected, "omegaA")) > 1e-12 || Math.abs(ob - number(expected, "omegaB")) > 1e-12)
        throw new IllegalStateException("EOS coefficient audit failed");
      audit.addObject().put("name", c.getComponentName()).put("tcK", c.getTC()).put("pcBar", c.getPC())
          .put("omegaA", oa).put("omegaB", ob).put("alphaAtT", c.alpha(t));
    }
    return result;
  }

  /**
   * Convert formula-unit salt feed to ionic EOS moles and recover formula-unit phase results.
   *
   * @param input explicit model and calculation settings
   * @param point temperature in degrees Celsius and pressure in MPa
   * @return calculated values and diagnostics
   */
  public static ObjectNode calculateBrine(JsonNode input, JsonNode point) {
    JsonNode components = input.path("components");
    int salt = -1, water = -1;
    for (int i = 0; i < components.size(); i++) {
      JsonNode c = components.get(i);
      if (c.path("name").asText().equals("H2O"))
        water = i;
      if (c.has("saltIonNumber")) {
        if (salt >= 0 || !c.path("name").asText().equals("NaCl") || !c.path("templateName").asText().equals("NaCl")
            || !c.path("componentType").asText().equals("ion") || number(c, "saltIonNumber") != 2
            || number(c, "crystalWater") != 0)
          throw new IllegalArgumentException("Unverified salt specification");
        salt = i;
      }
    }
    double tc = number(point, "temperatureC"), pmpa = number(point, "pressureMPa");
    if (salt < 0 || water < 0 || !input.path("eos").asText().equals("PengRobinson78Peneloux")
        || !input.path("mixingRule").asText().equals("HV") || tc < 25 || tc > 150 || pmpa < 1 || pmpa > 20)
      throw new IllegalArgumentException("NaCl benchmark outside verified model/range");
    double expansion = 1 + number(components.get(salt), "z");
    ObjectNode ionic = input.deepCopy();
    for (int i = 0; i < components.size(); i++) {
      ObjectNode c = (ObjectNode) ionic.path("components").get(i);
      double nu = i == salt ? 2 : 1;
      c.put("z", number(c, "z") * nu / expansion);
      c.put("mwKgMol", number(c, "mwKgMol") / nu);
    }
    ObjectNode result = calculateEos(ionic, point, false);
    result.put("internalMoleBasis", "salt ions; returned phases use salt formula units");
    result.put("ionicMolesPerFormulaFeedMole", expansion);
    result.put("maxIonicMaterialBalanceError", result.path("maxMaterialBalanceError").asDouble());
    double balance = 0, outsideSalt = 0;
    double[] recovered = new double[components.size()];
    for (JsonNode phase : result.path("phases")) {
      ObjectNode row = (ObjectNode) phase;
      ArrayNode x = (ArrayNode) row.path("composition");
      double reduction = 1 - x.get(salt).asDouble() / 2;
      row.put("moleFraction", number(row, "moleFraction") * reduction * expansion);
      row.put("molarMassKgMol", number(row, "molarMassKgMol") / reduction);
      row.put("z", number(row, "z") / reduction);
      row.put("untranslatedZBasis", "ionic EOS mole basis");
      row.put("mixtureABBasis", "ionic EOS mole basis");
      for (int i = 0; i < components.size(); i++) {
        x.set(i, JSON.getNodeFactory().numberNode(x.get(i).asDouble() / reduction / (i == salt ? 2 : 1)));
        recovered[i] += number(row, "moleFraction") * x.get(i).asDouble();
      }
      if (!row.path("type").asText().equals("aqueous")) {
        outsideSalt += number(row, "moleFraction") * x.get(salt).asDouble();
        continue;
      }
      if (!(x.get(water).asDouble() > 0))
        throw new IllegalStateException("Salt phase contains no water");
      double ratio = x.get(salt).asDouble() / x.get(water).asDouble();
      double tk = tc + 273.15;
      double limit = tk < 382.98 ? 0.07986 + 0.0001048 * tk : 0.01506 + 0.000274 * tk;
      ObjectNode audit = row.putObject("saltAudit");
      audit.put("saltWaterMoleRatio", ratio);
      audit.put("solubilityLimit", limit);
      if (ratio > limit) {
        result.put("status", "failed");
        result.put("failureCode", "SALT_SUPERSATURATED");
        result.put("error", "Aqueous NaCl/water mole ratio " + ratio + " exceeds solubility limit " + limit);
      }
      // Rebuild the salt-free liquid with identical EOS parameters at its fixed
      // dissolved-gas composition. No second phase equilibrium calculation.
      ObjectNode free = input.deepCopy();
      ArrayNode freeComponents = free.putArray("components"), pairs = free.putArray("interactions");
      double nonsalt = 1 - x.get(salt).asDouble();
      for (int i = 0; i < components.size(); i++)
        if (i != salt) {
          ObjectNode c = components.get(i).deepCopy();
          c.put("z", x.get(i).asDouble() / nonsalt);
          freeComponents.add(c);
        }
      for (JsonNode pair : input.path("interactions")) {
        int i = pair.path("i").asInt(), j = pair.path("j").asInt();
        if (i == salt || j == salt)
          continue;
        ObjectNode q = pair.deepCopy();
        q.put("i", i > salt ? i - 1 : i);
        q.put("j", j > salt ? j - 1 : j);
        pairs.add(q);
      }
      ObjectNode fixed = calculateEos(free, point, true);
      double rho0 = number(fixed.path("phases").get(0), "densityKgM3");
      double mass = 0;
      for (int i = 0; i < components.size(); i++)
        mass += x.get(i).asDouble() * number(components.get(i), "mwKgMol");
      double ws = x.get(salt).asDouble() * number(components.get(salt), "mwKgMol") / mass;
      double vs = (ws + 1.0166 + 0.014624 * tc)
          / ((-0.00433 * ws + 0.06471) * Math.exp(1e-6 * Math.pow(tc + 3315.6, 2)));
      double rho = 1 / ((1 - ws) / rho0 + ws * vs);
      if (!Double.isFinite(rho) || rho <= 0)
        throw new IllegalStateException("Invalid brine density");
      audit.put("saltFreeEosDensityKgM3", rho0);
      audit.put("saltMassFraction", ws);
      audit.put("partialSpecificVolumeM3Kg", vs);
      audit.put("saltFreeCompositionError", fixed.path("maxMaterialBalanceError").asDouble());
      row.put("densityKgM3", rho);
      row.put("z", mass / rho * pmpa * 1e6 / (R * tk));
    }
    for (int i = 0; i < components.size(); i++)
      balance = Math.max(balance, Math.abs(recovered[i] - number(components.get(i), "z")));
    if (!Double.isFinite(balance) || balance > 1e-7 || outsideSalt > 1e-12)
      throw new IllegalStateException("Salt formula-unit balance/aqueous retention failed");
    result.put("maxMaterialBalanceError", balance);
    result.put("saltOutsideAqueousPerFeedMole", outsideSalt);
    return result;
  }
}
