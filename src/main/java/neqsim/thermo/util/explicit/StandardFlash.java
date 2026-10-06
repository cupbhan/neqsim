package neqsim.thermo.util.explicit;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.ObjectNode;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Database-component PR/SRK flash with constant symmetric interactions and audited phase properties.
 *
 * @author NeqSim development team
 * @version 1.0
 */
public final class StandardFlash {
  private static final ObjectMapper JSON = new ObjectMapper();

  /**
   * Enforce a calculation invariant.
   *
   * @param ok condition that must hold
   * @param message failure description
   */
  private static void require(boolean ok, String message) {
    if (!ok)
      throw new IllegalStateException(message);
  }

  /**
   * Read the cubic EOS mixing rule.
   *
   * @param phase phase type and mole-fraction composition
   * @return evaluated value
   */
  private static EosMixingRulesInterface mixing(PhaseInterface phase) {
    return ((PhaseEosInterface) phase).getEosMixingRule();
  }

  /**
   * Evaluate the requested model and return its diagnostics.
   *
   * @param input explicit model and calculation settings
   * @return calculated values and diagnostics
   */
  public static ObjectNode calculate(JsonNode input) {
    ObjectNode response = JSON.createObjectNode();
    try {
      double temperature = input.path("temperatureC").asDouble() + 273.15;
      double pressure = input.path("pressureMPa").asDouble() * 10.0;
      SystemInterface fluid = input.path("eos").asText().equals("PR") ? new SystemPrEos(temperature, pressure)
          : new SystemSrkEos(temperature, pressure);
      JsonNode components = input.path("components");
      for (JsonNode c : components)
        fluid.addComponent(c.path("nativeName").asText(), c.path("moleFraction").asDouble());
      fluid.setMixingRule(2); // CLASSIC: constant symmetric k_ij, never temperature-dependent defaults.
      fluid.useVolumeCorrection(false);
      fluid.setMultiPhaseCheck(false); // This adapter accepts dry, nonpolar gas/oil systems only.
      fluid.init(0);
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        for (int i = 0; i < components.size(); i++)
          fluid.getPhase(p).getComponent(i).setUseIdealGasEnthalpyOfFormation(false);
        EosMixingRulesInterface rule = mixing(fluid.getPhase(p));
        for (JsonNode pair : input.path("interactions"))
          rule.setBinaryInteractionParameter(pair.path("i").asInt(), pair.path("j").asInt(),
              pair.path("value").asDouble());
      }
      new ThermodynamicOperations(fluid).TPflash();
      fluid.init(3);
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        fluid.getPhase(p).getPhysicalProperties().setViscosityModel("PFCT");
        fluid.getPhase(p).getPhysicalProperties().setConductivityModel("PFCT");
      }
      fluid.initProperties();

      ArrayNode bindings = response.putArray("components");
      for (int i = 0; i < components.size(); i++) {
        ObjectNode binding = bindings.addObject();
        binding.put("id", components.get(i).path("id").asText());
        binding.put("nativeName", fluid.getPhase(0).getComponent(i).getComponentName());
        binding.put("inputMolarMass", components.get(i).path("molarMass").asDouble());
        binding.put("engineMolarMass", fluid.getPhase(0).getComponent(i).getMolarMass() * 1000.0);
        binding.put("criticalTemperatureK", fluid.getPhase(0).getComponent(i).getTC());
        binding.put("criticalPressureBar", fluid.getPhase(0).getComponent(i).getPC());
        binding.put("acentricFactor", fluid.getPhase(0).getComponent(i).getAcentricFactor());
        binding.put("criticalVolumeM3Mol", fluid.getPhase(0).getComponent(i).getCriticalVolume() * 1.0e-6);
        binding.put("parachor", fluid.getPhase(0).getComponent(i).getParachorParameter());
        ArrayNode cp = binding.putArray("idealGasCpCoefficients");
        cp.add(fluid.getPhase(0).getComponent(i).getCpA());
        cp.add(fluid.getPhase(0).getComponent(i).getCpB());
        cp.add(fluid.getPhase(0).getComponent(i).getCpC());
        cp.add(fluid.getPhase(0).getComponent(i).getCpD());
        cp.add(fluid.getPhase(0).getComponent(i).getCpE());
      }
      ArrayNode phases = response.putArray("phases");
      double[] balance = new double[components.size()];
      double maxFugacityResidual = 0.0;
      double betaTotal = 0.0;
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        String type = phase.getType().toString().toLowerCase(java.util.Locale.ROOT);
        if (type.equals("liquid"))
          type = "oil";
        require(type.equals("gas") || type.equals("oil"), "Unsupported phase: " + type);
        require(!phase.useVolumeCorrection(), "Volume correction unexpectedly enabled");
        for (JsonNode pair : input.path("interactions")) {
          int i = pair.path("i").asInt(), j = pair.path("j").asInt();
          double k = pair.path("value").asDouble();
          require(
              Math.abs(mixing(phase).getBinaryInteractionParameter(i, j) - k) < 1e-14
                  && Math.abs(mixing(phase).getBinaryInteractionParameter(j, i) - k) < 1e-14,
              "BIP not applied symmetrically");
        }
        double beta = fluid.getBeta(p);
        require(Double.isFinite(beta) && beta >= 0 && beta <= 1, "Invalid phase fraction");
        betaTotal += beta;
        ObjectNode result = phases.addObject();
        result.put("phase", type);
        result.put("moleFraction", beta);
        double density = phase.getDensity("kg/m3");
        double viscosity = phase.getViscosity("cP"); // cP == mPa.s
        double enthalpy = phase.getEnthalpy("kJ/kg");
        require(Double.isFinite(density) && density > 0 && Double.isFinite(viscosity) && viscosity > 0
            && Double.isFinite(enthalpy), "Non-finite or nonphysical properties");
        result.put("density", density);
        result.put("viscosity", viscosity);
        result.put("enthalpy", enthalpy);
        double cp = phase.getCp("J/kgK");
        double z = phase.getZ();
        if (Double.isFinite(cp) && cp > 0)
          result.put("heatCapacity", cp);
        else
          result.putNull("heatCapacity");
        if (Double.isFinite(z) && z > 0)
          result.put("compressibilityFactor", z);
        else
          result.putNull("compressibilityFactor");
        result.putArray("unavailableProperties");
        ArrayNode composition = result.putArray("composition");
        double xTotal = 0.0;
        for (int i = 0; i < components.size(); i++) {
          double x = phase.getComponent(i).getx();
          require(!phase.getComponent(i).isUsingIdealGasEnthalpyOfFormation(),
              "Unexpected formation enthalpy reference");
          require(Double.isFinite(x) && x >= 0 && x <= 1, "Invalid phase composition");
          composition.addObject().put("id", components.get(i).path("id").asText()).put("moleFraction", x);
          xTotal += x;
          balance[i] += beta * x;
          if (p > 0 && components.get(i).path("moleFraction").asDouble() > 1e-12) {
            double f0 = fluid.getPhase(0).getComponent(i).getx()
                * fluid.getPhase(0).getComponent(i).getFugacityCoefficient();
            double fp = x * phase.getComponent(i).getFugacityCoefficient();
            double residual = Math.abs(Math.log(fp / f0));
            require(Double.isFinite(residual), "Invalid fugacity residual");
            maxFugacityResidual = Math.max(maxFugacityResidual, residual);
          }
        }
        require(Math.abs(xTotal - 1) < 1e-8, "Phase composition does not sum to one");
      }
      double maxBalanceResidual = 0;
      double totalInput = 0;
      for (JsonNode c : components)
        totalInput += c.path("moleFraction").asDouble();
      for (int i = 0; i < components.size(); i++)
        maxBalanceResidual = Math.max(maxBalanceResidual,
            Math.abs(balance[i] - components.get(i).path("moleFraction").asDouble() / totalInput));
      require(Math.abs(betaTotal - 1) < 1e-9 && maxBalanceResidual < 1e-7 && maxFugacityResidual < 1e-6,
          "Flash did not pass material balance / fugacity checks");
      response.putObject("diagnostics").put("maxMaterialBalanceResidual", maxBalanceResidual)
          .put("maxLogFugacityResidual", maxFugacityResidual).put("fugacityCheckApplicable", phases.size() > 1)
          .put("inputMoleFractionTotal", totalInput).put("volumeCorrection", false).put("mixingRule", "CLASSIC (2)")
          .put("viscosityModel", "PFCT").put("conductivityModel", "PFCT")
          .put("javaVersion", System.getProperty("java.version"));
      response.set("interactions", input.path("interactions"));
      response.put("ok", true);
    } catch (Exception e) {
      response.removeAll();
      response.put("ok", false).put("code", "ENGINE_CALCULATION_FAILED").put("message",
          e.getClass().getSimpleName() + ": " + e.getMessage());
    }
    return response;
  }
}
