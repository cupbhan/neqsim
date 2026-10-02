package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import com.google.gson.JsonObject;
import java.io.Serializable;
import java.util.EnumMap;
import java.util.Map;
import neqsim.thermo.mixingrule.EosMixingRulesInterface;
import neqsim.thermo.phase.PhaseEosInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;

/**
 * Opt-in model profile for physically labelled hydrocarbon-water boundary families.
 *
 * <p>
 * The profile clones the caller's fluid and applies a water binary-interaction multiplier only to the requested
 * boundary family. The caller's fluid and the default NeqSim envelope path are never modified. Family-specific
 * multipliers are intended for explicitly governed compatibility studies; they are not a universal EOS calibration and
 * must not be promoted without independent physical validation.
 *
 *
 * @author NeqSim contributors
 * @version 1.0
 */
public final class HydrocarbonWaterBoundaryModelProfile implements Serializable {
  private static final long serialVersionUID = 1000L;

  private final String name;
  private final Map<BoundaryFamily, Double> waterKijMultipliers = new EnumMap<BoundaryFamily, Double>(
      BoundaryFamily.class);

  /**
   * Creates an identity profile with a caller-supplied audit name.
   *
   * @param name non-empty profile name
   */
  public HydrocarbonWaterBoundaryModelProfile(String name) {
    if (name == null || name.trim().isEmpty()) {
      throw new IllegalArgumentException("a non-empty boundary model profile name is required");
    }
    this.name = name.trim();
  }

  /**
   * Sets the water-Kij multiplier for one physical boundary family.
   *
   * @param family physical hydrocarbon-water boundary family
   * @param multiplier positive finite multiplier
   * @return this profile
   */
  public HydrocarbonWaterBoundaryModelProfile setWaterKijMultiplier(BoundaryFamily family, double multiplier) {
    if (family == null) {
      throw new IllegalArgumentException("a boundary family is required");
    }
    if (Double.isNaN(multiplier) || Double.isInfinite(multiplier) || multiplier <= 0.0) {
      throw new IllegalArgumentException("the water-Kij multiplier must be positive and finite");
    }
    waterKijMultipliers.put(family, multiplier);
    return this;
  }

  /**
   * Returns the configured multiplier, or unity when the family is not overridden.
   *
   * @param family physical boundary family
   * @return water-Kij multiplier
   */
  public double getWaterKijMultiplier(BoundaryFamily family) {
    if (family == null) {
      throw new IllegalArgumentException("a boundary family is required");
    }
    Double multiplier = waterKijMultipliers.get(family);
    return multiplier == null ? 1.0 : multiplier.doubleValue();
  }

  /**
   * Returns a cloned family-specific fluid without modifying the supplied baseline.
   *
   * @param baseline fully configured EOS fluid containing water
   * @param family requested physical boundary family
   * @return cloned fluid with the configured water-Kij multiplier applied
   */
  public SystemInterface createTemplate(SystemInterface baseline, BoundaryFamily family) {
    if (baseline == null) {
      throw new IllegalArgumentException("a baseline fluid is required");
    }
    if (!baseline.hasComponent("water")) {
      throw new IllegalArgumentException("the baseline fluid must contain water");
    }
    double multiplier = getWaterKijMultiplier(family);
    SystemInterface candidate = baseline.clone();
    if (Math.abs(multiplier - 1.0) <= 1.0e-15) {
      return candidate;
    }

    int waterIndex = baseline.getPhase(0).getComponent("water").getComponentNumber();
    for (int phaseIndex = 0; phaseIndex < candidate.getMaxNumberOfPhases(); phaseIndex++) {
      PhaseInterface sourcePhase = baseline.getPhase(phaseIndex);
      PhaseInterface targetPhase = candidate.getPhase(phaseIndex);
      if (!(sourcePhase instanceof PhaseEosInterface) || !(targetPhase instanceof PhaseEosInterface)) {
        continue;
      }
      EosMixingRulesInterface sourceRule = ((PhaseEosInterface) sourcePhase).getEosMixingRule();
      EosMixingRulesInterface targetRule = ((PhaseEosInterface) targetPhase).getEosMixingRule();
      for (int componentIndex = 0; componentIndex < sourcePhase.getNumberOfComponents(); componentIndex++) {
        if (componentIndex == waterIndex) {
          continue;
        }
        double baselineKij = sourceRule.getBinaryInteractionParameter(waterIndex, componentIndex);
        double scaledKij = baselineKij * multiplier;
        if (!Double.isFinite(scaledKij)) {
          throw new IllegalArgumentException("scaled water binary interaction must remain finite");
        }
        targetRule.setBinaryInteractionParameter(waterIndex, componentIndex, scaledKij);
      }
    }
    candidate.init(0);
    return candidate;
  }

  /**
   * Describes the active family-specific profile for audit output.
   *
   * @param family requested physical boundary family
   * @return JSON description
   */
  public JsonObject toJson(BoundaryFamily family) {
    JsonObject output = new JsonObject();
    output.addProperty("name", name);
    output.addProperty("boundaryFamily", family.name());
    output.addProperty("waterKijMultiplier", getWaterKijMultiplier(family));
    output.addProperty("identity", Math.abs(getWaterKijMultiplier(family) - 1.0) <= 1.0e-15);
    output.addProperty("defaultModelModified", false);
    output.addProperty("experimentalCompatibilityProfile", true);
    output.addProperty("requiresIndependentPhysicalValidation", true);
    output.addProperty("productionPromotionAllowed", false);
    output.addProperty("calibrationBasis", "user-supplied research multiplier; not experimental validation");
    return output;
  }
}
