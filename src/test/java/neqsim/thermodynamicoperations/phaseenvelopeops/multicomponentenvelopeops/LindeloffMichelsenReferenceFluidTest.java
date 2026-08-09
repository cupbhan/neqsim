package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Tag;
import neqsim.NeqSimTest;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos1978;
import neqsim.thermo.system.SystemSrkPenelouxEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Provisional reference-fluid fixtures from the unfinished Issue #530 branch.
 *
 * <p>
 * The compositions are source-traceable but are not yet accepted numerical truth. Expected PT boundary coordinates will
 * only be frozen after the Lindeloff-Michelsen full text and its tables/figures have been independently verified.
 * </p>
 */
class LindeloffMichelsenReferenceFluidTest extends NeqSimTest {

  @Test
  void fluidOneInitializesThreeConsecutiveLowPressureDewPointsWithSrkHv() throws Exception {
    SystemInterface reference = fluidOne(false);
    assertEquals("Huron-Vidal", reference.getMixingRuleName());
    for (int pressureIndex = 1; pressureIndex <= 3; pressureIndex++) {
      SystemInterface point = reference.clone();
      point.setPressure(pressureIndex * 1.01325);
      new ThermodynamicOperations(point).dewPointTemperatureFlash(false);
      assertTrue(
          Double.isFinite(point.getTemperature()) && point.getTemperature() > 100.0 && point.getTemperature() < 1200.0,
          "Fluid 1 low-pressure dew point must be finite");
    }
  }

  @Test
  void fluidTwoCanRunThreePhaseTpFlashWithBothFrozenEosFamilies() {
    for (boolean pengRobinson : new boolean[] { false, true }) {
      SystemInterface point = fluidTwo(pengRobinson);
      point.setTemperature(273.15 + 25.0);
      point.setPressure(10.0);
      new ThermodynamicOperations(point).TPflash();
      point.initProperties();
      assertTrue(point.getNumberOfPhases() >= 1 && point.getNumberOfPhases() <= 3);
      assertTrue(Double.isFinite(point.getDensity("kg/m3")) && point.getDensity("kg/m3") > 0.0);
      assertEquals("Huron-Vidal", point.getMixingRuleName());
    }
  }

  @Test
  @Tag("slow")
  void fluidOneProducesAContinuousConventionalBoundaryBeforeThreePhaseExtension() {
    for (boolean pengRobinson : new boolean[] { false, true }) {
      TwoHydrocarbonPhaseEnvelopeSolver.Result result = new TwoHydrocarbonPhaseEnvelopeSolver(fluidOne(pengRobinson))
          .setPressureRange(0.1, 1000.0).setMaximumSteps(5.0, 5.0).setMaximumContinuationIterations(1200)
          .setTopologyFallbackEnabled(false).setStabilityAnalysisEnabled(true).setThreePhaseRefinement(16, 0.02, 0.02)
          .solve();
      long failedStabilitySamples = result.getSecondaryStabilitySamples().stream()
          .filter(sample -> sample.getFailureMessage() != null).count();
      long validPrimarySamples = result.getSecondaryStabilitySamples().stream()
          .filter(PTPhaseEnvelopeMichelsen.SecondaryStabilitySample::isPrimaryBoundaryValid).count();
      long stableAqueousSamples = result.getSecondaryStabilitySamples().stream()
          .filter(PTPhaseEnvelopeMichelsen.SecondaryStabilitySample::isPrimaryBoundaryValid)
          .filter(sample -> Double.isFinite(sample.getStabilityFunction()) && sample.getStabilityFunction() >= 0.0)
          .count();
      double minimumAqueousStability = result.getSecondaryStabilitySamples().stream()
          .mapToDouble(PTPhaseEnvelopeMichelsen.SecondaryStabilitySample::getStabilityFunction).filter(Double::isFinite)
          .min().orElse(Double.NaN);
      double maximumAqueousStability = result.getSecondaryStabilitySamples().stream()
          .mapToDouble(PTPhaseEnvelopeMichelsen.SecondaryStabilitySample::getStabilityFunction).filter(Double::isFinite)
          .max().orElse(Double.NaN);
      System.out.printf(
          "provisional Fluid 1 %s envelope: method=%s points=%d segments=%d closed=%s critical=%d " + "failure=%s%n",
          pengRobinson ? "PR78-HV" : "SRK-Peneloux-HV", result.getMethod(), result.getPointCount(),
          result.getSegments().size(), result.isEnvelopeClosed(), result.getCriticalPointCount(),
          result.getFailureMessage());
      System.out.printf("  secondary stability: samples=%d unstable=%d failures=%d brackets=%d%n",
          result.getSecondaryStabilitySamples().size(), result.getThreePhaseStabilityPointCount(),
          failedStabilitySamples, result.getSecondaryStabilityBrackets().size());
      System.out.printf("  primary-valid=%d aqueous-stable=%d aqueous-S-range=[%.8g, %.8g]%n", validPrimarySamples,
          stableAqueousSamples, minimumAqueousStability, maximumAqueousStability);
      for (PTPhaseEnvelopeMichelsen.SecondaryStabilityBracket bracket : result.getSecondaryStabilityBrackets()) {
        PTPhaseEnvelopeMichelsen.SecondaryStabilitySample twoPhase = bracket.getTwoPhaseSide();
        PTPhaseEnvelopeMichelsen.SecondaryStabilitySample threePhase = bracket.getThreePhaseSide();
        System.out.printf("    bracket %s: 2P=(%.3f K, %.3f bara) 3P=(%.3f K, %.3f bara, %s)%n", twoPhase.getBranch(),
            twoPhase.getTemperature(), twoPhase.getPressure(), threePhase.getTemperature(), threePhase.getPressure(),
            threePhase.getIncipientPhase());
      }
      for (PTPhaseEnvelopeMichelsen.ThreePhasePointCandidate candidate : result.getThreePhasePointCandidates()) {
        System.out.printf(
            "    refined: T=%.5f K P=%.5f bara widths=(%.5g K, %.5g bara) iterations=%d "
                + "phase=%s converged=%s failure=%s%n",
            candidate.getTemperature(), candidate.getPressure(), candidate.getTemperatureWidth(),
            candidate.getPressureWidth(), candidate.getIterations(), candidate.getIncipientPhase(),
            candidate.isConverged(), candidate.getFailureMessage());
      }
      assertTrue(result.hasFinitePhysicalValues());
      assertTrue(result.getPointCount() >= 6,
          "Fluid 1 must produce a traceable conventional gas-oil boundary before water-topology switching");
      assertTrue(result.getSecondaryStabilitySamples().size() >= 6,
          "A closed conventional boundary must expose enough points for secondary-phase bracketing");
      assertEquals(result.getSecondaryStabilityBrackets().size(), result.getThreePhasePointCandidates().size());
      assertTrue(validPrimarySamples >= 6,
          "Secondary TPD diagnostics are only meaningful where the conventional hydrocarbon boundary remains valid");
    }
  }

  static SystemInterface fluidOne(boolean pengRobinson) {
    SystemInterface fluid = baseSystem(pengRobinson);
    fluid.addComponent("nitrogen", 0.34);
    fluid.addComponent("CO2", 0.84);
    fluid.addComponent("methane", 89.95);
    fluid.addComponent("ethane", 5.17);
    fluid.addComponent("propane", 2.04);
    fluid.addComponent("i-butane", 0.36);
    fluid.addComponent("n-butane", 0.55);
    fluid.addComponent("i-pentane", 0.14);
    fluid.addComponent("n-pentane", 0.10);
    fluid.addComponent("n-hexane", 0.01);
    fluid.addComponent("water", 500.0e-4);
    return finish(fluid);
  }

  static SystemInterface fluidTwo(boolean pengRobinson) {
    SystemInterface fluid = baseSystem(pengRobinson);
    fluid.addComponent("CO2", 2.79);
    fluid.addComponent("methane", 71.51);
    fluid.addComponent("ethane", 5.77);
    fluid.addComponent("propane", 4.10);
    fluid.addComponent("i-butane", 1.32);
    fluid.addComponent("n-butane", 1.60);
    fluid.addComponent("i-pentane", 0.82);
    fluid.addComponent("n-pentane", 0.64);
    fluid.addComponent("n-hexane", 1.05);
    fluid.addTBPfraction("C7", 10.40, 0.191, 0.82);
    fluid.addComponent("water", 8.0);
    return finish(fluid);
  }

  private static SystemInterface baseSystem(boolean pengRobinson) {
    if (pengRobinson) {
      return new SystemPrEos1978(298.15, 10.0);
    }
    return new SystemSrkPenelouxEos(298.15, 10.0);
  }

  private static SystemInterface finish(SystemInterface fluid) {
    fluid.createDatabase(true);
    fluid.setMixingRule("HV", "NRTL");
    fluid.setMultiPhaseCheck(true);
    fluid.init(0);
    return fluid;
  }
}
