package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Small cubic-EOS fixtures cross-checked against the ordinary TP flash path.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
final class ResearchEquilibriumFixtures {
  /** Prevents construction of the fixture utility. */
  private ResearchEquilibriumFixtures() {
  }

  /**
   * Creates a two-phase hydrocarbon or three-phase hydrocarbon-water equilibrium.
   *
   * @param wet whether to include a water-rich phase
   * @return independently flashed reference state
   */
  static SystemInterface equilibrium(boolean wet) {
    SystemInterface fluid = new SystemSrkEos(350.0, 10.0);
    fluid.addComponent("methane", wet ? 0.20 : 0.60);
    fluid.addComponent("nC10", wet ? 0.15 : 0.40);
    if (wet) {
      fluid.addComponent("water", 0.65);
    }
    fluid.createDatabase(true);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    fluid.setMaxNumberOfPhases(3);
    new ThermodynamicOperations(fluid).TPflash();
    fluid.init(1);
    assertEquals(wet ? 3 : 2, fluid.getNumberOfPhases(), "reference phase count");
    return fluid;
  }

  /**
   * Extracts a defensive array of phase compositions.
   *
   * @param fluid equilibrated system
   * @return phase-by-component mole fractions
   */
  static double[][] compositions(SystemInterface fluid) {
    double[][] values = new double[fluid.getNumberOfPhases()][fluid.getNumberOfComponents()];
    for (int phase = 0; phase < values.length; phase++) {
      for (int component = 0; component < values[phase].length; component++) {
        values[phase][component] = fluid.getPhase(phase).getComponent(component).getx();
      }
    }
    return values;
  }

  /**
   * Extracts reference phase amounts.
   *
   * @param fluid equilibrated system
   * @return phase mole fractions
   */
  static double[] fractions(SystemInterface fluid) {
    double[] values = new double[fluid.getNumberOfPhases()];
    for (int phase = 0; phase < values.length; phase++) {
      values[phase] = fluid.getBeta(phase);
    }
    return values;
  }

  /**
   * Reads the actual reference phase types.
   *
   * @param fluid equilibrated system
   * @return matching physical phase families
   */
  static CandidatePhase[] slots(SystemInterface fluid) {
    CandidatePhase[] values = new CandidatePhase[fluid.getNumberOfPhases()];
    for (int phase = 0; phase < values.length; phase++) {
      values[phase] = CandidatePhase.valueOf(fluid.getPhase(phase).getType().name());
    }
    return values;
  }
}
