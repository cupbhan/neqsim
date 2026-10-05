package neqsim.thermo.mixingrule;

import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;
import neqsim.thermo.ThermodynamicConstantsInterface;
import neqsim.thermo.component.ComponentEosInterface;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;

/**
 * Verifies fixed-state Huron-Vidal calculations independently of initialization and component order.
 *
 * @author NeqSim contributors
 * @version 1.0
 */
class HuronVidalInitializationConsistencyTest {
  private static final String[] NAMES = {"methane", "water", "n-heptane"};
  private static final double[] AMOUNTS = {0.65, 0.10, 0.25};
  private static final int[][] ORDERS = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}};

  /** Verifies the direct modified-NRTL excess energy at all supported derivative levels. */
  @Test
  void mixedPairsMatchIndependentExcessEnergyAtEveryInitializationLevel() {
    for (boolean pr : new boolean[] {false, true}) {
      for (int[] order : ORDERS) {
        SystemInterface fluid = createFluid(pr, order, true);
        for (int level : new int[] {1, 2, 3, 1}) {
          fluid.init(level);
          PhaseEos phase = (PhaseEos) fluid.getPhase(0);
          double expected = independentMolarAttraction(phase, order);
          double actual = phase.getA() / Math.pow(phase.getNumberOfMolesInPhase(), 2.0);
          assertEquals(expected, actual, Math.abs(expected) * 1.0e-11,
              "Independent HV attraction for PR=" + pr + ", init=" + level);
        }
      }
    }
  }

  /** Verifies fugacity coefficients and Z across initialization levels and component permutations. */
  @Test
  void fixedStateFugacitiesDoNotDependOnInitializationOrComponentOrder() {
    for (boolean pr : new boolean[] {false, true}) {
      double[] reference = null;
      for (int[] order : ORDERS) {
        SystemInterface fluid = createFluid(pr, order, true);
        for (int level : new int[] {1, 2, 3, 1}) {
          fluid.init(level);
          PhaseEos phase = (PhaseEos) fluid.getPhase(0);
          double[] values = new double[NAMES.length + 1];
          values[0] = phase.getZ();
          for (int i = 0; i < order.length; i++) {
            values[order[i] + 1] = phase.getComponent(i).getLogFugacityCoefficient();
          }
          if (reference == null) {
            reference = values;
          } else {
            for (int i = 0; i < values.length; i++) {
              assertEquals(reference[i], values[i], 1.0e-11,
                  "Fixed-state property " + i + " for PR=" + pr + ", init=" + level);
            }
          }
        }
      }
    }
  }

  /** Verifies that an HV rule with exclusively classical pairs reduces to the cubic mixing rule. */
  @Test
  void classicalLimitMatchesQuadraticMixing() {
    for (boolean pr : new boolean[] {false, true}) {
      for (int[] order : ORDERS) {
        SystemInterface fluid = createFluid(pr, order, false);
        for (int level : new int[] {1, 2, 3}) {
          fluid.init(level);
          PhaseEos phase = (PhaseEos) fluid.getPhase(0);
          double expected = 0.0;
          for (int i = 0; i < order.length; i++) {
            for (int j = 0; j < order.length; j++) {
              expected += phase.getComponent(i).getx() * phase.getComponent(j).getx()
                  * Math.sqrt(((ComponentEosInterface) phase.getComponent(i)).aT(phase.getTemperature())
                      * ((ComponentEosInterface) phase.getComponent(j)).aT(phase.getTemperature()))
                  * (1.0 - binaryInteraction(order[i], order[j]));
            }
          }
          assertEquals(expected, phase.getA() / Math.pow(phase.getNumberOfMolesInPhase(), 2.0),
              Math.abs(expected) * 1.0e-11);
        }
      }
    }
  }

  /**
   * Creates an explicit synthetic model with no dependence on database interaction parameters.
   *
   * @param pr whether to use PR instead of SRK
   * @param order permutation of canonical component indices
   * @param useHv whether the methane-water pair uses HV
   * @return initialized fixed-composition fluid at 320 K and 30 bara
   */
  private SystemInterface createFluid(boolean pr, int[] order, boolean useHv) {
    SystemInterface fluid = pr ? new SystemPrEos(320.0, 30.0) : new SystemSrkEos(320.0, 30.0);
    for (int index : order) {
      fluid.addComponent(NAMES[index], AMOUNTS[index]);
    }
    fluid.setMixingRule(4);
    for (int phaseIndex = 0; phaseIndex < fluid.getNumberOfPhases(); phaseIndex++) {
      PhaseEos phase = (PhaseEos) fluid.getPhase(phaseIndex);
      HVMixingRulesInterface rule = (HVMixingRulesInterface) phase.getEosMixingRule();
      for (int i = 0; i < order.length; i++) {
        for (int j = 0; j < order.length; j++) {
          boolean hv = useHv && isHvPair(order[i], order[j]);
          rule.setClassicOrHV(i, j, hv ? "HV" : "Classic");
          rule.setBinaryInteractionParameter(i, j, binaryInteraction(order[i], order[j]));
          rule.setHValphaParameter(i, j, hv ? 0.35 : 0.0);
          rule.setHVDijParameter(i, j, hv ? (order[i] == 0 ? 1200.0 : -1000.0) : 0.0);
          rule.setHVDijTParameter(i, j, hv ? (order[i] == 0 ? 1.0 : -0.5) : 0.0);
        }
      }
    }
    fluid.init(0);
    return fluid;
  }

  /**
   * Identifies the only nonclassical pair in the synthetic model.
   *
   * @param first canonical component index
   * @param second canonical component index
   * @return whether the indices form the methane-water pair
   */
  private boolean isHvPair(int first, int second) {
    return (first == 0 && second == 1) || (first == 1 && second == 0);
  }

  /**
   * Supplies the synthetic symmetric classical binary interaction coefficient.
   *
   * @param first canonical component index
   * @param second canonical component index
   * @return dimensionless interaction coefficient
   */
  private double binaryInteraction(int first, int second) {
    return first == second ? 0.0 : (first == 1 || second == 1 ? 0.30 : 0.02);
  }

  /**
   * Computes HV attraction from the scalar modified-NRTL excess-energy expression, without activity coefficients.
   *
   * @param phase fixed-state cubic phase
   * @param order permutation of canonical component indices
   * @return molar attraction in the phase's internal EOS units
   */
  private double independentMolarAttraction(PhaseEos phase, int[] order) {
    double temperature = phase.getTemperature();
    double rt = ThermodynamicConstantsInterface.R * temperature;
    double[] delta = ((ComponentEosInterface) phase.getComponent(0)).getDeltaEosParameters();
    double lambda = Math.log((1.0 + delta[1]) / (1.0 + delta[0])) / (delta[1] - delta[0]);
    double mixtureB = 0.0;
    double pureSum = 0.0;
    double excess = 0.0;
    for (int i = 0; i < order.length; i++) {
      ComponentEosInterface ci = (ComponentEosInterface) phase.getComponent(i);
      mixtureB += ci.getx() * ci.getb();
      pureSum += ci.getx() * ci.aT(temperature) / ci.getb();
      double numerator = 0.0;
      double denominator = 0.0;
      for (int j = 0; j < order.length; j++) {
        ComponentEosInterface cj = (ComponentEosInterface) phase.getComponent(j);
        double tau;
        double alpha = 0.0;
        if (isHvPair(order[j], order[i])) {
          tau = (order[j] == 0 ? 1200.0 : -1000.0) / temperature + (order[j] == 0 ? 1.0 : -0.5);
          alpha = 0.35;
        } else {
          double gii = -lambda * ci.aT(temperature) / ci.getb();
          double gjj = -lambda * cj.aT(temperature) / cj.getb();
          double gji = -2.0 * Math.sqrt(ci.getb() * cj.getb() * gii * gjj) / (ci.getb() + cj.getb())
              * (1.0 - binaryInteraction(order[j], order[i]));
          tau = (gji - gii) / rt;
        }
        double weight = cj.getx() * cj.getb() * Math.exp(-alpha * tau);
        numerator += weight * tau;
        denominator += weight;
      }
      excess += ci.getx() * numerator / denominator;
    }
    return mixtureB * (pureSum - rt * excess / lambda);
  }
}
