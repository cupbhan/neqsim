package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import Jama.Matrix;

/**
 * Checks continuation history interpolation independently of expensive EOS evaluations.
 *
 * @author NeqSim developers
 * @version 1.0
 */
class EnvelopeHistoryPredictionTest {
  /** Last accepted specification observed in the multicomponent UMR-PRU regression. */
  private static final double SPECIFICATION_ORIGIN = 49.62464976852;

  /** Extrapolation distance measured in history intervals. */
  private static final double NEXT_POSITION = 0.5;

  /** Shifting the specification origin must not change the predicted thermodynamic variables. */
  @Test
  void predictionDoesNotDependOnSpecificationOffset() {
    double[] reference = predict(0.0, 1.0e-4);
    double[] shifted = predict(SPECIFICATION_ORIGIN, 1.0e-4);
    for (int variable = 1; variable < reference.length; variable++) {
      assertEquals(reference[variable], shifted[variable], 1.0e-8,
          "An arbitrary specification offset must not change continuation variable " + variable);
    }
  }

  /** Both forward and reverse affine specification scales must reproduce the same cubic history. */
  @Test
  void predictionReproducesCubicUnderAffineSpecificationChanges() {
    for (double spacing : new double[] {1.0e-4, 1.0e-2, -3.0e-4}) {
      double[] predicted = predict(SPECIFICATION_ORIGIN, spacing);
      assertEquals(SPECIFICATION_ORIGIN + spacing * NEXT_POSITION, predicted[0], 1.0e-10,
          "The specification itself must advance by the requested signed step");
      for (int variable = 1; variable < predicted.length; variable++) {
        assertEquals(historyValue(variable, NEXT_POSITION), predicted[variable], 1.0e-8,
            "A cubic predictor must reproduce cubic variable " + variable + " with spacing " + spacing);
      }
      double temperature = Math.exp(predicted[2]);
      double pressure = Math.exp(predicted[3]);
      assertTrue(Double.isFinite(temperature) && temperature > 100.0 && temperature < 150.0,
          "The synthetic logarithmic temperature must remain on its local physical scale");
      assertTrue(Double.isFinite(pressure) && pressure > 3.0 && pressure < 4.0,
          "The synthetic logarithmic pressure must remain on its local physical scale");
    }
  }

  /**
   * Predicts an exactly specified cubic history without invoking an EOS or weakening any corrector gate.
   *
   * @param origin arbitrary specification origin
   * @param spacing signed distance between consecutive history specification values
   * @return predicted log-state vector
   */
  private static double[] predict(double origin, double spacing) {
    SysNewtonRhapsonPhaseEnvelope solver = new SysNewtonRhapsonPhaseEnvelope();
    solver.numberOfComponents = 2;
    solver.neq = 4;
    solver.speceq = 0;
    solver.Xgij = new Matrix(4, 4);
    solver.u = new Matrix(4, 1);
    solver.dxds = new Matrix(4, 1);
    solver.dxds.set(0, 0, 1.0);
    solver.ds = spacing * NEXT_POSITION;
    for (int point = 0; point < 4; point++) {
      double position = point - 3.0;
      double specification = origin + spacing * position;
      solver.s.set(0, point, specification);
      solver.Xgij.set(0, point, specification);
      for (int variable = 1; variable < 4; variable++) {
        solver.Xgij.set(variable, point, historyValue(variable, position));
      }
    }
    solver.u = solver.Xgij.getMatrix(0, 3, 3, 3);
    solver.calcInc2(5);
    return solver.u.getColumnPackedCopy();
  }

  /**
   * Evaluates synthetic cubic state histories on realistic log-temperature and log-pressure scales.
   *
   * @param variable state variable index: one for log K, two for log T, three for log P
   * @param position history coordinate measured from the last accepted point
   * @return exact synthetic logarithmic state variable
   */
  private static double historyValue(int variable, double position) {
    if (variable == 1) {
      return -0.8 + position * (0.04 + position * (0.01 + position * 0.003));
    }
    if (variable == 2) {
      return Math.log(127.118982) + position * (0.001 + position * (0.0002 + position * 0.00003));
    }
    return Math.log(3.259684) + position * (0.002 + position * (0.0004 + position * 0.00005));
  }
}
