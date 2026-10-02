package neqsim.thermodynamicoperations.flashops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import org.junit.jupiter.api.Test;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermo.util.readwrite.JsonFluidReadWrite;

/** Tests explicit diagnostics from the multiphase phase-fraction solver. */
class TPmultiflashSolveStatusTest {

  /** Verifies that the personal status gate accepts an upstream recovery only at a balanced equilibrium. */
  @Test
  void acceptsValidatedColdWaterEndpointRecovery() {
    SystemInterface system = new SystemSrkEos(232.0, 92.0);
    String[] names = {"CO2", "methane", "ethane", "nC10", "water"};
    double[] composition = {0.74, 0.15, 0.05, 0.01, 0.05};
    for (int index = 0; index < names.length; index++) {
      system.addComponent(names[index], composition[index]);
    }
    system.setMixingRule(2);
    system.setMultiPhaseCheck(true);
    TPflash operation = new TPflash(system);
    operation.run();
    system.init(1);

    assertTrue(operation.isLastMultiphaseSolveAccepted(), operation.getLastMultiphaseSolveMessage());
    assertEquals(2, system.getNumberOfPhases());
    assertTrue(system.hasPhaseType(PhaseType.OIL));
    assertTrue(system.hasPhaseType(PhaseType.AQUEOUS));
    assertTrue(massBalanceResidual(system) < 1.0e-10);
    assertTrue(fugacityResidual(system) < 1.0e-8);
  }

  @Test
  void propagatesValidatedMultiphaseStatusThroughPublicTpFlash() throws Exception {
    SystemInterface system = buildFieldSystem();
    TPflash operation = new TPflash(system);

    operation.run();

    assertEquals(3, system.getNumberOfPhases());
    assertTrue(operation.getLastMultiphaseSolveStatus().isConverged(), operation.getLastMultiphaseSolveMessage());
    assertTrue(operation.isLastMultiphaseSolveAccepted());
    assertTrue(operation.getLastMultiphaseSolveAttemptCount() > 0);
    assertTrue(operation.getLastMultiphaseMassBalanceResidual() <= 1.0e-9);
  }

  @Test
  void solvesPhaseFractionsAfterLateAqueousPhaseSeeding() throws Exception {
    SystemInterface system = buildFieldSystem(315.0, 280.0);
    TPflash operation = new TPflash(system);

    operation.run();

    assertTrue(operation.isLastMultiphaseSolveAccepted(), operation.getLastMultiphaseSolveMessage());
    assertTrue(operation.getLastMultiphaseSolveAttemptCount() > 1, operation.getLastMultiphaseSolveMessage());
    assertTrue(massBalanceResidual(system) <= 1.0e-9);
  }

  @Test
  void repairsWaterRichSplitBeforeStabilityAnalysis() throws Exception {
    double[][] states = {{215.0, 280.0}, {220.0, 300.0}, {225.0, 300.0}, {230.0, 300.0}, {235.0, 300.0}, {245.0, 300.0},
        {250.0, 300.0}};

    for (double[] state : states) {
      SystemInterface system = buildFieldSystem(state[0], state[1]);
      TPflash operation = new TPflash(system);
      operation.run();

      String diagnostic = state[1] + " bara / " + state[0] + " C: " + operation.getLastMultiphaseSolveMessage();
      assertTrue(operation.isLastMultiphaseSolveAccepted(), diagnostic);
      assertTrue(operation.getLastMultiphaseSolveAttemptCount() > 0, diagnostic);
      assertTrue(system.hasPhaseType(PhaseType.AQUEOUS), diagnostic);
      assertTrue(massBalanceResidual(system) <= 1.0e-9, diagnostic);
      assertTrue(fugacityResidual(system) <= 1.0e-7, diagnostic);
    }
  }

  @Test
  void recoversLateInvalidBetaVectorForNeutralOilWaterSplit() throws Exception {
    SystemInterface system = buildH2O65RegressionSystem();
    TPflash operation = new TPflash(system);

    operation.run();

    String diagnostic = operation.getLastMultiphaseSolveMessage();
    assertEquals(2, system.getNumberOfPhases(), diagnostic);
    assertTrue(system.hasPhaseType(PhaseType.AQUEOUS), diagnostic);
    assertTrue(system.hasPhaseType(PhaseType.OIL), diagnostic);
    assertTrue(operation.isLastMultiphaseSolveAccepted(), diagnostic);
    assertTrue(massBalanceResidual(system) <= 1.0e-9, diagnostic);
    assertTrue(fugacityResidual(system) <= 1.0e-7, diagnostic);
    assertEquals(1.0, system.getBeta(0) + system.getBeta(1), 1.0e-12, diagnostic);
  }

  @Test
  void preservesWaterSubcriticalPressureStabilityPath() throws Exception {
    for (double temperatureC : new double[] {320.0, 325.0}) {
      SystemInterface system = buildFieldSystem(temperatureC, 200.0);
      TPflash operation = new TPflash(system);
      operation.run();

      String diagnostic = "200 bara / " + temperatureC + " C: " + operation.getLastMultiphaseSolveMessage();
      assertTrue(operation.isLastMultiphaseSolveAccepted(), diagnostic);
      assertTrue(massBalanceResidual(system) <= 1.0e-9, diagnostic);
      assertTrue(fugacityResidual(system) <= 1.0e-7, diagnostic);
    }
  }

  @Test
  void repairsPreviouslySpeculativeBetaStall() throws Exception {
    SystemInterface system = buildFieldSystem();
    system.setMultiPhaseCheck(false);
    system.setMaxNumberOfPhases(2);
    new TPflash(system).run();
    system.setMultiPhaseCheck(true);
    system.setMaxNumberOfPhases(3);

    TPmultiflash operation = new TPmultiflash(system, false);
    operation.run();

    String diagnostic = operation.getSolveStatus() + ": " + operation.getSolveStatusMessage();
    assertTrue(operation.getLastSolveBetaStatus().isConverged(), diagnostic);
    assertTrue(operation.getLastSolveBetaIterations() > 0, diagnostic);
    assertTrue(operation.getSolveStatus().isConverged(), diagnostic);
    assertTrue(!operation.isPhaseCleanupSkipped(), diagnostic);
    assertTrue(operation.getSolveBetaAttemptCount() > 0, diagnostic);
    assertTrue(massBalanceResidual(system) <= 1.0e-9, diagnostic);
    assertTrue(fugacityResidual(system) <= 1.0e-7, diagnostic);
  }

  @Test
  void reportsConvergedStatusForValidatedThreePhaseContinuation() throws Exception {
    SystemInterface system = buildFieldSystem();
    new TPflash(system).run();
    assertEquals(3, system.getNumberOfPhases());

    TPmultiflash operation = new TPmultiflash(system, false);
    setBooleanField(operation, "doStabilityAnalysis", false);
    setBooleanField(operation, "multiPhaseTest", true);
    operation.run();

    assertEquals(3, system.getNumberOfPhases());
    assertTrue(system.hasPhaseType(PhaseType.GAS));
    assertTrue(system.hasPhaseType(PhaseType.AQUEOUS));
    assertTrue(system.hasPhaseType(PhaseType.OIL));
    assertTrue(massBalanceResidual(system) <= 1.0e-9);
    assertTrue(operation.getSolveStatus().isConverged(), operation.getSolveStatusMessage());
    assertTrue(operation.getSolveBetaAttemptCount() > 0);
    assertTrue(operation.getFinalMassBalanceResidual() <= 1.0e-9);
  }

  @Test
  void reportsSingularStatusWhenNewtonMatrixCannotBeSolved() {
    SystemInterface system = new SystemSrkEos(298.15, 10.0);
    system.addComponent("methane", 0.8);
    system.addComponent("ethane", 0.2);
    system.setMixingRule("classic");
    system.setNumberOfPhases(2);
    system.setBeta(0, 0.5);
    system.setBeta(1, 0.5);
    system.init(0);
    system.init(1);

    TPmultiflash operation = new SingularMatrixTPmultiflash(system);
    operation.setDoubleArrays();
    operation.solveBeta();

    assertEquals(TPmultiflash.SolveStatus.SINGULAR_OR_ILL_CONDITIONED, operation.getLastSolveBetaStatus());
    assertTrue(!operation.getLastSolveBetaStatus().isConverged());
  }

  private static final class SingularMatrixTPmultiflash extends TPmultiflash {
    private static final long serialVersionUID = 1L;

    private SingularMatrixTPmultiflash(SystemInterface system) {
      super(system, false);
    }

    @Override
    public double calcQ() {
      for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
        dQdbeta[phase][0] = 1.0;
        for (int otherPhase = 0; otherPhase < system.getNumberOfPhases(); otherPhase++) {
          Qmatrix[phase][otherPhase] = 0.0;
        }
      }
      return 0.0;
    }
  }

  private static SystemInterface buildFieldSystem() throws Exception {
    return buildFieldSystem(215.0, 240.0);
  }

  private static SystemInterface buildFieldSystem(double temperatureC, double pressureBara) throws Exception {
    JsonObject payload;
    InputStream input = TPmultiflashSolveStatusTest.class
        .getResourceAsStream("/neqsim/mcp/runners/field_wet_classic_payload.json");
    assertTrue(input != null);
    try (InputStreamReader reader = new InputStreamReader(input, StandardCharsets.UTF_8)) {
      payload = JsonParser.parseReader(reader).getAsJsonObject();
    }

    Class<?> runner = Class.forName("neqsim.mcp.runners.FieldFluidRunner");
    Method build = runner.getDeclaredMethod("buildSystem", JsonObject.class, JsonObject.class, String.class,
        boolean.class, double.class, double.class);
    build.setAccessible(true);
    return (SystemInterface) build.invoke(null, payload.getAsJsonObject("components"),
        payload.getAsJsonObject("fluidDefinition"), "SRK", false, temperatureC + 273.15, pressureBara);
  }

  private static SystemInterface buildH2O65RegressionSystem() throws Exception {
    InputStream input = TPmultiflashSolveStatusTest.class
        .getResourceAsStream("/neqsim/thermodynamicoperations/flashops/h2o65_beta_recovery.json");
    assertTrue(input != null);
    StringBuilder json = new StringBuilder();
    try (InputStreamReader reader = new InputStreamReader(input, StandardCharsets.UTF_8)) {
      char[] buffer = new char[4096];
      int count;
      while ((count = reader.read(buffer)) >= 0) {
        json.append(buffer, 0, count);
      }
    }
    SystemInterface system = JsonFluidReadWrite.readString(json.toString());
    system.setTemperature(250.0, "C");
    system.setPressure(66.0, "bara");
    system.setMultiPhaseCheck(true);
    return system;
  }

  private static void setBooleanField(TPmultiflash operation, String name, boolean value) throws Exception {
    Field field = TPmultiflash.class.getDeclaredField(name);
    field.setAccessible(true);
    field.setBoolean(operation, value);
  }

  private static double massBalanceResidual(SystemInterface system) {
    double maximumResidual = 0.0;
    for (int component = 0; component < system.getPhase(0).getNumberOfComponents(); component++) {
      double reconstructed = 0.0;
      for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
        reconstructed += system.getBeta(phase) * system.getPhase(phase).getComponent(component).getx();
      }
      maximumResidual = Math.max(maximumResidual,
          Math.abs(reconstructed - system.getPhase(0).getComponent(component).getz()));
    }
    return maximumResidual;
  }

  private static double fugacityResidual(SystemInterface system) {
    if (system.getNumberOfPhases() < 2) {
      return 0.0;
    }
    double maximumResidual = 0.0;
    for (int component = 0; component < system.getPhase(0).getNumberOfComponents(); component++) {
      if (system.getPhase(0).getComponent(component).getz() < 1.0e-12) {
        continue;
      }
      double minimum = Double.POSITIVE_INFINITY;
      double maximum = Double.NEGATIVE_INFINITY;
      for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
        double moleFraction = Math.max(system.getPhase(phase).getComponent(component).getx(), 1.0e-100);
        double logFugacity = Math.log(moleFraction)
            + system.getPhase(phase).getComponent(component).getLogFugacityCoefficient();
        minimum = Math.min(minimum, logFugacity);
        maximum = Math.max(maximum, logFugacity);
      }
      maximumResidual = Math.max(maximumResidual, maximum - minimum);
    }
    return maximumResidual;
  }
}
