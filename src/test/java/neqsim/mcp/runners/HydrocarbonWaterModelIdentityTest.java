package neqsim.mcp.runners;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import org.junit.jupiter.api.Test;
import com.google.gson.JsonObject;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

class HydrocarbonWaterModelIdentityTest extends neqsim.NeqSimTest {
  @Test
  void liveIdentityIsDeterministicAndSeparatesCompositionFromParameters() {
    SystemInterface first = fluid(0.8, 0.2);
    SystemInterface second = fluid(0.8, 0.2);
    JsonObject firstIdentity = HydrocarbonWaterModelIdentity.describe(first, "SRK", "Classic", null);
    JsonObject secondIdentity = HydrocarbonWaterModelIdentity.describe(second, "SRK", "Classic", null);

    assertEquals("resolved", firstIdentity.get("status").getAsString());
    assertEquals(64, firstIdentity.get("compositionFingerprint").getAsString().length());
    assertEquals(firstIdentity.get("fullParameterFingerprint"), secondIdentity.get("fullParameterFingerprint"));
    assertEquals(firstIdentity.get("compositionFingerprint"), secondIdentity.get("compositionFingerprint"));

    SystemInterface changedComposition = fluid(0.7, 0.3);
    JsonObject compositionIdentity = HydrocarbonWaterModelIdentity.describe(changedComposition, "SRK", "Classic", null);
    assertNotEquals(firstIdentity.get("compositionFingerprint"), compositionIdentity.get("compositionFingerprint"));
    assertEquals(firstIdentity.get("fullParameterFingerprint"), compositionIdentity.get("fullParameterFingerprint"));

    second.setBinaryInteractionParameter("methane", "water", 0.33);
    JsonObject changedParameters = HydrocarbonWaterModelIdentity.describe(second, "SRK", "Classic", null);
    assertNotEquals(firstIdentity.get("fullParameterFingerprint"), changedParameters.get("fullParameterFingerprint"));
  }

  private static SystemInterface fluid(double methane, double water) {
    SystemInterface system = new SystemSrkEos(350.0, 10.0);
    system.addComponent("methane", methane);
    system.addComponent("water", water);
    system.setMixingRule(2);
    system.init(0);
    return system;
  }
}
