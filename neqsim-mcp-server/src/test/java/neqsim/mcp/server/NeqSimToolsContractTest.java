package neqsim.mcp.server;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import io.quarkiverse.mcp.server.Tool;
import java.lang.reflect.Method;
import org.junit.jupiter.api.Test;

class NeqSimToolsContractTest {
  private final NeqSimTools tools = new NeqSimTools();

  @Test
  void runsTpFlashThroughThePublishedToolContract() {
    JsonObject result = JsonParser
        .parseString(tools.runFlash("{\"methane\":0.9,\"ethane\":0.1}", 30.0, "C", 50.0, "bara", "PR", "TP"))
        .getAsJsonObject();

    assertFalse(result.has("error"), result.toString());
    assertTrue(result.has("flash"), result.toString());
    assertTrue(result.has("fluid"), result.toString());
  }

  @Test
  void exposesWaterIf97WithFiniteStateAndSaturationCurve() {
    JsonObject result = JsonParser.parseString(tools.runWaterIF97(473.15, 15.0, 373.15, 473.15, 5)).getAsJsonObject();

    assertEquals("success", result.get("status").getAsString(), result.toString());
    assertTrue(result.has("fluid"), result.toString());
    assertEquals(11, result.getAsJsonArray("saturationCurve").size(), result.toString());
  }

  @Test
  void rejectsMalformedFlashCompositionAtTheServiceBoundary() {
    JsonObject result = JsonParser.parseString(tools.runFlash("[]", 30.0, "C", 50.0, "bara", "PR", "TP"))
        .getAsJsonObject();

    assertTrue(result.has("error") || "error".equals(result.has("status") ? result.get("status").getAsString() : ""),
        result.toString());
  }

  @Test
  void exposedToolParametersAvoidNondeterministicPrimitiveInvokerGeneration() {
    for (Method method : NeqSimTools.class.getDeclaredMethods()) {
      if (!method.isAnnotationPresent(Tool.class)) {
        continue;
      }
      for (Class<?> parameterType : method.getParameterTypes()) {
        assertFalse(parameterType.isPrimitive(),
            () -> method.getName() + " exposes primitive parameter " + parameterType.getName());
      }
    }
  }
}
