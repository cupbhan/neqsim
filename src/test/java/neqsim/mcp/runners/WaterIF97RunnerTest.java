package neqsim.mcp.runners;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

/** Tests for {@link WaterIF97Runner}. */
class WaterIF97RunnerTest {
  @Test
  void calculatesStateAndMonotonicSaturationCurve() {
    JsonObject result = JsonParser.parseString(WaterIF97Runner.run(473.15, 15.0, 293.15, 623.15, 31)).getAsJsonObject();

    assertEquals("success", result.get("status").getAsString());
    assertEquals("IAPWS-IF97", result.get("model").getAsString());
    assertTrue(result.getAsJsonObject("flash").get("numberOfPhases").getAsInt() >= 1);
    assertTrue(result.getAsJsonObject("fluid").has("properties"));
    JsonArray curve = result.getAsJsonArray("envelope");
    assertEquals(31, curve.size());
    for (int index = 1; index < curve.size(); index++) {
      assertTrue(curve.get(index).getAsJsonObject().get("pressure_bara").getAsDouble() > curve.get(index - 1)
          .getAsJsonObject().get("pressure_bara").getAsDouble());
    }
  }
}
