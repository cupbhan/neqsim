package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import java.util.Arrays;
import org.junit.jupiter.api.Test;
import neqsim.NeqSimTest;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.BoundaryDefinition;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Phase;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.Region;

/** Tests for the physical gas-oil-aqueous topology contract. */
class HydrocarbonWaterPhaseTopologyTest extends NeqSimTest {

  @Test
  void classifiesEverySupportedStablePhaseRegion() {
    assertEquals(Region.GAS, Region.fromPhaseTypeNames(Arrays.asList("gas")));
    assertEquals(Region.OIL, Region.fromPhaseTypeNames(Arrays.asList("oil")));
    assertEquals(Region.AQUEOUS, Region.fromPhaseTypeNames(Arrays.asList("aqueous")));
    assertEquals(Region.GAS_OIL, Region.fromPhaseTypeNames(Arrays.asList("oil", "gas")));
    assertEquals(Region.GAS_AQUEOUS, Region.fromPhaseTypeNames(Arrays.asList("aqueous", "gas")));
    assertEquals(Region.OIL_AQUEOUS, Region.fromPhaseTypeNames(Arrays.asList("water", "liquid")));
    assertEquals(Region.GAS_OIL_AQUEOUS, Region.fromPhaseTypeNames(Arrays.asList("gas", "oil", "aqueous")));
  }

  @Test
  void permitsOnlyOnePhaseTopologyTransitions() {
    assertTrue(Region.GAS.isAdjacentTo(Region.GAS_OIL));
    assertTrue(Region.GAS_OIL.isAdjacentTo(Region.GAS_OIL_AQUEOUS));
    assertFalse(Region.GAS.isAdjacentTo(Region.OIL));
    assertFalse(Region.GAS.isAdjacentTo(Region.GAS_OIL_AQUEOUS));

    BoundaryDefinition waterFromGasOil = new BoundaryDefinition(Region.GAS_OIL, Region.GAS_OIL_AQUEOUS);
    assertEquals(Phase.AQUEOUS, waterFromGasOil.getIncipientPhase());
    assertEquals("GO->GOW:AQUEOUS", waterFromGasOil.getCode());

    assertThrows(IllegalArgumentException.class, () -> new BoundaryDefinition(Region.GAS, Region.GAS_OIL_AQUEOUS));
    assertThrows(IllegalArgumentException.class, () -> new BoundaryDefinition(Region.GAS_OIL, Region.GAS));
  }
}
