"""Typed, validated process definitions via pydantic.

The Java ``ProcessSystem.fromJsonAndRun`` accepts a declarative JSON flowsheet.
Hand-writing that JSON as raw dicts gives no autocomplete and no validation until
Java fails. These pydantic models mirror the JSON schema so you get editor
autocomplete, type checking, and clear validation errors *before* the JVM runs.

pydantic is an optional dependency. Install it with::

    pip install "neqsim[schema]"    # or: pip install pydantic

Example
-------
    >>> from neqsim.process.schema import ProcessModel, Fluid, Unit
    >>> model = ProcessModel(
    ...     fluid=Fluid(eos="srk", components={"methane": 0.9, "ethane": 0.1}),
    ...     process=[
    ...         Unit(type="Stream", name="feed",
    ...              properties={"flowRate": [50000.0, "kg/hr"],
    ...                          "temperature": [30.0, "C"],
    ...                          "pressure": [50.0, "bara"]}),
    ...         Unit(type="Separator", name="HP Sep", inlet="feed"),
    ...         Unit(type="Compressor", name="Comp", inlet="HP Sep.gasOut",
    ...              properties={"outletPressure": [100.0, "bara"]}),
    ...     ],
    ... )
    >>> result = model.run()          # doctest: +SKIP
    >>> process = result.getProcessSystem()   # doctest: +SKIP
"""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

try:
    from pydantic import BaseModel, ConfigDict, field_validator
except Exception as exc:  # pragma: no cover - exercised only without pydantic
    raise ImportError(
        "neqsim.process.schema requires pydantic. "
        'Install it with: pip install "neqsim[schema]"'
    ) from exc

__all__ = ["Fluid", "Unit", "ProcessModel"]


class Fluid(BaseModel):
    """A fluid definition for the JSON process builder.

    Attributes:
        eos: Equation of state (e.g. ``"srk"``, ``"pr"``, ``"cpa"``).
        mixingRule: Mixing rule name (e.g. ``"classic"``).
        components: Mapping of component name to mole fraction (or flow).
    """

    model_config = ConfigDict(extra="allow")

    eos: str = "srk"
    mixingRule: str = "classic"
    components: Dict[str, float] = {}


class Unit(BaseModel):
    """A single process unit (equipment or stream) definition.

    Attributes:
        type: Equipment type, e.g. ``"Stream"``, ``"Separator"``,
            ``"Compressor"`` (matched case-insensitively by the Java builder).
        name: Unique unit name used for stream wiring by reference.
        inlet: Name (or ``"unit.port"``) of a single inlet stream, if any.
        inlets: Names of multiple inlet streams, if any (e.g. for mixers).
        fluid: Name of a named fluid to use for a stream (optional).
        properties: Free-form property map. Values are commonly
            ``[value, "unit"]`` pairs (e.g. ``{"flowRate": [50000.0, "kg/hr"]}``).
    """

    model_config = ConfigDict(extra="allow")

    type: str
    name: str
    inlet: Optional[str] = None
    inlets: Optional[List[str]] = None
    fluid: Optional[str] = None
    properties: Optional[Dict[str, Any]] = None

    @field_validator("type", "name")
    @classmethod
    def _non_empty(cls, value: str) -> str:
        """Ensure ``type`` and ``name`` are non-empty.

        Args:
            value: The field value.

        Returns:
            The validated value.

        Raises:
            ValueError: If the value is blank.
        """
        if not value or not value.strip():
            raise ValueError("must be a non-empty string")
        return value


class ProcessModel(BaseModel):
    """A complete declarative flowsheet.

    Attributes:
        fluid: The default fluid (used by streams without an explicit fluid).
        fluids: Named fluids map, for multi-fluid flowsheets.
        process: The ordered list of units.
    """

    model_config = ConfigDict(extra="allow")

    fluid: Optional[Fluid] = None
    fluids: Optional[Dict[str, Fluid]] = None
    process: List[Unit] = []

    @field_validator("process")
    @classmethod
    def _unique_names(cls, units: List[Unit]) -> List[Unit]:
        """Ensure all unit names are unique.

        Args:
            units: The list of units.

        Returns:
            The validated list.

        Raises:
            ValueError: If a duplicate name is found.
        """
        seen: set = set()
        for unit in units:
            if unit.name in seen:
                raise ValueError(f"duplicate unit name: {unit.name!r}")
            seen.add(unit.name)
        return units

    def to_json(self) -> str:
        """Serialize the model to the JSON string the Java builder expects.

        Returns:
            A JSON string.
        """
        data = self.model_dump(exclude_none=True)
        return json.dumps(data)

    def run(self, fluid: Any = None) -> Any:
        """Build and run the process via ``ProcessSystem.fromJsonAndRun``.

        Args:
            fluid: Optional pre-built Java fluid to use as the default fluid
                (e.g. one imported from an Eclipse E300 file). When provided, the
                model's ``fluid`` section is ignored by Java.

        Returns:
            The Java ``SimulationResult``. Use ``.isError()`` /
            ``.getProcessSystem()`` to inspect it.
        """
        from neqsim import jneqsim

        ProcessSystem = jneqsim.process.processmodel.ProcessSystem
        if fluid is not None:
            return ProcessSystem.fromJsonAndRun(self.to_json(), fluid)
        return ProcessSystem.fromJsonAndRun(self.to_json())
