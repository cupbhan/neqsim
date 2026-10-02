"""Tests for the pydantic-based process schema (neqsim.process.schema)."""

import pytest

pydantic = pytest.importorskip("pydantic")

from neqsim.process.schema import Fluid, ProcessModel, Unit  # noqa: E402


def _model():
    return ProcessModel(
        fluid=Fluid(eos="srk", components={"methane": 0.9, "ethane": 0.1}),
        process=[
            Unit(
                type="Stream",
                name="feed",
                properties={
                    "flowRate": [50000.0, "kg/hr"],
                    "temperature": [30.0, "C"],
                    "pressure": [50.0, "bara"],
                },
            ),
            Unit(type="Separator", name="HP Sep", inlet="feed"),
        ],
    )


def test_to_json_contains_fluid_and_process():
    data = _model().to_json()
    assert '"fluid"' in data
    assert '"process"' in data
    assert '"HP Sep"' in data


def test_duplicate_names_rejected():
    with pytest.raises(Exception):
        ProcessModel(
            process=[
                Unit(type="Stream", name="dup"),
                Unit(type="Separator", name="dup"),
            ]
        )


def test_empty_type_rejected():
    with pytest.raises(Exception):
        Unit(type="", name="x")


def test_run_builds_process():
    result = _model().run()
    assert not result.isError()
    names = list(result.getProcessSystem().getAllUnitNames())
    assert "feed" in names
    assert "HP Sep" in names
