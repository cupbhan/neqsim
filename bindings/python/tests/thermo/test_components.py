"""Tests for neqsim.thermo.components helpers."""

from neqsim.thermo.components import (
    find_components,
    is_valid_component,
    list_components,
    suggest_component,
)


def test_list_components_nonempty():
    comps = list_components()
    assert len(comps) > 50
    assert "methane" in comps


def test_is_valid_component_case_insensitive():
    assert is_valid_component("methane")
    assert is_valid_component("METHANE")
    assert not is_valid_component("definitely_not_a_component")


def test_find_components_keyword():
    matches = find_components("methan")
    assert "methane" in matches
    assert "methanol" in matches


def test_suggest_component_for_typo():
    suggestions = suggest_component("methan")
    assert "methane" in suggestions
