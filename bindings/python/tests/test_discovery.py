"""Tests for the neqsim.discovery API-introspection helpers."""

from neqsim import discovery


def test_list_classes_finds_many():
    classes = discovery.list_classes("neqsim")
    assert len(classes) > 100
    assert all(name.startswith("neqsim.") for name in classes)


def test_list_equipment_includes_known_classes():
    equipment = discovery.list_equipment()
    assert any(name.endswith(".Compressor") for name in equipment)
    assert any(name.endswith(".Separator") for name in equipment)


def test_list_classes_non_recursive_excludes_subpackages():
    direct = discovery.list_classes("neqsim.process.equipment", recursive=False)
    assert all(
        name.count(".") == len("neqsim.process.equipment".split(".")) for name in direct
    )
    # Compressor lives in a sub-package, so it must NOT be in the direct listing.
    assert not any(name.endswith(".compressor.Compressor") for name in direct)


def test_list_packages_returns_subpackages():
    packages = discovery.list_packages("process.equipment")
    assert "neqsim.process.equipment.compressor" in packages
    assert "neqsim.process.equipment.separator" in packages


def test_find_classes_keyword_search():
    matches = discovery.find_classes("scrubber")
    assert any(name.endswith(".GasScrubber") for name in matches)


def test_get_class_simple_and_full_name():
    from_simple = discovery.get_class("Compressor")
    from_full = discovery.get_class("neqsim.process.equipment.compressor.Compressor")
    assert str(from_simple) == str(from_full)


def test_get_class_unknown_raises():
    import pytest

    with pytest.raises(ValueError):
        discovery.get_class("ThisClassDoesNotExist12345")


def test_describe_lists_constructors_and_methods():
    text = discovery.describe("Compressor")
    assert "Constructors" in text
    assert "Public methods" in text
    assert "getPower" in text


def test_manifest_offline_describe_when_present():
    # When the packaged manifest exists, describe works without live reflection.
    manifest = discovery._load_manifest()
    if manifest is None:
        import pytest

        pytest.skip("no offline manifest packaged")
    assert "neqsim.process.equipment.compressor.Compressor" in manifest["classes"]
    text = discovery.describe("neqsim.process.equipment.compressor.Compressor")
    assert "getPower" in text
