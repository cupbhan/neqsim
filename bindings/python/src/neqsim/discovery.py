"""Discoverability helpers for the full NeqSim Java API.

NeqSim exposes hundreds of Java classes (equipment, mechanical design, safety,
field development, automation, thermodynamics). Only a small, curated subset has
hand-written Python wrappers. Every Java class is, however, always reachable
through direct access::

    from neqsim import jneqsim
    sep = jneqsim.process.equipment.separator.ThreePhaseSeparator("sep", stream)

The problem with the direct path is *discoverability* - there is no autocomplete
for a ``jpype.JPackage`` and no easy way to see which classes exist under, e.g.,
``neqsim.process.equipment``. This module solves that by scanning the JVM
classpath (the packaged ``lib/*.jar`` in a normal install, or ``target/classes``
in a development checkout) and letting you list, search, and describe the whole
API from Python.

Examples
--------
List every equipment class::

    from neqsim.discovery import list_equipment
    for name in list_equipment():
        print(name)

Search the whole API for anything matching a keyword::

    from neqsim.discovery import find_classes
    find_classes("compressor")

Inspect the constructors and methods of a class without leaving Python::

    from neqsim.discovery import describe
    print(describe("neqsim.process.equipment.compressor.Compressor"))

Get a JClass handle by (simple or fully-qualified) name::

    from neqsim.discovery import get_class
    Compressor = get_class("Compressor")
"""

from __future__ import annotations

import gzip
import json
import os
import zipfile
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

import jpype

from neqsim import jneqsim  # noqa: F401  (ensures the JVM is started)

__all__ = [
    "list_classes",
    "list_packages",
    "list_equipment",
    "find_classes",
    "describe",
    "get_class",
    "refresh",
]

_ROOT_PACKAGE = "neqsim"
_DATA_DIR = Path(__file__).resolve().parent / "data"
_MANIFEST_GZ = _DATA_DIR / "neqsim_api.json.gz"
_MANIFEST_JSON = _DATA_DIR / "neqsim_api.json"


@lru_cache(maxsize=1)
def _load_manifest() -> Optional[Dict]:
    """Load the packaged offline API manifest, if present.

    The manifest (``neqsim/data/neqsim_api.json.gz``, or an uncompressed
    ``neqsim_api.json``) is produced at release time by
    ``scripts/generate_api_manifest.py``. When present it makes class listing,
    search, and :func:`describe` work offline and instantly, without scanning a
    jar or reflecting at runtime. When absent, discovery falls back to live
    reflection.

    Returns:
        The parsed manifest dict, or None if no manifest is available or it
        cannot be read.
    """
    try:
        if _MANIFEST_GZ.is_file():
            with gzip.open(_MANIFEST_GZ, "rt", encoding="utf-8") as handle:
                data = json.load(handle)
        elif _MANIFEST_JSON.is_file():
            with open(_MANIFEST_JSON, "r", encoding="utf-8") as handle:
                data = json.load(handle)
        else:
            return None
        if isinstance(data, dict) and isinstance(data.get("classes"), dict):
            return data
    except Exception:
        pass
    return None


def _code_source_entries() -> List[str]:
    """Return the jar/dir a loaded NeqSim class was actually loaded from.

    This is the most reliable source: NeqSim's jar is added to a JPype dynamic
    classloader *after* JVM start, so it does not appear on
    ``java.class.path``. The protection domain's code source, however, always
    points at the real jar (packaged install) or class directory (dev checkout).

    Returns:
        A list with the code-source path(s), or an empty list on failure.
    """
    entries: List[str] = []
    for probe in (
        "neqsim.thermo.system.SystemSrkEos",
        "neqsim.process.processmodel.ProcessSystem",
    ):
        try:
            jclass = jpype.JClass(probe)
            location = jclass.class_.getProtectionDomain().getCodeSource().getLocation()
            path = str(jpype.JClass("java.io.File")(location.toURI()).getPath())
            if path:
                entries.append(path)
        except Exception:
            continue
    return entries


def _classpath_entries() -> List[str]:
    """Return every classpath location worth scanning for NeqSim classes.

    Combines (in priority order) the code source of a loaded NeqSim class, the
    packaged ``lib/*.jar`` folder, and the JVM ``java.class.path`` entries. This
    is correct for both a packaged install (``lib/*.jar``) and a development
    checkout (``target/classes``).

    Returns:
        A de-duplicated list of filesystem paths (jars and/or directories).
    """
    entries: List[str] = list(_code_source_entries())

    lib_dir = Path(__file__).resolve().parent / "lib"
    if lib_dir.is_dir():
        entries.extend(str(p) for p in lib_dir.glob("*.jar"))

    if jpype.isJVMStarted():
        try:
            System = jpype.JClass("java.lang.System")
            raw = str(System.getProperty("java.class.path") or "")
            entries.extend(p for p in raw.split(os.pathsep) if p)
        except Exception:
            pass

    seen: set = set()
    unique: List[str] = []
    for entry in entries:
        if entry and entry not in seen:
            seen.add(entry)
            unique.append(entry)
    return unique


def _iter_class_names_from_jar(jar_path: str) -> List[str]:
    """Enumerate fully-qualified class names inside a jar under ``neqsim``.

    Args:
        jar_path: Path to a ``.jar`` file.

    Returns:
        A list of dotted, fully-qualified class names (no nested/anonymous
        classes).
    """
    names: List[str] = []
    try:
        with zipfile.ZipFile(jar_path) as zf:
            for entry in zf.namelist():
                if not entry.endswith(".class"):
                    continue
                if "$" in entry:  # skip nested / anonymous classes
                    continue
                dotted = entry[:-6].replace("/", ".")
                if dotted.startswith(_ROOT_PACKAGE + "."):
                    names.append(dotted)
    except (zipfile.BadZipFile, FileNotFoundError, OSError):
        pass
    return names


def _iter_class_names_from_dir(dir_path: str) -> List[str]:
    """Enumerate fully-qualified class names under a directory (``neqsim/...``).

    Args:
        dir_path: Path to a classpath directory (e.g. ``target/classes``).

    Returns:
        A list of dotted, fully-qualified class names (no nested/anonymous
        classes).
    """
    names: List[str] = []
    root = Path(dir_path)
    base = root / _ROOT_PACKAGE
    if not base.is_dir():
        return names
    for path in base.rglob("*.class"):
        if "$" in path.name:  # skip nested / anonymous classes
            continue
        rel = path.relative_to(root).with_suffix("")
        dotted = ".".join(rel.parts)
        names.append(dotted)
    return names


@lru_cache(maxsize=1)
def _all_class_names() -> List[str]:
    """Return the sorted, de-duplicated list of all ``neqsim`` class names.

    Prefers the packaged offline manifest; falls back to scanning the classpath
    (jar or ``target/classes``). The result is cached. Call :func:`refresh` to
    rebuild it (e.g. after adding a jar to the classpath at runtime).

    Returns:
        A sorted list of fully-qualified class names.
    """
    manifest = _load_manifest()
    if manifest is not None:
        return sorted(manifest["classes"].keys())

    seen: set = set()
    for entry in _classpath_entries():
        if entry.lower().endswith(".jar"):
            names = _iter_class_names_from_jar(entry)
        elif os.path.isdir(entry):
            names = _iter_class_names_from_dir(entry)
        else:
            names = []
        seen.update(names)
    return sorted(seen)


def refresh() -> int:
    """Clear the cached class list so it is rebuilt on next use.

    Useful after adding a jar to the classpath at runtime, or after (re)placing
    the offline API manifest.

    Returns:
        The number of classes discovered after refreshing.
    """
    _all_class_names.cache_clear()
    _load_manifest.cache_clear()
    return len(_all_class_names())


def list_classes(
    package_prefix: str = _ROOT_PACKAGE, recursive: bool = True
) -> List[str]:
    """List Java class names under a package.

    Args:
        package_prefix: Package to list, e.g. ``"neqsim.process.equipment"``.
            The leading ``"neqsim."`` may be omitted (``"process.equipment"``).
        recursive: If True (default), include classes in sub-packages. If
            False, only classes directly in ``package_prefix``.

    Returns:
        A sorted list of fully-qualified class names.

    Example:
        >>> from neqsim.discovery import list_classes
        >>> list_classes("process.equipment.compressor")  # doctest: +SKIP
        ['neqsim.process.equipment.compressor.Compressor', ...]
    """
    prefix = _normalize_package(package_prefix)
    result: List[str] = []
    for name in _all_class_names():
        if not name.startswith(prefix + "."):
            continue
        if not recursive:
            remainder = name[len(prefix) + 1 :]
            if "." in remainder:  # class lives in a sub-package
                continue
        result.append(name)
    return result


def list_packages(package_prefix: str = _ROOT_PACKAGE) -> List[str]:
    """List the immediate sub-packages of a package.

    Args:
        package_prefix: Package to inspect, e.g. ``"neqsim.process"``. The
            leading ``"neqsim."`` may be omitted.

    Returns:
        A sorted list of fully-qualified sub-package names.

    Example:
        >>> from neqsim.discovery import list_packages
        >>> list_packages("process.equipment")  # doctest: +SKIP
        ['neqsim.process.equipment.compressor', ...]
    """
    prefix = _normalize_package(package_prefix)
    subpackages: set = set()
    for name in _all_class_names():
        if not name.startswith(prefix + "."):
            continue
        remainder = name[len(prefix) + 1 :]
        if "." in remainder:
            subpackages.add(prefix + "." + remainder.split(".", 1)[0])
    return sorted(subpackages)


def list_equipment(recursive: bool = True) -> List[str]:
    """List every process-equipment class.

    Convenience shortcut for ``list_classes("neqsim.process.equipment")``.

    Args:
        recursive: If True (default), include sub-packages (compressor,
            separator, valve, ...).

    Returns:
        A sorted list of fully-qualified equipment class names.
    """
    return list_classes("neqsim.process.equipment", recursive=recursive)


def find_classes(keyword: str, package_prefix: str = _ROOT_PACKAGE) -> List[str]:
    """Search the API for classes whose name contains a keyword.

    Args:
        keyword: Case-insensitive substring to search for in the class name.
        package_prefix: Restrict the search to this package. Defaults to the
            whole ``neqsim`` tree.

    Returns:
        A sorted list of matching fully-qualified class names.

    Example:
        >>> from neqsim.discovery import find_classes
        >>> find_classes("scrubber")  # doctest: +SKIP
        ['neqsim.process.equipment.separator.GasScrubber', ...]
    """
    prefix = _normalize_package(package_prefix)
    needle = keyword.lower()
    return [
        name
        for name in _all_class_names()
        if name.startswith(prefix + ".") and needle in name.rsplit(".", 1)[-1].lower()
    ]


def get_class(name: str):
    """Return a ``jpype`` class handle from a simple or fully-qualified name.

    Args:
        name: Either a fully-qualified name
            (``"neqsim.process.equipment.compressor.Compressor"``) or a simple
            class name (``"Compressor"``). If a simple name is ambiguous, the
            first match is returned and the alternatives are noted in the raised
            error only when there is no match.

    Returns:
        The ``jpype`` class object, ready to instantiate.

    Raises:
        ValueError: If no class matches ``name``.

    Example:
        >>> from neqsim.discovery import get_class
        >>> Compressor = get_class("Compressor")  # doctest: +SKIP
    """
    if "." in name:
        return jpype.JClass(name)
    matches = [n for n in _all_class_names() if n.rsplit(".", 1)[-1] == name]
    if not matches:
        raise ValueError(
            f"No NeqSim class named '{name}' found. "
            f"Use find_classes('{name}') to search for similar names."
        )
    return jpype.JClass(matches[0])


def describe(cls_or_name) -> str:
    """Describe a class's constructors and public methods.

    Uses the offline API manifest when available (fast, JVM-free); otherwise
    falls back to live Java reflection.

    Args:
        cls_or_name: A fully-qualified/simple class name, or a ``jpype`` class
            handle.

    Returns:
        A human-readable, multi-line description of the constructors and public
        methods.

    Example:
        >>> from neqsim.discovery import describe
        >>> print(describe("Compressor"))  # doctest: +SKIP
    """
    manifest = _load_manifest()
    if manifest is not None and isinstance(cls_or_name, str):
        fqcn = cls_or_name
        if "." not in fqcn:
            candidates = [
                n for n in manifest["classes"] if n.rsplit(".", 1)[-1] == fqcn
            ]
            fqcn = candidates[0] if candidates else fqcn
        entry = manifest["classes"].get(fqcn)
        if entry is not None:
            return _describe_from_manifest(fqcn, entry)

    if isinstance(cls_or_name, str):
        jclass = get_class(cls_or_name)
    else:
        jclass = cls_or_name
    java_class = jclass.class_

    lines: List[str] = [f"class {java_class.getName()}", ""]

    constructors = sorted(
        (str(c) for c in java_class.getConstructors()),
    )
    lines.append(f"Constructors ({len(constructors)}):")
    for ctor in constructors:
        lines.append(f"  {_simplify_signature(ctor)}")
    lines.append("")

    methods = sorted(
        {_simplify_signature(str(m)) for m in java_class.getMethods()},
    )
    lines.append(f"Public methods ({len(methods)}):")
    for method in methods:
        lines.append(f"  {method}")

    return "\n".join(lines)


def _describe_from_manifest(fqcn: str, entry: Dict) -> str:
    """Render a :func:`describe` string from a manifest class entry.

    Args:
        fqcn: The fully-qualified class name.
        entry: The manifest entry with ``constructors`` and ``methods`` lists.

    Returns:
        A human-readable, multi-line description.
    """
    constructors = list(entry.get("constructors", []))
    methods = list(entry.get("methods", []))
    lines: List[str] = [f"class {fqcn}", ""]
    lines.append(f"Constructors ({len(constructors)}):")
    for ctor in constructors:
        lines.append(f"  {ctor}")
    lines.append("")
    lines.append(f"Public methods ({len(methods)}):")
    for method in sorted(methods):
        lines.append(f"  {method}")
    return "\n".join(lines)


def _normalize_package(package_prefix: str) -> str:
    """Prepend ``neqsim.`` to a package prefix when it is omitted.

    Args:
        package_prefix: A package name with or without the leading
            ``"neqsim."``.

    Returns:
        The fully-qualified package prefix.
    """
    prefix = package_prefix.strip(". ")
    if prefix == _ROOT_PACKAGE or prefix.startswith(_ROOT_PACKAGE + "."):
        return prefix
    return _ROOT_PACKAGE + "." + prefix


def _simplify_signature(signature: str) -> str:
    """Strip package qualifiers from a reflected member signature for reading.

    Args:
        signature: The raw ``toString()`` of a Java ``Constructor`` or
            ``Method``.

    Returns:
        A shorter signature with fully-qualified type names reduced to their
        simple names.
    """
    out: List[str] = []
    token = ""
    for ch in signature:
        if ch.isalnum() or ch in "_$.":
            token += ch
        else:
            out.append(_shorten_token(token))
            token = ""
            out.append(ch)
    out.append(_shorten_token(token))
    return "".join(out)


def _shorten_token(token: str) -> str:
    """Reduce a dotted type token to its simple name (keeps ``neqsim`` tokens).

    Args:
        token: A candidate type token (possibly dotted, possibly empty).

    Returns:
        The simple type name if the token looks like a qualified type,
        otherwise the token unchanged.
    """
    if "." in token and token[0].islower():
        return token.rsplit(".", 1)[-1]
    return token
