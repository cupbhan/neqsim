"""Component-name helpers to reduce typos when building fluids.

NeqSim ships a database of ~2000 component names. ``addComponent("methan")``
fails at the Java layer with an unhelpful error. These helpers let you validate,
search, and get suggestions for component names from Python before touching Java.

Examples
--------
    >>> from neqsim.thermo.components import list_components, find_components
    >>> "methane" in list_components()            # doctest: +SKIP
    True
    >>> find_components("glycol")                  # doctest: +SKIP
    ['MEG', 'DEG', 'TEG', ...]
    >>> from neqsim.thermo.components import suggest_component
    >>> suggest_component("methan")                # doctest: +SKIP
    ['methane', 'methanol', ...]
"""

from __future__ import annotations

import difflib
from functools import lru_cache
from typing import List

import jpype

__all__ = [
    "list_components",
    "is_valid_component",
    "find_components",
    "suggest_component",
]

_RESOURCE = "/neqsim_component_names.txt"


@lru_cache(maxsize=1)
def list_components() -> List[str]:
    """Return every component name known to NeqSim.

    The list is read once from the ``neqsim_component_names.txt`` resource inside
    the NeqSim jar and cached.

    Returns:
        A sorted list of component names.
    """
    names: List[str] = []
    try:
        cls = jpype.JClass("neqsim.util.database.NeqSimDataBase").class_
        stream = cls.getResourceAsStream(_RESOURCE)
        if stream is None:
            loader = cls.getClassLoader()
            if loader is not None:
                stream = loader.getResourceAsStream(_RESOURCE.lstrip("/"))
        if stream is not None:
            reader = jpype.JClass("java.io.BufferedReader")(
                jpype.JClass("java.io.InputStreamReader")(stream, "UTF-8")
            )
            line = reader.readLine()
            while line is not None:
                text = str(line).strip()
                if text and text.upper() != "NAME":
                    names.append(text)
                line = reader.readLine()
            reader.close()
    except Exception:
        names = []
    return sorted(set(names))


def is_valid_component(name: str) -> bool:
    """Return whether a name is a known NeqSim component (case-insensitive).

    Args:
        name: Candidate component name.

    Returns:
        True if the name matches a known component.
    """
    lowered = name.lower()
    return any(lowered == c.lower() for c in list_components())


def find_components(keyword: str) -> List[str]:
    """Return components whose name contains a keyword (case-insensitive).

    Args:
        keyword: Substring to search for.

    Returns:
        A sorted list of matching component names.
    """
    needle = keyword.lower()
    return [c for c in list_components() if needle in c.lower()]


def suggest_component(name: str, n: int = 5) -> List[str]:
    """Return the closest component names to a (possibly misspelled) input.

    Args:
        name: The candidate/misspelled component name.
        n: Maximum number of suggestions to return.

    Returns:
        A list of up to ``n`` close matches, best first.
    """
    return difflib.get_close_matches(name, list_components(), n=n, cutoff=0.6)
