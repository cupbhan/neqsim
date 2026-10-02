"""Convert NeqSim process results to pandas DataFrames.

A single :func:`stream_table` (aliased :func:`to_dataframe`) covers dozens of
equipment types, because every piece of equipment exposes its connected streams
through ``getInletStreams()`` / ``getOutletStreams()``. Instead of writing a
result wrapper per equipment class, this walks the whole flowsheet and builds a
tidy stream table (one row per stream: flow, temperature, pressure, phase count,
density, molar mass).

Examples
--------
    >>> from neqsim.process import stream_table, runProcess
    >>> runProcess()                       # doctest: +SKIP
    >>> df = stream_table()                # global process  # doctest: +SKIP
    >>> df = stream_table(my_process)      # explicit process  # doctest: +SKIP

    >>> from neqsim.process import equipment_table
    >>> equipment_table(my_process)        # doctest: +SKIP
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import jpype
import pandas as pd

__all__ = ["stream_table", "to_dataframe", "equipment_table"]


def _identity_hash(obj: Any) -> int:
    """Return a stable identity key for a Java object.

    Args:
        obj: Any object (typically a Java stream).

    Returns:
        ``System.identityHashCode`` for Java objects, or Python ``id`` as a
        fallback.
    """
    try:
        return int(jpype.JClass("java.lang.System").identityHashCode(obj))
    except Exception:
        return id(obj)


def _resolve_process(process: Any = None) -> Any:
    """Resolve the target ``ProcessSystem`` from various input types.

    Args:
        process: A ``ProcessSystem``, a ``ProcessContext``/``ProcessBuilder``
            wrapper exposing ``process``/``getProcess``, or None to use the
            global process.

    Returns:
        The underlying Java ``ProcessSystem`` object.
    """
    if process is None:
        from neqsim.process.processTools import getProcess

        return getProcess()
    for attr in ("process", "_process"):
        inner = getattr(process, attr, None)
        if inner is not None:
            return inner
    getter = getattr(process, "getProcess", None)
    if callable(getter):
        return getter()
    return process


def _safe(func, *args) -> Optional[float]:
    """Call a Java getter and return a float, or None on any failure.

    Args:
        func: A callable (Java method handle) or None.
        *args: Arguments to pass to ``func``.

    Returns:
        The float result, or None if the call fails or ``func`` is None.
    """
    if func is None:
        return None
    try:
        value = func(*args)
        return float(value)
    except Exception:
        return None


def _collect_streams(system: Any) -> List[Any]:
    """Collect all unique streams referenced by a flowsheet.

    Walks every unit's inlet and outlet streams, de-duplicating by Java object
    identity so a stream shared between two units appears once.

    Args:
        system: A Java ``ProcessSystem``.

    Returns:
        A list of unique stream objects, in first-seen order.
    """
    seen: Dict[int, Any] = {}
    try:
        units = list(system.getUnitOperations())
    except Exception:
        units = []

    for unit in units:
        for getter_name in ("getInletStreams", "getOutletStreams"):
            getter = getattr(unit, getter_name, None)
            if getter is None:
                continue
            try:
                for stream in getter():
                    if stream is not None:
                        seen.setdefault(_identity_hash(stream), stream)
            except Exception:
                continue
    return list(seen.values())


def _stream_row(stream: Any) -> Dict[str, Any]:
    """Build one stream-table row from a stream object.

    Args:
        stream: A Java stream object.

    Returns:
        A dict of column name -> value for the stream.
    """
    try:
        name = str(stream.getName())
    except Exception:
        name = None

    getflow = getattr(stream, "getFlowRate", None)
    gettemp = getattr(stream, "getTemperature", None)
    getpres = getattr(stream, "getPressure", None)

    row: Dict[str, Any] = {
        "name": name,
        "flow_kg_hr": _safe(getflow, "kg/hr"),
        "flow_m3_hr": _safe(getflow, "m3/hr"),
        "temperature_C": _safe(gettemp, "C"),
        "pressure_bara": _safe(getpres, "bara"),
    }

    try:
        fluid = stream.getFluid()
    except Exception:
        fluid = None

    if fluid is not None:
        row["n_phases"] = _safe(getattr(fluid, "getNumberOfPhases", None))
        row["molar_mass_kg_mol"] = _safe(getattr(fluid, "getMolarMass", None))
        row["density_kg_m3"] = _safe(getattr(fluid, "getDensity", None), "kg/m3")
    return row


def stream_table(process: Any = None) -> pd.DataFrame:
    """Return a stream table for a process as a pandas DataFrame.

    One row per unique stream in the flowsheet, with flow, temperature,
    pressure, phase count, molar mass, and density. Missing/uncomputed values
    appear as ``NaN``.

    Args:
        process: A ``ProcessSystem``, ``ProcessContext``, ``ProcessBuilder``, or
            None to use the global process.

    Returns:
        A ``pandas.DataFrame`` of stream properties.

    Example:
        >>> from neqsim.process import stream_table, runProcess
        >>> runProcess()                    # doctest: +SKIP
        >>> stream_table()                  # doctest: +SKIP
    """
    system = _resolve_process(process)
    rows = [_stream_row(s) for s in _collect_streams(system)]
    return pd.DataFrame(rows)


def equipment_table(process: Any = None) -> pd.DataFrame:
    """Return a summary of equipment in a process as a pandas DataFrame.

    One row per unit operation, with its name, Java class (type), and the number
    of inlet and outlet streams.

    Args:
        process: A ``ProcessSystem``, ``ProcessContext``, ``ProcessBuilder``, or
            None to use the global process.

    Returns:
        A ``pandas.DataFrame`` with columns ``name``, ``type``, ``n_inlets``,
        ``n_outlets``.
    """
    system = _resolve_process(process)
    try:
        units = list(system.getUnitOperations())
    except Exception:
        units = []

    rows: List[Dict[str, Any]] = []
    for unit in units:
        name = None
        try:
            name = str(unit.getName())
        except Exception:
            pass
        type_name = None
        try:
            type_name = str(unit.getClass().getSimpleName())
        except Exception:
            pass
        n_in = _count_streams(getattr(unit, "getInletStreams", None))
        n_out = _count_streams(getattr(unit, "getOutletStreams", None))
        rows.append(
            {"name": name, "type": type_name, "n_inlets": n_in, "n_outlets": n_out}
        )
    return pd.DataFrame(rows)


def _count_streams(getter) -> Optional[int]:
    """Return the number of streams returned by a getter, or None.

    Args:
        getter: A callable returning a list of streams, or None.

    Returns:
        The stream count, or None if the getter is missing or fails.
    """
    if getter is None:
        return None
    try:
        return int(len(list(getter())))
    except Exception:
        return None


# Primary, discoverable name for the stream-table helper.
to_dataframe = stream_table
