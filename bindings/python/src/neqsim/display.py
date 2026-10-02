"""Rich Jupyter display for NeqSim streams and processes.

Importing this module attaches ``_repr_html_`` to the relevant Java classes via
JPype class customizers, so a stream or a ``ProcessSystem`` renders as a
formatted HTML table in a notebook instead of a Java ``toString()``.

It is imported automatically by ``neqsim`` (best-effort). If the JVM is not
started, or JPype customizers are unavailable, importing is a silent no-op so it
never breaks a normal ``import neqsim``.

    >>> from neqsim.process import stream, runProcess   # doctest: +SKIP
    >>> s = stream("feed", fluid)                        # doctest: +SKIP
    >>> s        # in Jupyter, shows an HTML table        # doctest: +SKIP
"""

from __future__ import annotations

from typing import Any, Optional

__all__ = ["register"]

_registered = False


def _fmt(value: Optional[float], digits: int = 3) -> str:
    """Format a float for display, or return an empty string for None/NaN.

    Args:
        value: The value to format.
        digits: Number of decimal places.

    Returns:
        A formatted string, or ``""`` if the value is missing.
    """
    try:
        if value is None:
            return ""
        f = float(value)
        if f != f:  # NaN
            return ""
        return f"{f:.{digits}f}"
    except Exception:
        return ""


def _safe(func, *args) -> Optional[float]:
    """Call a Java getter and return a float, or None on failure.

    Args:
        func: A callable or None.
        *args: Arguments to pass.

    Returns:
        The float result, or None.
    """
    if func is None:
        return None
    try:
        return float(func(*args))
    except Exception:
        return None


def _stream_html(self: Any) -> str:
    """Return an HTML table describing a stream.

    Args:
        self: The Java stream object (bound as a method).

    Returns:
        An HTML string.
    """
    try:
        name = str(self.getName())
    except Exception:
        name = "stream"
    rows = [
        ("Flow [kg/hr]", _fmt(_safe(getattr(self, "getFlowRate", None), "kg/hr"))),
        ("Temperature [C]", _fmt(_safe(getattr(self, "getTemperature", None), "C"))),
        ("Pressure [bara]", _fmt(_safe(getattr(self, "getPressure", None), "bara"))),
    ]
    try:
        fluid = self.getFluid()
    except Exception:
        fluid = None
    if fluid is not None:
        rows.append(
            ("Phases", _fmt(_safe(getattr(fluid, "getNumberOfPhases", None)), 0))
        )
        rows.append(
            (
                "Density [kg/m3]",
                _fmt(_safe(getattr(fluid, "getDensity", None), "kg/m3")),
            )
        )
    body = "".join(
        f"<tr><td style='text-align:left'>{k}</td>"
        f"<td style='text-align:right'>{v}</td></tr>"
        for k, v in rows
    )
    return (
        f"<table><caption style='text-align:left'><b>Stream:</b> {name}</caption>"
        f"{body}</table>"
    )


def _process_html(self: Any) -> str:
    """Return an HTML table summarising a process's equipment.

    Args:
        self: The Java ``ProcessSystem`` object (bound as a method).

    Returns:
        An HTML string.
    """
    try:
        units = list(self.getUnitOperations())
    except Exception:
        units = []
    header = (
        "<tr><th style='text-align:left'>Name</th>"
        "<th style='text-align:left'>Type</th></tr>"
    )
    rows = []
    for unit in units:
        try:
            name = str(unit.getName())
        except Exception:
            name = ""
        try:
            type_name = str(unit.getClass().getSimpleName())
        except Exception:
            type_name = ""
        rows.append(
            f"<tr><td style='text-align:left'>{name}</td>"
            f"<td style='text-align:left'>{type_name}</td></tr>"
        )
    try:
        title = str(self.getName())
    except Exception:
        title = "process"
    return (
        f"<table><caption style='text-align:left'><b>Process:</b> {title} "
        f"({len(units)} units)</caption>{header}{''.join(rows)}</table>"
    )


def register() -> bool:
    """Attach ``_repr_html_`` to NeqSim stream and process classes.

    Safe to call multiple times and safe to call when the JVM is not started
    (returns False without raising).

    Returns:
        True if the customizers were registered, False otherwise.
    """
    global _registered
    if _registered:
        return True
    try:
        import jpype

        if not jpype.isJVMStarted():
            return False

        @jpype.JImplementationFor("neqsim.process.equipment.stream.StreamInterface")
        class _StreamDisplay:  # noqa: D401 - JPype customizer
            def _repr_html_(self):
                return _stream_html(self)

        @jpype.JImplementationFor("neqsim.process.processmodel.ProcessSystem")
        class _ProcessDisplay:  # noqa: D401 - JPype customizer
            def _repr_html_(self):
                return _process_html(self)

        _registered = True
        return True
    except Exception:
        return False
