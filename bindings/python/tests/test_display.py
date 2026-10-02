"""Tests for the rich Jupyter display (neqsim.display)."""

from neqsim.process import clearProcess, getProcess, runProcess, separator, stream
from neqsim.thermo import fluid


def _build():
    clearProcess()
    f = fluid("srk")
    f.addComponent("methane", 0.9)
    f.addComponent("ethane", 0.1)
    f.setTemperature(30.0, "C")
    f.setPressure(50.0, "bara")
    s = stream("feed", f)
    separator("sep", s)
    runProcess()
    return s


def test_stream_repr_html():
    s = _build()
    html = s._repr_html_()
    assert "<table" in html
    assert "feed" in html
    assert "Pressure" in html


def test_process_repr_html():
    _build()
    html = getProcess()._repr_html_()
    assert "<table" in html
    assert "sep" in html
