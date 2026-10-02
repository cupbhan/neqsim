"""Tests for the neqsim.process pandas DataFrame helpers."""

import pandas as pd

from neqsim.process import (
    clearProcess,
    compressor,
    cooler,
    equipment_table,
    runProcess,
    stream,
    stream_table,
    to_dataframe,
)
from neqsim.thermo import fluid


def _build_process():
    clearProcess()
    f = fluid("srk")
    f.addComponent("methane", 0.9)
    f.addComponent("ethane", 0.1)
    f.setTemperature(30.0, "C")
    f.setPressure(10.0, "bara")
    inlet = stream("inlet", f)
    comp = compressor("comp1", inlet, pres=50.0)
    cool = cooler("cool1", comp.getOutletStream())
    cool.setOutTemperature(298.15)
    runProcess()


def test_stream_table_columns_and_rows():
    _build_process()
    df = stream_table()
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    for col in ("name", "flow_kg_hr", "temperature_C", "pressure_bara"):
        assert col in df.columns
    # inlet + compressor out + cooler out
    assert len(df) >= 3


def test_to_dataframe_is_stream_table_alias():
    assert to_dataframe is stream_table


def test_equipment_table_lists_units():
    _build_process()
    df = equipment_table()
    assert isinstance(df, pd.DataFrame)
    types = set(df["type"].tolist())
    assert "Compressor" in types
    assert "Cooler" in types


def test_stream_table_values_reasonable():
    _build_process()
    df = stream_table().set_index("name")
    # Compressor discharge pressure was set to 50 bara.
    assert df.loc["comp1 out stream", "pressure_bara"] == 50.0
