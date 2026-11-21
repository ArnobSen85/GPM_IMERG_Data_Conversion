"""Tests with small synthetic IMERG-like HDF5 files. Run from the repo root:  pytest -q"""
import os
import sys
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import pytest

h5py = pytest.importorskip("h5py")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import imerg_utils as iu  # noqa: E402

LAT = np.round(np.arange(20.05, 27.0, 0.1), 2)   # small 0.1-degree grid over Bangladesh
LON = np.round(np.arange(87.05, 93.0, 0.1), 2)


def rate(t):
    """Synthetic rain rate (mm/hr) that varies in space and time."""
    i = np.arange(len(LON))[:, None]
    j = np.arange(len(LAT))[None, :]
    return 0.01 * i + 0.001 * j + t.hour


def write_file(folder, t, data):
    e = t + timedelta(minutes=29, seconds=59)
    name = f"3B-HHR.MS.MRG.3IMERG.{t:%Y%m%d}-S{t:%H%M%S}-E{e:%H%M%S}.0000.V07B.HDF5"
    with h5py.File(os.path.join(folder, name), "w") as h:
        g = h.create_group("Grid")
        g["lat"], g["lon"] = LAT, LON
        g["precipitation"] = data[None].astype("float32")   # (time, lon, lat) like IMERG


@pytest.fixture
def setup(tmp_path):
    t = datetime(2024, 7, 19, 0)
    while t < datetime(2024, 7, 19, 12):
        write_file(tmp_path, t, rate(t))
        t += timedelta(minutes=30)
    st = pd.DataFrame({"station_id": [1, 2], "name": ["A", "B"], "lat": [23.75, 22.35], "lon": [90.38, 91.82]})
    return tmp_path, st


def expected(st_row, end, hours, window):
    i = int(np.abs(LON - st_row.lon).argmin()); j = int(np.abs(LAT - st_row.lat).argmin()); h = window // 2
    total = sum(rate(end - timedelta(minutes=30 * k)) * 0.5 for k in range(1, hours * 2 + 1))
    return total[i - h:i + h + 1, j - h:j + h + 1].mean()


@pytest.mark.parametrize("hours,window", [(6, 1), (6, 3), (3, 1), (3, 3)])
def test_accumulation(setup, hours, window):
    folder, st = setup
    files = iu.index_imerg_files(folder)
    assert len(files) == 24
    bins = iu.make_bin_ends("2024-07-19 00:00", "2024-07-19 12:00", hours)
    table, missing = iu.station_accumulation(files, st, bins, hours, window, decimals=6)
    assert missing == []
    assert len(table) == len(st) * len(bins)
    for r in table.itertuples():
        s = st[st.station_id == r.station_id].iloc[0]
        assert r.imerg_mm == pytest.approx(expected(s, r.datetime_end_utc.to_pydatetime(), hours, window), abs=1e-5)


def test_missing_file_gives_nan(setup):
    folder, st = setup
    os.remove(next(p for p in folder.iterdir() if "S013000" in p.name))
    files = iu.index_imerg_files(folder)
    bins = iu.make_bin_ends("2024-07-19 00:00", "2024-07-19 12:00", 6)
    table, missing = iu.station_accumulation(files, st, bins, 6, 3)
    assert missing == [datetime(2024, 7, 19, 1, 30)]
    first = table[table.datetime_end_utc == pd.Timestamp("2024-07-19 06:00")]
    assert first.imerg_mm.isna().all()
    assert table[table.datetime_end_utc == pd.Timestamp("2024-07-19 12:00")].imerg_mm.notna().all()


def test_bins_and_window():
    bins = iu.make_bin_ends("2024-07-19 06:00", "2024-07-20 06:00", 6, utc_offset_hours=6)
    assert bins[0] == pd.Timestamp("2024-07-19 06:00") and len(bins) == 4
    assert iu.download_window(iu.make_bin_ends("2024-07-19 00:00", "2024-07-19 12:00", 6), 6) == (
        pd.Timestamp("2024-07-19 00:00"), pd.Timestamp("2024-07-19 12:00"))
    with pytest.raises(ValueError):
        iu.make_bin_ends("2024-07-19 03:00", "2024-07-19 12:00", 6)


def test_merge(setup):
    folder, st = setup
    files = iu.index_imerg_files(folder)
    bins = iu.make_bin_ends("2024-07-19 00:00", "2024-07-19 12:00", 6)
    table, _ = iu.station_accumulation(files, st, bins, 6, 3)
    obs = pd.DataFrame({"sid": [1, 2, 1], "time": ["2024-07-19 06:00", "2024-07-19 12:00", "2024-07-20 00:00"],
                        "observed": [1.0, 2.0, 3.0], "other": [0, 0, 0]})
    m = iu.merge_into_table(obs, table, "sid", "time", after_col="observed")
    assert list(m.columns) == ["sid", "time", "observed", "IMERG", "other"]
    assert m.IMERG.iloc[:2].notna().all() and np.isnan(m.IMERG.iloc[2])


def test_scores():
    s = iu.scores([0, 2, 5, 0], [0.5, 3, 4, 2])
    assert (s["hits"], s["misses"], s["false_alarms"]) == (2, 0, 1)
    assert s["csi"] == pytest.approx(2 / 3, abs=1e-3)
