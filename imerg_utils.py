"""
imerg_utils – convert GPM IMERG half-hourly files into 3- or 6-hourly
accumulated rainfall (mm) at station locations.

Used by the notebooks in ``notebooks/``. All times are UTC unless a
``utc_offset_hours`` is given.

Conventions
-----------
* An IMERG file is named by the START of its 30-minute period,
  e.g. ``...3IMERG.20240719-S053000-E055959...`` covers 05:30-06:00 UTC.
* An accumulation row is labelled by the END of its period:
  the 6-hour total labelled 06:00 = rain from 00:00 to 06:00
  = the 12 files starting 00:00, 00:30, ..., 05:30.
* IMERG ``precipitation`` is a rate in mm/hr, so each file adds rate x 0.5 h.
"""

from __future__ import annotations

import glob
import os
import re
from datetime import datetime, timedelta

import numpy as np
import pandas as pd

__all__ = [
    "PRODUCTS",
    "make_bin_ends",
    "required_file_times",
    "download_window",
    "index_imerg_files",
    "load_stations",
    "ImergReader",
    "station_accumulation",
    "merge_into_table",
    "scores",
    "download_imerg",
]

# earthaccess short names of the half-hourly IMERG V07 runs
PRODUCTS = {
    "final": "GPM_3IMERGHH",   # research quality, gauge-adjusted, ~3.5 months latency
    "late": "GPM_3IMERGHHL",   # ~14 hours latency
    "early": "GPM_3IMERGHHE",  # ~4 hours latency
}

_NAME_RE = re.compile(r"3IMERG\.(\d{8})-S(\d{6})")


# --------------------------------------------------------------------------- times
def make_bin_ends(start, end, hours: int = 6, utc_offset_hours: float = 0) -> pd.DatetimeIndex:
    """End times (UTC) of the accumulation periods that tile ``[start, end)``.

    ``start`` and ``end`` are given in local time when ``utc_offset_hours`` is
    set (e.g. 6 for Bangladesh), otherwise in UTC.
    """
    if hours <= 0 or (hours * 2) != int(hours * 2):
        raise ValueError("hours must be a positive multiple of 0.5")
    offset = pd.Timedelta(hours=utc_offset_hours)
    start = pd.Timestamp(start) - offset
    end = pd.Timestamp(end) - offset
    step = pd.Timedelta(hours=hours)
    if end <= start:
        raise ValueError("end must be after start")
    if (end - start) % step != pd.Timedelta(0):
        raise ValueError(f"the period {start} -> {end} (UTC) is not a whole number of {hours}-hour bins")
    return pd.date_range(start + step, end, freq=step)


def required_file_times(bin_ends, hours: int = 6) -> list[datetime]:
    """Start times of every half-hourly file needed for the given bins."""
    need = set()
    n = int(hours * 2)
    for e in pd.DatetimeIndex(bin_ends):
        e = e.to_pydatetime()
        need.update(e - timedelta(minutes=30 * k) for k in range(1, n + 1))
    return sorted(need)


def download_window(bin_ends, hours: int = 6) -> tuple[pd.Timestamp, pd.Timestamp]:
    """(first file start, last file end) to download for the given bins."""
    b = pd.DatetimeIndex(bin_ends)
    return b.min() - pd.Timedelta(hours=hours), b.max()


# --------------------------------------------------------------------------- files
def index_imerg_files(folder: str) -> dict[datetime, str]:
    """Map file start time (UTC) -> path for all IMERG half-hourly files in ``folder``.

    Works for Final (``.HDF5``) and Early/Late (``.RT-H5``) file names.
    """
    files = {}
    for path in glob.glob(os.path.join(folder, "*3IMERG*")):
        if not path.upper().endswith((".HDF5", ".RT-H5", ".H5")):
            continue
        m = _NAME_RE.search(os.path.basename(path))
        if m:
            files[datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S")] = path
    return dict(sorted(files.items()))


def load_stations(path: str, id_col="station_id", name_col="name", lat_col="lat", lon_col="lon") -> pd.DataFrame:
    """Read a station list and return columns station_id, name, lat, lon."""
    df = pd.read_csv(path)
    missing = [c for c in (id_col, lat_col, lon_col) if c not in df.columns]
    if missing:
        raise KeyError(f"{path}: missing column(s) {missing}; found {list(df.columns)}")
    out = pd.DataFrame({
        "station_id": df[id_col],
        "name": df[name_col] if name_col in df.columns else df[id_col].astype(str),
        "lat": df[lat_col].astype(float),
        "lon": df[lon_col].astype(float),
    })
    return out.drop_duplicates("station_id").reset_index(drop=True)


# --------------------------------------------------------------------------- reading
def _h5_read_grid(path):
    import h5py
    with h5py.File(path, "r") as h:
        return h["Grid/lat"][:], h["Grid/lon"][:]


def _h5_read_precip(path, lon_slice, lat_slice):
    """mm/hr array (lon, lat) for the requested box."""
    import h5py
    with h5py.File(path, "r") as h:
        name = "precipitation" if "Grid/precipitation" in h else "precipitationCal"  # V07 / V06
        return h[f"Grid/{name}"][0, lon_slice, lat_slice].astype("float64")


class ImergReader:
    """Reads only a small box around the stations from each file, and caches it."""

    def __init__(self, files: dict, stations: pd.DataFrame, window: int = 3, margin_deg: float = 0.5):
        if window < 1 or window % 2 == 0:
            raise ValueError("window must be an odd number: 1 (nearest cell), 3 (3x3), 5 (5x5) ...")
        if not files:
            raise FileNotFoundError("no IMERG files given")
        self.files = files
        self.half = window // 2
        lat, lon = _h5_read_grid(next(iter(files.values())))
        self.lat, self.lon = np.asarray(lat), np.asarray(lon)

        lo = lambda axis, v: max(int(np.searchsorted(axis, v)), 0)
        hi = lambda axis, v: min(int(np.searchsorted(axis, v)), len(axis))
        self.lon_slice = slice(lo(self.lon, stations.lon.min() - margin_deg), hi(self.lon, stations.lon.max() + margin_deg))
        self.lat_slice = slice(lo(self.lat, stations.lat.min() - margin_deg), hi(self.lat, stations.lat.max() + margin_deg))

        # nearest cell centre of each station, relative to the box
        self.ix = np.array([np.abs(self.lon - v).argmin() for v in stations.lon]) - self.lon_slice.start
        self.iy = np.array([np.abs(self.lat - v).argmin() for v in stations.lat]) - self.lat_slice.start
        nx = self.lon_slice.stop - self.lon_slice.start
        ny = self.lat_slice.stop - self.lat_slice.start
        if ((self.ix < self.half) | (self.ix >= nx - self.half) | (self.iy < self.half) | (self.iy >= ny - self.half)).any():
            raise ValueError("a station is too close to the edge of the IMERG grid for this window size")
        self._cache = {}

    def halfhour_mm(self, t: datetime):
        """Rain (mm) in the 30 minutes starting at ``t`` over the box, or None if the file is missing."""
        if t not in self._cache:
            path = self.files.get(t)
            if path is None:
                self._cache[t] = None
            else:
                p = _h5_read_precip(path, self.lon_slice, self.lat_slice)
                p[p < 0] = np.nan            # fill value -9999.9
                self._cache[t] = p * 0.5     # mm/hr x 0.5 h
        return self._cache[t]

    def station_values(self, total: np.ndarray) -> np.ndarray:
        """Mean of the window x window cells around each station."""
        h = self.half
        if h == 0:
            return total[self.ix, self.iy]
        return np.array([np.nanmean(total[x - h:x + h + 1, y - h:y + h + 1]) for x, y in zip(self.ix, self.iy)])


# --------------------------------------------------------------------------- main step
def station_accumulation(files: dict, stations: pd.DataFrame, bin_ends, hours: int = 6,
                         window: int = 3, utc_offset_hours: float = 0, decimals: int = 3):
    """Accumulated IMERG rainfall (mm) per station and period.

    Returns ``(table, missing)``: a long table with one row per station and
    period, and the sorted list of file start times that were missing.
    A period with any missing file is left empty (NaN).
    """
    reader = ImergReader(files, stations, window=window)
    n = int(hours * 2)
    rows, missing = [], set()
    for end in pd.DatetimeIndex(bin_ends):
        end = end.to_pydatetime()
        total, ok = None, True
        for k in range(n, 0, -1):
            t = end - timedelta(minutes=30 * k)
            b = reader.halfhour_mm(t)
            if b is None:
                missing.add(t)
                ok = False
                continue
            total = b.copy() if total is None else total + b
        vals = reader.station_values(total) if ok else np.full(len(stations), np.nan)
        rows.append(pd.DataFrame({
            "station_id": stations.station_id.values,
            "name": stations.name.values,
            "lat": stations.lat.values,
            "lon": stations.lon.values,
            "datetime_end_utc": end,
            "imerg_mm": np.round(vals, decimals),
        }))
    table = pd.concat(rows, ignore_index=True)
    if utc_offset_hours:
        table.insert(5, "datetime_end_local", table.datetime_end_utc + pd.Timedelta(hours=utc_offset_hours))
    table = table.sort_values(["station_id", "datetime_end_utc"]).reset_index(drop=True)
    return table, sorted(missing)


def merge_into_table(obs: pd.DataFrame, imerg: pd.DataFrame, obs_station_col: str, obs_time_col: str,
                     obs_time_format: str | None = None, after_col: str | None = None,
                     new_col: str = "IMERG", utc_offset_hours: float = 0) -> pd.DataFrame:
    """Add IMERG values to your own table (e.g. station observations).

    Rows are matched on station ID and period end time. ``obs_time_col``
    must be the END of the accumulation period, in UTC unless
    ``utc_offset_hours`` says it is local time. The new column is placed
    right after ``after_col`` (or at the end).
    """
    out = obs.copy()
    t = pd.to_datetime(out[obs_time_col], format=obs_time_format) - pd.Timedelta(hours=utc_offset_hours)
    key = pd.DataFrame({"station_id": out[obs_station_col].astype(str).values, "datetime_end_utc": t.values})
    right = imerg[["station_id", "datetime_end_utc", "imerg_mm"]].copy()
    right["station_id"] = right["station_id"].astype(str)
    vals = key.merge(right, on=["station_id", "datetime_end_utc"], how="left")["imerg_mm"].values
    if new_col in out.columns:
        out = out.drop(columns=new_col)
    pos = out.columns.get_loc(after_col) + 1 if after_col else len(out.columns)
    out.insert(pos, new_col, vals)
    return out


# --------------------------------------------------------------------------- scores
def scores(obs, est, wet_mm: float = 1.0) -> dict:
    """Continuous and rain/no-rain scores of ``est`` against ``obs``.

    POD = hits / (hits + misses), FAR = false alarms / (hits + false alarms),
    CSI = hits / (hits + misses + false alarms); "rain" means >= ``wet_mm``.
    """
    o, e = np.asarray(obs, float), np.asarray(est, float)
    ok = ~np.isnan(o) & ~np.isnan(e)
    o, e = o[ok], e[ok]
    hit = int(((o >= wet_mm) & (e >= wet_mm)).sum())
    miss = int(((o >= wet_mm) & (e < wet_mm)).sum())
    fa = int(((o < wet_mm) & (e >= wet_mm)).sum())
    div = lambda a, b: round(a / b, 3) if b else np.nan
    return {
        "n": int(ok.sum()),
        "pearson_r": round(float(np.corrcoef(o, e)[0, 1]), 3) if len(o) > 1 and o.std() > 0 and e.std() > 0 else np.nan,
        "rmse_mm": round(float(np.sqrt(np.mean((e - o) ** 2))), 3) if len(o) else np.nan,
        "mae_mm": round(float(np.mean(np.abs(e - o))), 3) if len(o) else np.nan,
        "bias_ratio": div(e.sum(), o.sum()),
        "hits": hit, "misses": miss, "false_alarms": fa,
        "pod": div(hit, hit + miss), "far": div(fa, hit + fa), "csi": div(hit, hit + miss + fa),
    }


# --------------------------------------------------------------------------- download
def download_imerg(start, end, out_dir: str, product: str = "final", version: str = "07"):
    """Download every half-hourly file whose period starts in ``[start, end)`` (UTC).

    Needs a free NASA Earthdata account and ``earthaccess``. Files already
    in ``out_dir`` are skipped. Returns the list of local paths.
    """
    import earthaccess

    start, end = pd.Timestamp(start), pd.Timestamp(end)
    os.makedirs(out_dir, exist_ok=True)
    granules = earthaccess.search_data(
        short_name=PRODUCTS[product],
        version=version,
        temporal=(start.strftime("%Y-%m-%d %H:%M:%S"),
                  (end - pd.Timedelta(seconds=1)).strftime("%Y-%m-%d %H:%M:%S")),
    )
    expected = int((end - start) / pd.Timedelta(minutes=30))
    print(f"{start} -> {end}: {len(granules)} files found, {expected} expected")
    return earthaccess.download(granules, out_dir)
