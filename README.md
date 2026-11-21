# GPM IMERG → station rainfall (3- or 6-hourly CSV)

Turn NASA **GPM IMERG** half-hourly satellite rainfall files (HDF5) into **3- or 6-hourly
accumulated rainfall (mm) at your station locations**, saved as a simple CSV that you can
compare with rain-gauge data or model forecasts.

```
station_id,name,lat,lon,datetime_end_utc,imerg_mm
<one row per station and period: the IMERG rain total in mm for the period ending at datetime_end_utc>
```

You only need a list of stations (ID, name, latitude, longitude) and a free NASA Earthdata account.
No observation data is required.

## Features

- Downloads exactly the IMERG V07 files you need (Final, Late or Early run) with `earthaccess`, including the
  6 hours before your first period, so that no row is left empty.
- 3- or 6-hourly totals (any multiple of 30 minutes works), aligned to UTC or to local time.
- Value at each station from the **nearest 0.1° cell** or the **mean of the 3 × 3 (or 5 × 5) cells** around it.
- Reads only a small box around your stations from each file, so it is fast even for thousands of files.
- Reports every missing file and leaves the affected periods empty instead of giving wrong totals.
- Optional: adds an `IMERG` column to your own observation table and prints verification scores
  (Pearson r, RMSE, MAE, bias, POD, FAR, CSI).

## Repository layout

```
GPM_IMERG_Data_Conversion/
├── imerg_utils.py                      # all the processing functions (used by the notebooks)
├── notebooks/
│   ├── 01_download_imerg.ipynb         # download the half-hourly HDF5 files
│   ├── 02_imerg_to_station_csv.ipynb   # half-hourly → 3/6-hourly CSV at stations
│   └── 03_merge_with_observations.ipynb  # optional: add IMERG to your gauge table and score it
├── examples/
│   ├── stations_bangladesh.csv         # 41 Bangladesh Meteorological Department stations (ID, name, lat, lon, elevation)
│   ├── periods.csv                     # periods to process (a 30-hour demo period)
│   └── observations_template.csv       # layout for notebook 3
├── tests/test_imerg_utils.py
├── requirements.txt  /  environment.yml
└── LICENSE
```

## Quick start

```bash
git clone https://github.com/ArnobSen85/GPM_IMERG_Data_Conversion.git
cd GPM_IMERG_Data_Conversion
pip install -r requirements.txt          # or: conda env create -f environment.yml && conda activate imerg
jupyter notebook notebooks/
```

1. Create a free NASA Earthdata account at <https://urs.earthdata.nasa.gov/users/new>.
2. Open **`01_download_imerg.ipynb`** and run all cells. You are asked for your Earthdata username and
   password once. The demo period (2024-07-19 00:00 to 2024-07-20 06:00 UTC) needs 60 files, about 0.5 GB.
3. Open **`02_imerg_to_station_csv.ipynb`** and run all cells. The result is written to
   `output/imerg_6h_demo_2024_07_19.csv`.

## Using your own stations and periods

**Stations**: a CSV with at least these columns (extra columns are ignored):

```
station_id,name,lat,lon
41923,Dhaka,23.7667,90.3833
```

**Periods**: one row per period, in UTC (or local time, see below). The length of each
period must be a whole number of accumulation steps.

```
name,start,end
monsoon_event_1,2024-07-19 00:00,2024-07-20 06:00
```

Then set these at the top of notebooks 1 and 2 (use the same values in both):

| Setting | Meaning |
|---|---|
| `HOURS` | accumulation length, e.g. `3` or `6` |
| `UTC_OFFSET_HOURS` | `0` if your periods are in UTC; e.g. `6` for Bangladesh time (UTC+6) |
| `PRODUCT` | `'final'` (gauge-adjusted, about 3.5 months behind), `'late'` (about 14 h) or `'early'` (about 4 h) |
| `WINDOW` | `1` = nearest cell, `3` = mean of 3 × 3 cells (default), `5` = 5 × 5 |

## How the values are computed

- IMERG files are named by the **start** of their 30-minute period: `...20240719-S053000-E055959...`
  covers 05:30–06:00 UTC.
- Each output row is labelled by the **end** of its accumulation period. The 6-hour total at
  `2024-07-19 06:00` is the rain from 00:00 to 06:00 UTC, which is the sum of the 12 files starting
  00:00, 00:30, …, 05:30.
- The IMERG `precipitation` variable is a rate in mm/hr, so each file contributes `rate × 0.5 h`.
  Fill values (−9999.9) are treated as missing.
- The 6-hour total is first computed for every grid cell. The station value is then either the nearest
  cell or the mean of the `WINDOW × WINDOW` cells around it.
- If any file of a period is missing, that period is left empty and the file is listed.

**Why a 3 × 3 mean by default?** A rain gauge measures one point, while one IMERG cell is about
10 km × 10 km. Averaging a few cells reduces the effect of small storms falling just inside or
outside the cell, and usually matches gauge totals a little better. Use `WINDOW = 1` if you need
the raw cell value.

## Optional: compare with your own observations

`03_merge_with_observations.ipynb` reads your own table (station ID, period end time, observed rain)
and adds an `IMERG` column right after the column you choose. Time formats like
`2024:07:19-06:00:00` are supported through `OBS_TIME_FORMAT`. It then prints:

| Score | Meaning |
|---|---|
| Pearson r, RMSE, MAE | agreement of the amounts |
| Bias | IMERG total ÷ observed total |
| POD, FAR, CSI | detection of rain ≥ `WET_MM` (default 1 mm) |

## Using the functions directly

```python
import imerg_utils as iu

stations = iu.load_stations('examples/stations_bangladesh.csv')
files    = iu.index_imerg_files('data/imerg')
bins     = iu.make_bin_ends('2024-07-19 00:00', '2024-07-20 06:00', hours=6)
table, missing = iu.station_accumulation(files, stations, bins, hours=6, window=3)
table.to_csv('imerg_6h.csv', index=False)
```

## Tests

```bash
pytest -q
```

The tests write small IMERG-like HDF5 files and check the totals, the window mean, missing-file
handling, time alignment and the merge.

## Data and citation

IMERG data are produced by NASA and are free to use. If you use them, cite the dataset:

> Huffman, G. J., E. F. Stocker, D. T. Bolvin, E. J. Nelkin, and J. Tan (2023). *GPM IMERG Final Precipitation
> L3 Half Hourly 0.1 degree x 0.1 degree V07.* Greenbelt, MD, Goddard Earth Sciences Data and Information
> Services Center (GES DISC). https://doi.org/10.5067/GPM/IMERG/3B-HH/07

Dataset page: <https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGHH_07/summary>

The station coordinates in `examples/stations_bangladesh.csv` are for weather stations of the
Bangladesh Meteorological Department (BMD). No observation data is included in this repository.

## License

MIT. See [LICENSE](LICENSE).
