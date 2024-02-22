![downloads](https://img.shields.io/pypi/dw/wapordl) [![version](https://img.shields.io/pypi/v/wapordl)](https://anaconda.org/conda-forge/wapordl)

# WaPORDL

Download data from the WaPOR3 dataset as spatially aggregated timeseries or as spatial data clipped to a bounding-box or shapefile.

## Installation

### Conda (recommended)
`conda install -c conda-forge wapordl`

### Pip (make sure GDAL is already installed)
`pip install wapordl`

## Usage

To download a timeseries for a certain region:

```python
from wapordl import wapor_ts

region = "path/to/some/shape.geojson"
variable = "L2-AETI-D"
period = ["2021-01-01", "2021-07-01"]
overview = 3 # set to "NONE" to use native resolution data.

df = wapor_ts(region, variable, period, overview)

df

>>>     minimum  maximum    mean start_date   end_date number_of_days
>>> 0       0.0      7.3  0.5962 2021-01-01 2021-01-10        10 days
>>> 1       0.0      6.7  0.3923 2021-01-11 2021-01-20        10 days
>>> 2       0.0      6.7  0.3158 2021-01-21 2021-01-31        11 days
>>> ...
>>> 16      0.0      4.8  0.3192 2021-06-11 2021-06-20        10 days
>>> 17      0.0      4.9  0.4197 2021-06-21 2021-06-30        10 days
>>> 18      0.0      5.0  0.5727 2021-07-01 2021-07-10        10 days

df.attrs

>>> {'long_name': 'Actual EvapoTranspiration and Interception',
>>> 'units': 'mm/day',
>>> 'overview': 3}
```

To download a timerseries and convert its unit provide the `unit_conversion` keyword:

```python
unit = "dekad" # or choose "day", "month", "year", "none" (default).

df = wapor_ts(region, variable, period, overview, unit_conversion = unit)

df

>>>     minimum  maximum    mean start_date   end_date number_of_days
>>> 0       0.0     73.0  5.9617 2021-01-01 2021-01-10        10 days
>>> 1       0.0     67.0  3.9235 2021-01-11 2021-01-20        10 days
>>> 2       0.0     73.7  3.4740 2021-01-21 2021-01-31        11 days
>>> ...
>>> 16      0.0     48.0  3.1919 2021-06-11 2021-06-20        10 days
>>> 17      0.0     49.0  4.1972 2021-06-21 2021-06-30        10 days
>>> 18      0.0     50.0  5.7273 2021-07-01 2021-07-10        10 days

df.attrs

>>> {'long_name': 'Actual EvapoTranspiration and Interception',
>>> 'units': 'mm/dekad',
>>> 'overview': 3,
>>> 'original_units': 'mm/day'}
```

To download a geotiff for a certain region and period of time:

```python
from wapordl import wapor_map

region = "path/to/some/my_region.geojson"
folder = "path/to/some/output/folder"
variable = "L2-AETI-D"
period = ["2021-01-01", "2021-07-01"]

fp = wapor_map(region, variable, period, folder)

fp

>>> 'path/to/some/output/folder/my_region_L2-AETI-D_NONE_none.tif'
```

To download a timeseries and a netcdf for a bounding-box:

```python
region = [35.75, 33.70, 35.82, 33.75] # [xmin, ymin, xmax, ymax]
folder = "path/to/some/output/folder"
variable = "L3-AETI-D"
period = ["2021-01-01", "2021-07-01"]
overview = 3

df = wapor_ts(region, variable, period, overview)
fp = wapor_map(region, variable, period, folder, extension = ".nc")
```

When working with level-3 data, an entire L3 region can be downloaded by specifying a three letter region code:
    
```python
region = "BKA"
folder = "path/to/some/output/folder"
variable = "L3-T-D"
period = ["2021-01-01", "2021-07-01"]
overview = 3

df = wapor_ts(region, variable, period, overview)
fp = wapor_map(region, variable, period, folder, unit_conversion = "year")
```

## Upcoming

- Automatic overview selection based on the size of the shape.
- Support for variables with daily resolution (i.e. `L1-PCP-E` and `L1-RET-E`).
- Easily download a lower level variable for a level-3 region.
- Support for agERA5 variables.
- ~~Determine `l3_region` automatically from `region`.~~ ✅
- ~~Select unit of output.~~ ✅
- ~~Download a region from a bounding-box (i.e. without a shape).~~ ✅
- ~~A progress bar.~~ ✅
- ~~A warning if the given shape doesnt cover an area for which data is available.~~ ✅
- ~~Support for other output formats besides geotiff (e.g. netcdf).~~ ✅
- ~~Installation with conda.~~ ✅
- ~~More metadata in the output files.~~ ✅
- ~~More log information.~~ ✅
- ~~Option to select region for Level-3 data.~~ ✅

Got a feature-request or a question? Don't hesitate to contact me at bert.coerver@fao.org.