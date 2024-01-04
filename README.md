![downloads](https://img.shields.io/pypi/dw/wapordl) [![version](https://img.shields.io/pypi/v/wapordl)](https://anaconda.org/conda-forge/wapordl)

# WaPORDL

This package allows users to download data from the WaPOR3 dataset as spatially aggregated timeseries or as spatial data clipped to a bounding-box or shapefile.

## Installation

### Conda (recommended)
`conda install -c conda-forge wapordl`

### Pip (make sure GDAL is already installed)
`pip install wapordl`

## Usage

To download a timeseries for a certain region:

    region = "path/to/some/shape.geojson"
    variable = "L2-AETI-D"
    period = ["2021-01-01", "2021-07-01"]
    overview = 3

    df = wapor_ts(region, variable, period, overview)
    
    df

    >>>     minimum  maximum   mean        date
    >>> 0       0.1     4.40  2.312  2021-01-01
    >>> 1       0.0     4.80  2.388  2021-01-11
    >>> 2       0.1     4.40  2.132  2021-01-21
    >>> ...
    >>> 16      0.0     4.50  2.299  2021-06-11
    >>> 17      0.0     4.60  2.364  2021-06-21
    >>> 18      0.0     5.40  2.325  2021-07-01

    df.attrs

    >>> {'long_name': 'Actual EvapoTranspiration and Interception',
    >>> 'units': 'mm/day',
    >>> 'overview': 3}

To download a geotiff for a certain region and period of time:

    region = "path/to/some/my_region.geojson"
    folder = "path/to/some/output/folder"
    variable = "L2-AETI-D"
    period = ["2021-01-01", "2021-07-01"]

    fp = wapor_map(region, variable, period, folder)

    fp

    >>> 'path/to/some/output/folder/my_region_L2-AETI-D_NONE.tif'

To download a timeseries and a netcdf for a bounding-box:

    bb = [25, -17, 26, -16]
    folder = "path/to/some/output/folder"
    variable = "L2-AETI-D"
    period = ["2021-01-01", "2021-07-01"]
    overview = 3

    df = wapor_ts(bb, variable, period, overview)
    fp = wapor_map(bb, variable, period, folder, extenstion = ".nc")

When working with level-3 data, a region-code needs to be specified and the region argument can (optionally) be set to `None`:
    
    region = "path/to/some/my_region.geojson"
    folder = "path/to/some/output/folder"
    variable = "L3-T-D"
    period = ["2021-01-01", "2021-07-01"]
    overview = 3
    l3_region = "BKA"

    df = wapor_ts(region, variable, period, overview, l3_region = l3_region)
    fp = wapor_map(None, variable, period, folder, l3_region = l3_region)

## Upcoming

- Automatic overview selection based on the size of the shape.
- ~~Download a region from a bounding-box (i.e. without a shape).~~ ✅
- ~~A progress bar.~~ ✅
- ~~A warning if the given shape doesnt cover an area for which data is available.~~ ✅
- ~~Support for other output formats besides geotiff (e.g. netcdf).~~ ✅
- ~~Installation with conda.~~ ✅
- ~~More metadata in the output files.~~ ✅
- ~~More log information.~~ ✅
- ~~Option to select region for Level-3 data.~~ ✅

Got a feature-request or a question? Don't hesitate to contact me at bert.coerver@fao.org.