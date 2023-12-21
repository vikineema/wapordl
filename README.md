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

    df = wapor_ts(region, variable, period, overview = overview)
    
    df

    >>>     minimum  maximum    mean  stddev       date
    >>> 0       1.0     44.0  23.127   5.269 2021-01-01
    >>> 1       0.0     48.0  23.881   6.102 2021-01-11
    >>> 2       1.0     44.0  21.329   5.249 2021-01-21
    >>> ...
    >>> 16      0.0     45.0  22.991   7.253 2021-06-11
    >>> 17      0.0     46.0  23.644   8.104 2021-06-21
    >>> 18      0.0     54.0  23.253   8.932 2021-07-01


To download a geotiff for a certain region and period of time:

    region = "path/to/some/my_region.geojson"
    folder = "path/to/some/output/folder"
    variable = "L2-AETI-D"
    period = ["2021-01-01", "2021-07-01"]

    fp = wapor_map(region, variable, folder, period)

    fp

    >>> 'path/to/some/output/folder/my_region_L2-AETI-D_NONE.tif'

## Upcoming

- Download a region from a bounding-box (i.e. without a shape).
- A progress bar.
- A warning if the given shape doesnt cover an area for which data is available.
- Automatic overview selection based on the size of the shape.
- Support for other input formats besides geojson (e.g. shapefiles, geopackages etc.)
- Support for other output formats besides geotiff (e.g. netcdf).
- ~~Installation with conda.~~ ✅
- More metadata in the output files.
- More log information.
- Option to select region for Level-3 data.