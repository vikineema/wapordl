import os
import requests
import logging
import shapely
import numpy as np
import pandas as pd
from osgeo import gdal
gdal.UseExceptions()

def collect_responses(url, info = ["code"]):
    data = {"links": [{"rel": "next", "href": url}]}
    output = list()
    while "next" in [x["rel"] for x in data["links"]]:
        url_ = [x["href"] for x in data["links"] if x["rel"] == "next"][0]
        response = requests.get(url_)
        response.raise_for_status()
        data = response.json()["response"]
        if isinstance(info, list):
            output += [tuple(x.get(y) for y in info) for x in data["items"]]
        else:
            output += data["items"]
    if isinstance(info, list):
        output = sorted(output)
    return output

def generate_urls_v3(variable):
    if "L1" in variable:
        base_url = f"https://data.apps.fao.org/gismgr/api/v2/catalog/workspaces/WAPOR-3/mapsets"
    elif "L2" in variable:
        base_url = f"https://data.apps.fao.org/gismgr/api/v2/catalog/workspaces/WAPOR-3/mapsets"
    elif "L3" in variable:
        base_url = f"https://data.apps.fao.org/gismgr/api/v2/catalog/workspaces/WAPOR-3/mosaicsets"
    else:
        raise ValueError
    mapset_url = f"{base_url}/{variable}/rasters"
    urls = [x[0] for x in collect_responses(mapset_url, info = ["downloadUrl"])]
    return tuple(sorted(urls))

def wapor_dl(region, variable,
             l3_region = None,
             period = ["2021-01-01", "2022-01-01"], 
             overview = "NONE", 
             req_stats = ["minimum", "maximum", "mean"],
             folder = None):
    """_summary_

    Parameters
    ----------
    region : str, list
        Path to a geojson file, or a list of floats specifying a bounding-box [<xmin> <ymin> <xmax> <ymax>].
    variable : str
        Name of the variable to download.
    period : list, optional
        List of a start and end date, by default ["2021-01-01", "2022-01-01"]
    overview : str, int, optional
        Which overview to use, specify "NONE" to not use an overview, 0 uses the first overview, etc., by default "NONE"
    req_stats : list, optional
        Specify which statistics to export, by default ["minimum", "maximum", "mean"]
    folder : str, optional
        Folder to store output files, by default None

    Returns
    -------
    str, pd.Dataframe
        If `req_stats` is not None, returns a pd.Dataframe. Otherwise a path to file is returned.
    """
    ## Retrieve info from variable name.
    level, var_code, tres = variable.split("-")

    ## Check l3_region code.
    if level == "L3" and isinstance(l3_region, type(None)):
        raise ValueError(f"Please specify a `l3_region` for a {level} variable ({variable})")
    if not isinstance(l3_region, type(None)):
        valid_l3_regions = ['BKA', 'ERB', 'GAR', 'GEZ', 'JAF', 'JEN', 'JVA', 'KAI', 'ODN',
                                'PAL', 'SED']
        if not bool(np.isin(l3_region, valid_l3_regions)):
            raise ValueError(f"Invalid `l3_region`, please specify one of {valid_l3_regions}.")
        if level != "L3":
            logging.warning("l3_region will be ignored (no l3 variable selected).")
    
    ## Check if region is valid.
    if isinstance(region, str):
        if not os.path.isfile(region) or os.path.splitext(region)[-1] != ".geojson":
            raise ValueError(f"Geojson file not found.")
        else:
            region_code = os.path.split(region)[-1].replace(".geojson", "")
            region_shape = shapely.from_geojson(open(region,'r').read())
    elif isinstance(region, list):
        if not all([region[2] > region[0], region[3] > region[1]]):
            raise ValueError(f"Invalid bounding box.")
        else:
            region_code = "bb"
            region_shape = shapely.Polygon([(region[0], region[1]), 
                                            (region[2], region[1]), 
                                            (region[2], region[3]), 
                                            (region[0], region[3]), 
                                            (region[0], region[1])])
    elif isinstance(region, type(None)):
        region_code = l3_region
        region_shape = None
        if level != "L3":
            raise ValueError("Only level-3 variables can be processed without specifying a region.")

    ## Parse the dates in period.
    if not isinstance(period, type(None)):
        period = [pd.Timestamp(x) for x in period]
        if period[0] > period[1]:
            raise ValueError(f"Invalid period.")

    ## Define function to get pd.Timestamp from urls.
    if tres == "D":
        def date_func(url):
            year, month, dekad = os.path.split(url)[-1].split(".")[-2].split("-")
            day = {'D1': '01', 'D2': '11', 'D3': '21'}[dekad]
            return pd.Timestamp(f"{year}-{month}-{day}")
    elif tres == "M":
        def date_func(url):
            year, month = os.path.split(url)[-1].split(".")[-2].split("-")
            return pd.Timestamp(f"{year}-{month}-01")
    elif tres == "A":
        def date_func(url):
            year = os.path.split(url)[-1].split(".")[-2]
            return pd.Timestamp(f"{year}-01-01")
    else:
        raise ValueError("Invalid temporal resolution.")

    ## Collect urls for requested variable.
    public_urls = generate_urls_v3(variable)

    ## TODO: Check if these filters can be applied directly in the API Query.
    ## Filter L3 region.
    if level == "L3" and not isinstance(l3_region, type(None)):
        region_urls = [url for url in public_urls if l3_region in url]
    else:
        region_urls = public_urls

    ## Determine date for each url and filter on requested period.
    date_urls = [(date_func(url), url) for url in region_urls]
    if not isinstance(period, type(None)):
        date_urls = [x for x in date_urls if (x[0] >= period[0]) & (x[0] <= period[1])]
        ## Print overview statement.
        print(f"Found {len(date_urls)} files for {variable} between {period[0]} and {period[1]}.")
    else:
        print(f"Found {len(date_urls)} files for {variable}.")

    ## Determine required output resolution.
    # NOTE maybe move this to external function (assumes info the same for all urls)
    info = gdal.Info("/vsicurl/" + date_urls[0][1], format = "json")
    overview_ = -1 if overview == "NONE" else overview
    xres, yres = info["geoTransform"][1::4]

    ## Check if region overlaps with datasets bounding-box.
    if not isinstance(region_shape, type(None)):
        data_bb = shapely.Polygon(np.array(info["wgs84Extent"]["coordinates"])[0])
        if not data_bb.intersects(region_shape):
            info_lbl1 = region_code if region_code != "bb" else str(region)
            info_lbl2 = variable if isinstance(l3_region, type(None)) else f"{variable}.{l3_region}"
            raise ValueError(f"Selected region ({info_lbl1}) has no overlap with the datasets ({info_lbl2}) bounding-box.")

    ## Get scale and offset factor.
    scale = info["bands"][0]["scale"]
    offset = info["bands"][0]["offset"]

    ## Check offset factor.
    if offset != 0:
        logging.warning("Offset factor is not zero, statistics might be wrong.")

    if folder:
        if not os.path.isdir(folder):
            os.makedirs(folder)
        vrt_fn = os.path.join(folder, f"{region_code}_{variable}_{overview}.vrt")
        warp_fn = os.path.join(folder, f"{region_code}_{variable}_{overview}.tif")
    else:
        vrt_fn = f"/vsimem/{region_code}_{variable}_{overview}.vrt"
        warp_fn = f"/vsimem/{region_code}_{variable}_{overview}.tif"

    ## Build VRT with all the required data.
    vrt_options = gdal.BuildVRTOptions(
        separate=True,
    )
    vrt = gdal.BuildVRT(vrt_fn, ["/vsicurl/" + x[1] for x in date_urls], options = vrt_options)
    vrt.FlushCache()

    if isinstance(region, list):
        region_option = {"outputBounds": region,
                         "outputBoundsSRS": "epsg:4326",
                         }
    elif isinstance(region, str):
        region_option = {"cutlineDSName": region}
    else:
        region_option = {}

    ## Download the data.
    warp_options = gdal.WarpOptions(
        cropToCutline = True,
        overviewLevel = overview,
        multithread = True,
        targetAlignedPixels = True,
        xRes = abs(xres) * 2**(overview_ + 1),
        yRes = abs(yres) * 2**(overview_ + 1),
        creationOptions = ["COMPRESS=LZW"],
        **region_option
    )
    warp = gdal.Warp(warp_fn, vrt_fn, options = warp_options)
    warp.FlushCache()

    ## Collect the stats into a pd.Dataframe if necessary.
    if not isinstance(req_stats, type(None)):
        stats = gdal.Info(warp_fn, format = "json", stats = True)
        data = {statistic: [x["stats"][statistic] for x in stats["stac"]["raster:bands"]] for statistic in req_stats}
        data = pd.DataFrame(data) * scale
        data["date"] = [x[0] for x in date_urls]
    else:
        data = warp_fn

    ## Unlink memory files.
    if "/vsimem/" in vrt_fn:
        _ = gdal.Unlink(vrt_fn)
    if "/vsimem/" in warp_fn:
        _ = gdal.Unlink(warp_fn)

    return data

def wapor_map(region, variable, period, folder, l3_region = None, overview = "NONE"):

    ## Check if raw-data will be downloaded.
    if overview != "NONE":
        logging.warning("Downloading an overview instead of original data.")

    ## Check if a valid path to download into has been defined.
    if not os.path.isdir(folder):
        os.makedirs(folder)

    ## Call wapor_dl to create a GeoTIFF.
    fp = wapor_dl(region, variable,
                  l3_region = l3_region,
                  folder = folder, 
                  period = period,
                  overview = overview,
                  req_stats = None,
                  )
    return fp

def wapor_ts(region, variable, period, overview, l3_region = None,
             req_stats = ["minimum", "maximum", "mean"]):
    
    ## Check if valid statistics have been selected.
    if isinstance(req_stats, type(None)):
        raise ValueError("Please specify a list of required statistics.")
    valid_stats = np.isin(req_stats, ["minimum", "maximum", "mean"])
    req_stats = np.array(req_stats)[valid_stats].tolist()
    if len(req_stats) == 0:
        raise ValueError(f"Please select at least one valid statistic from {valid_stats}.")
    if False in valid_stats:
        print(f"Invalid statistics detected, continuing with `{', '.join(req_stats)}`.")

    ## Call wapor_dl to create a timeseries.
    df = wapor_dl(
            region, variable, 
            l3_region = l3_region,
            period = period, 
            overview = overview, 
            req_stats = req_stats,
            folder = None,
             )
    
    return df

def l3_codes(variable):
    public_urls = generate_urls_v3(variable)
    valids = np.unique([os.path.split(x)[-1].split(".")[2] for x in public_urls])
    return valids.tolist()

if __name__ == "__main__":

    import glob

    regions = glob.glob(r"/Users/hmcoerver/Local/wapor_validation/SHAPES/BASINS/*.geojson")
    region = regions[0]

    bb = [25, -17, 26, -16]

    overview = 3

    period = ["2021-01-01", "2021-03-01"]
    req_stats = ["minimum", "maximum", "mean"]
    variable = "L3-T-D"

    folder = r"/Users/hmcoerver/Local/test"

    l3_region = "BKA"
    region = [35.75,33.70,35.82,33.75]

    # df1 = wapor_ts(region, "L2-AETI-D", period, overview)

    # out = wapor_dl(bb, variable, l3_region, period = period, folder = folder)

