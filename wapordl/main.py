import os
import requests
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
             period = ["2021-01-01", "2022-01-01"], 
             overview = "auto", 
             req_stats = ["minimum", "maximum", "mean"],
             folder = None):

    ## Check if region is valid.
    if not os.path.isfile(region) or os.path.splitext(region)[-1] != ".geojson":
        raise ValueError
    else:
        region_code = os.path.split(region)[-1].replace(".geojson", "")

    ## Collect urls for requested variable.
    public_urls = generate_urls_v3(variable)

    ## Parse the dates in period.
    if not isinstance(period, type(None)):
        period = [pd.Timestamp(x) for x in period]
        if period[0] > period[1]:
            raise ValueError

    ## Retrieve info from variable name.
    level, var_code, tres = variable.split("-")

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
        raise ValueError

    ## Determine date for each url and filter on requested period.
    date_urls = [(date_func(url), url) for url in public_urls]
    if not isinstance(period, type(None)):
        date_urls = [x for x in date_urls if (x[0] >= period[0]) & (x[0] <= period[1])]

    ## Print overview statement.
    print(f"Found {len(date_urls)} files for {variable} between {period[0]} and {period[1]}.")

    ## Determine required output resolution.
    # NOTE maybe move this to external function (assumes info the same for all urls)
    info = gdal.Info("/vsicurl/" + date_urls[0][1], format = "json")
    overview_ = -1 if overview == "NONE" else overview
    xres, yres = info["geoTransform"][1::4]

    ## Get scale and offset factor.
    scale = info["bands"][0]["scale"]
    offset = info["bands"][0]["offset"]

    ## Check offset factor.
    if offset != 0:
        print("WARNING: offset factor is not zero, statistics might be wrong.")

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

    ## Download the data.
    warp_options = gdal.WarpOptions(
        cutlineDSName = region,
        cropToCutline = True,
        overviewLevel = overview,
        multithread = True,
        targetAlignedPixels = True,
        xRes = abs(xres) * 2**(overview_ + 1),
        yRes = abs(yres) * 2**(overview_ + 1),
        creationOptions = ["COMPRESS=LZW"],
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

def wapor_map(region, variable, folder, period, overview = "NONE"):

    ## Check if raw-data will be downloaded.
    if overview != "NONE":
        print("Downloading an overview instead of original data.")

    ## Check if a valid path to download into has been defined.
    if not os.path.isdir(folder):
        os.makedirs(folder)

    ## Call wapor_dl to create a GeoTIFF.
    fp = wapor_dl(region, variable, 
                  folder = folder, 
                  period = period,
                  overview = overview,
                  req_stats = None,
                  )
    return fp

def wapor_ts(region, variable, period, overview, 
             req_stats = ["minimum", "maximum", "mean"]):
    
    ## Check if valid statistics have been selected.
    if isinstance(req_stats, type(None)):
        raise ValueError
    valid_stats = np.isin(req_stats, ["minimum", "maximum", "mean"])
    req_stats = np.array(req_stats)[valid_stats].tolist()
    if len(req_stats) == 0:
        raise ValueError
    if False in valid_stats:
        print(f"Invalid statistics detected, continuing with `{', '.join(req_stats)}`.")

    ## Call wapor_dl to create a timeseries.
    df = wapor_dl(
            region, variable, 
            period = period, 
            overview = overview, 
            req_stats = req_stats,
            folder = None,
             )
    
    return df

if __name__ == "__main__":

    import glob

    regions = glob.glob(r"/Users/hmcoerver/Local/wapor_validation/SHAPES/BASINS/*.geojson")
    region = regions[0]

    period1 = ["2021-01-01", "2021-07-01"]
    period2 = ["2021-01-01", "2022-02-01"]
    req_stats = ["minimum", "maximum", "mean"]
    variable = "L2-AETI-D"

    df1 = wapor_ts(region, "L2-AETI-D", period = period1, overview = 6)
    df2 = wapor_ts(region, "L2-AETI-M", period = period2, overview = 5)
    df3 = wapor_ts(region, "L2-AETI-A", period = period2, overview = 5)

    fp = wapor_map(region, "L2-AETI-D", period = period1, folder = r"/Users/hmcoerver/Local/test")