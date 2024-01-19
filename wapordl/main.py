import os
import requests
import logging
import shapely
import numpy as np
import pandas as pd
from tqdm import tqdm
from osgeo import gdal
from osgeo_utils import gdal_calc
from string import ascii_lowercase, ascii_uppercase
gdal.UseExceptions()
logging.basicConfig(encoding='utf-8', level=logging.INFO, format='%(levelname)s: %(message)s')

L3_BBS = {
    'AWA': [[39.1751869, 8.9148245], [39.1749088, 8.3098793], [40.0254231, 8.3085969], [40.0270531, 8.9134473], [39.1751869, 8.9148245]], 
    'BKA': [[35.7339813, 34.0450172], [35.7204902, 33.6205171], [36.2189392, 33.6085397], [36.2348925, 34.0328476], [35.7339813, 34.0450172]], 
    'BUS': [[34.0888682, 0.7782889], [34.0887832, 0.3004971], [34.4348724, 0.3004569], [34.4349843, 0.7781846], [34.0888682, 0.7782889]], 
    'ERB': [[43.4914673, 36.4664115], [43.5044881, 35.7830724], [44.1182727, 35.7891403], [44.1105935, 36.4726323], [43.4914673, 36.4664115]], 
    'GAR': [[45.76927, 32.4318047], [45.7627701, 31.6513125], [46.1095137, 31.6487693], [46.1189671, 32.4291837], [45.76927, 32.4318047]], 
    'GEZ': [[33.1560646, 14.4532485], [33.1558745, 14.1778375], [33.5854923, 14.1771733], [33.5862063, 14.4525708], [33.1560646, 14.4532485]], 
    'JAF': [[35.4431064, 30.5318429], [35.4299372, 29.9988137], [36.2127498, 29.9820272], [36.2301458, 30.5146957], [35.4431064, 30.5318429]], 
    'JEN': [[8.4492905, 36.6560571], [8.4511646, 36.3912042], [9.2256932, 36.3922519], [9.2264639, 36.657115], [8.4492905, 36.6560571]], 
    'JVA': [[35.5548028, 32.6838097], [35.5452241, 32.3447739], [35.6558413, 32.3424921], [35.6658352, 32.6814981], [35.5548028, 32.6838097]], 
    'KAI': [[9.6792133, 35.7356369], [9.6772114, 35.4985191], [10.0485184, 35.4958636], [10.0516175, 35.7329582], [9.6792133, 35.7356369]], 
    'KOG': [[36.9808785, 11.5497162], [36.9826246, 11.302976], [37.2648864, 11.3047648], [37.2633842, 11.551545], [36.9808785, 11.5497162]], 
    'LAM': [[34.0513688, -19.1827337], [34.0530979, -19.4529053], [34.4506668, -19.4501609], [34.4482854, -19.1800303], [34.0513688, -19.1827337]], 
    'LCE': [[14.2008944, 26.7493197], [14.2026255, 26.5001416], [14.5460084, 26.5016471], [14.5450227, 26.7508416], [14.2008944, 26.7493197]], 
    'LDA': [[16.1209927, 29.178662], [16.1187095, 28.9674972], [16.4144337, 28.9647125], [16.41732, 29.1758532], [16.1209927, 29.178662]], 
    'MAL': [[80.1829979, 8.6817104], [80.1841408, 8.1356067], [80.5350546, 8.1361599], [80.5344033, 8.6823012], [80.1829979, 8.6817104]], 
    'MIT': [[2.3861237, 36.804671], [2.3893986, 36.3901832], [3.403877, 36.3910611], [3.4060432, 36.8055622], [2.3861237, 36.804671]], 
    'MUV': [[30.2315244, -1.0474462], [30.2323749, -1.6835342], [30.474523, -1.6831148], [30.4736091, -1.0471853], [30.2315244, -1.0474462]], 
    'ODN': [[-6.2554174, 14.5579879], [-6.2650165, 13.7577987], [-5.8657445, 13.7530413], [-5.8547493, 14.5529423], [-6.2554174, 14.5579879]], 
    'PAL': [[35.3712698, 32.1997283], [35.359621, 31.7457763], [35.5587674, 31.7419309], [35.5713968, 32.1958148], [35.3712698, 32.1997283]], 
    'SAN': [[43.9733487, 15.607452], [43.9752772, 15.2142363], [44.4893856, 15.2159917], [44.4884245, 15.609255], [43.9733487, 15.607452]], 
    'SED': [[-16.2610895, 16.4920046], [-16.2593842, 16.2262901], [-15.8685096, 16.22825], [-15.8696858, 16.4939983], [-16.2610895, 16.4920046]], 
    'YAN': [[29.9186281, -1.7554126], [29.9189336, -1.9422648], [30.0389016, -1.9420519], [30.0385836, -1.7552203], [29.9186281, -1.7554126]],
    }

def guess_l3_region(region_shape):

    checks = {x: shapely.Polygon(np.array(bb)).intersects(region_shape) for x, bb in L3_BBS.items()}
    number_of_results = sum(checks.values())
    if number_of_results == 0:
        raise ValueError(f"`region` can't be linked to any L3 region.") # NOTE: TESTED
    
    l3_regions = [k for k, v in checks.items() if v]
    l3_region = l3_regions[0]
    if number_of_results > 1:
        logging.warning(f"`region` intersects with multiple L3 regions ({l3_regions}), continuing with {l3_region} only.")
    else:
        logging.info(f"Given `region` matches with `{l3_region}` L3 region.")
    
    return l3_region      

def collect_responses(url, info = ["code"]):
    data = {"links": [{"rel": "next", "href": url}]}
    output = list()
    while "next" in [x["rel"] for x in data["links"]]:
        url_ = [x["href"] for x in data["links"] if x["rel"] == "next"][0]
        response = requests.get(url_)
        response.raise_for_status()
        data = response.json()["response"]
        if isinstance(info, list) and "items" in data.keys():
            output += [tuple(x.get(y) for y in info) for x in data["items"]]
        elif "items" in data.keys():
            output += data["items"]
        else:
            output.append(data)
    if isinstance(info, list):
        output = sorted(output)
    return output

def date_func(url, tres):
    if tres == "D":
        year, month, dekad = os.path.split(url)[-1].split(".")[-2].split("-")
        start_day = {'D1': '01', 'D2': '11', 'D3': '21'}[dekad]
        start_date = f"{year}-{month}-{start_day}"
        end_day = {'D1': '10', 'D2': '20', 'D3': pd.Timestamp(start_date).daysinmonth}[dekad]
        end_date = f"{year}-{month}-{end_day}"
    elif tres == "M":
        year, month = os.path.split(url)[-1].split(".")[-2].split("-")
        start_date = f"{year}-{month}-01"
        end_date = f"{year}-{month}-{pd.Timestamp(start_date).days_in_month}"
    elif tres == "A":
        year = os.path.split(url)[-1].split(".")[-2]
        start_date = f"{year}-01-01"
        end_date = f"{year}-12-31"
    else:
        raise ValueError("Invalid temporal resolution.") # NOTE: TESTED
    number_of_days = (pd.Timestamp(end_date) - pd.Timestamp(start_date) + pd.Timedelta(1, "D")).days
    return {"start_date": start_date, "end_date": end_date, "number_of_days": number_of_days}

def collect_metadata(variable):
    if "L1" in variable:
        base_url = f"https://data.apps.fao.org/gismgr/api/v2/catalog/workspaces/WAPOR-3/mapsets"
    elif "L2" in variable:
        base_url = f"https://data.apps.fao.org/gismgr/api/v2/catalog/workspaces/WAPOR-3/mapsets"
    elif "L3" in variable:
        base_url = f"https://data.apps.fao.org/gismgr/api/v2/catalog/workspaces/WAPOR-3/mosaicsets"
    else:
        raise ValueError(f"Invalid variable name {variable}.") # NOTE: TESTED
    info = ["code", "measureCaption", "measureUnit"]
    var_codes = {x[0]: {"long_name": x[1], "units": x[2]} for x in collect_responses(base_url, info = info)}
    return var_codes[variable]

def generate_urls_v3(variable, l3_region = None, period = None):
    
    level, _, tres = variable.split("-")

    if (level == "L1") or (level == "L2"):
        base_url = f"https://data.apps.fao.org/gismgr/api/v2/catalog/workspaces/WAPOR-3/mapsets"
    elif level == "L3":
        base_url = f"https://data.apps.fao.org/gismgr/api/v2/catalog/workspaces/WAPOR-3/mosaicsets"
    else:
        raise ValueError(f"Invalid level {level}.") # NOTE: TESTED
    
    mapset_url = f"{base_url}/{variable}/rasters?filter="
    if not isinstance(l3_region, type(None)):
        mapset_url += f"code:CONTAINS:{l3_region};"
    if not isinstance(period, type(None)):
        tres_translator = {"D": "dekad", "M": "month", "A": "year"}
        mapset_url += f"{tres_translator[tres]}:OVERLAPS:{period[0]}:{period[1]};"
    
    urls = [x[0] for x in collect_responses(mapset_url, info = ["downloadUrl"])]
    
    return tuple(sorted(urls))

def cog_dl(urls, out_fn, overview = "NONE", warp_kwargs = {}, vrt_options = {"separate": True}, unit_conversion = "none"):

    out_ext = os.path.splitext(out_fn)[-1]
    valid_ext = {".nc": "netCDF", ".tif": "GTiff"}
    valid_cos = {".nc": ["COMPRESS=DEFLATE", "FORMAT=NC4C"], ".tif": ["COMPRESS=LZW"]}
    if not bool(np.isin(out_ext, list(valid_ext.keys()))):
        raise ValueError(f"Please use one of {list(valid_ext.keys())} as extension for `out_fn`, not {out_ext}") # NOTE: TESTED
    vrt_fn = out_fn.replace(out_ext, ".vrt")

    ## Build VRT with all the required data.
    vrt_options_ = gdal.BuildVRTOptions(
        **vrt_options
    )
    vrt = gdal.BuildVRT(vrt_fn, ["/vsicurl/" + x[1] for x in urls], options = vrt_options_)
    vrt.FlushCache()

    # Create waitbar.
    waitbar = tqdm(desc = f"Downloading {len(urls)} COGs", leave = False, total = 100, bar_format='{l_bar}{bar}|')
    # Define callback function for waitbar progress.
    def _callback_func(info, *args):
        waitbar.update(info * 100 - waitbar.n)

    ## Download the data.
    warp_options = gdal.WarpOptions(
        format = valid_ext[out_ext],
        cropToCutline = True,
        overviewLevel = overview,
        multithread = True,
        targetAlignedPixels = True,
        creationOptions = valid_cos[out_ext],
        callback = _callback_func,
        **warp_kwargs,
    )
    warp = gdal.Warp(out_fn, vrt_fn, options = warp_options)
    waitbar.close()
    nbands = warp.RasterCount
    
    if nbands == len(urls) and unit_conversion != "none":
        input_files = dict()
        input_bands = dict()
        calc = list()
        for i, (md, _) in enumerate(urls):
            band_number = i+1
            letters = ascii_lowercase[i] # TODO make sure letters is longer than RasterCount
            input_files[letters] = out_fn
            input_bands[f"{letters}_band"] = band_number
            number_of_days = md.get("number_of_days", "unknown") # TODO handle "unknown"
            days_in_month = pd.Timestamp(md.get("start_date", "unknown")).daysinmonth # TODO handle "unknown"
            source_unit = md.get("units", "unknown") # TODO handle "unknown"
            source_unit_q, source_unit_time = source_unit.split("/") # TODO what if there is no "/" in the string
            conversion = {
                ("day", "day"): 1,
                ("day", "dekad"): number_of_days,
                ("day", "month"): days_in_month,
                ("day", "year"): 365,
                ("month", "day"): 1/days_in_month,
                ("month", "dekad"): 1/3,
                ("month", "month"): 1,
                ("month", "year"): 12,
                ("year", "dekad"): 1/36,
                ("year", "day"): 1/365,
                ("year", "month"): 1/12,
                ("year", "year"): 1,
            }[(source_unit_time, unit_conversion)]
            calc.append(f"{letters}*{conversion}")
            md["units"] = f"{source_unit_q}/{unit_conversion}"
            md["unit_conversion"] = f"[{source_unit}] --*{conversion:.3f}--> [{source_unit_q}/{unit_conversion}]"

        logging.debug(f"input_files: {input_files}\ninput_bands: {input_bands}\ncalc: {calc}")
        logging.info(f"Convertin units: [{source_unit}] --> [{source_unit_q}/{unit_conversion}]")

        warp.FlushCache()
        warp = gdal_calc.Calc(
            calc = calc,
            outfile = out_fn,
            overwrite = True,
            quiet = True,
            **input_files,
            **input_bands,
            )

    if nbands == len(urls):
        for i, (md, _) in enumerate(urls):
            if not isinstance(md, type(None)):
                band = warp.GetRasterBand(i + 1)
                band.SetMetadata(md)

    warp.FlushCache()

    if os.path.isfile(vrt_fn):
        try:
            os.remove(vrt_fn)
        except PermissionError:
            ...

    return out_fn

def wapor_dl(region, variable,
             period = ["2021-01-01", "2022-01-01"], 
             overview = "NONE",
             unit_conversion = "none", 
             req_stats = ["minimum", "maximum", "mean"],
             extension = ".tif",
             folder = None):
    """_summary_

    Parameters
    ----------
    region : str, list, None
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

    ## Check if region is valid.
    if all([isinstance(region, str), 
            # not os.path.isfile(region),
            region in list(L3_BBS.keys()), 
            len(region) == 3]):
        l3_region = region[:]
        region = None
        region_code = l3_region[:]
        region_shape = None
    elif isinstance(region, str):
        if not os.path.isfile(region) or os.path.splitext(region)[-1] != ".geojson":
            raise ValueError(f"Geojson file not found.") # NOTE: TESTED
        else:
            region_code = os.path.split(region)[-1].replace(".geojson", "")
            region_shape = shapely.from_geojson(open(region,'r').read())
        l3_region = None
    elif isinstance(region, list):
        if not all([region[2] > region[0], region[3] > region[1]]):
            raise ValueError(f"Invalid bounding box.") # NOTE: TESTED
        else:
            region_code = "bb"
            region_shape = shapely.Polygon([(region[0], region[1]), 
                                            (region[2], region[1]), 
                                            (region[2], region[3]), 
                                            (region[0], region[3]), 
                                            (region[0], region[1])])
        l3_region = None
    else:
        raise ValueError(f"Invalid value for region ({region}).") # NOTE: TESTED

    ## Check l3_region code.
    if level == "L3" and isinstance(l3_region, type(None)):
        l3_region = guess_l3_region(region_shape)
        region_code += f".{l3_region}"

    ## Check the dates in period.
    if not isinstance(period, type(None)):
        period = [pd.Timestamp(x) for x in period]
        if period[0] > period[1]:
            raise ValueError(f"Invalid period.") # NOTE: TESTED
        period = [x.strftime("%Y-%m-%d") for x in period]

    ## Collect urls for requested variable.
    urls = generate_urls_v3(variable, l3_region = l3_region, period = period)

    if len(urls) == 0:
        raise ValueError("No files found for selected region, variable and period.")  # NOTE: TESTED

    ## Determine date for each url.
    md = collect_metadata(variable)
    md["overview"] = overview
    md_urls = [({**date_func(url, tres), **md}, url) for url in urls]

    logging.info(f"Found {len(md_urls)} files for {variable}.")

    ## Determine required output resolution.
    # NOTE maybe move this to external function (assumes info the same for all urls)
    info = gdal.Info("/vsicurl/" + md_urls[0][1], format = "json")
    overview_ = -1 if overview == "NONE" else overview
    xres, yres = info["geoTransform"][1::4]
    warp_kwargs = {
        "xRes": abs(xres) * 2**(overview_ + 1),
        "yRes": abs(yres) * 2**(overview_ + 1),
    }

    if isinstance(region, list):
        warp_kwargs["outputBounds"] = region
        warp_kwargs["outputBoundsSRS"] = "epsg:4326"
    elif isinstance(region, str):
        warp_kwargs["cutlineDSName"] = region
    else:
        ...

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
        warp_fn = os.path.join(folder, f"{region_code}_{variable}_{overview}{extension}")
    else:
        warp_fn = f"/vsimem/{region_code}_{variable}_{overview}{extension}"

    warp_fn = cog_dl(md_urls, warp_fn, overview = overview_, warp_kwargs = warp_kwargs)

    ## Collect the stats into a pd.Dataframe if necessary.
    if not isinstance(req_stats, type(None)):
        stats = gdal.Info(warp_fn, format = "json", stats = True)
        data = {statistic: [x["stats"][statistic] for x in stats["stac"]["raster:bands"]] for statistic in req_stats}
        data = pd.DataFrame(data) * scale
        data["start_date"] = [pd.Timestamp(x[0]["start_date"]) for x in md_urls]
        data["end_date"] = [pd.Timestamp(x[0]["end_date"]) for x in md_urls]
        data["number_of_days"] = [pd.Timedelta(x[0]["number_of_days"], "days") for x in md_urls]
        data.attrs = md
    else:
        data = warp_fn

    ## Unlink memory files.
    if "/vsimem/" in warp_fn.replace(".tif", ".vrt"):
        _ = gdal.Unlink(warp_fn.replace(".tif", ".vrt"))
    if "/vsimem/" in warp_fn:
        _ = gdal.Unlink(warp_fn)

    return data

def wapor_map(region, variable, period, folder, 
              unit_conversion = "none",
              overview = "NONE", extension = ".tif"):

    ## Check if raw-data will be downloaded.
    if overview != "NONE":
        logging.warning("Downloading an overview instead of original data.")

    ## Check if a valid path to download into has been defined.
    if not os.path.isdir(folder):
        os.makedirs(folder)

    if not unit_conversion in ["none", "match", "day", "month", "year"]:
        raise ValueError

    ## Call wapor_dl to create a GeoTIFF.
    fp = wapor_dl(region, variable,
                  folder = folder, 
                  period = period,
                  overview = overview,
                  extension = extension,
                  unit_conversion = unit_conversion,
                  req_stats = None,
                  )
    return fp

def wapor_ts(region, variable, period, overview,
             req_stats = ["minimum", "maximum", "mean"]):
    
    ## Check if valid statistics have been selected.
    if not isinstance(req_stats, list):
        raise ValueError("Please specify a list of required statistics.") # NOTE: TESTED
    valid_stats = np.isin(req_stats, ["minimum", "maximum", "mean"])
    req_stats = np.array(req_stats)[valid_stats].tolist()
    if len(req_stats) == 0:
        raise ValueError(f"Please select at least one valid statistic from {valid_stats}.") # NOTE: TESTED
    if False in valid_stats:
        logging.warning(f"Invalid statistics detected, continuing with `{', '.join(req_stats)}`.")

    ## Call wapor_dl to create a timeseries.
    df = wapor_dl(
            region, variable, 
            period = period, 
            overview = overview, 
            req_stats = req_stats,
            folder = None,
             )
    
    return df

def __l3_codes__(variable = "L3-T-A"):
    public_urls = generate_urls_v3(variable, period = ["2019-01-01", "2019-02-01"])
    valids = np.unique([os.path.split(x)[-1].split(".")[2] for x in public_urls])
    return valids.tolist()

def __l3_bounding_boxes__(variable = "L3-T-A"):
    urls = generate_urls_v3(variable, period = ["2019-01-01", "2019-02-01"])
    l3_bbs = {}
    for region_code, url in zip([os.path.split(x)[-1].split(".")[-3] for x in urls], urls):
        info = gdal.Info("/vsicurl/" + url, format = "json")
        bb = info["wgs84Extent"]["coordinates"][0]
        l3_bbs[region_code] = bb
    return l3_bbs

if __name__ == "__main__":

    import glob

    regions = glob.glob(r"/Users/hmcoerver/Local/wapor_validation/SHAPES/BASINS/*.geojson")
    region = regions[0]

    bb = [25, -17, 26, -16]

    overview = 3

    period = ["2021-01-01", "2021-03-01"]
    req_stats = ["minimum", "maximum", "mean"]
    variable = "L2-T-D"
    extension = ".tif"

    folder = r"/Users/hmcoerver/Local/test"
