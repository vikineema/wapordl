import os
import requests
import logging
import shapely
import numpy as np
import pandas as pd
from tqdm import tqdm
from osgeo import gdal
from osgeo_utils import gdal_calc
from osgeo import gdalconst
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

L2_BB = """
{
"type": "FeatureCollection",
"name": "L2_BB",
"crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:OGC:1.3:CRS84" } },
"features": [
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 78.606357023900003, 9.4349173999 ], [ 79.791447473800005, 10.843036789099999 ], [ 81.532302093699997, 10.0567161562 ], [ 82.443434906199997, 5.7509010313 ], [ 79.385666656300003, 5.3514869688 ], [ 78.606357023900003, 9.4349173999 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ -18.988762401500001, 28.401098830700001 ], [ -18.111502625699998, 29.879377722600001 ], [ -14.972341179, 29.111181072800001 ], [ -13.5009765625, 30.4189453125 ], [ -12.115385095, 29.2105313884 ], [ -6.2942422389, 36.778086065899998 ], [ -2.0569908152, 36.191945442600002 ], [ 8.9326171875, 38.591796875 ], [ 11.915295971200001, 37.680847568700003 ], [ 12.123803220699999, 34.088306982699997 ], [ 18.887130778300001, 31.573339802900001 ], [ 21.70703125, 33.9912109375 ], [ 33.623480691499999, 32.209304431200003 ], [ 36.036889218100001, 37.743060045500002 ], [ 43.200531648800002, 38.467604023200003 ], [ 44.6142578125, 40.8779296875 ], [ 48.342594953499997, 40.749464548100001 ], [ 51.677209209600001, 38.2333984375 ], [ 57.2255859375, 39.3798828125 ], [ 61.646432345500003, 37.621966253499998 ], [ 62.848104302099998, 30.483598847700001 ], [ 68.276899746599995, 32.879153056200003 ], [ 71.915217546500003, 37.752394235600001 ], [ 75.641735145300004, 37.980093310199997 ], [ 78.8806588315, 35.588582621 ], [ 75.510902048600002, 33.486328125 ], [ 75.621966253500005, 30.584036404500001 ], [ 71.2128556694, 26.7056396043 ], [ 71.74453125, 23.559375 ], [ 68.123046875, 22.64453125 ], [ 65.723022931, 24.3057800218 ], [ 59.372680356399997, 24.366895656899999 ], [ 60.862370486300001, 22.1314994502 ], [ 58.581130218699997, 18.2548072813 ], [ 45.678914320799997, 11.911479749 ], [ 53.5263671875, 13.7646484375 ], [ 55.576947894, 12.637410746 ], [ 52.237831130700002, 11.0995614591 ], [ 48.783104245700002, 3.8200205213 ], [ 40.864452432900002, -4.0321830045 ], [ 41.822601286599998, -15.183992647 ], [ 35.994355375600001, -20.3394637892 ], [ 36.497511423100001, -24.435368391 ], [ 34.020526029899997, -26.139145276 ], [ 33.222311489100001, -29.182952948299999 ], [ 27.520455783, -34.490130316 ], [ 20.130870672099999, -35.877043501700001 ], [ 17.441510309, -34.747018283 ], [ 10.837449514399999, -18.424626203500001 ], [ 12.7608733379, -11.181594454700001 ], [ 11.099264440800001, -5.565783305 ], [ 7.7085432644, -1.021305891 ], [ 8.245721747599999, 0.905887832 ], [ 6.9114297389, -0.9548438784 ], [ 5.5171370144, -0.314274641 ], [ 7.3658002468, 3.3406550233 ], [ 3.849663333, 5.3854824289 ], [ -7.9152135632, 3.3315597774 ], [ -13.8096458437, 6.8866432188 ], [ -17.376047613200001, 10.7335511196 ], [ -18.577719391799999, 14.676933354399999 ], [ -17.08984375, 17.852845460499999 ], [ -17.980250570700001, 22.0637266942 ], [ -15.377485955, 26.7154968436 ], [ -18.455026095499999, 26.706158746500002 ], [ -18.988762401500001, 28.401098830700001 ] ], [ [ 44.629896785500001, 11.4972001985 ], [ 44.663639299899998, 11.474891229200001 ], [ 45.135579571100003, 11.7055316784 ], [ 44.5260602432, 11.6065297549 ], [ 44.629896785500001, 11.4972001985 ] ], [ [ 40.431145729500003, 17.560889204399999 ], [ 40.010343883600001, 18.758173957 ], [ 38.301348287499998, 20.544170626300001 ], [ 38.3173828125, 19.795202525899999 ], [ 40.431145729500003, 17.560889204399999 ] ], [ [ 37.946878988500004, 22.382604237500001 ], [ 37.570236265200002, 23.127954984700001 ], [ 37.207601586499997, 23.262849484299998 ], [ 37.918475282499998, 22.289412114 ], [ 37.946878988500004, 22.382604237500001 ] ], [ [ 8.369399142600001, 2.206189029 ], [ 8.4789315626, 1.941754375 ], [ 8.5003244219, 1.8222512998 ], [ 8.680302471899999, 2.158203125 ], [ 8.369399142600001, 2.206189029 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 54.246341378799997, -20.6914705895 ], [ 57.494246281, -18.938556559799999 ], [ 58.837177414099997, -20.015080286300002 ], [ 55.854863018700001, -22.4009517595 ], [ 54.246341378799997, -20.6914705895 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 55.547732956399997, -5.8544921875 ], [ 54.165323943799997, -4.3351411438 ], [ 55.814533088600001, -3.2443975153 ], [ 56.9052734375, -4.37890625 ], [ 55.547732956399997, -5.8544921875 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 42.1806640625, -22.156880701399999 ], [ 43.539959982900001, -15.6689001241 ], [ 49.877942383700002, -11.0966606432 ], [ 51.502644806200003, -15.6004057312 ], [ 47.895154819699997, -25.679988982899999 ], [ 43.773424803799998, -26.1425971693 ], [ 42.1806640625, -22.156880701399999 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 42.1796875, -11.731138914500001 ], [ 43.752490629100002, -10.386913009900001 ], [ 46.275806562600003, -13.043316875 ], [ 44.753368745899998, -13.970508865099999 ], [ 42.1796875, -11.731138914500001 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ -18.0589314855, 33.465806658299996 ], [ -16.139828643800001, 34.127644806200003 ], [ -15.3116387962, 32.731022068100003 ], [ -17.288508030799999, 31.647540866300002 ], [ -18.0589314855, 33.465806658299996 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ -26.381789271700001, 16.800942259399999 ], [ -25.370465312499999, 18.2142831251 ], [ -22.619878676900001, 17.861146597699999 ], [ -21.6431387499, 15.8209409375 ], [ -24.924631760099999, 13.773325574299999 ], [ -26.381789271700001, 16.800942259399999 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ -80.053586131, 1.8210000352 ], [ -78.253525059400005, 9.285827860099999 ], [ -75.354440158, 12.0243100035 ], [ -70.584289327400001, 13.147079746899999 ], [ -71.2841796875, 8.0822342226 ], [ -66.700119781300003, 6.9336692812 ], [ -65.803995818800004, 1.4978666687 ], [ -68.5805109251, -0.2643027361 ], [ -69.126309970400001, -4.8989734622 ], [ -80.053586131, 1.8210000352 ] ] ] } },
{ "type": "Feature", "properties": { }, "geometry": { "type": "Polygon", "coordinates": [ [ [ -82.784960882099995, 12.512207031299999 ], [ -82.005273437499994, 14.236523437500001 ], [ -80.474590330699996, 13.952161133700001 ], [ -81.09921875, 11.638476562499999 ], [ -82.784960882099995, 12.512207031299999 ] ] ] } }
]
}
"""

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

def __make_band_names__(length):
    letters = [x for x in ascii_lowercase + ascii_uppercase]
    i = 2
    while len(letters) < length:
        for letter in letters[:52]:
            letters.append(letter*i)
        i += 1
    return letters[:length]

def unit_convertor(urls, out_fn, unit_conversion, warp):

    input_files = dict()
    input_bands = dict()
    calc = list()
    should_convert = list()
    letters = __make_band_names__(len(urls))

    for i, (md, _) in enumerate(urls):
        band_number = i+1
        letter = letters[i]
        input_files[letter] = out_fn
        input_bands[f"{letter}_band"] = band_number
        number_of_days = md.get("number_of_days", "unknown")
        days_in_month = pd.Timestamp(md.get("start_date", "nat")).daysinmonth
        source_unit = md.get("units", "unknown")
        source_unit_split = source_unit.split("/")
        source_unit_q = "/".join(source_unit_split[:-1])
        source_unit_time = source_unit_split[-1]
        if any([
                source_unit_time not in ["day", "month", "year"],
                number_of_days == "unknown",
                source_unit == "unknown",
                pd.isnull(days_in_month)
            ]):
            calc.append(f"{letter}.astype(numpy.float64)")
            md["units"] = source_unit
            md["units_conversion_factor"] = "N/A"
            md["original_units"] = "N/A"
            should_convert.append(False)
        else:
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
            calc.append(f"{letter}.astype(numpy.float64)*{conversion}")
            should_convert.append(True)
            md["units"] = f"{source_unit_q}/{unit_conversion}"
            md["units_conversion_factor"] = conversion
            md["original_units"] = source_unit

    logging.debug(f"\ninput_files: {input_files}\ninput_bands: {input_bands}\ncalc: {calc}")

    conversion_is_one = [x["units_conversion_factor"] == 1.0 for x, _ in urls]

    # NOTE See todo just below.
    scales = [warp.GetRasterBand(i+1).GetScale() for i in range(warp.RasterCount)]
    offsets = [warp.GetRasterBand(i+1).GetOffset() for i in range(warp.RasterCount)]

    logging.debug(f"\nSCALES: {scales}\nOFFSETS: {offsets}")

    if all(should_convert) and not all(conversion_is_one):
        logging.info(f"Converting units from [{source_unit}] to [{source_unit_q}/{unit_conversion}].")
        ndv = warp.GetRasterBand(1).GetNoDataValue()
        warp = gdal_calc.Calc(
            calc = calc,
            outfile = out_fn,
            overwrite = True,
            quiet = True,
            type = gdalconst.GDT_Int32, # NOTE unit conversion can increase the DN's, 
                                        # causing the data to not fit inside Int16 anymore...
                                        # so for now just moving up to Int32. Especially necessary
                                        # for NPP (which has a scale-factor of 0.001).
            NoDataValue = ndv,
            **input_files,
            **input_bands,
            )
        # TODO make bug report on GDAL for gdal_calc removing scale/offset factors
        for i, (scale, offset) in enumerate(zip(scales, offsets)):
            warp.GetRasterBand(i+1).SetScale(scale)
            warp.GetRasterBand(i+1).SetOffset(offset)

        warp.FlushCache()
    else:
        if all(conversion_is_one):
            logging.info(f"Units are already as requested, no conversion needed.")
        else:
            logging.warning(f"Couldn't succesfully determine unit conversion factors, keeping original units.")
        for i, (md, _) in enumerate(urls):
            if md["units_conversion_factor"] != "N/A":
                md["units"] = md["original_units"]
                md["units_conversion_factor"] = f"N/A"
                md["original_units"] = "N/A"

    return warp

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

    n_urls = len(urls)

    # Create waitbar.
    waitbar = tqdm(desc = f"Downloading {n_urls} COGs", leave = False, total = 100, bar_format='{l_bar}{bar}|')
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
    warp.FlushCache() # NOTE do not remove this.
    waitbar.close()
    nbands = warp.RasterCount
    
    if nbands == n_urls and unit_conversion != "none":
        warp = unit_convertor(urls, out_fn, unit_conversion, warp)

    if nbands == n_urls:
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
    ## NOTE This should be corrected in the GISMGR2.0 API instead of here:
    if all([
        "NPP" in var_code, 
        "year" not in md["units"], 
        "month" not in md["units"], 
        "day" not in md["units"]
        ]):
        md["units"] += {"D": "/day", "M": "/month", "A": "/year"}[tres]
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
        if level == "L2":
            data_bb = shapely.from_geojson(L2_BB)
        else:
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
        warp_fn = os.path.join(folder, f"{region_code}_{variable}_{overview}_{unit_conversion}.tif")
    else:
        warp_fn = f"/vsimem/{pd.Timestamp.now()}_{region_code}_{variable}_{overview}_{unit_conversion}.tif"

    warp_fn = cog_dl(md_urls, warp_fn, overview = overview_, warp_kwargs = warp_kwargs, unit_conversion = unit_conversion)

    ## Collect the stats into a pd.Dataframe if necessary.
    if not isinstance(req_stats, type(None)):
        stats = gdal.Info(warp_fn, format = "json", stats = True)
        data = {statistic: [x.get(statistic, np.nan) for x in stats["bands"]] for statistic in req_stats}
        data = pd.DataFrame(data) * scale
        data["start_date"] = [pd.Timestamp(x.get("metadata", {}).get("", {}).get("start_date", "nat")) for x in stats["bands"]]
        data["end_date"] = [pd.Timestamp(x.get("metadata", {}).get("", {}).get("end_date", "nat")) for x in stats["bands"]]
        data["number_of_days"] = [pd.Timedelta(float(x.get("metadata", {}).get("", {}).get("number_of_days", np.nan)), "days") for x in stats["bands"]]
        out_md = {k: v for k, v in md_urls[0][0].items() if k in ['long_name', 'units', 'overview', 'original_units']}
        data.attrs = out_md
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

    valid_units = ["none", "dekad", "day", "month", "year"]
    if not unit_conversion in valid_units:
        raise ValueError(f"Please select one of {valid_units} instead of {unit_conversion}.") # NOTE: TESTED

    ## Call wapor_dl to create a GeoTIFF.
    fp = wapor_dl(region, variable,
                  folder = folder, 
                  period = period,
                  overview = overview,
                  unit_conversion = unit_conversion,
                  req_stats = None,
                  )
    
    if extension != ".tif":
        logging.info(f"Converting from `.tif` to `{extension}`.")
        toptions = {".nc": {"creationOptions": ["COMPRESS=DEFLATE", "FORMAT=NC4C"]}}
        options = gdal.TranslateOptions(
            **toptions.get(extension, {})
            )
        new_fp = fp.replace(".tif", extension)
        ds = gdal.Translate(new_fp, fp, options = options)
        ds.FlushCache()
        try:
            os.remove(fp)
        except PermissionError:
            ...
        return new_fp
    else:
        return fp

def wapor_ts(region, variable, period, overview,
             unit_conversion = "none",
             req_stats = ["minimum", "maximum", "mean"]):
    
    valid_units = ["none", "dekad", "day", "month", "year"]
    if not unit_conversion in valid_units:
        raise ValueError(f"Please select one of {valid_units} instead of {unit_conversion}.")  # NOTE: TESTED
    
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
            unit_conversion = unit_conversion,
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

    region = "/Users/hmcoerver/Library/Mobile Documents/com~apple~CloudDocs/GitHub/wapordl/wapordl/test_data/1237500.geojson"
    variable = "L1-AETI-D"
    period = ["2021-01-12", "2021-01-28"]
    overview = 3
    unit_conversion = "none"
    req_stats = ["minimum", "maximum", "mean"]
    extension = ".nc"
    unit_conversion = "dekad"
    folder = r"/Users/hmcoerver/Local/test"

    map1 = wapor_map(region, variable, period, folder, unit_conversion=unit_conversion)
    map2 = wapor_map(region, variable, period, folder, unit_conversion=unit_conversion, extension=extension)

