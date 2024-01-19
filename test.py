import wapordl
module_path = wapordl.__path__[0]
assert "conda" not in module_path
from wapordl import wapor_ts, wapor_map, wapor_dl, date_func, collect_metadata, generate_urls_v3
import pandas as pd
from osgeo import gdal, osr
import numpy as np
import pathlib
import xarray as xr
import matplotlib.pyplot as plt

test_data_folder = pathlib.Path(module_path).parent / "test_data"
region = str(list(test_data_folder.glob("[0-9]*.geojson"))[0])
l3_regions = list(test_data_folder.glob("test_*.geojson"))
folder = r"/Users/hmcoerver/Local/test"

bb = [25, -17, 26, -16]

overview = 3

period = ["2021-01-12", "2021-01-25"]
nodata_period = ["2015-01-01", "2016-01-01"]
req_stats = ["minimum", "maximum", "mean"]
variable = "L2-AETI-D"
l3_region = "BKA"

####
# START TESTS
####

df1 = wapor_ts(region, "L2-AETI-D", period, overview)
assert df1.iloc[0].start_date == pd.Timestamp('2021-01-11 00:00:00')
assert df1.attrs == {
    'long_name': 'Actual EvapoTranspiration and Interception',
    'units': 'mm/day',
    'overview': overview
    }

df2 = wapor_ts(region, "L2-AETI-M", period, overview)
assert df2.attrs["units"] == "mm/month"

df3 = wapor_ts(region, "L2-AETI-A", period, overview)
assert df3.attrs["units"] == "mm/year"

df4 = wapor_ts(region, "L1-AETI-D", period, overview)
df5 = wapor_ts(region, "L1-AETI-M", period, overview)
df6 = wapor_ts(region, "L1-AETI-A", period, overview)

df7 = wapor_ts(bb, "L2-AETI-D", period, overview)
df8 = wapor_ts(bb, "L2-AETI-M", period, overview)
df9 = wapor_ts(bb, "L2-AETI-A", period, overview)

df10 = wapor_ts(bb, "L1-AETI-D", period, overview)
df11 = wapor_ts(bb, "L1-AETI-M", period, overview)
df12 = wapor_ts(bb, "L1-AETI-A", period, overview)

fp1 = wapor_map(region, "L2-AETI-D", period, folder)
ds = gdal.Open(fp1)
assert ds.RasterCount == 2
band = ds.GetRasterBand(1)
ndv = band.GetNoDataValue()
scale = band.GetScale()
array = band.ReadAsArray() * scale
array[array == ndv*scale] = np.nan
mean = np.nanmean(array)
md = band.GetMetadata()
assert md == {'end_date': '2021-01-20',
 'long_name': 'Actual EvapoTranspiration and Interception',
 'number_of_days': '10',
 'overview': 'NONE',
 'start_date': '2021-01-11',
 'units': 'mm/day'}
assert mean > 0.0
assert mean < 15.0
proj = osr.SpatialReference(wkt=ds.GetProjection())
assert proj.GetAttrValue('AUTHORITY',1) == "4326"
ds = ds.FlushCache()

fp2 = wapor_map(region, "L2-AETI-M", period, folder)
ds = gdal.Open(fp2)
band = ds.GetRasterBand(1)
md = band.GetMetadata()
assert md["units"] == "mm/month"
ds = ds.FlushCache()

fp3 = wapor_map(region, "L2-AETI-A", period, folder)
ds = gdal.Open(fp3)
band = ds.GetRasterBand(1)
md = band.GetMetadata()
assert md["units"] == "mm/year"
ds = ds.FlushCache()

fp4 = wapor_map(region, "L1-AETI-D", period, folder)
fp5 = wapor_map(region, "L1-AETI-M", period, folder)
fp6 = wapor_map(region, "L1-AETI-A", period, folder)

fp7 = wapor_map(bb, "L2-AETI-D", period, folder)
ds = gdal.Open(fp7)
geot = ds.GetGeoTransform()
nx = ds.RasterXSize
ny = ds.RasterYSize
assert bb[0] == geot[0]
assert bb[2] == geot[0] + nx * geot[1]
assert bb[1] == geot[3] + ny * geot[5]
assert bb[3] == geot[3]
band = ds.GetRasterBand(1)
md = band.GetMetadata()
scale = band.GetScale()
array = band.ReadAsArray() * scale
ndv = band.GetNoDataValue()
array[array == ndv*scale] = np.nan
mean = np.nanmean(array)
assert mean > 0.0
assert mean < 15.0
proj = osr.SpatialReference(wkt=ds.GetProjection())
assert proj.GetAttrValue('AUTHORITY',1) == "4326"
ds = ds.FlushCache()

fp8 = wapor_map(bb, "L2-AETI-M", period, folder)
fp9 = wapor_map(bb, "L2-AETI-A", period, folder)

fp10 = wapor_map(bb, "L1-AETI-D", period, folder)
fp11 = wapor_map(bb, "L1-AETI-M", period, folder)
fp12 = wapor_map(bb, "L1-AETI-A", period, folder)

fp13 = wapor_map(region, "L1-AETI-D", period, folder, extension=".nc")
info = gdal.Info(fp13, format = "json")
assert len(info["metadata"]["SUBDATASETS"]) == 4
ds = gdal.Open(info["metadata"]["SUBDATASETS"]["SUBDATASET_1_NAME"])
band = ds.GetRasterBand(1)
md = band.GetMetadata()
scale = band.GetScale()
array = band.ReadAsArray() * scale
ndv = band.GetNoDataValue()
array[array == ndv*scale] = np.nan
mean = np.nanmean(array)
assert mean > 0.0
assert mean < 15.0
proj = osr.SpatialReference(wkt=ds.GetProjection())
assert proj.GetAttrValue('AUTHORITY',1) == "4326"
ds = ds.FlushCache()

fp14 = wapor_map("BKA", "L3-T-D", period, folder)
ds = gdal.Open(fp14)
band = ds.GetRasterBand(1)
md = band.GetMetadata()
proj = osr.SpatialReference(wkt=ds.GetProjection())
assert proj.GetAttrValue('AUTHORITY',1) == "32636"

bbBKA = [35.75,33.70,35.82,33.75]
fp15 = wapor_map(bbBKA, "L3-T-D", period, folder)
assert "bb.BKA" in fp15
ds = gdal.Open(fp15)
band = ds.GetRasterBand(1)
md = band.GetMetadata()
scale = band.GetScale()
array = band.ReadAsArray() * scale
ndv = band.GetNoDataValue()
array[array == ndv*scale] = np.nan
mean = np.nanmean(array)
assert mean > 0.0
assert mean < 15.0
proj = osr.SpatialReference(wkt=ds.GetProjection())
assert proj.GetAttrValue('AUTHORITY',1) == "32636"
ds = ds.FlushCache()

region_GEZ = [str(x) for x in l3_regions if "GEZ" in str(x)][0]
fp16 = wapor_map(region_GEZ, "L3-T-D", period, folder)

region_MUV = [str(x) for x in l3_regions if "MUV" in str(x)][0]
fp17 = wapor_map(region_MUV, "L3-T-D", period, folder)

region_MULTIPLE = [str(x) for x in l3_regions if "MULTIPLE" in str(x)][0]
fp18 = wapor_map(region_MULTIPLE, "L3-T-D", period, folder)

try: ##
    region_FAIL = [str(x) for x in l3_regions if "FAIL" in str(x)][0]
    fp19 = wapor_map(region_FAIL, "L3-T-D", period, folder)
except ValueError as e:
    if "`region` can't be linked to any L3 region." in str(e):
        print("succes")
    else:
        raise e

try: ##
    _ = wapor_map(region, "L1-AETI-M", period, folder, extension=".vrt")
except ValueError as e:
    if "Please use one of " in str(e):
        print("succes")
    else:
        raise e
    
try: ##
    _ = date_func("https://storage.googleapis.com/fao-gismgr-wapor-3-data/DATA/WAPOR-3/MOSAICSET/L3-T-D/WAPOR-3.L3-T-D.BKA.2021-01-D1.tif", "W")
except ValueError as e:
    if "Invalid temporal resolution." in str(e):
        print("succes")
    else:
        raise e
    
try: ##
    _ = collect_metadata("L4-AETI-D")
except ValueError as e:
    if "Invalid variable name" in str(e):
        print("succes")
    else:
        raise e
    
try: ##
    _ = generate_urls_v3("L4-AETI-D", l3_region = None, period = None)
except ValueError as e:
    if "Invalid level " in str(e):
        print("succes")
    else:
        raise e

try: ##
    _ = wapor_ts(region.replace(".geojson", "blabla.geojson"), "L1-AETI-A", period, overview)
except ValueError as e:
    if "Geojson file not found." in str(e):
        print("succes")
    else:
        raise e

try: ##
    _ = wapor_ts([25, -17, 24, -16], "L1-AETI-A", period, overview)
except ValueError as e:
    if "Invalid bounding box." in str(e):
        print("succes")
    else:
        raise e

try: ##
    _ = wapor_ts((500.0, "lala"), "L3-AETI-A", period, overview)
except ValueError as e:
    if "Invalid value for" in str(e):
        print("succes")
    else:
        raise e
    
try: ##
    _ = wapor_ts(region, "L1-AETI-A", ["2021-01-15", "2021-01-01"], overview)
except ValueError as e:
    if "Invalid period." in str(e):
        print("succes")
    else:
        raise e

try: ##
    _ = wapor_ts(region, "L1-AETI-A", ["2021-01-15", "2021-01-01"], overview, req_stats = None)
except ValueError as e:
    if "Please specify a list of required statistics." in str(e):
        print("succes")
    else:
        raise e

try: ##
    _ = wapor_ts(region, "L1-AETI-A", ["2021-01-15", "2021-01-01"], overview, req_stats = ["std"])
except ValueError as e:
    if "Please select at least one valid statistic from" in str(e):
        print("succes")
    else:
        raise e
    
try: ##
    _ = wapor_ts(region, "L2-AETI-D", nodata_period, overview)
except ValueError as e:
    if "No files found for selected region, variable and period." in str(e):
        print("succes")
    else:
        raise e    

# try:
#     bb_south_america = [-68.203125,-18.979026,-55.371094,-6.839170]
#     _ = wapor_map(bb_south_america, "L2-T-D", period, folder)
# except ValueError as e:
#     if "has no overlap with the datasets" in str(e):
#         print("succes")
#     else:
#         raise e
    
#####
# SUMMATION CHECK
#####
    
bb = [30.2, 28.6, 31.3, 30.5]
folder = "/Users/hmcoerver/Local/wapordl_test"
variable = "L2-AETI-D"
period = ["2018-01-01", "2018-12-31"]

fp_a_nc = wapordl.wapor_map(bb, "L2-AETI-A", period, folder, extension = ".nc")
fp_d_nc = wapordl.wapor_map(bb, "L2-AETI-D", period, folder, extension = ".nc")

ds_d = xr.open_dataset(fp_d_nc, decode_coords = "all")
coords = [np.datetime64(da.attrs["start_date"], "ns") for da in ds_d.data_vars.values()]
da_d = ds_d.to_array("time").assign_coords({"time": coords})

length = xr.where(da_d["time"].dt.day != 21, 10, da_d["time"].dt.daysinmonth - 20)
da_d = (da_d * length).sum(dim = "time")

ds_a = xr.open_dataset(fp_a_nc, decode_coords = "all")
da_a = ds_a["Band1"]

assert abs((da_a - da_d).mean().values) < 0.00001
    
# fig, axs = plt.subplots(1, 2, figsize = (15, 5))
# da_d.plot(ax = axs[0])
# axs[0].set_title("summed dekad")
# da_a.plot(ax = axs[1])
# axs[1].set_title("annual")