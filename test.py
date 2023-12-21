import glob
from wapordl import wapor_ts, wapor_map

regions = glob.glob(r"/Users/hmcoerver/Local/wapor_validation/SHAPES/BASINS/*.geojson")
folder = r"/Users/hmcoerver/Local/test"
region = regions[0]

bb = [25, -17, 26, -16]

overview = 3

period = ["2021-01-01", "2021-03-01"]
req_stats = ["minimum", "maximum", "mean"]
variable = "L2-AETI-D"

df1 = wapor_ts(region, "L2-AETI-D", period, overview)
df2 = wapor_ts(region, "L2-AETI-M", period, overview)
df3 = wapor_ts(region, "L2-AETI-A", period, overview)

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
fp2 = wapor_map(region, "L2-AETI-M", period, folder)
fp3 = wapor_map(region, "L2-AETI-A", period, folder)

fp4 = wapor_map(region, "L1-AETI-D", period, folder)
fp5 = wapor_map(region, "L1-AETI-M", period, folder)
fp6 = wapor_map(region, "L1-AETI-A", period, folder)

fp7 = wapor_map(bb, "L2-AETI-D", period, folder)
fp8 = wapor_map(bb, "L2-AETI-M", period, folder)
fp9 = wapor_map(bb, "L2-AETI-A", period, folder)

fp10 = wapor_map(bb, "L1-AETI-D", period, folder)
fp11 = wapor_map(bb, "L1-AETI-M", period, folder)
fp12 = wapor_map(bb, "L1-AETI-A", period, folder)

try:
    _ = wapor_ts(bb, "L3-AETI-A", period, overview)
except ValueError as e:
    if "Level-3 data will be available soon." in str(e):
        print("succes")
    else:
        raise e
    
try:
    _ = wapor_ts([25, -17, 24, -16], "L1-AETI-A", period, overview)
except ValueError as e:
    if "Invalid bounding box." in str(e):
        print("succes")
    else:
        raise e

try:
    _ = wapor_ts(region.replace(".geojson", "blabla.geojson"), "L1-AETI-A", period, overview)
except ValueError as e:
    if "Geojson file not found." in str(e):
        print("succes")
    else:
        raise e
    
try:
    _ = wapor_ts(None, "L1-AETI-A", period, overview)
except ValueError as e:
    if "Only level-3 variables can be processed without specifying a region." in str(e):
        print("succes")
    else:
        raise e

try:
    _ = wapor_ts(region, "L1-AETI-A", ["2021-01-15", "2021-01-01"], overview)
except ValueError as e:
    if "Invalid period." in str(e):
        print("succes")
    else:
        raise e
    
try:
    _ = wapor_ts(bb, "L1-AETI-W", period, overview)
except ValueError as e:
    if "Invalid temporal resolution." in str(e):
        print("succes")
    else:
        raise e
