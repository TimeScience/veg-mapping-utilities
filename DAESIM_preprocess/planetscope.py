# +
# Planet CLI refernce is here: https://planet-sdk-for-python.readthedocs.io/en/latest/cli/reference.html
# Planet band info is here: https://developers.planet.com/docs/apis/data/sensors/

# +
# Assumes tiff files have already been downloaded using John's bash script

# +
# Standard library
import os
import sys
import json
import pickle
from datetime import datetime

# Dependencies
import pandas as pd
import numpy as np
import xarray as xr
import rioxarray as rxr


# -


sys.path.insert(1, '../Tools/')
from dea_tools.plotting import rgb

# +
orderid = ""
outpath = ""
stub = "MILG"
base_dir = "/g/data/xe2/datasets/Planet/Farms"

order_ids = "5534c295-c405-4e10-8d51-92c7c7c02427", "570f6f33-9471-4dcb-a553-2c97928c53aa"


# +
def find_timestamps(base_dir, stub, order_id):
    """Find all the timestamps in a planetscope order folder"""
    timestamps = set()
    dir_name = os.path.join(base_dir, stub, order_id, "PSScene")
    for s in os.listdir(dir_name):
        timestamp = s[:23]
        timestamps.add(timestamp)
    if 'PSScene_collection.json' in timestamps:
        timestamps.remove('PSScene_collection.json')
    timestamps = sorted(timestamps)
    return timestamps

timestamps = find_timestamps(base_dir, stub, order_ids[0])
timestamps[:5]


# +
def load_image(base_dir, stub, order_id, timestamp, bands_tiff_prefix="_3B_Visual_clip"):
    """Load a single planetscope image into an xarray"""
    
    filename = os.path.join(base_dir, stub, order_id, "PSScene", f"{timestamp}{bands_tiff_prefix}.tif")
    da = rxr.open_rasterio(filename)
    
    # In the 3band tiff, there is actually 4 bands, with the 4th being a cloud mask. Here we apply the cloud mask to the other bands and then drop it.
    if bands_tiff_prefix == "_3B_Visual_clip":
        cloud_mask = da.sel(band=4)
        da_masked = da.where(cloud_mask != 0, other=np.nan)
        da_3band = da_masked.sel(band=slice(1, 3))

    # Extract the bands into their own variables to match sentinel
    ds = da_3band.to_dataset(dim='band')
    ds_named = ds.rename({1: 'nbart_red', 2: 'nbart_blue', 3: 'nbart_green'})


    return ds_named

ds = load_image(base_dir, stub, order_ids[0], timestamps[0])
print(ds)
rgb(ds)


# +
# %%time
def load_directory(base_dir, stub, order_id, bands3_tiff_prefix="_3B_Visual_clip", limit=None):
    """Load all the in a single planetscope order folder into a single xarray"""
    
    timestamps = find_timestamps(base_dir, stub, order_id)
    
    if limit:
        timestamps = timestamps[:limit]
    
    dss = []
    for timestamp in timestamps:
        ds = load_image(base_dir, stub, order_id, timestamp, bands3_tiff_prefix)
        time = pd.to_datetime(ds.attrs['TIFFTAG_DATETIME'], format='%Y:%m:%d %H:%M:%S')
        ds_timed = ds.expand_dims(time=[time])
        dss.append(ds_timed)
        
    combined_ds = xr.concat(dss, dim='time')
    return combined_ds

ds = load_directory(base_dir, stub, order_ids[1], limit=3)
print(ds)
rgb(ds, col="time")
# +
def load_directories(base_dir, stub, order_ids, bands_tiff_prefix="_3B_Visual_clip", limit=None):
    """Load images from multiple directories"""
    dss = []
    for order_id in order_ids:
        ds = load_directory(base_dir, stub, order_id, bands_tiff_prefix, limit)
        dss.append(ds)    
    combined_ds = xr.concat(dss, dim='time')
    return combined_ds

ds = load_directories(base_dir, stub, order_ids, limit=3)
print(ds)
rgb(ds, col='time')
# -

if __name__ == '__main__':

    base_directory = args.indir
    order_id = args.orderid
    outpath = args.outpath

    print(f"{datetime.now()} Starting planet_xarray_3b in {base_directory} for {order_id}")

    xarray = load_single_order_3band(base_directory, order_id)

    with open(outpath, 'wb') as handle:
        pickle.dump(xarray, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"{datetime.now()} Finished planet_xarray_3b and exported to {outpath}")

