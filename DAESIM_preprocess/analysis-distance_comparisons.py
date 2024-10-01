# +
# Aim of this notebook is to combine all the datasources into a single xarray

# +
# Standard library
import os
import pickle

# Dependencies
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
from shapely.geometry import box, Polygon
from rasterio.enums import Resampling
from rasterio import features
import scipy.ndimage
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator

# Local imports
os.chdir(os.path.join(os.path.expanduser('~'), "Projects/PaddockTS"))
from DAESIM_preprocess.util import gdata_dir, scratch_dir, transform_bbox
from DAESIM_preprocess.topography import pysheds_accumulation, calculate_slope
from DAESIM_preprocess.silo_daily import merge_ozwald_silo, resample_weekly, visualise_water, visualise_temp
# -

stubs = {
    "MULL": "Mulloon",
    "CRGM": "Craig Moritz Farm",
    "MILG": "Milgadara",
    "ARBO": "Arboreturm",
    "KOWN": "Kowen Forest",
    "ADAM": "Canowindra",
    "LCHV": "Lachlan Valley"
}

# Filepaths
outdir = os.path.join(gdata_dir, "Data/PadSeg/")
stub = "MILG"

# %%time
# Sentinel imagery
filename = os.path.join(outdir, f"{stub}_ds2.pkl")
with open(filename, 'rb') as file:
    ds = pickle.load(file)

# +
# Weather data
filename_silo = os.path.join(outdir, f"{stub}_silo_daily.nc")
ds_silo = xr.open_dataset(filename_silo)

variable = 'max_temp'
df_drought = pd.DataFrame(ds_silo.isel(lat=0, lon=0)[variable], index=ds_silo.isel(lat=0, lon=0).time, columns=[variable])


# -

def add_tiff_band(ds, variable, resampling_method, outdir, stub):
    """Add a new band to the xarray from a tiff file using the given resampling method"""
    filename = os.path.join(outdir, f"{stub}_{variable}.tif")
    array = rxr.open_rasterio(filename)
    reprojected = array.rio.reproject_match(ds, resampling=resampling_method)
    ds[variable] = reprojected.isel(band=0).drop_vars('band')
    return ds


ds = add_tiff_band(ds, "canopy_height", Resampling.max, outdir, stub)

# The resampling often messes up the boundary, so we trim the outside pixels after adding all the resampled bounds
ds_trimmed = ds.isel(
    y=slice(1, -1),
    x=slice(1, -1) 
)

# +
# ds_trimmed['canopy_height'].plot()

# +
# ds_trimmed['NDVI'].isel(time=0).plot()
# -

ds_original = ds.copy()
ds = ds_trimmed

# +
# %%time
# Shelterscore showing the number of trees 
pixel_size = 10  # metres
# distance = 20  # This corresponds to a 200m radius if the pixel size is 10m
distances = 5, 10, 20, 30, 50, 100

# Classify anything with a height greater than 1 as a tree
tree_threshold = 1
tree_mask = ds['canopy_height'] >= tree_threshold

# # Find the pixels adjacent to trees
structuring_element = np.ones((3, 3))  # This defines adjacency (including diagonals)
adjacent_mask = scipy.ndimage.binary_dilation(tree_mask, structure=structuring_element)

for distance in distances:
    # Calculate the number of trees within a given distance for each pixel
    y, x = np.ogrid[-distance:distance+1, -distance:distance+1]
    kernel = x**2 + y**2 <= distance**2
    kernel = kernel.astype(float)
    shelter_score = scipy.ndimage.convolve(tree_mask.astype(float), kernel, mode='constant', cval=0.0)
    
    # Mask out trees and adjacent pixels
    shelter_score[np.where(adjacent_mask)] = np.nan
    
    # Add the shelter_score to the xarray
    shelter_score_da = xr.DataArray(
        shelter_score, 
        dims=("y", "x"),  
        coords={"y": ds.coords["y"], "x": ds.coords["x"]}, 
        name="shelter_score" 
    )
    layer_name = f"num_trees_{pixel_size * distance}m"
    ds[layer_name] = shelter_score_da
    print("Added layer: {layer_name}")


# +
layer_name = "num_trees_200m"
filename = os.path.join(scratch_dir, f'{stub}_{layer_name}.tif')
ds[layer_name].rio.to_raster(filename)
print(filename)

ds[layer_name].plot()

# +
# Example shelter vs productivity score
time = '2020-01-01'
ndvi = ds.sel(time=time, method='nearest')['NDVI']
productivity_score1 = ndvi.where(~adjacent_mask)
s = ds['num_trees_200m'].values

# Flatten the arrays for plotting
y = productivity_score1.values.flatten()
y_values = y[~np.isnan(y)]   # Remove all pixels that are trees or adjacent to trees
x = s.flatten()
x_values = x[~np.isnan(y)]   # Match the shape of the x_values
len(y_values)
# -


# Example 2d histogram
plt.hist2d(x_values, y_values, bins=100, norm=mcolors.PowerNorm(0.1))
plt.ylabel('NDVI', fontsize=12)
pixel_size = 10
plt.xlabel(f'Number of tree pixels within {distance * pixel_size}m', fontsize=12)
plt.title(stub + ": " + str(time)[:10], fontsize=14)
plt.show()

# +
# Example linear regression
res = stats.linregress(x_values, y_values)
print(f"Sample size: {len(x_values)}")
print(f"R-squared: {res.rvalue**2:.6f}")
print(f"Slope: {res.slope:.6f}")
plt.plot(x_values, y_values, 'o', label='original data')
plt.plot(x_values, res.intercept + res.slope*x_values, 'r', label='fitted line')

plt.ylabel('NDVI', fontsize=12)
pixel_size = 10
plt.xlabel(f'Number of tree pixels within {distance * pixel_size}m', fontsize=12)
plt.title(stub + ": " + str(time)[:10], fontsize=14)

plt.legend()
plt.show()
# -

# Example sheltered vs unsheltered threshold
distance = 20
shelter_threshold = 0.1  # Percentage tree cover
num_trees_threshold = ((distance * 2) ** 2) * shelter_threshold # Number of tree pixels
num_trees_threshold

# +
# Example box plot
sheltered = y_values[np.where(x_values >= num_trees_threshold)]
unsheltered = y_values[np.where(x_values < num_trees_threshold)]

plt.boxplot([unsheltered, sheltered], labels=['Unsheltered', 'Sheltered'], showfliers=False)

plt.title(stub + ": " + str(time)[:10], fontsize=14)
plt.ylabel('NDVI', fontsize=12)

print(f"Shelter threshold = {int(num_trees_threshold)} tree pixels within {distance * pixel_size}m")
print("Number of sheltered pixels: ", len(sheltered))
print("Number of unsheltered pixels: ", len(unsheltered))

plt.show()
# -

# This is a thing
standardized_mean_difference = (np.mean(sheltered) - np.mean(unsheltered))/np.std(y_values)
standardized_mean_difference

# Is this a thing?
standardized_median_difference = (np.median(sheltered) - np.median(unsheltered)) / stats.median_abs_deviation(y_values)
standardized_median_difference

# +
# Find the max_temp for each date with satellite imagery 
xarray_times = ds['time'].values
nearest_times_indices = df_drought.index.get_indexer(xarray_times, method='nearest')
nearest_df_times = df_drought.index[nearest_times_indices]

# Find the timepoints where the temp is greater than 25 degrees
temp_threshold = 25
selected_times = nearest_df_times[df_drought['max_temp'].iloc[nearest_times_indices] > temp_threshold]
ds_drought = ds.sel(time=selected_times, method='nearest')
# -

len(ds.time.values)

len(ds_drought.time.values)

# +
# %%time
# Calculate the productivity score

shelter_thresholds = 0.025, 0.05, 0.075, 0.1, 0.2, 0.3   # Percentage tree cover

total_benefits = []
benefit_scores_dict = {}

for i, shelter_threshold in enumerate(shelter_thresholds):
    print(f"Working on {i}/{len(shelter_thresholds)}", shelter_threshold)
    num_trees_threshold = ((distance * 2) ** 2) * shelter_threshold # Number of tree pixels

    benefit_scores = []

    # I should just loop through the hot days, instead of looping through everything and filtering later
    for i, time in enumerate(ds_drought.time.values):
        ndvi = ds.sel(time=time, method='nearest')['NDVI']
        productivity_score1 = ndvi.where(~adjacent_mask)
        s = ds['num_trees_200m'].values
        
        # Flatten the arrays for plotting
        y = productivity_score1.values.flatten()
        y_values = y[~np.isnan(y)]   # Remove all pixels that are trees or adjacent to trees
        x = s.flatten()
        x_values = x[~np.isnan(y)]   # Match the shape of the x_values
    
        # Select sheltered pixels and calculate z scores for NDVI at each pixel
        sheltered = y_values[np.where(x_values >= num_trees_threshold)]
        unsheltered = y_values[np.where(x_values < num_trees_threshold)]

        median_diff = np.median(sheltered) - np.median(unsheltered)
        median_diff_standard = (np.median(sheltered) - np.median(unsheltered)) / stats.median_abs_deviation(y_values)
        mean_diff_standard = (np.mean(sheltered) - np.mean(unsheltered))/np.std(y_values)

        benefit_score = {
            "time":time,
            "median_diff":median_diff,
            "median_diff_standard": median_diff_standard,
            "mean_diff_standard": mean_diff_standard,
        }
        benefit_score = benefit_score
        benefit_scores.append(benefit_score)

    benefit_scores_dict[shelter_threshold] = benefit_scores
    df_shelter = pd.DataFrame(benefit_scores)
    df_shelter = df_shelter.set_index('time')
    
    df_merged = df_shelter

    overall_median_diff = df_merged['median_diff'].median()
    overall_median_diff_standard = df_merged['median_diff_standard'].median()
    overall_mean_diff_standard = df_merged['mean_diff_standard'].median()

    total_benefit = {
            "shelter_threshold":shelter_threshold,
            "sheltered_pixels":len(sheltered),
            "unsheltered_pixels":len(y_values) - len(sheltered),
            "median_difference": overall_median_diff,
            "standardized_median_difference": overall_median_diff_standard,
            "standardized_mean_difference": overall_mean_diff_standard,
        }
    total_benefits.append(total_benefit)

print(len(total_benefits))
pd.DataFrame(total_benefits)
# -

# Plotting benefits over time at different productivity percentiles
df = pd.DataFrame(benefit_scores_dict[0.1])
df = df.set_index('time')
df = df.astype(float)
df.index = pd.to_datetime(df.index)
df = df[['p10', 'p25', 'p50', 'p75', 'p90']]
df.plot(figsize=(50,20))
ax = plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(100))
ax.axhline(0, color='black', linestyle='--', linewidth=1)
plt.xticks(rotation=45)
plt.show()

# +
# Plot a comparison of max temp and shelter benefit
df = df_merged
fig, ax1 = plt.subplots(figsize=(50, 20))  # Create the base figure and axis

# Maximum temperature
ax1.plot(df.index, df['max_temp'], color='tab:red', label='Maximum temperature')
ax1.set_xlabel('Date')
ax1.set_ylabel('Maximum temperature (°C)', color='tab:red')  # Set the y-axis label
ax1.tick_params(axis='y', labelcolor='tab:red')  # Change tick colors to match the line

# Shelter benefit
ax2 = ax1.twinx()  
ax2.plot(df.index, df['p50'], color='tab:green', label='Shelter benefit')
ax2.set_ylabel('Shelter benefit', color='tab:green')  # Set the second y-axis label
ax2.tick_params(axis='y', labelcolor='tab:green')  # Change tick colors to match the line

plt.xticks(rotation=45)
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
plt.show()
