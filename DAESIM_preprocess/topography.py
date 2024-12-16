# +
# Pysheds documentation is here: https://mattbartos.com/pysheds/

import os

# Dependencies
import numpy as np
from scipy import ndimage
from pysheds.grid import Grid
import matplotlib.pyplot as plt
from matplotlib import colors
import rasterio

# Local imports
os.chdir(os.path.join(os.path.expanduser('~'), "Projects/PaddockTS"))
from DAESIM_preprocess.util import gdata_dir, scratch_dir, plot_categorical

dirmap = (64, 128, 1, 2, 4, 8, 16, 32)

def pysheds_accumulation(tiff_file):
    """Read in the grid and dem and calculate the water flow direction and accumulation"""
    
    # Load both the dem (basically a numpy array), and the grid (all the metadata like the extent)
    grid = Grid.from_raster(tiff_file)
    dem = grid.read_raster(tiff_file)

    # Hydrologically enforce the DEM so water can flow downhill to the edge and not get stuck
    pit_filled_dem = grid.fill_pits(dem)
    flooded_dem = grid.fill_depressions(pit_filled_dem)
    inflated_dem = grid.resolve_flats(flooded_dem)

    # Calculate the direction and accumulation of water
    fdir = grid.flowdir(inflated_dem)
    acc = grid.accumulation(fdir)

    return grid, dem, fdir, acc

def find_segment_above(acc, coord, branches_np):
    """Look for a segment upstream with the highest accumulation"""
    segment_above = None
    acc_above = -1
    for i, branch in enumerate(branches_np):
        if branch[-1] == coord:
            branch_acc = acc[branch[-2][0], branch[-2][1]] 
            if branch_acc > acc_above:
                segment_above = i
                acc_above = branch_acc
    return segment_above

def catchment_gullies(grid, fdir, acc, num_catchments=10):
    """Find the largest gullies"""

    # Extract the branches
    branches = grid.extract_river_network(fdir, acc > np.max(acc)/(num_catchments*10))

    # Convert the branches to numpy coordinates 
    branches_np = []
    for i, feature in enumerate(branches["features"]):
        line_coords = feature['geometry']['coordinates']
        branch_np = []
        for coord in line_coords:
            col, row = ~grid.affine * (coord[0], coord[1])
            row, col = int(round(row)), int(round(col))
            branch_np.append([row,col])
        branches_np.append(branch_np)

    # Repeatedly find the main segments to the branch with the highest accumulation. 
    full_branches = []
    for i in range(num_catchments):
        
        # Using the second last pixel before it's merged with another branch.
        branch_accs = [acc[branch[-2][0], branch[-2][1]] for branch in branches_np]
        largest_branch = np.argmax(branch_accs)

        # Follow the stream all the way up this branch
        branch_segment_ids = []
        while largest_branch != None:
            upper_coord = branches_np[largest_branch][0]
            branch_segment_ids.append(largest_branch)
            largest_branch = find_segment_above(acc, upper_coord, branches_np)

        # Combine the segments in this branch
        branch_segments = [branches_np[i] for i in sorted(branch_segment_ids)]
        branch_combined = [item for sublist in branch_segments for item in sublist]
        full_branches.append(branch_combined)

        # Remove all the segments from that branch and start again
        branch_segments_sorted = sorted(branch_segment_ids, reverse=True)
        for i in branch_segments_sorted:
            del branches_np[i]

    # Extract the gullies
    gullies = np.zeros(acc.shape, dtype=bool)
    for branch in full_branches:
        for x, y in branch:
            gullies[x, y] = True

    return gullies, full_branches


def catchment_ridges(grid, fdir, acc, full_branches):
    """Finds the ridges/catchment boundaries corresponding to those gullies"""

    # Progressively delineate each catchment
    catchment_id = 1
    all_catchments = np.zeros(acc.shape, dtype=int)
    for branch in full_branches:
        
        # Find the coordinate with second highest accumulation
        coords = branch[-2]

        # Convert from numpy coordinate to geographic coordinate
        x, y = grid.affine * (coords[1], coords[0])

        # Generate the catchment above that pixel
        catch = grid.catchment(x=x, y=y, fdir=fdir, 
                            xytype='coordinate')

        # Override relevant pixels in all_catchments with this new catchment_id
        all_catchments[catch] = catchment_id
        catchment_id += 1

    # Find the edges of the catchments
    sobel_x = ndimage.sobel(all_catchments, axis=0)
    sobel_y = ndimage.sobel(all_catchments, axis=1)  
    edges = np.hypot(sobel_x, sobel_y) 
    ridges = edges > 0

    return ridges

def show_acc(acc, outdir=scratch_dir, stub="Test", reference_tiff=None):
    """Visualization of water accumulation"""
    # Save as GeoTIFF
    if reference_tiff:
        geotiff_filename = os.path.join(outdir, f"{stub}_topographic_index.tif")
        save_geotiff(acc, geotiff_filename, reference_tiff)
        print("Saved:", geotiff_filename)
    
    # Original PNG visualization
    fig, ax = plt.subplots(figsize=(8,6))
    fig.patch.set_alpha(0)
    plt.grid('on', zorder=0)
    im = ax.imshow(acc, zorder=2,
                   cmap='cubehelix',
                   norm=colors.LogNorm(1, acc.max()),
                   interpolation='bilinear')
    plt.colorbar(im, ax=ax, label='Upstream Cells')
    plt.title('Topographic Index', size=14)
    plt.tight_layout()
    png_filename = os.path.join(outdir, f"{stub}_topographic_index.png")
    plt.savefig(png_filename)
    print("Saved:", png_filename)
    plt.close()

def show_ridge_gullies(dem, ridges, gullies, outdir=scratch_dir, stub="Test", reference_tiff=None):
    """Visualization of ridges and gullies"""
    # Save as GeoTIFF
    if reference_tiff:
        # Save ridges
        ridge_filename = os.path.join(outdir, f"{stub}_ridges.tif")
        save_geotiff(ridges.astype('uint8'), ridge_filename, reference_tiff)
        print("Saved:", ridge_filename)
        
        # Save gullies
        gully_filename = os.path.join(outdir, f"{stub}_gullies.tif")
        save_geotiff(gullies.astype('uint8'), gully_filename, reference_tiff)
        print("Saved:", gully_filename)
    
    # Original PNG visualization
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)
    im = ax.imshow(dem, cmap='terrain', zorder=1, interpolation='bilinear')
    plt.colorbar(im, ax=ax, label='Elevation (m)')
    ax.contour(ridges, levels=[0.5], colors='red', linewidths=1.5, zorder=2)
    ax.contour(gullies, levels=[0.5], colors='blue', linewidths=1.5, zorder=3)
    ax.contour(dem, colors='black', linewidths=0.5, zorder=4, alpha=0.5)
    plt.title('Ridges and Gullies', size=14)
    plt.tight_layout()
    png_filename = os.path.join(outdir, f"{stub}_ridge_gullies.png")
    plt.savefig(png_filename)
    print("Saved:", png_filename)
    plt.close()

def show_aspect(fdir, outdir=scratch_dir, stub="Test", reference_tiff=None):
    """Visualization of aspect"""
    # Save as GeoTIFF
    if reference_tiff:
        geotiff_filename = os.path.join(outdir, f"{stub}_aspect.tif")
        save_geotiff(fdir, geotiff_filename, reference_tiff)
        print("Saved:", geotiff_filename)

    # Original visualization code...
    directions = {
        64: "North", 128: "Northeast", 1: "East", 2: "Southeast",
        4: "South", 8: "Southwest", 16: "West", 32: "Northwest",
        -1: "East", -2: "East"
    }
    colour_dict = {
        "East": '#EE82EE', "Southeast": '#00008B', "South": '#ADD8E6',
        "Southwest": '#006400', "West": '#90EE90', "Northwest": '#FFFF00',
        "North": '#FFA500', "Northeast": '#DC143C',
    }
    fdir_categorized = np.vectorize(directions.get)(fdir)
    png_filename = os.path.join(outdir, f"{stub}_aspect.png")
    plot_categorical(fdir_categorized, colour_dict, "Aspect", png_filename)

def show_slope(slope, outdir=scratch_dir, stub="Test", reference_tiff=None):
    """Visualization of slope"""
    # Save as GeoTIFF
    if reference_tiff:
        geotiff_filename = os.path.join(outdir, f"{stub}_slope.tif")
        save_geotiff(slope, geotiff_filename, reference_tiff)
        print("Saved:", geotiff_filename)

    # Original PNG visualization
    bin_edges = np.arange(0, 16, 1)
    slope_categories = np.digitize(slope, bin_edges, right=True)
    colours = plt.cm.viridis(np.linspace(0, 1, len(bin_edges) - 1))
    cmap = colors.ListedColormap(colours)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(slope_categories, cmap=cmap)
    labels = [f'{bin_edges[i]}°' for i in range(len(bin_edges))]
    cbar = plt.colorbar(im, ticks=np.arange(len(bin_edges)))
    cbar.ax.set_yticklabels(labels)
    plt.title('Slope', size=14)
    plt.tight_layout()
    png_filename = os.path.join(outdir, f"{stub}_slope.png")
    plt.savefig(png_filename)
    print("Saved:", png_filename)
    plt.close()

def visualise_canopy_height(filename, outdir=scratch_dir, stub="Test"):
    """Visualization of canopy height"""
    with rasterio.open(filename) as src:
        canopy_height = src.read(1)
        
    # Save as GeoTIFF
    geotiff_filename = os.path.join(outdir, f"{stub}_canopy_height.tif")
    save_geotiff(canopy_height, geotiff_filename, filename)
    print("Saved:", geotiff_filename)

    # Original PNG visualization
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(canopy_height, cmap='viridis')
    plt.colorbar(im, ax=ax, label='Height (m)')
    plt.title('Canopy Height', size=14)
    plt.tight_layout()
    png_filename = os.path.join(outdir, f"{stub}_canopy_height.png")
    plt.savefig(png_filename)
    print("Saved:", png_filename)
    plt.close()

def visualise_soil_texture(outdir_tb, outdir=scratch_dir, stub="Test"):
    """Visualization of soil texture"""
    filename = os.path.join(outdir_tb, f"{stub}_soil_texture.tif")
    with rasterio.open(filename) as src:
        soil_texture = src.read(1)
        
    # Save as GeoTIFF
    geotiff_filename = os.path.join(outdir, f"{stub}_soil_texture.tif")
    save_geotiff(soil_texture, geotiff_filename, filename)
    print("Saved:", geotiff_filename)

    # Original PNG visualization code...
    # (Add your existing visualization code here)

def visualise_soil_pH(outdir_tb, outdir=scratch_dir, stub="Test"):
    """Visualization of soil pH"""
    filename = os.path.join(outdir_tb, f"{stub}_soil_pH.tif")
    with rasterio.open(filename) as src:
        soil_pH = src.read(1)
        
    # Save as GeoTIFF
    geotiff_filename = os.path.join(outdir, f"{stub}_soil_pH.tif")
    save_geotiff(soil_pH, geotiff_filename, filename)
    print("Saved:", geotiff_filename)

    # Original PNG visualization code...
    # (Add your existing visualization code here)

def calculate_slope(tiff_file):
    """Calculate the slope of a DEM"""
    with rasterio.open(tiff_file) as src:
        dem = src.read(1)  
        transform = src.transform 
    gradient_y, gradient_x = np.gradient(dem, transform[4], transform[0])
    slope = np.arctan(np.sqrt(gradient_x**2 + gradient_y**2)) * (180 / np.pi)
    return slope


def save_geotiff(data, output_path, reference_tiff):
    """
    Save data as a GeoTIFF using the same spatial reference as reference_tiff
    """
    with rasterio.open(reference_tiff) as src:
        transform = src.transform
        crs = src.crs
        
    with rasterio.open(
        output_path,
        'w',
        driver='GTiff',
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
    ) as dst:
        dst.write(data, 1)


if __name__ == '__main__':
    filepath = "/g/data/xe2/cb8590/Data/PadSeg/MILG_terrain.tif"
    grid, dem, fdir, acc = pysheds_accumulation(filepath)
    # show_acc(acc)
    show_aspect(fdir)
    
    # num_catchments = 10
    # gullies, full_branches = catchment_gullies(grid, fdir, acc, num_catchments)
    # ridges = catchment_ridges(grid, fdir, acc, full_branches)
    # show_ridge_gullies(dem, ridges, gullies)

    # slope = calculate_slope(filepath)
    # show_slope(slope)


