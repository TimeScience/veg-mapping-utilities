"""
Description:
This script downloads SENTINEL2 data (using DEA) as an xarray dataset, saves it as a pickle, and also prepares a 3-band Fourier Transform of the NDWI time series that can be used as input in SAMGeo for paddock segmentation. 

Usage:
The script is designed to be executed from the command line, where the user can specify the stub name for file naming and the directory for output files, and other parameters:

Requirements:
- Runs using a python environment designed for DEA Sandbox use on NCI. 
module use /g/data/v10/public/modules/modulefiles
module load dea/20231204
- A fair amount of memory for downloading large regions of data. 

Inputs:
- stub name
- coordinates
- buffer (degrees)
- start/end date

Outputs:
- A pickle dump of the xarray (ds2) containing sentinel data and MANY band incides, and metadata, for the ROI and TOI
- a tif image of the Fourier Transform reprojection of the time-stack into 3 axes. This is a 3-band image normalised between 1-255 to be compatible with SAM model.
- a pickle containg dict of the query used to generate ds2
"""
import argparse
import os
import sys
import logging
import pickle
import numpy as np
import xarray as xr
import rioxarray
import datacube
from dea_tools.temporal import xr_phenology, temporal_statistics
from dea_tools.datahandling import load_ard
from dea_tools.bandindices import calculate_indices
from dea_tools.plotting import display_map, rgb
from dea_tools.dask import create_local_dask_cluster
import hdstats

# Adjust logging configuration for the script
logging.basicConfig(level=logging.INFO)

# This line of code determines if the data is downloaded from the NCI project or from the AWS version of the same data
# Set AWS environment variable for public access:
os.environ["AWS_NO_SIGN_REQUEST"] = "YES" # Comment this out if you want to download data from NCI rather than AWS

def str2bool(v): #required because passing booleans to python from bash seems to be impossible without hours of suffering
    """Convert string to boolean"""
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Download and save Sentinel data, prepare input image for SAMGeo.
        Example usage:
        python3 Code/01_pre-segment.py --stub test1 --outdir /g/data/xe2/John/Data/PadSeg/ --lat -34.3890 --lon 148.4695 --buffer 0.01 --start_time '2020-01-01' --end_time '2020-03-31'""",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--stub", type=str, required=True, help="Stub name for file naming")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory for saved files")
    parser.add_argument("--lat", type=float, required=True, help="Latitude of the center of the area of interest")
    parser.add_argument("--lon", type=float, required=True, help="Longitude of the center of the area of interest")
    parser.add_argument("--buffer", type=float, required=True, help="Buffer in degrees to define the area around the center point")
    parser.add_argument("--start_time", type=str, required=True, help="Start time for the data query (YYYY-MM-DD)")
    parser.add_argument("--end_time", type=str, required=True, help="End time for the data query (YYYY-MM-DD)")
    parser.add_argument("--use_geojson", type=str2bool, required=True,
                       help="Whether to use geojson (true) or lat/lon/buffer (false)")   
    parser.add_argument("--geojson_path", type=str, required=False, help="Path to geojson file containing bounding box")

    args = parser.parse_args()
    # Validate arguments based on use_geojson
    logging.info(f"args =  {args}")
    if args.use_geojson == True:
        logging.info(f"args.use_geojson =  {args}")
        if not args.geojson_path:
            parser.error("--geojson_path is required when use_geojson=1")
    else:
        logging.info(f"else part 2: args.use_geojson =  {args.use_geojson}")
        if None in (args.lat, args.lon, args.buffer):
            parser.error("--lat, --lon, and --buffer are required when use_geojson=0")
            
    return parser.parse_args()

def define_query_lat_lon(lat, lon, buffer, time_range):
    lat_range = (lat-buffer, lat+buffer)
    lon_range = (lon-buffer, lon+buffer)
    query = {
        'centre': (lat, lon),
        'y': lat_range,
        'x': lon_range,
        'time': time_range,
        'resolution': (-10, 10),
        'output_crs': 'epsg:6933',
        'group_by': 'solar_day'
    }
    # note that centre is not recognized as query option in load_arc, but we want to output it as a record.
    return query

#Use this if we are passing a bounding box
def define_query_bbox(x_range, y_range, time_range):
    """Define query using bounding box coordinates"""
    query = {
        'x': x_range,
        'y': y_range,
        'time': time_range,
        'resolution': (-10, 10),
        'output_crs': 'epsg:6933',
        'group_by': 'solar_day'
    }
    return query
    
def get_bbox_coords(geojson_path):
    """Extract bounding box coordinates from geojson file"""
    import json
    
    with open(geojson_path) as f:
        geojson = json.load(f)
    
    # Access the first feature's geometry
    coords = geojson['features'][0]['geometry']['coordinates'][0]
    
    # Get min/max x and y coordinates
    x_coords = [coord[0] for coord in coords]
    y_coords = [coord[1] for coord in coords]
    
    x_range = (min(x_coords), max(x_coords))
    y_range = (min(y_coords), max(y_coords))
    
    print(f"Extracted bounding box coordinates:")
    print(f"X range: {x_range}")
    print(f"Y range: {y_range}")
    
    return x_range, y_range


def load_and_process_data(dc, query):
    query = {k: v for k, v in query.items() if k != 'centre'} # this takes centre out of the query	
    ds = load_ard(
        dc=dc,
        products=['ga_s2am_ard_3', 'ga_s2bm_ard_3'],
        cloud_mask='s2cloudless',
        min_gooddata=0.9,
        **query
    )
    ds = calculate_indices(ds, 
                           #index=['NDVI', 'kNDVI', 'EVI', 'LAI', 'SAVI', 'MSAVI', 'NDMI', 'NDWI', 'MNDWI', 'NBR', 'NDCI', 'NDTI', 'BSI'],
			   index=['NDVI', 'LAI', 'SAVI'],
                           collection='ga_s2_3')
    return ds
    
def transform(ds):
	keep_vars = ['nbart_red','nbart_green','nbart_blue','nbart_nir_1']
	data = ds[keep_vars].to_array().transpose('y', 'x','variable', 'time').values.astype(np.float32)
	data[data == 0] = np.nan
	data /= 10000.
	ndwi_obs = (data[:,:,1,:]-data[:,:,3,:])/(data[:,:,1,:]+data[:,:,3,:]) # w = water. (g-nir)/(g+nir)
	ndwi = hdstats.completion(ndwi_obs)
	f2 = hdstats.fourier_mean(ndwi)
	return f2
	
def rescale(im):
    '''rescale raster (im) to between 0 and 255.
    Attempts to rescale each band separately, then join them back together to achieve exact same shape as input.
    Note. Assumes multiple bands, otherwise breaks'''
    n_bands = im.shape[2]
    _im = np.empty(im.shape)
    for n in range(0,n_bands):
        matrix = im[:,:,n]
        scaled_matrix = (255*(matrix - np.min(matrix))/np.ptp(matrix)).astype(int)
        _im[:,:,n] = scaled_matrix
    print('output shape equals input:', im.shape == im.shape)
    return(_im)

def export_for_segmentation(ds, inp, out_stub):
    '''prepares a 3-band image for SAMgeo. 
    First rescale bands in the image. Then convert to xarray with original geo info. Then save geotif'''
    if inp.shape[2] == 3:
        image = rescale(inp) # 3d array 
        lat = list(ds.y.values) # latitude is the same size as the first axis
        lon = list(ds.x.values) # longitude is the same size as second axis
        bands = list(range(1,image.shape[2]+1)) # band is the 3rd axis
        crs = ds.rio.crs
        # create xarray object
        data_xr = xr.DataArray(image, 
                       coords={'y': lat,'x': lon,'band': bands}, 
                       dims=["y", "x", "band"])
        data_xr.rio.write_crs(crs, inplace=True)
        # save as geotif:
        data_xr.transpose('band', 'y', 'x').rio.to_raster(out_stub + '.tif')
    else:
        print("Input image is wrong shape! No action taken")
        
def main(args):
    client = create_local_dask_cluster(return_client=True)
    dc = datacube.Datacube(app='Vegetation_phenology')
    logging.info(f"args.use_geojson =  {args.use_geojson}")

    if args.use_geojson:
        logging.info(f"in use_geojson should be TRUE | use_geojson =  {args.use_geojson}")
        # Get bounding box coordinates from geojson
        x_range, y_range = get_bbox_coords(args.geojson_path) 
        query = define_query_bbox(x_range, y_range, (args.start_time, args.end_time))
    else:
        logging.info(f"running lat lon query")
        query = define_query_lat_lon(args.lat, args.lon, args.buffer, (args.start_time, args.end_time))
    
            
    ds = load_and_process_data(dc, query)
    client.close()
    f2 = transform(ds)
    im = rescale(f2)
    export_for_segmentation(ds, im, args.outdir+args.stub)
    
    # save ds for later
    outnamepath = os.path.join(args.outdir + '/sentinel-pkl', args.stub + '_ds2.pkl')
    with open(outnamepath, 'wb') as handle:
        pickle.dump(ds, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    logging.info(f"Data saved successfully to {outnamepath}")

# save query for record keeping
    outnamepath = os.path.join(args.outdir + '/sentinel-pkl', args.stub + '_ds2_query.pkl')
    with open(outnamepath, 'wb') as f:
        pickle.dump(query, f)
    logging.info(f"Query saved successfully to {outnamepath}")

if __name__ == "__main__":
    args = parse_arguments()
    main(args)