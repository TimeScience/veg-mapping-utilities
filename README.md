Project forked from John Burley and Chris Bradley's [PaddocksTS](https://github.com/johnburley3000/PaddockTS/) project

Main changes thus far:
* Bash files and Pre-segment code adapted to pull Sentinel-2 data from DEA datasets on AWS rather than on NCI
* 01_pre-segment.py now can take either a geojson file or lat/lon/buffer
* If using a geojson file, update the filepath to your geojson folder in the _main.sh file before running it and change the STUB name accordingly

#### Workflow to download new sentinel datasets:
1. Open one of the _main.sh files and save it with the name of your new dataset

1. Create a new geojson file with the area you want to download
2. Name the geojson by the area being downloaded.
3. Make a copy of one of the _main.sh files and save it with the name of your new dataset
4. In the _main.sh files, change the STUB name to be the same as the filename and the geojson area
5. Include a timestamp in the name, for example:
	geojson: ACT-NE-v1.geojson

  One year:
    * stub: ACT-NE-v1_2021
    * main.sh file: ACT-NE-v1_2021_main.s
	
	Multiple years:
    * stub: ACT-NE-v1_2020-2024
    * main.sh file: ACT-NE-v1_2020-2_main.sh


# PaddockTS notebooks have been adapted for ACT Veg mapping using various ACT datasets.
* New Jupyter notebooks for processing this data are now in the [act-environment repo](https://github.com/TimeScience/act-environment) (at least until I have time to organise them better)

### Note for JupyterLab setup on GADI/ARE at NCI
When launching a new JupyterLab the following must be set:

* Walltime: _Whatever you are likely to use_
* Queue: _expressbw_
* Compute Size: *cpus=28 mem=252G*  <- this can be smaller if you aren't loading huge time-series pkl files
* Project: _xe2_
* Storage: _gdata/xe2+scratch/xe2+gdata/v10+gdata/ub8_
* Module directories: _/g/data/v10/public/modules/modulefiles/_
* Modules: _ffmpeg/4.3.1 tensorflow/2.15.0 gdal/3.6.4_




