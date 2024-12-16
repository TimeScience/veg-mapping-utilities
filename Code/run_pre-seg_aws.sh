#!/bin/bash
#PBS -N PreSegment
#PBS -l mem=186GB
#PBS -l ncpus=24
#PBS -l jobfs=24GB
#PBS -P xe2
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/xe2+gdata/v10+gdata/ka08
#PBS -q normal


# Print out input variables to the error log
echo "Running with the following input variables:"
echo "stub: $stub"
echo "outdir: $dir"
echo "latitude: $lat"
echo "longitude: $lon"
echo "buffer: $buffer"
echo "start date: $start_time"
echo "end date: $end_time"
echo "use_geojson: $use_geojson"
echo "use_geojson: $geojson_path"
echo "-------------------"

# Requirements:
# Needs access to project v10 to load the dea modules
# (also ka08 and xe2)

cd $wd

# Setup DEA environment modules for running the Python script
module use /g/data/v10/public/modules/modulefiles
module load dea/20231204

# Run the Python script with required parameters (provide these as command-line arguments or modify the script to set defaults)
if [[ $use_geojson == true ]]
then
	lat=0
	lon=0
	buffer=0
    python3 Code/01_pre-segment.py --stub $stub --outdir $dir --lat $lat --lon $lon --buffer $buffer --start_time $start_time --end_time $end_time --use_geojson true --geojson_path $geojson_path
else
    echo "in pre_seg_aws use_geojson=0, passing lat, lon, buffer"
    echo $use_geojson
    python3 Code/01_pre-segment.py --stub $stub --outdir $dir --lat $lat --lon $lon --buffer $buffer --start_time $start_time --end_time $end_time --use_geojson false
fi


