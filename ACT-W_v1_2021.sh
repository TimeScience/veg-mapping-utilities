#!/bin/bash

## Config:
# specify working directory:
wd=/home/801/tbb801/Projects/PaddockTS
dir=/g/data/xe2/tbb801/Data/PadSeg
tmpdir=/scratch/xe2/tbb801/tmp  
# params


#lat=-35.187472
#lon=149.145069
#buffer=0.02
stub=ACT-W_v1_2021
start=2021-01-01
end_=2021-12-31

stub=test-9

# NOTE: When "use_geojson=1", lat/lon/buffer are ignored

use_geojson=true # MUST be true or false (lowercase), not 0,1 or True or False due to stupid issues with how python parses booleans when called by bash

geojson_path=/g/data/xe2/tbb801/Data/PadSeg/geojson/ACT-NW-01_v1.geojson

# params for paddock filtering
min_area_ha=10
max_area_ha=1500
max_perim_area_ratio=40

if [[ $use_geojson == true ]]
then
	lat=0
	lon=0
	buffer=0
	echo "use_geojson = TRUE so using geojson file: $geojson_path and setting lat, lon and buffer=0"

else
	echo "use_geojson= false, passing lat, lon, buffer"
	geojson_path=""

fi

## Run first job
#job_id1=$(qsub -v wd=$wd,stub=$stub,dir=$dir,lat=$lat,lon=$lon,buffer=$buffer,start_time=$start,end_time=$end_ Code/run_pre-seg_aws.sh)

## This version downloads the Sentinel data from AWS rather than NCI
job_id1=$(qsub -v wd=$wd,stub=$stub,dir=$dir,lat=$lat,lon=$lon,buffer=$buffer,start_time=$start,end_time=$end_,use_geojson=$use_geojson,geojson_path=$geojson_path Code/run_pre-seg_aws.sh)

echo "First job submitted running run_pre-seg_aws - with ID $job_id1"

## Run second job (if job ID was produced, and when job complete)
if [[ -z "$job_id1" ]]
then
    echo "Failed to submit the first job."
    exit 1
else
    echo "NOT Submitting SAMGeo job - downloading data only"
    #qsub -W depend=afterok:$job_id1 -v wd=$wd,stub=$stub,dir=$dir,min_area_ha=$min_area_ha,max_area_ha=$max_area_ha,max_perim_area_ratio=$max_perim_area_ratio Code/run_SAMGeo_paddocks-ts.sh
    #qsub -v stub=$stub,dir=$dir,min_area_ha=$min_area_ha,max_area_ha=$max_area_ha,max_perim_area_ratio=$max_perim_area_ratio Code/run_SAMGeo_paddocks-ts.sh
fi

