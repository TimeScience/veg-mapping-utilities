#!/bin/bash

## Config:
# specify working directory:
wd=/home/801/tbb801/Projects/PaddockTS
dir=/g/data/xe2/tbb801/Data/PadSeg
tmpdir=/scratch/xe2/tbb801/tmp  
# params

lat=-35.558408
lon=149.079151
buffer=0.032616
stub="riverview-1"
start='2017-01-01'
end_='2024-12-31'


# params for paddock filtering
min_area_ha=10
max_area_ha=1500
max_perim_area_ratio=40


## Run first job
job_id1=$(qsub -v wd=$wd,stub=$stub,dir=$dir,lat=$lat,lon=$lon,buffer=$buffer,start_time=$start,end_time=$end_ Code/run_pre-seg_aws.sh)
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

