#!/bin/bash

## Config:
# specify working directory:
wd=/home/801/tbb801/Projects/PaddockTS
# params
stub=CRGU
dir=/g/data/xe2/tbb801/Data/PadSeg
lat=-35.217844
lon=149.140914
buffer=0.05
start='2015-01-01'
end_='2024-12-31'
# params for paddock filtering
min_area_ha=10
max_area_ha=1500
max_perim_area_ratio=40


## Run first job
job_id1=$(qsub -v wd=$wd,stub=$stub,dir=$dir,lat=$lat,lon=$lon,buffer=$buffer,start_time=$start,end_time=$end_ Code/run_pre-seg.sh)
echo "First job submitted with ID $job_id1"

## Run second job (if job ID was produced, and when job complete)
if [[ -z "$job_id1" ]]
then
    echo "Failed to submit the first job."
    exit 1
else
    echo "Submitting second job, dependent on the completion of the first."
    qsub -W depend=afterok:$job_id1 -v wd=$wd,stub=$stub,dir=$dir,min_area_ha=$min_area_ha,max_area_ha=$max_area_ha,max_perim_area_ratio=$max_perim_area_ratio Code/run_SAMGeo_paddocks-ts.sh
    #qsub -v stub=$stub,dir=$dir,min_area_ha=$min_area_ha,max_area_ha=$max_area_ha,max_perim_area_ratio=$max_perim_area_ratio Code/run_SAMGeo_paddocks-ts.sh
fi
