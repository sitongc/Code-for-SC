#!/bin/bash

path="/Users/lyang/OneDrive/Hansen_lab/spatial_data/spatial/batch2"

#src="/mnt/c/Users/zhuzh/OneDrive - University of Toronto/Lin_scripts/projects/spatial"

mkdir $path/2_scanpy_process
for File in $(cat ${path}/sam_list.txt);

do
  sample=$(basename "$File" .h5ad)
  mkdir $path/2_scanpy_process/$sample
  python ./scanpy_process.py -i $path/scanpy_object/$File -o $path/2_scanpy_process/ \
         -s ${sample} -m $path/singleR_labels/${sample}.csv
done
