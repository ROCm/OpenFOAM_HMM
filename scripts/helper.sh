#!/bin/bash

#assume server has 8 GPUs
let NGPUS=8

#set visible device for each process  based on local to each server rank ID and total number of visible devices

let RANKS_PER_GPU=$((($OMPI_COMM_WORLD_LOCAL_SIZE+$NGPUS-1)/$NGPUS))
export ROCR_VISIBLE_DEVICES=$(($OMPI_COMM_WORLD_LOCAL_RANK/$RANKS_PER_GPU+(0)))
echo "deviceID: " $ROCR_VISIBLE_DEVICES "  rank = " $OMPI_COMM_WORLD_LOCAL_RANK  "RANKS_PER_GPU = " $RANKS_PER_GPU

#---------------#
#to use rocprof we need to create names for directories and files for rocprof output
pid="$$"
outdir="openfoam.${OMPI_COMM_WORLD_LOCAL_RANK}.${pid}"
outfile="${pid}_${OMPI_COMM_WORLD_LOCAL_RANK}.csv"

#uncomment the following line in order to use rocprof
#eval "rocprof -d ${outdir} -o ${outdir}/${outfile}  $*" 


#--------------#

$@
