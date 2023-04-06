#!/bin/bash

if [ -z "${RANK_STRIDE}" ]; then
    let RANK_STRIDE=128/${OMPI_COMM_WORLD_LOCAL_SIZE}
fi

if [ -z "${OMP_STRIDE}" ]; then
    let OMP_STRIDE=1
fi

if [ -z "${NUM_GPUS}" ]; then
    let NUM_GPUS=4
fi

if [ -z "${CPU_SHIFT}" ]; then
    let CPU_SHIFT=0
fi

if [ -z "${GPU_START}" ]; then
    let GPU_START=0
fi

if [ -z "${GPU_STRIDE}" ]; then
    let GPU_STRIDE=1
fi

let ranks_per_gpu=$(((${OMPI_COMM_WORLD_LOCAL_SIZE}+${NUM_GPUS}-1)/${NUM_GPUS}))
echo $ranks_per_gpu
let my_gpu=$(($OMPI_COMM_WORLD_LOCAL_RANK*$GPU_STRIDE/$ranks_per_gpu))+${GPU_START}

let cpu_start=$(($RANK_STRIDE*$OMPI_COMM_WORLD_LOCAL_RANK))+${GPU_START}+${CPU_SHIFT}
let cpu_stop=$(($cpu_start+$OMP_NUM_THREADS*$OMP_STRIDE-1))
#export GOMP_CPU_AFFINITY=$cpu_start-$cpu_stop:$OMP_STRIDE
export OMP_PLACES="{$cpu_start:$OMP_NUM_THREADS:$OMP_STRIDE}"

export ROCR_VISIBLE_DEVICES=$my_gpu

#echo "rank_local= " $OMPI_COMM_WORLD_LOCAL_RANK "  GOMP_CPU_AFFINITY= " $GOMP_CPU_AFFINITY "  ROCR_VISIBLE_DEVICES= " $ROCR_VISIBLE_DEVICES
echo "rank_local= " $OMPI_COMM_WORLD_LOCAL_RANK "  OMP_PLACES= " $OMP_PLACES "  ROCR_VISIBLE_DEVICES= " $ROCR_VISIBLE_DEVICES

prof_name=motorBike_piso_E

#eval "${ROCM_PATH}/bin/rocprof --hsa-trace  --roctx-trace   -d ${prof_name}.${OMPI_COMM_WORLD_RANK} -o ${prof_name}.${OMPI_COMM_WORLD_RANK}/${prof_name}.${OMPI_COMM_WORLD_RANK}.csv $*"

$@

