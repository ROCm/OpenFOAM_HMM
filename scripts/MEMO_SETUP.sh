export HSA_XNACK=1
export OMP_NUM_THREADS=2

export ROCM4FOAM=/global/software/rocm-afar001-732
export MPI4FOAM=/global/software/openmpi/rocm-afar001/openmpi-5.0.0rc2


export PATH=${MPI4FOAM}/build/bin:${MPI4FOAM}/build/lib:$PATH
export LIBRARY_PATH=${MPI4FOAM}/build/lib:$LIBRARY_PATH

export LD_LIBRARY_PATH=${ROCM4FOAM}/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=${ROCM4FOAM}/lib:$LIBRARY_PATH

#for profiler
export PATH=${ROCM4FOAM}/bin:$PATH

source etc/bashrc




