

export HSA_XNACK=1
export OMP_NUM_THREADS=1
export OMPX_DISABLE_USM_MAPS=1
export GPU_MAX_HW_QUEUES=4
#export HSA_FORCE_STORESC1=0 

#export ROCM4FOAM=${HOME}/OPENMP_COMPILER_UVM/rocm-5.5.0
#export ROCM4FOAM=/opt/rocm-5.7.0-1207
export ROCM4FOAM=/opt/rocm-5.7.0-1330
#export ROCM4FOAM=/home/rlieberm/rocm-5.7.0-12374
#export MPI4FOAM=${THERA_OMPI_DIR}

#sh5/MI300A
export MPI4FOAM=/opt/ompi

#sh5/MI300A
#export MPI4FOAM=${HOME}/openmpi-5.0.0rc12

export OMPI_CXX=clang++
export OMPI_CC=clang

export PATH=${MPI4FOAM}/bin:${MPI4FOAM}/lib:${MPI4FOAM}/include::$PATH
export LIBRARY_PATH=${MPI4FOAM}/lib:$LIBRARY_PATH

export PATH=${ROCM4FOAM}/llvm/bin:$PATH
export LD_LIBRARY_PATH=${ROCM4FOAM}/llvm/lib:${ROCM4FOAM}/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=${ROCM4FOAM}/llvm/lib:${ROCM4FOAM}/lib:$LIBRARY_PATH



#for profiler
export PATH=${ROCM4FOAM}/bin:$PATH



source etc/bashrc




