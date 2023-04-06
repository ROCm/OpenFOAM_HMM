module load openmpi/4.0.3-ucx1.13.0-rocm5.4.0


export HSA_XNACK=1
export OMP_NUM_THREADS=4
export OMPX_DISABLE_USM_MAPS=1
export GPU_MAX_HW_QUEUES=4


#export ROCM4FOAM=${HOME}/OPENMP_COMPILER_UVM/rocm-5.5.0
export ROCM4FOAM=${HOME}/OPENMP_COMPILER_UVM/rocm-afar-001-afar
#export MPI4FOAM=${THERA_OMPI_DIR}

export MPI4FOAM=/home/lgrinber/OPENMP_COMPILER_UVM/openmpi-5.0.0rc2
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




