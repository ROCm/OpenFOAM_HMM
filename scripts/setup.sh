

export HSA_XNACK=1
export OMP_NUM_THREADS=1
export OMPX_DISABLE_USM_MAPS=1
export GPU_MAX_HW_QUEUES=4

#set ROCm and MPI paths
export ROCM4FOAM=${ROCM_PATH:-/opt/rocm}
export MPI4FOAM=${MPI_PATH:-/opt/ompi}
export UMPIRE4FOAM=${UMPIRE_PATH:-/opt/umpire-6.0.0}

#set OMPI compilers
export OMPI_CXX=clang++
export OMPI_CC=clang

#add OMPI and ROCm to PATH and LD_LIBRARY_PATH
export PATH=${MPI4FOAM}/bin:$PATH
export LD_LIBRARY_PATH=${MPI4FOAM}/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=${MPI4FOAM}/lib:$LIBRARY_PATH

export PATH=${ROCM4FOAM}/bin:${ROCM4FOAM}/llvm/bin:$PATH
export LD_LIBRARY_PATH=${ROCM4FOAM}/llvm/lib:${ROCM4FOAM}/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=${ROCM4FOAM}/llvm/lib:${ROCM4FOAM}/lib:$LIBRARY_PATH

#source OpenFOAM environment
source ./etc/bashrc




