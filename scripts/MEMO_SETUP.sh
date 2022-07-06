export HSA_XNACK=1
export OMP_NUM_THREADS=2

source etc/bashrc


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/lgrinber/rocm-afar001-460/llvm/lib
export PATH=/work/lgrinber/rocm-afar001-460/llvm/bin:$PATH

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/lgrinber/rocm-afar001-432/llvm/lib
#export PATH=/work/lgrinber/rocm-afar001-432/llvm/bin:$PATH
#export PATH=/work/lgrinber/MPI_CLANG/openmpi-5.0.0rc2/build/bin:$PATH

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/lgrinber/OPENMP_COMPILER_UVM/aomp_15.0-1/lib
#export PATH=/work/lgrinber/OPENMP_COMPILER_UVM/aomp_15.0-1/bin:/work/lgrinber/OPENMP_COMPILER_UVM/aomp_15.0-1/include:$PATH
#export PATH=/work/lgrinber/MPI_CLANG/openmpi-5.0.0rc2/build/bin:$PATH


export PATH=$PATH:/work/lgrinber/MPI_CLANG/openmpi-5.0.0rc2/build/bin:/work/lgrinber/MPI_CLANG/openmpi-5.0.0rc2/build/lib
export LIBRARY_PATH=/work/lgrinber/MPI_CLANG/openmpi-5.0.0rc2/build/lib:$LIBRARY_PATH

#if fails to execute try:
#export LD_LIBRARY_PATH=/opt/rocm-5.0.2/hsa/lib:$LD_LIBRARY_PATH


