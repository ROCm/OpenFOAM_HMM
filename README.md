# OpenFOAM with HMM and OpenMP offloading
---
```
#------------------------------------------------------------------------------
# =========               |
# \\      /  F ield       | OpenFOAM: The Open Source CFD Toolbox  
#  \\    /   O peration   |
#   \\  /    A nd         |www.openfoam.com
#    \\/     M anipulation|
#------------------------------------------------------------------------------
```
To run [OpenFOAM](https://www.openfoam.com) on A+A systems, we demonstrate the use of OpenMP offloading with HMM. Follow the steps below to build OpenFOAM with HMM.

## Requirements

on the ACP cloud all the requirements have been already satisfied, so no need to install any additional software.

In general, OpenFOAM has the following dependencies. The installation has been tested
with the mentioned versions of the libraries, which are therefore recommended.

1. gcc-8.3.1
2. MPI (openmpi, etc.): openmpi/4.0.3 built with ucx1.8.0 support
3. boost/1.75.0
4. cmake/3.18.4
5. ROCm
6. BLAS (openblas, etc.): openblas/0.3.15
4. Kokkos and Kokkos-Kernels (for PETSc)

## Build OpenFOAM with HMM
```
instractions for compiling OpenFOAM on the ACP cloud

1. after clonnning OpenFOAM load the following modules:

module load rocm-afar/001-732 openmpi/5.0.0-rocm-afar001-732

2. change directory to OpenFOAM-v2112 and execute:

source ./scripts/MEMO_SETUP.sh

3. compile :
   ./Allwmake -l -q -j 32

#------------------------------------------------------------------------------

```
Example of simulation using 8 MPI ranks

mpirun -np 8 ./helper.sh pisoFoam -parallel

Here the role of the helper.sh script is to assign a GPU to each MPI rank.
make sure you set 
chmod +x  helper.sh
```


profiling with rocprof
1. make sure to uncomment the following line in the helper.sh file
    #eval "rocprof -d ${outdir} -o ${outdir}/${outfile}  $*" 

2. execute the code using :
   mpirun -np 8 ./helper.sh --roctx-trace --hsa-trace  pisoFoam -parallel
   
3. at the end of the execution you will see a directory created for each MPI rank
   time lines of execution are saved in *json file inside each directory
   you can visualize time line corrponding to each MPI rank separately using 
   for example   https://ui.perfetto.dev/
   it is recomended to 
   gzip *json file and load the compressed file into Perfetto.  Perfetto will gunzip it automatically.
   
 4.It is possible to merge the output from rocprof tracing from all MPI ranks into a single merged output
   to do so :
   create a new directory , for example   "mkdir MERGE", then execute the following command:
   
   /global/software/rocm-afar001-732/libexec/rocprofiler/merge_traces.sh -o MERGED <list-of-files>
   
   for <the-of-files>  you can list space separated list of all the files you want to merge , or you can use 
   something like:
   
   /global/software/rocm-afar001-732/libexec/rocprofiler/merge_traces.sh -o MERGED openfoam.*.12*
    
   if you want to execute without profiling - do not forget to comment  the following line in the helper.sh 
   
    eval "rocprof -d ${outdir} -o ${outdir}/${outfile}  $*"
   
   ```
   
   
    




@author	: Suyash Tandon, Leopold Grinberg<br>
@date	: June 20, 2022<br>
<span style="color:blue">For AMD internal use only</span>
