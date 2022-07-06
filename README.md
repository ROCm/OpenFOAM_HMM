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
OpenFOAM has the following dependencies. The installation has been tested
with the mentioned versions of the libraries, which are therefore recommended.

1. gcc-8.3.1
2. MPI (openmpi, etc.): openmpi/4.0.3 built with ucx1.8.0 support
3. boost/1.75.0
4. cmake/3.18.4
5. ROCm/4.4.0-7272
6. BLAS (openblas, etc.): openblas/0.3.15
4. Kokkos and Kokkos-Kernels (for PETSc)

## Build OpenFOAM with HMM
```
TODO Add build instructions
```

@author	: Suyash Tandon<br>
@date	: June 20, 2022<br>
<span style="color:blue">For AMD internal use only</span>
