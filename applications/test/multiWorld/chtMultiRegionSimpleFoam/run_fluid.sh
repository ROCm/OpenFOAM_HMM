#!/bin/bash
#. /home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/etc/bashrc 
. $WM_PROJECT_DIR/etc/bashrc
#cd /home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/applications/test/multiWorld/chtMultiRegionSimpleFoam
chtMultiRegionSimpleFoam -case ./fluid -world fluid 2>&1 | tee run_fluid.log
read dummy
