#!/bin/bash
. $WM_PROJECT_DIR/etc/bashrc
#. /home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/etc/bashrc 
#cd /home/mattijs/OpenFOAM/OpenFOAM-plus.feature-localWorld/applications/test/multiWorld/chtMultiRegionSimpleFoam
chtMultiRegionSimpleFoam -case ./solid -world solid 2>&1 | tee run_solid.log
read dummy
