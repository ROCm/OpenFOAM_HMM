#!/bin/bash

cd ${0%/*} || exit 1    # run from this directory

#cp rawSurfaces/MRF.stl MRF.stl
#cp rawSurfaces/baffles.stl baffles.stl

cp rawSurfaces/vessel.stl vessel.stl
cp rawSurfaces/sparger.stl sparger.stl
cp rawSurfaces/shaftRotating.stl shaftRotating.stl
cp rawSurfaces/shaftStatic.stl shaftStatic.stl
cp rawSurfaces/gasInlet.stl gasInlet.stl
cp rawSurfaces/stirrer.stl stirrer.stl

# Vessel surface
surfaceAutoPatch vessel.stl vessel.stl 120

# Sparger
surfaceCheck sparger.stl
surfaceAdd gasInlet.stl sparger_0.obj spargerInlet.stl
surfaceConvert sparger_1.obj spargerShaft.stl
surfaceOrient -inside spargerInlet.stl "(-100 -50 500)" spargerInlet.stl
surfaceOrient -inside spargerShaft.stl "(-50 -20 -100)" spargerShaft.stl

# Rotating shaft
surfaceOrient -inside shaftRotating.stl "(-100 -50 500)" shaftRotating.stl

# Static shaft
surfaceOrient -inside shaftStatic.stl "(15 -200 1000)" shaftStatic.stl

# Stirrer
surfaceSplitByTopology stirrer.stl stirrer.stl
surfaceOrient -inside stirrer.stl "(15 -200 1000)" stirrer.stl


# ----------------------------------------------------------------- end-of-file
