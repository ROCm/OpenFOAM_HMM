cmake_minimum_required(VERSION 2.8)

project(runTimePostProcessing)

# Set policy for CMP0002 needed for cmake > 3
cmake_policy(SET CMP0002 OLD)

if (EXISTS "$ENV{VTK_DIR}")
    message("Building with VTK from $ENV{VTK_DIR}")
    find_package(VTK REQUIRED HINTS $ENV{VTK_DIR})
    include(${VTK_USE_FILE})
elseif (EXISTS "$ENV{ParaView_DIR}")
    message("Building with Paraview from $ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED)
    include(${VTK_USE_FILE})
    set(
        VTK_VERSION
        "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION}"
    )
else ()
    message (FATAL_ERROR "VTK not found using VTK_DIR or ParaView_DIR")
endif ()

include(CMakeLists-Common.txt)


#-----------------------------------------------------------------------------
