#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0002 NEW)   # Policy CMP0002 required for for cmake >= 3

# Fail if not building out-of-source
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds disallowed. Use a separate build directory")
endif()

#-----------------------------------------------------------------------------
# Simple discovery and sanity checks

if (EXISTS "$ENV{VTK_DIR}")
    message("Building with VTK from $ENV{VTK_DIR}")
    find_package(VTK REQUIRED HINTS $ENV{VTK_DIR})
    include(${VTK_USE_FILE})
elseif (EXISTS "$ENV{ParaView_DIR}")
    message("Building with Paraview from $ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED HINTS $ENV{ParaView_DIR})
    include(${VTK_USE_FILE})
    set(
        VTK_VERSION
        "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION}"
    )
else()
    message(FATAL_ERROR "VTK not found using VTK_DIR or ParaView_DIR")
endif()

#-----------------------------------------------------------------------------

include(CMakeLists-OpenFOAM.txt)
include(CMakeLists-Project.txt)

#-----------------------------------------------------------------------------
