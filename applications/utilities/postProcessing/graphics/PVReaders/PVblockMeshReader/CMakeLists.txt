#-----------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0002 NEW)   # Policy CMP0002 required for for cmake >= 3

# Fail if not building out-of-source
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds disallowed. Use a separate build directory")
endif()

#-----------------------------------------------------------------------------
# Simple discovery and sanity checks

find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})

#-----------------------------------------------------------------------------

include(CMakeLists-OpenFOAM.txt)
include(CMakeLists-Project.txt)

#-----------------------------------------------------------------------------
