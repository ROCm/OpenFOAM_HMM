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

# Qualify OpenFOAM/ParaView libraries with -pvMAJ.MIN
# The pv_api can also be used for paraview-MAJ.MIN
if (DEFINED ENV{PARAVIEW_API})
    set(pv_api "$ENV{PARAVIEW_API}")
else()
    if (${ParaView_VERSION} GREATER_EQUAL 5.7)
        string(REGEX MATCH [0-9]+[.][0-9]+ pv_api ${ParaView_VERSION})
        message("The PARAVIEW_API env was unset, "
            "extracted value from version="
            "${ParaView_VERSION} is ${pv_api}")
    else()
        set(pv_api "${PARAVIEW_VERSION_MAJOR}.${ParaView_VERSION_MINOR}")
        message("The PARAVIEW_API env was unset, "
            "assembled from maj/min="
            "${ParaView_VERSION_MAJOR}/${ParaView_VERSION_MINOR} is ${pv_api}")
    endif()
endif()

set(foamPvLibQualifier "-pv${pv_api}")
message("OpenFOAM libraries qualified with ${foamPvLibQualifier}")

#-----------------------------------------------------------------------------

include(CMakeLists-OpenFOAM.txt)

if (${ParaView_VERSION} GREATER_EQUAL 5.7)
  include(CMakeLists-Project57.txt)
else()
  include(${PARAVIEW_USE_FILE})
  include(CMakeLists-Project.txt)
endif()

#-----------------------------------------------------------------------------
