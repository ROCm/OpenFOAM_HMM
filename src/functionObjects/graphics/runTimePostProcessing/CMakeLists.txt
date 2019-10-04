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

unset(VTK_VERSION)
unset(VTK_VERSION_MAJOR)

if (EXISTS "$ENV{VTK_DIR}")
    message("Building with VTK from $ENV{VTK_DIR}")
    find_package(VTK REQUIRED HINTS $ENV{VTK_DIR})
elseif (EXISTS "$ENV{ParaView_DIR}")
    message("Building with Paraview from $ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED HINTS $ENV{ParaView_DIR})
else()
    message(FATAL_ERROR "VTK not found using VTK_DIR or ParaView_DIR")
endif()

if (VTK_USE_FILE)
    message("VTK from ${VTK_USE_FILE}")
    include(${VTK_USE_FILE})
endif()

#-----------------------------------------------------------------------------
# Handle name changes (eg, VTK_MAJOR_VERSION to VTK_VERSION_MAJOR etc.)

if (VTK_MAJOR_VERSION AND NOT VTK_VERSION_MAJOR)
    message("Found older VTK version naming")
    set(VTK_VERSION_MAJOR ${VTK_MAJOR_VERSION})
    set(VTK_VERSION_MINOR ${VTK_MINOR_VERSION})
    set(VTK_VERSION_PATCH ${VTK_BUILD_VERSION})

    if (NOT VTK_VERSION)
        set(
            VTK_VERSION
            "${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}.${VTK_VERSION_PATCH}"
        )
        message("Synthesized VTK version: " ${VTK_VERSION})
    endif()
endif()

#-----------------------------------------------------------------------------

include(CMakeLists-OpenFOAM.txt)
include(CMakeLists-Project.txt)

#-----------------------------------------------------------------------------
