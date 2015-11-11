cmake_minimum_required(VERSION 2.8)

project(runTimePostProcessing)

if (EXISTS $ENV{VTK_DIR})
    message("Building with VTK from $ENV{VTK_DIR}")
    find_package(VTK REQUIRED HINTS $ENV{VTK_DIR})
    include(${VTK_USE_FILE})

else (EXISTS $ENV{ParaView_DIR})
    message("Building with Paraview from $ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED)
    include(${VTK_USE_FILE})

    set(
        VTK_VERSION
        "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION}"
    )
endif (EXISTS $ENV{VTK_DIR})

include(CMakeLists-Common.txt)

#-----------------------------------------------------------------------------
