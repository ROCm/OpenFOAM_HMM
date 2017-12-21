# Create a plugin to add a reader to the ParaView GUI

cmake_minimum_required(VERSION 2.8)

find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})

link_directories(
    $ENV{FOAM_LIBBIN}
)

include_directories(
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/OSspecific/$ENV{WM_OSTYPE}/lnInclude
    $ENV{WM_PROJECT_DIR}/src/conversion/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteArea/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    ${PROJECT_SOURCE_DIR}/../foamPv
    ${PROJECT_SOURCE_DIR}/../vtkPVFoam
)

add_definitions(
    -std=c++11
    -DWM_$ENV{WM_PRECISION_OPTION}
    -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
)

# Set output library destination to plugin directory
set(
    LIBRARY_OUTPUT_PATH $ENV{PV_PLUGIN_PATH}
    CACHE INTERNAL
    "Single output directory for building all libraries."
)


if (PARAVIEW_QT_VERSION VERSION_GREATER "4")
  qt5_wrap_cpp(MOC_SRCS
      pqFoamReaderControls.h
  )
else()
  qt4_wrap_cpp(MOC_SRCS
      pqFoamReaderControls.h
  )
endif()

add_paraview_property_group_widget(IFACES0 IFACES0_SRCS
    TYPE "openfoam_reader_general_controls"
    CLASS_NAME pqFoamReaderControls
)

add_paraview_plugin(
    PVFoamReader_SM "1.0"
    SERVER_MANAGER_XML PVFoamReader_SM.xml
    SERVER_MANAGER_SOURCES vtkPVFoamReader.cxx
    GUI_INTERFACES
        ${IFACES0}
    SOURCES
        ${IFACES0_SRCS}
        ${MOC_SRCS}
        pqFoamReaderControls.cxx
)


target_link_libraries(
    PVFoamReader_SM
    LINK_PUBLIC
    vtkPVFoam-pv${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}
    foamPv-pv${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}
    conversion
    finiteArea
    finiteVolume
    OpenFOAM
)

#-----------------------------------------------------------------------------
