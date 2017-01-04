# Create a plugin to add a reader to the ParaView GUI

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

LINK_DIRECTORIES(
    $ENV{FOAM_LIBBIN}
    $ENV{FOAM_EXT_LIBBIN}
)

INCLUDE_DIRECTORIES(
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/OSspecific/$ENV{WM_OSTYPE}/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshing/blockMesh/lnInclude
    ${PROJECT_SOURCE_DIR}/../vtkPVblockMesh
)

ADD_DEFINITIONS(
    -std=c++11
    -DWM_$ENV{WM_PRECISION_OPTION}
    -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
)

# Set output library destination to plugin directory
SET(
    LIBRARY_OUTPUT_PATH $ENV{PV_PLUGIN_PATH}
    CACHE INTERNAL
    "Single output directory for building all libraries."
)


IF (PARAVIEW_QT_VERSION VERSION_GREATER "4")
  QT5_WRAP_CPP(MOC_SRCS
      pqRefreshProperty.h
      pqShowPointNumbersProperty.h
  )
ELSE()
  QT4_WRAP_CPP(MOC_SRCS
      pqRefreshProperty.h
      pqShowPointNumbersProperty.h
  )
ENDIF()

ADD_PARAVIEW_PROPERTY_WIDGET(IFACES0 IFACES0_SRCS
    TYPE "openfoam_refresh_button"
    CLASS_NAME pqRefreshProperty
)

ADD_PARAVIEW_PROPERTY_WIDGET(IFACES1 IFACES1_SRCS
    TYPE "openfoam_show_point_numbers"
    CLASS_NAME pqShowPointNumbersProperty
)

ADD_PARAVIEW_PLUGIN(
    PVblockMeshReader_SM "1.0"
    SERVER_MANAGER_XML PVblockMeshReader_SM.xml
    SERVER_MANAGER_SOURCES  vtkPVblockMeshReader.cxx
    GUI_INTERFACES
        ${IFACES0}
        ${IFACES1}
    SOURCES
        ${IFACES0_SRCS}
        ${IFACES1_SRCS}
        ${MOC_SRCS}
        pqRefreshProperty.cxx
        pqShowPointNumbersProperty.cxx
)


TARGET_LINK_LIBRARIES(
    PVblockMeshReader_SM
    LINK_PUBLIC
    vtkPVblockMesh
    blockMesh
    OpenFOAM
)

#-----------------------------------------------------------------------------
