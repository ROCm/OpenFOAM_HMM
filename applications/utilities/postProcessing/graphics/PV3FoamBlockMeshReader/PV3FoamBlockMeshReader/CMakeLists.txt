# create a plugin that adds a reader to the ParaView GUI
# it is added in the file dialog when doing opens/saves.

# The qrc file is processed by Qt's resource compiler (rcc)
# the qrc file must have a resource prefix of "/ParaViewResources"
# and ParaView will read anything contained under that prefix
# the pqReader.xml file contains xml defining readers with their
# file extensions and descriptions.

cmake_minimum_required(VERSION 2.4)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

LINK_DIRECTORIES(
    $ENV{FOAM_LIBBIN}
)

INCLUDE_DIRECTORIES(
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshing/blockMesh/lnInclude
    ${PROJECT_SOURCE_DIR}/../vtkPV3FoamBlockMesh
)

ADD_DEFINITIONS(
    -DWM_$ENV{WM_PRECISION_OPTION}
)

# Set output library destination to plugin folder
SET(
    LIBRARY_OUTPUT_PATH $ENV{PV_PLUGIN_PATH}
    CACHE INTERNAL
    "Single output directory for building all libraries."
)

# Build the server-side plugin
ADD_PARAVIEW_PLUGIN(
    PV3FoamBlockMeshReader_SM "1.0"
    SERVER_MANAGER_XML PV3FoamBlockMeshReader_SM.xml
    SERVER_MANAGER_SOURCES vtkPV3FoamBlockMeshReader.cxx
)

# Build the client-side plugin
ADD_PARAVIEW_PLUGIN(
    PV3FoamBlockMeshReader
    "1.0"
    GUI_RESOURCES PV3FoamBlockMeshReader.qrc
)

TARGET_LINK_LIBRARIES(
    PV3FoamBlockMeshReader_SM
    OpenFOAM
    blockMesh
    vtkPV3FoamBlockMesh
)
#-----------------------------------------------------------------------------
