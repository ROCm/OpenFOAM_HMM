PROJECT(fastDualOctree)

cmake_minimum_required(VERSION 2.8)

#---------------------------------------------------
#----- Dependencies
#---------------------------------------------------
FIND_PACKAGE(OpenGL REQUIRED)
if ( OPENGL_FOUND )
   message( "-- OpenGL found. include=${OPENGL_INCLUDE_DIR}" )
   include_directories( ${OPENGL_INCLUDE_DIR} )
endif ( OPENGL_FOUND )

FIND_PACKAGE(ZLIB REQUIRED)
if ( ZLIB_FOUND )
   message( "-- zlib found. include=${ZLIB_INCLUDE_DIR}" )
   include_directories( ${ZLIB_INCLUDE_DIR} )
endif ( ZLIB_FOUND )

FIND_PACKAGE(PNG REQUIRED)
if ( PNG_FOUND )
   message( "-- PNG found. include=${PNG_INCLUDE_DIR}" )
   include_directories( ${PNG_INCLUDE_DIR} )
endif ( PNG_FOUND )

##-- GLUI dependency 
#FIND_PATH(GLUI_INCLUDE_DIR GL/glui.h ~/pub/system/glui-2.36/src/include)
#FIND_LIBRARY(GLUI_LIBRARIES
#      NAMES GLUI glui
#            GLUT glut
#      PATHS /usr/local/lib/
#            /usr/lib/
#            /opt/graphics/OpenGL/lib
#            /usr/openwin/lib
#            /usr/shlib
#            /usr/X11R6/lib
#            ~/pub/system/glui-2.36/src/lib
#    )
#if ( GLUI_INCLUDE_DIR )
#   message( "-- GL/glui.h found. include=${GLUI_INCLUDE_DIR}" )
#   include_directories( ${GLUI_INCLUDE_DIR} )
#else ( GLUI_INCLUDE_DIR )
#   message( "-- GL/glui.h not found." )
#endif ( GLUI_INCLUDE_DIR )
#if ( GLUI_LIBRARIES )
#   message( "-- libglui found. include=${GLUI_LIBRARIES}" )
#else( GLUI_LIBRARIES )
#   message( "-- libglui not found." )
#endif ( GLUI_LIBRARIES )

#---------------------------------------------------
# -- Environement variables
#---------------------------------------------------
if (UNIX)
  add_definitions(-DUNIX)
endif (UNIX)
if (WIN32)
  add_definitions(-DWIN32)
endif (WIN32)
if (APPLE)
  add_definitions(-DAPPLE)
endif (APPLE)

#---------------------------------------------------
# -- SRCs
#---------------------------------------------------

SET(Perf_SRC ${Perf_SRC}
  data_access.cpp  
  hash_octree.cpp    
  mc_draw.cpp     
  opt_octree.cpp  
  fparser.cpp      
  leaf_octree.cpp    
  mem_octree.cpp  
  fpoptimizer.cpp  
  MarchingCubes.cpp  
  morton.cpp      
  ptr_octree.cpp
)

#SET(VIZ_SRC ${VIZ_SRC}
# viz_glui_cmdline.cpp  
# viz_glui_controls.cpp  
# viz_glui_draws.cpp  
# viz_glui_export.cpp  
# viz_glui_mouse.cpp 
#)


add_executable(perf_main  perf_main ${Perf_SRC})
target_link_libraries(perf_main ${OPENGL_LIBRARIES} ${ZLIB_LIBRARIES})


#add_executable(viz_glui_main viz_glui_main ${Perf_SRC} ${VIZ_SRC} )
#target_link_libraries(viz_glui_main ${OPENGL_LIBRARIES} ${ZLIB_LIBRARIES} ${PNG_LIBRARIES} ${GLUI_LIBRARIES})
