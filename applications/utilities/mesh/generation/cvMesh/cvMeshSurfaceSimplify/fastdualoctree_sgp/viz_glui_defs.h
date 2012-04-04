/**
 * @file    glui_defs.h
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.1
 * @date    30/05/2006
 *
 * @brief   Octree graphical interface
 */
//________________________________________________


#pragma once

#if !defined(WIN32) || defined(__CYGWIN__)
#pragma interface
#endif // WIN32


#include <GL/glui.h>      // openGL user interface


#define OCTREE_SWITCH 4

//_____________________________________________________________________________
// octree switch
#ifdef OCTREE_SWITCH
# if   OCTREE_SWITCH==0
#  define OCTREE_PTR   1
# elif OCTREE_SWITCH==1
#  define OCTREE_HASH  1
# elif OCTREE_SWITCH==2
#  define OCTREE_OPT   1
# elif OCTREE_SWITCH==3
#  define OCTREE_LEAF  1
# elif OCTREE_SWITCH==4
#  define OCTREE_MEM   1
# endif // OCTREE_SWITCH
#endif // OCTREE_SWITCH


// switch values between pointer and hash octree
#if !OCTREE_PTR && !OCTREE_HASH && !OCTREE_OPT && !OCTREE_LEAF && !OCTREE_MEM
// # define OCTREE_PTR  1
// # define OCTREE_HASH 1
# define OCTREE_OPT  1
// # define OCTREE_LEAF 1
// # define OCTREE_MEM  1
#endif // !OCTREE_PTR && !OCTREE_HASH && !OCTREE_OPT && !OCTREE_LEAF && !OCTREE_MEM


#ifdef OCTREE_PTR
# undef OCTREE_HASH
# undef OCTREE_OPT  
# undef OCTREE_LEAF 
# undef OCTREE_MEM  
# define OCTREE_STRING "pointer octree"
# include "ptr_octree.h"
/// octree type
typedef PtrOctree Octree ;
#endif //OCTREE_PTR


#ifdef OCTREE_HASH
# undef OCTREE_PTR
# undef OCTREE_OPT  
# undef OCTREE_LEAF 
# undef OCTREE_MEM  
# define OCTREE_STRING "hash octree"
# include "hash_octree.h"
/// octree type
typedef HashOctree Octree ;
#endif //OCTREE_HASH


#ifdef OCTREE_OPT
# undef OCTREE_PTR
# undef OCTREE_HASH
# undef OCTREE_LEAF 
# undef OCTREE_MEM  
# define OCTREE_STRING "optimized octree"
# include "opt_octree.h"
/// octree type
typedef OptOctree Octree ;
#endif //OCTREE_OPT


#ifdef OCTREE_LEAF
# undef OCTREE_PTR
# undef OCTREE_HASH
# undef OCTREE_OPT  
# undef OCTREE_MEM  
# define OCTREE_STRING "leaf octree"
# include "leaf_octree.h"
/// octree type
typedef LeafOctree Octree ;
#endif //OCTREE_LEAF


#ifdef OCTREE_MEM
# undef OCTREE_PTR
# undef OCTREE_HASH
# undef OCTREE_OPT  
# undef OCTREE_LEAF 
# define OCTREE_STRING "mem octree"
# include "mem_octree.h"
/// octree type
typedef MemOctree Octree ;
#endif //OCTREE_MEM



#ifdef _DEBUG
#define PRINT_GL_DEBUG  { if( ::glGetError() != GL_NO_ERROR ) printf( "openGL watch at line %d: %s\n", __LINE__, ::gluErrorString( ::glGetError() ) ) ; }
#else  // _DEBUG
#define PRINT_GL_DEBUG  {}
#endif // _DEBUG



//_____________________________________________________________________________
// Main objects

  /// main object: octree
  extern Octree  octree  ;

  /// switch between implicit funciont or volumetric data
  extern int        impl_data   ;

  /// switch values between implicit funciont or volumetric data
  enum { SWITCH_IMPL = 0, SWITCH_DATA = 1 } ;

  /// octree implicit data wrapper
  extern data_func *dat_func ;

  /// octree data file wrapper
  extern data_file *dat_data  ;

  /// maximal level of the tree
  extern int   max_level ;

  /// isovalue
  extern float iso_val ;


//-----------------------------------------------------------------------------
// Display elements

  /// display element switch: octree wireframe
  extern int   show_octree        ;
  /// display element switch: octree cells' centers
  extern int   show_nodes         ;
  /// display element switch: octree dual
  extern int   show_dual          ;  
  /// display element switch: isosurface of the field
  extern int   show_iso           ;  
  /// display element switch: direct drawing isosurface of the field
  extern int   show_direct_iso    ;
  /// display element switch: direct drawing wire isosurface of the field
  extern int   show_direct_iso_wire ;  


//-----------------------------------------------------------------------------
// GLUI windows

  /// main window id
  extern int  main_window ;

  /// main glui class: (right) side panel
  extern GLUI *glui_side   ;

  /// bottom panel
  extern GLUI *glui_bottom ;


/// create side panel
void create_side_panel() ;

/// create bottom panel
void create_bottom_panel() ;

/// control events callback
void control_cb( int control ) ;

/// parse command line
bool parse_command_line(int argc, char* argv[]) ;


//-----------------------------------------------------------------------------
// Lights

  /// enable blue light
  extern int   light0_enabled   ;

  /// enable orange light
  extern int   light1_enabled   ;

  /// blue light intensity
  extern int   light0_intensity ;

  /// orange light intensity
  extern int   light1_intensity ;

  /// blue light intensity (specular)
  extern int   light0_intensity2 ;

  /// orange light intensity (specular)
  extern int   light1_intensity2 ;

  /// blue light diffuse color
  extern float light0_diffuse[4] ;

  /// orange light diffuse color
  extern float light1_diffuse[4] ;

  /// blue light position
  extern float light0_rotation[16] ;

  /// orange light position
  extern float light1_rotation[16] ;


//-----------------------------------------------------------------------------
// mouse and object movements

  /// motion type (-1 -> no motion, 0 -> rotate, 1 -> zoom, 2 -> translate)
  extern int  motion_type ;

  /// window trackball
  extern GLUI_Rotation     mouse_rot   ;
  /// panel trackball
  extern GLUI_Rotation    *objects_rot ;
  /// window translation
  extern GLUI_Translation  mouse_mv    ;
  /// panel translation
  extern GLUI_Translation *objects_mv  ;
  /// window zoom
  extern GLUI_Translation  mouse_zm    ;
  /// panel zoom
  extern GLUI_Translation *objects_zm  ;

  /// number of calls for updating the GLUI control
  extern int ncalls ;

  /// export movie switch
  extern int export_movie ;

  /// automatic rotation
  extern int auto_rotate  ;

  /// automatic x translation
  extern int auto_translate_x ;

  /// automatic y translation
  extern int auto_translate_y ;

  /// automatic zoom
  extern int auto_zoom ;

  /// automatic plane
  extern int auto_translate_d ;

  /// mouse x position
  extern int mouse_x ;

  /// mouse y position
  extern int mouse_y ;

/// init mouse and window controls
void init_trackballs() ;

/// mouse events tracking
void mouse(int button, int button_state, int x, int y ) ;

/// mouse motion tracking
void motion(int x, int y ) ;

/// keyboard events
void keyboard(unsigned char key, int x, int y);

//-----------------------------------------------------------------------------
// implicit function

  /// implicit function formula
  extern GLUI_String impl_fun  ;


//-----------------------------------------------------------------------------
// i/o filenames

  /// iso file
  extern GLUI_String filename  ;

  /// file browser
  extern GLUI_FileBrowser *file_browser ;



//-----------------------------------------------------------------------------
// drawing parameters

  /// orthographic / perspective projection switch
  extern int ortho ;
  /// object rotation
  extern float view_rotate[16] ;
  /// object translation
  extern float obj_pos    [3 ] ;

  /// point size
  extern int point_size;

/// 3D projection matrix
void project( int o, float sx, float sy, int &tx, int &ty, int &tw, int &th ) ;

/// window resizing
void reshape( int x, int y ) ;

/// main drawing function
void display() ;

/// print screen
void export_png() ;

/// load viewport
void load_viewport() ;
/// save viewport
void save_viewport() ;


//_____________________________________________________________________________




//_____________________________________________________________________________
/// Callback ids
enum
{
  SET_IMPL_DATA_ID     ,
  SET_VAL_ID           ,
  REFINE_ID            ,
  ADAPT_ID             ,
  ISO_BUILD_ID         ,
  SEARCH_ID            ,
  FUN_ID               ,

  SAVE_VIEWPORT_ID     ,
  LOAD_VIEWPORT_ID     ,

  RESET_ROTATION_ID    ,
  RESET_TRANSLATION_ID ,
  RESET_ZOOM_ID        ,

  REDRAW_ID            ,
  FILENAME_ID          ,
  PROJ_ID              ,  

  EXIT_ID              
};
//_____________________________________________________________________________

