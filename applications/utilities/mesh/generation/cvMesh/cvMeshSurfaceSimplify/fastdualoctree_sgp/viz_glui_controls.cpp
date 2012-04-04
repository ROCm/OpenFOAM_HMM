/**
 * @file    glui_controls.cpp
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.1
 * @date    30/05/2006
 *
 * @brief   Octree graphical interface: interface controls
 */
//________________________________________________


#if !defined(WIN32) || defined(__CYGWIN__)
#pragma implementation "viz_glui_defs.h"
#endif // WIN32


#include "viz_glui_defs.h"
#include "implfuns.h"


//_____________________________________________________________________________
// declarations of this file

// main object: octree
Octree  octree  ;

// switch between implicit funciont or volumetric data
int impl_data = SWITCH_DATA ;

// octree implicit data wrapper
data_func *dat_func = NULL ;

// octree data file wrapper
data_file  *dat_data  = NULL ;


// maximal level of the tree
int   max_level  = 4;

// isovalue
float iso_val    = 0.0 ;

// curvature threshold
float curv_thres = 0.5 ;

// selected implicit function
int curr_fun = 0 ;


//-----------------------------------------------------------------------------

// main glui class: (right) side panel
GLUI *glui_side   = NULL ;

// bottom panel
GLUI *glui_bottom = NULL ;

// name of the import file
GLUI_FileBrowser *file_browser  ;

// iso file
GLUI_String filename  ;

// formula of the implicit function
GLUI_String impl_fun     ;

//-----------------------------------------------------------------------------

// control events callback
void control_cb( int control ) ;

// create side panel
void create_side_panel() ;

// create bottom panel
void create_bottom_panel() ;

/// set file extension
int set_ext( const char ext[4] ) ;

//_____________________________________________________________________________




//_____________________________________________________________________________
// control events callback
void control_cb( int control )
//-----------------------------------------------------------------------------
{
  data_access *ref = NULL ;
  switch( impl_data )
  {
    case SWITCH_IMPL   : ref = dat_func   ; break ;
    case SWITCH_DATA   : ref = dat_data   ; break ;
  }
  
  switch( control )
  {
    case SET_IMPL_DATA_ID :
      switch( impl_data )
      {
      case SWITCH_IMPL   : delete dat_func   ;  dat_func   = new data_func      ( max_level, iso_val, impl_fun ) ;  break ;
      case SWITCH_DATA   : delete dat_data   ;  dat_data   = new data_file      ( max_level, iso_val, filename.c_str() ) ;  break ;
      } 
      break ;

    case SET_VAL_ID :
      if( !ref ) break ;
      ref->_max_level = max_level ;
      ref->_iso_val = iso_val ;
      octree.set_impl( ref ) ;
      octree.check() ;
      break ;

    case REFINE_ID        :
      if( !ref ) break ;
      ref->_max_level = max_level ;
      ref->_iso_val = iso_val ;
      octree.clear() ;
      octree.refine( ref ) ;
      octree.set_impl( ref ) ;
      octree.check() ;
      break ;
      
    case ADAPT_ID        :
      if( !ref ) break ;
      ref->_max_level = max_level ;
      ref->_iso_val   = iso_val ;
      octree.adapt( ref ) ;
      octree.set_impl( ref ) ;
      octree.check() ;
      break ;

    case ISO_BUILD_ID :
      if( !ref ) break ;
      ref->_iso_val = iso_val ;
      octree.build_isosurface(ref) ;
      octree.mc().writeOFF("iso.off") ;
      octree.check() ;
      break ;
      
      // set implicit function
    case FUN_ID   :
      if( curr_fun == 0 ) break ;
      impl_fun.clear() ;
      impl_fun = fun_def[curr_fun] ;
      control_cb( SET_IMPL_DATA_ID ) ;
      control_cb( ADAPT_ID ) ;
      glui_bottom->sync_live() ;
      break ;
      
      
      // load/save viewpoint
    case SAVE_VIEWPORT_ID   : save_viewport() ;  break ;
      
    case LOAD_VIEWPORT_ID   : load_viewport() ;  glui_bottom->sync_live() ;  break ;
      
      // reset rotation
    case RESET_ROTATION_ID    :
      view_rotate[ 0] = view_rotate[ 5] = view_rotate[10] = view_rotate[15] = 1.0f ;
      view_rotate[ 1] = view_rotate[ 2] = view_rotate[ 3] = view_rotate[ 4] = 0.0f ;
      view_rotate[ 6] = view_rotate[ 7] = view_rotate[ 8] = view_rotate[ 9] = 0.0f ;
      view_rotate[11] = view_rotate[12] = view_rotate[13] = view_rotate[14] = 0.0f ;
      break ;
      
      // reset translation
    case RESET_TRANSLATION_ID :  obj_pos[0] = obj_pos[1] = 0.0f ;  break ;
      // reset zoom
    case RESET_ZOOM_ID        :  obj_pos[2] = 0.0f ;  break ;
      
      // changed a parameter
    case REDRAW_ID          :  break ;
      
      // changed a parameter
    case FILENAME_ID        :  filename = file_browser->current_dir + file_browser->get_file() ;
      
      // orthographic/perspective projection
    case PROJ_ID            :  reshape   (0,0) ;  break ;
      
    case EXIT_ID            : 
      delete dat_func   ;
      delete dat_data   ;
      octree.clear() ;
      exit(0) ;
      break ;
      
    default :  break ;
  }
  
  ::glutPostRedisplay();
}
//_____________________________________________________________________________




//_____________________________________________________________________________
//_____________________________________________________________________________



//_____________________________________________________________________________
// create side panel
void create_side_panel()
//-----------------------------------------------------------------------------
{
  glui_side = NULL ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// create bottom panel
void create_bottom_panel()
//-----------------------------------------------------------------------------
{
  GLUI_Rollout     *roll ;
  GLUI_EditText    *text ;
  GLUI_Rotation    *rot  ;
  GLUI_Translation *trans;
  GLUI_RadioGroup  *radio;
  GLUI_Listbox     *list ;
  GLUI_Scrollbar   *sb   ;
  
  glui_bottom = GLUI_Master.create_glui("controls") ; //create_glui_subwindow( main_window, GLUI_SUBWINDOW_BOTTOM );
  
  //--------------------------------------------------//

  file_browser = new GLUI_FileBrowser( glui_bottom, "iso file", false, FILENAME_ID, control_cb ) ;
  file_browser->current_dir =  "data" ;
  file_browser->fbreaddir( "data" ) ;

  filename = "data/engine.iso.gz";
  text = glui_bottom ->add_edittext( "iso file"  , filename ) ;
  text->set_w( 200 ) ;

  list = glui_bottom->add_listbox( "Implicit Functions:", &curr_fun, FUN_ID, control_cb );
  for( int i=0; i<NFUNS; i++ ) list->add_item( i, fun_list[i] );
  
  impl_fun = "";
  text = glui_bottom ->add_edittext( "impl fun"  , impl_fun ) ;
  text->set_w( 200 ) ;
  

  //--------------------------------------------------//
  // 
  roll  = glui_bottom ->add_rollout( "iso value", true );
  sb = new GLUI_Scrollbar( roll, "iso value in [-1,1]",GLUI_SCROLL_HORIZONTAL, &iso_val, ADAPT_ID,control_cb);
  sb->set_float_limits(-1,1);
  sb->set_speed( .005f );
  
  roll  = glui_bottom ->add_rollout( "level in [2,7]", true );
  sb = new GLUI_Scrollbar( roll, "level", GLUI_SCROLL_HORIZONTAL, &max_level, ADAPT_ID, control_cb);
  sb->set_int_limits(2,7);
  sb->set_speed( .005f );

  //--------------------------------------------------//  
  glui_bottom->add_column( true );
  radio = glui_bottom->add_radiogroup( &impl_data, SET_IMPL_DATA_ID, control_cb ) ;
  glui_bottom->add_radiobutton_to_group( radio, "Implicit function" ) ;
  glui_bottom->add_radiobutton_to_group( radio, "Direct data file" ) ;

  //--------------------------------------------------//
  
  
  glui_bottom ->add_button  ( "load"      , SET_IMPL_DATA_ID, control_cb ) ;
  glui_bottom ->add_button  ( "set values", SET_VAL_ID   , control_cb ) ;
  glui_bottom ->add_button  ( "refine"    , REFINE_ID, control_cb ) ;
  glui_bottom ->add_button  ( "adapt"     , ADAPT_ID , control_cb ) ;
  
  
  //--------------------------------------------------//
  // display element
  glui_bottom->add_checkbox( "ortho"    , &ortho               , PROJ_ID  , control_cb );
  glui_bottom->add_checkbox( "octree"   , &show_octree         , REDRAW_ID, control_cb );
  glui_bottom->add_checkbox( "nodes"    , &show_nodes          , REDRAW_ID, control_cb );
  glui_bottom->add_checkbox( "dual"     , &show_dual           , REDRAW_ID, control_cb );
  glui_bottom->add_checkbox( "surface"  , &show_iso            , REDRAW_ID, control_cb );
  glui_bottom->add_checkbox( "dir surf" , &show_direct_iso     , SET_VAL_ID, control_cb );
  glui_bottom->add_checkbox( "dir sur w", &show_direct_iso_wire, SET_VAL_ID, control_cb );
  

  glui_bottom ->add_button  ( "Isosurface", ISO_BUILD_ID, control_cb ) ;
  
//  glui_bottom->add_button( "Redraw", REDRAW_ID, control_cb ) ;
  
  glui_bottom ->add_button( "Open View", LOAD_VIEWPORT_ID, control_cb ) ;
  glui_bottom ->add_button( "Save View", SAVE_VIEWPORT_ID, control_cb ) ;
  glui_bottom ->add_button( "Quit", EXIT_ID, control_cb );
  glui_bottom->add_column( true );
  
  
  glui_bottom->add_column( true );
  
  //--------------------------------------------------//
  // position
  roll  = glui_bottom ->add_rollout( "3D position", true );
  objects_rot = glui_bottom->add_rotation_to_panel( roll, "Objects", view_rotate );
  objects_rot->set_spin( 1.0f );
  glui_bottom->add_button_to_panel( roll, "Reset", RESET_ROTATION_ID, control_cb ) ;
  glui_bottom->add_column_to_panel( roll, false );
  
  objects_mv = glui_bottom->add_translation_to_panel( roll, "Objects XY", GLUI_TRANSLATION_XY, obj_pos );
  objects_mv->set_speed( .005f );
  glui_bottom->add_button_to_panel( roll, "Reset", RESET_TRANSLATION_ID, control_cb ) ;
  glui_bottom->add_column_to_panel( roll, false );
  
  objects_zm = glui_bottom->add_translation_to_panel( roll, "Objects Z", GLUI_TRANSLATION_Z, &obj_pos[2] );
  objects_zm->set_speed( .005f );
  glui_bottom->add_button_to_panel( roll, "Reset", RESET_ZOOM_ID, control_cb ) ;
  
  
  //--------------------------------------------------//
  // Lights
  roll = glui_bottom->add_rollout( "Lights", true );
  rot  = glui_bottom->add_rotation_to_panel( roll, "Blue Light", light0_rotation );
  rot->set_spin( .82f );
  
  rot  = glui_bottom->add_rotation_to_panel( roll, "Orange Light", light1_rotation );
  rot->set_spin( .82f );
  
  
  glui_bottom->add_column_to_panel( roll, false );
  
  glui_bottom->add_checkbox_to_panel( roll, "On", &light0_enabled );
  text = glui_bottom->add_edittext_to_panel( roll, "Id (%)", GLUI_EDITTEXT_INT, &light0_intensity );
  text->set_int_limits( 0, 100 );
  text->set_w(4) ;
  text = glui_bottom->add_edittext_to_panel( roll, "Is (%)", GLUI_EDITTEXT_INT, &light0_intensity2 );
  text->set_int_limits( 0, 100 );
  text->set_w(4) ;
  
  glui_bottom->add_checkbox_to_panel( roll, "On", &light1_enabled );
  text = glui_bottom->add_edittext_to_panel( roll, "Id (%)", GLUI_EDITTEXT_INT, &light1_intensity );
  text->set_int_limits( 0, 100 );
  text->set_w(4) ;
  text = glui_bottom->add_edittext_to_panel( roll, "Is (%)", GLUI_EDITTEXT_INT, &light1_intensity2 );
  text->set_int_limits( 0, 100 );
  text->set_w(4) ;
}
//_____________________________________________________________________________
