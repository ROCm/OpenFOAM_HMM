/**
 * @file    glui_cmdline.cpp
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.1
 * @date    20/08/2007
 *
 * @brief   Example Graphical interface: command line parser
 */
//________________________________________________


#ifndef WIN32
#pragma implementation "viz_glui_defs.h"
#endif // WIN32

#include "viz_glui_defs.h"


//_____________________________________________________________________________
// parse command line
bool parse_command_line(int argc, char* argv[])
//-----------------------------------------------------------------------------
{
  bool quit = false ;
  
  for( int i = 1 ; i < argc ; ++i )
  {
    if     ( !strcmp( argv[i], "-lv" ) )
    {
      if( ++i != argc ) { max_level = atoi( argv[i] ) ; }
    }
    
    else if( !strcmp( argv[i], "+ortho" ) )
    {
      ortho = true ;
    }
    else if( !strcmp( argv[i], "-ortho" ) )
    {
      ortho = false ;
    }
    else if( !strcmp( argv[i], "-pos" ) )
    {
      if( ++i != argc ) { obj_pos[0] = atof( argv[i] ) ; }
      if( ++i != argc ) { obj_pos[1] = atof( argv[i] ) ; }
      if( ++i != argc ) { obj_pos[2] = atof( argv[i] ) ; }
    }
    else if( !strcmp( argv[i], "-rot" ) )
    {
      for( int j = 0 ; j < 16 ; ++j )
        if( ++i != argc ) { view_rotate[j] = atof( argv[i] ) ; }
    }
    else if( !strcmp( argv[i], "-view" ) )
    {
      control_cb( LOAD_VIEWPORT_ID ) ;
    }
    else if( !strcmp( argv[i], "-q" ) )
    {
      quit = true ;
    }
    else if( !strcmp( argv[i], "-h" ) )
    {
      printf( "usage %s [-in file.xyz] [-q]\n", argv[0] ) ;
    }
  }
  //  glui_side  ->sync_live() ;
  glui_bottom->sync_live() ;
  
  return quit ;
}
//_____________________________________________________________________________
