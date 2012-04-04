/**
 * @file    glui_draws.cpp
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.1
 * @date    30/05/2006
 *
 * @brief   Octree graphical interface: drawing commands
 */
//________________________________________________


#if !defined(WIN32) || defined(__CYGWIN__)
#pragma implementation "viz_glui_defs.h"
#endif // WIN32


#include "viz_glui_defs.h"


//_____________________________________________________________________________
// declarations of this file


// display element switches
int   show_octree         = 1 ;
int   show_nodes          = 0 ;
int   show_dual           = 0 ;
int   show_iso            = 0 ;
int   show_direct_iso     = 0 ;
int   show_direct_iso_wire = 0 ;

// orthographic / perspective projection
int   ortho =   0   ;

// object transformation
float view_rotate[16] = { 1.0f,0.0f,0.0f,0.0f, 0.0f,1.0f,0.0f,0.0f, 0.0f,0.0f,1.0f,0.0f, 0.0f,0.0f,0.0f,1.0f };
float obj_pos    [3 ] = { 0.0f, 0.0f, 0.0f };


//-----------------------------------------------------------------------------
// lights
int   light0_enabled        =   1 ;
int   light1_enabled        =   0 ;
int   light0_intensity      = 100 ;
int   light1_intensity      =  60 ;
int   light0_intensity2      = 50 ;
int   light1_intensity2      = 30 ;
GLfloat light0_diffuse [ 4] = {.7f, .7f, 1.0f, 1.0f};
GLfloat light1_diffuse [ 4] = {.9f, .7f, 0.2f, 1.0f};
GLfloat light0_position[ 4] = {.5f, .5f, 10.0f, 0.0f};
GLfloat light1_position[ 4] = {-1.0f, -1.0f, 10.0f, 0.0f};
GLfloat light0_rotation[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
GLfloat light1_rotation[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };

//-----------------------------------------------------------------------------

// 3D projection matrix
void project( int o, float sx, float sy, int &tx, int &ty, int &tw, int &th ) ;

// window resizing
void reshape( int x, int y ) ;

// main drawing function
void display() ;

//_____________________________________________________________________________



//_____________________________________________________________________________
// 3D projection matrix
void project( int o, float sx, float sy, int &tx, int &ty, int &tw, int &th )
//-----------------------------------------------------------------------------
{  
  // get viewport
  GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );

  tw = (int)( tw * sx );
  th = (int)( th * sy );
  ::glViewport( tx, ty, tw, th );

  // sets the projection matrix
  ::glMatrixMode(GL_PROJECTION);
  ::glLoadIdentity();

  // sets the viewport
  if( o )
  {
    if( th > tw )
      ::glOrtho( -1.2, 1.2, -1.2*(double)th/tw, 1.2*(double)th/tw, 0.0, 10.0  ) ;
    else
      ::glOrtho( -1.2*(double)tw/th, 1.2*(double)tw/th, -1.2, 1.2, 0.0, 10.0  ) ;
  }
  else
    ::gluPerspective( 45.0, th>0?(double)tw/th:1.0, 0.5, 10.0 );

  ::gluLookAt( 0.0f,0.0f,3.0f, 0.0f,0.0f,0.5f, 0.0f,1.0f,0.0f ) ;

  // switch to the modelview matrix
  ::glMatrixMode( GL_MODELVIEW );
  ::glLoadIdentity();
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// window resizing
void reshape( int x, int y )
//-----------------------------------------------------------------------------
{
  ::glutSetWindow(main_window);
  
  int tx, ty, tw, th;
  project( ortho, 1.0f, 1.0f, tx, ty, tw, th ) ;

  // sets the window trackball
  mouse_rot.set_w( tw ) ;
  mouse_rot.set_h( th ) ;
  mouse_mv .set_w( tw ) ;
  mouse_mv .set_h( th ) ;
  mouse_zm .set_w( tw ) ;
  mouse_zm .set_h( th ) ;

  // redisplay
  ::glutPostRedisplay();
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// main drawing function
void display()
//-----------------------------------------------------------------------------
{
  ::glutSetWindow(main_window);

  // clear screen
  ::glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  if( !true ) { ::glutSwapBuffers();  return ; }

  int tx, ty, tw, th ;
  project( ortho, 1.0f, 1.0f, tx, ty, tw, th ) ;
  
  ::glEnable( GL_DEPTH_TEST );
  ::glShadeModel(GL_SMOOTH);
  ::glEnable(GL_LIGHTING);

  //::glDisable (GL_LIGHTING);
  float val[3] ;
  if ( light0_enabled )
  {
    ::glEnable( GL_LIGHT0 ) ;

    val[0] = light0_diffuse[0] * light0_intensity / 100 ;
    val[1] = light0_diffuse[1] * light0_intensity / 100 ;
    val[2] = light0_diffuse[2] * light0_intensity / 100 ;
    ::glLightfv(GL_LIGHT0, GL_DIFFUSE, val );
    
    val[0] = light0_diffuse[0] * light0_intensity2 / 100 ;
    val[1] = light0_diffuse[1] * light0_intensity2 / 100 ;
    val[2] = light0_diffuse[2] * light0_intensity2 / 100 ;
    ::glLightfv(GL_LIGHT0, GL_SPECULAR, val );    
    const GLfloat light0_ambient[4] =  {0.1f, 0.1f, 0.3f, 1.0f};
    ::glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    

    ::glLoadIdentity();
    ::glMultMatrixf( light0_rotation );
    ::glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

  }
  else 
    ::glDisable( GL_LIGHT0 ); 
  
  if ( light1_enabled )
  {
    ::glEnable( GL_LIGHT1 );

    val[0] = light1_diffuse[0] * light1_intensity / 100 ;
    val[1] = light1_diffuse[1] * light1_intensity / 100 ;
    val[2] = light1_diffuse[2] * light1_intensity / 100 ;
    ::glLightfv(GL_LIGHT1, GL_DIFFUSE, val );

    val[0] = light1_diffuse[0] * light1_intensity2 / 100 ;
    val[1] = light1_diffuse[1] * light1_intensity2 / 100 ;
    val[2] = light1_diffuse[2] * light1_intensity2 / 100 ;
    ::glLightfv(GL_LIGHT1, GL_SPECULAR, val );
    const GLfloat light1_ambient[4] =  {0.1f, 0.1f, 0.3f, 1.0f};
    ::glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
    
    ::glLoadIdentity();
    ::glMultMatrixf( light1_rotation );
    ::glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

  }
  else 
    ::glDisable( GL_LIGHT1 );  

  // transformation matrix
  ::glLoadIdentity();
  ::glTranslatef( obj_pos[0], obj_pos[1], obj_pos[2] );
  if( ortho )
  {
    ::glScalef( 1.0+obj_pos[2], 1.0+obj_pos[2], 1.0+obj_pos[2] ) ;
     ::glEnable( GL_NORMALIZE );
  }
  else
  {
    ::glDisable( GL_NORMALIZE );
  }

  ::glMultMatrixf( view_rotate );
  ::glTranslatef( -0.5f, -0.5f, -0.5f );
  PRINT_GL_DEBUG ;

//-----------------------------------------------------------------------------
// 3D display

  data_access *ref = NULL ;
  switch( impl_data )
  {
    case SWITCH_IMPL   : ref = dat_func   ; break ;
    case SWITCH_DATA   : ref = dat_data   ; break ;
  }
  if( ref )
    ::glScalef( ref->_sx, ref->_sy, ref->_sz ) ;

  // octree wireframe
  if( show_octree         )
  {
    ::glDisable( GL_TEXTURE_1D ) ;
    ::glLineWidth( 0.5 ) ;
    ::glColor3f( 0.3f, 0.3f, 0.5f ) ;
    ::glBegin( GL_LINES ) ;
    {
      octree.draw_wire() ;
    }
    ::glEnd() ;  // GL_LINES
  }


  // cells' centers
  if( show_nodes          )
  {
    ::glEnable( GL_TEXTURE_1D ) ;
    ::glPointSize( 4.0 ) ;
    ::glColor3f( 0.0f, 0.8f, 0.2f ) ;
    ::glBegin( GL_POINTS ) ;
    {
      octree.draw_centers() ;
    }
    ::glEnd() ;  // GL_POINTS
  }

  // octree dual
  if( show_dual           )
  {
    ::glDisable( GL_TEXTURE_1D ) ;
    ::glDisable(GL_LIGHTING);
    octree.draw_dual() ;	 
    ::glEnable(GL_LIGHTING);  
  }
  
  // isosurface of the field
  if( show_iso            )
  {
    ::glDisable( GL_TEXTURE_1D ) ;
    ::glDisable(GL_LIGHTING);
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
    ::glColor3f( 0.3f, 0.3f, 0.5f ) ;    
    ::glBegin( GL_TRIANGLES ) ;
    octree.draw_iso() ;	 
    ::glEnd() ; // GL_TRIANGLES

    ::glEnable(GL_LIGHTING);  
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
    ::glColor3f( 0.6f, 0.6f, 0.7f ) ;    
    ::glBegin( GL_TRIANGLES ) ;
    octree.draw_iso() ;	 
    ::glEnd() ; // GL_TRIANGLES
  }
  
  if( ref && show_direct_iso_wire )
  {
    ::glDisable( GL_TEXTURE_1D ) ;
    
    ::glDisable(GL_LIGHTING);
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
    ::glColor3f( 0.3f, 0.3f, 0.5f ) ;    
    ::glBegin( GL_TRIANGLES ) ;
    octree.direct_draw_isosurface(ref) ;	 
    ::glEnd() ; // GL_TRIANGLES
  }

  if( ref && show_direct_iso      )
  {
    ::glDisable( GL_TEXTURE_1D ) ;

    ::glEnable(GL_LIGHTING);  
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
    ::glColor3f( 0.6f, 0.6f, 0.7f ) ;    
    ::glBegin( GL_TRIANGLES ) ;
    octree.direct_draw_isosurface(ref) ;	 
    ::glEnd() ; // GL_TRIANGLES
  }
  
  PRINT_GL_DEBUG ;

  // next frame
  ::glutSwapBuffers();
}
//_____________________________________________________________________________





