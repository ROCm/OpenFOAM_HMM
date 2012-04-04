/**
 * @file    glui_export.cpp
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.1
 * @date    20/08/2007
 *
 * @brief   Example Graphical interface: screen image exporatation
 */
//________________________________________________


#ifndef WIN32
#pragma implementation "viz_glui_defs.h"
#endif // WIN32

#include <png.h>
#include <stdio.h>
#include "viz_glui_defs.h"



//_____________________________________________________________________________
// declarations of this file

// set file extension of filename
int  set_ext( const char ext[3], char *fn ) ;

// PPM export
void export_ppm() ;

// PNG image export
void export_png() ;

// Viewport save
void save_viewport() ;

// Viewport load
void load_viewport() ;

//_____________________________________________________________________________



//_____________________________________________________________________________
// set file extension of filename
int set_ext( const char ext[3], char *fn )
//-----------------------------------------------------------------------------
{
  strcpy( fn, filename.c_str() ) ;
  int l = (int)strlen(fn) ;
  if( l == 0 ) return 0 ;
  if( fn[l-4] != '.' )
  {
    strcat( fn, "." ) ;
    strcat( fn, ext ) ;
  }
  else if( strcmp( fn + l-3, ext ) )
  {
    fn[l-3] = ext[0] ;
    fn[l-2] = ext[1] ;
    fn[l-1] = ext[2] ;
  }

  return l ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// PPM export
void export_ppm()
//-----------------------------------------------------------------------------
{
  char fn[1024] ;
  set_ext( "ppm", fn ) ;
  FILE *fp = fopen( fn, "w" ) ;
  if( !fp ) return ;

  int tx, ty, tw, th;
  glutSetWindow(main_window);
  GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );

  // read openGL buffer
  int buf_size = (tw * th * 3);
  unsigned char *buffer = new unsigned char[buf_size] ;
  memset(buffer, 0, buf_size * sizeof(unsigned char)) ;
  display() ;
  glFinish() ;
  glPixelStorei( GL_PACK_ALIGNMENT,  1 ); // Force 1-byte alignment
  glPixelStorei( GL_PACK_ROW_LENGTH, 0 );
  glPixelStorei( GL_PACK_SKIP_ROWS,  0 );
  glPixelStorei( GL_PACK_SKIP_PIXELS,0 );
  glReadPixels(tx, ty, tw, th, GL_RGB, GL_UNSIGNED_BYTE, buffer);

  // write the header
  fprintf(fp, "P6\n%d %d\n", tw, th);
  fprintf(fp, "# %s (Affine Surface by Thomas Lewiner)\n", fn ) ;
  fprintf(fp, "255\n");

  // write file
  for (int i = 0; i < buf_size; ++i)
    fputc( buffer[i], fp ) ;

  delete [] buffer ;

  fclose(fp);
  printf( "exported image to PPM file %s\n", fn ) ;

}
//_____________________________________________________________________________


//_____________________________________________________________________________
//
void export_tga()
//-----------------------------------------------------------------------------
{
  char fn[1024] ;
  set_ext("tga", fn) ;
  FILE *fp = fopen( fn, "w" ) ;
  if (!fp)
     return ;

  int tx, ty, tw, th;
  glutSetWindow(main_window);
  GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);

  int             buf_size  = 18 + (tw* th * 3);  // HEADER_SIZE  ==> 18
  unsigned char*  buffer    = new unsigned char[buf_size] ;
  memset(buffer, 0, buf_size * sizeof(unsigned char)) ;

  // Header
  buffer[2] = 2;     // Non-compressed
  buffer[12] = tw & 255;
  buffer[13] = (tw >> 8) & 255 ;
  buffer[14] = th & 255;
  buffer[15] = (th >> 8) & 255 ;
  buffer[16] = 24;    // 24 bits per pixel

  // Read openGL buffer
  display() ;
  glFinish() ;
  glPixelStorei(GL_PACK_ALIGNMENT, 1); // Force 1-byte alignment
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
  glReadPixels(tx, ty, tw, th, GL_RGB, GL_UNSIGNED_BYTE, buffer + 18);

  // Conversion from RGB to BGR
  for (int i = 18; i < buf_size; i += 3)
  {
     unsigned char  temp  = buffer[i];
     buffer[i] = buffer[i + 2];
     buffer[i + 2] = temp;
  }

  // Write file
  fwrite(buffer, sizeof(unsigned char), buf_size, fp);
  fclose(fp);
  delete[] buffer ;
  printf("exported image to TGA file %s\n", fn ) ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
//_____________________________________________________________________________




//_____________________________________________________________________________
//
void save_viewport()
//-----------------------------------------------------------------------------
{
  char fn[1024] ;
  set_ext( "view", fn ) ;
  FILE *fp = fopen( fn, "w" ) ;
  if( !fp ) return ;
  fprintf( fp, "rotate:\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\n",
          view_rotate[ 0], view_rotate[ 1], view_rotate[ 2], view_rotate[ 3],
          view_rotate[ 4], view_rotate[ 5], view_rotate[ 6], view_rotate[ 7],
          view_rotate[ 8], view_rotate[ 9], view_rotate[10], view_rotate[11],
          view_rotate[12], view_rotate[13], view_rotate[14], view_rotate[15] ) ;
  fprintf( fp, "translate:\t%f %f %f\n", obj_pos[0], obj_pos[1], obj_pos[2] ) ;

  fprintf( fp, "plane rotate:\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f\n\n",
          1.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,1.0,0.0, 0.0,0.0,0.0,1.0, 0.0 ) ;
  
  fprintf( fp, "corner size:\t%f %f\n", 0.0,0.0 ) ;
  fprintf( fp, "corner trans:\t%f %f %f\n", 0.0,0.0,1.0 ) ;
  
  
  fprintf( fp, "light 0:\n\t%d %d %d\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\n\n",
          light0_enabled, light0_intensity, light0_intensity2, 
          light0_diffuse [ 0], light0_diffuse [ 1], light0_diffuse [ 2], light0_diffuse [ 3],
          light0_rotation[ 0], light0_rotation[ 1], light0_rotation[ 2], light0_rotation[ 3],
          light0_rotation[ 4], light0_rotation[ 5], light0_rotation[ 6], light0_rotation[ 7],
          light0_rotation[ 8], light0_rotation[ 9], light0_rotation[10], light0_rotation[11],
          light0_rotation[12], light0_rotation[13], light0_rotation[14], light0_rotation[15] ) ;
  
  fprintf( fp, "light 1:\n\t%d %d %d\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\n\n",
          light1_enabled, light1_intensity, light1_intensity2, 
          light1_diffuse [ 0], light1_diffuse [ 1], light1_diffuse [ 2], light1_diffuse [ 3],
          light1_rotation[ 0], light1_rotation[ 1], light1_rotation[ 2], light1_rotation[ 3],
          light1_rotation[ 4], light1_rotation[ 5], light1_rotation[ 6], light1_rotation[ 7],
          light1_rotation[ 8], light1_rotation[ 9], light1_rotation[10], light1_rotation[11],
          light1_rotation[12], light1_rotation[13], light1_rotation[14], light1_rotation[15] ) ;
  
  fclose( fp ) ;

  printf("saved viewport to file %s\n", fn ) ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
//
void load_viewport()
//-----------------------------------------------------------------------------
{
  char fn[1024] ;
  set_ext("view",fn) ;
  FILE* fp  = fopen( fn, "r" ) ;
  if( !fp ) return ;
  fscanf( fp, "rotate: %f %f %f %f  %f %f %f %f  %f %f %f %f  %f %f %f %f ",
         view_rotate +  0, view_rotate +  1, view_rotate +  2, view_rotate +  3,
         view_rotate +  4, view_rotate +  5, view_rotate +  6, view_rotate +  7,
         view_rotate +  8, view_rotate +  9, view_rotate + 10, view_rotate + 11,
         view_rotate + 12, view_rotate + 13, view_rotate + 14, view_rotate + 15 ) ;
  fscanf( fp, "translate: %f %f %f ", obj_pos + 0, obj_pos + 1, obj_pos + 2 ) ;

  fscanf ( fp, "plane rotate:\n\t%*f %*f %*f %*f\n\t%*f %*f %*f %*f\n\t%*f %*f %*f %*f\n\t%*f %*f %*f %*f\n\t%*f\n\n" ) ;
  
  fscanf ( fp, "corner size:\t%*f %*f\n" ) ;
  fscanf ( fp, "corner trans:\t%*f %*f %*f\n" ) ;
  
  
  fscanf ( fp, "light 0:\n\t%d %d %d\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\n\n",
          &light0_enabled, &light0_intensity, &light0_intensity2, 
          light0_diffuse  +  0, light0_diffuse  +  1, light0_diffuse  +  2, light0_diffuse  +  3,
          light0_rotation +  0, light0_rotation +  1, light0_rotation +  2, light0_rotation +  3,
          light0_rotation +  4, light0_rotation +  5, light0_rotation +  6, light0_rotation +  7,
          light0_rotation +  8, light0_rotation +  9, light0_rotation + 10, light0_rotation + 11,
          light0_rotation + 12, light0_rotation + 13, light0_rotation + 14, light0_rotation + 15 ) ;
  
  fscanf ( fp, "light 1:\n\t%d %d %d\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\t%f %f %f %f\n\n\n",
          &light1_enabled, &light1_intensity, &light1_intensity2, 
          light1_diffuse  +  0, light1_diffuse  +  1, light1_diffuse  +  2, light1_diffuse  +  3,
          light1_rotation +  0, light1_rotation +  1, light1_rotation +  2, light1_rotation +  3,
          light1_rotation +  4, light1_rotation +  5, light1_rotation +  6, light1_rotation +  7,
          light1_rotation +  8, light1_rotation +  9, light1_rotation + 10, light1_rotation + 11,
          light1_rotation + 12, light1_rotation + 13, light1_rotation + 14, light1_rotation + 15 ) ;

  fclose( fp ) ;
  
  printf("loaded viewport from file %s\n", fn ) ;
}
//_____________________________________________________________________________






//_____________________________________________________________________________
// PNG export
void export_png()
//-----------------------------------------------------------------------------
{
  char fn[1024] ;
  set_ext( "png", fn ) ;
  FILE *fp = fopen( fn, "wb" ) ;
  if( !fp ) return ;
  
  int tx, ty, tw, th;
  glutSetWindow(main_window);
  GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
  
  // read openGL buffer
  int buf_size = (tw * th * 3);
  unsigned char *buffer = new unsigned char[buf_size] ;
  memset(buffer, 0, buf_size * sizeof(unsigned char)) ;
  display() ;
  glFinish() ;
  glPixelStorei( GL_PACK_ALIGNMENT,  1 ); // Force 1-byte alignment
  glPixelStorei( GL_PACK_ROW_LENGTH, 0 );
  glPixelStorei( GL_PACK_SKIP_ROWS,  0 );
  glPixelStorei( GL_PACK_SKIP_PIXELS,0 );
  glReadPixels(tx, ty, tw, th, GL_RGB, GL_UNSIGNED_BYTE, buffer);
  

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL, NULL, NULL);
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (setjmp(png_jmpbuf(png_ptr)))
  {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    return;
  }
  png_init_io(png_ptr, fp);
  png_set_IHDR(png_ptr, info_ptr, tw, th, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);
  for (int i=0; i<th; i++)
    png_write_row(png_ptr, &(buffer[3*tw*((th-1) - i)]));
  png_write_end(png_ptr, NULL);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  
  delete [] buffer ;
  
  fclose(fp);
  printf( "exported image to PNG file %s\n", fn ) ;
}
//_____________________________________________________________________________


