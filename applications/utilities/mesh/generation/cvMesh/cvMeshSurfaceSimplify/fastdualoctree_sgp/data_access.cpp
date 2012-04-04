/**
 * @file    data_access.cpp
 * \author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
 * \author  Math Dept, PUC-Rio
 * \author  Matmidia Labs
 * \date    14/02/2006
 *
 * @brief   data accesss
 */
//________________________________________________



#ifndef WIN32
#pragma implementation
#endif // WIN32

#include "data_access.h"


//_____________________________________________________________________________
//_____________________________________________________________________________
// data_func

//_____________________________________________________________________________
// crossing data access: parser
void data_func::parse()
//-----------------------------------------------------------------------------
{
  if( _formula.empty() )
    return ;
  std::string vars("x,y,z") ;
  int i = _parser.Parse( _formula, vars ) ;
  if( _parser.GetParseErrorType() != FunctionParser::FP_NO_ERROR )
    printf( "%s\n%s\n% *d\n", _parser.ErrorMsg(), _formula.c_str(), i, 1 ) ;
  else
    _parser.Optimize() ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// crossing data access
bool data_func::need_refine( const Cube &c )
//-----------------------------------------------------------------------------
{
    if( _parser.GetParseErrorType() != FunctionParser::FP_NO_ERROR ) return false ;
  int l = c.lv() ;
  if( l < 2 ) return true ;
  if( l >= _max_level ) return false ;
  
  real cx = c.cx() ;
  real cy = c.cy() ;
  real cz = c.cz() ;
  _v[X] = cx ;
  _v[Y] = cy ;
  _v[Z] = cz ;
  
  real ref_val = _parser.Eval(_v) - _iso_val ;
  if( fabs(ref_val) < R_EPSILON ) return true ;
  
  real sz = c.sz() ;
  for( _v[X] = cx-sz ; _v[X] <= cx+sz ; _v[X] += sz )
  {
    for( _v[Y] = cy-sz ; _v[Y] <= cy+sz ; _v[Y] += sz )
    {
      for( _v[Z] = cz-sz ; _v[Z] <= cz+sz ; _v[Z] += sz )
      {
        real val = _parser.Eval(_v) - _iso_val ;
        if( ref_val * val < 0 ) return true ;
      }
    }
  }
  
  return false ;
}
//_____________________________________________________________________________


#ifndef BISSEC_ITER
#define BISSEC_ITER 3
#endif // BISSEC_ITER

//_____________________________________________________________________________
// interpolation function
void data_func::interpolate( real ci, real cj, const Point &pi, const Point &pj, Point &p )
//-----------------------------------------------------------------------------
{
    real u ;
  Point qi = pi ;
  Point qj = pj ;
  for( int i = 0 ; i < BISSEC_ITER ; ++i )
  {
    u = ( ci ) / ( ci - cj ) ;
    p = ((1-u)*qi) + (u*qj) ;
    
    real c = value_at(p) ;
    if( c*ci < 0 )
    {
      qj = p ;
      cj = c ;
    }
    else if( c*cj < 0 )
    {
      qi = p ;
      ci = c ;
    }
    else break ;
  }
  u = ( ci ) / ( ci - cj ) ;
  p = ((1-u)*qi) + (u*qj) ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
//_____________________________________________________________________________
// data_file


#include <errno.h>
#include <string.h>
#include <zlib.h>

//_____________________________________________________________________________
// data data access: load file
bool data_file::load_data( const char *fn  )
//-----------------------------------------------------------------------------
{
  errno = 0;
  gzFile isofile = gzopen( fn, "rb" ) ;
  if( !isofile )
  {
    printf( "data_file::load_data couldn't open file %s: %s\n", fn, strerror (errno));
    return false ;
  }
  
  // read size
  unsigned char buf[sizeof(float)] ;
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  _size_x = * (int*)buf ;
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  _size_y = * (int*)buf ;
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  _size_z = * (int*)buf ;
  
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  float xmin = * (float*)buf ;
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  float xmax = * (float*)buf ;
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  float ymin = * (float*)buf ;
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  float ymax = * (float*)buf ;
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  float zmin = * (float*)buf ;
  if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
  float zmax = * (float*)buf ;
  
  // allocate data
  delete [] _data ;
  _size_yz = _size_z * _size_y ;
  uint total_size = _size_x * _size_yz ;
  if( total_size == 0 ) return false ;
  
  // cube fitting
  _max_size = _size_x ;
  if( _max_size < _size_y ) _max_size = _size_y ;
  if( _max_size < _size_z ) _max_size = _size_z ;
  Level p = (Level)floor( log(_max_size)/log(2) ) ;
  if( _max_size != (1<<p) ) _max_size = 1 << (p+1) ;
  
  _data = new float[ total_size ] ;
  if( _data == (float*)NULL ) return false ;
  
  float min_field =  FLT_MAX ;
  float max_field = -FLT_MAX ;
  float *data_ptr = _data ;
  for( uint i = 0 ; i < total_size ; ++i, ++data_ptr )
  {
    if( gzread( isofile, buf, sizeof(float) ) < 1 ) return false ;
    float v = * (float*) buf ;
    *data_ptr = v ;
    if( v < min_field ) min_field = v ;
    if( v > max_field ) max_field = v ;
  }
  
  gzclose( isofile ) ;
  
  printf( "data %s: %dx%dx%d, field originally in [ %f , %f ]\n", fn, (int)_size_x, (int)_size_y, (int)_size_z, min_field, max_field ) ;
  
  // rescale data in [-1.1]
  float v_scale = 2.0 / (max_field - min_field) ;
  data_ptr = _data ;
  for( uint i = 0 ; i < total_size ; ++i, ++data_ptr )
  {
    *data_ptr = v_scale * ( (*data_ptr) - min_field ) - 1.0 ;
  }
  
  // cube scaling
  _sx = (xmax - xmin)/_size_x ;
  _sy = (ymax - ymin)/_size_y ;
  _sz = (zmax - zmin)/_size_z ;
  float s = _sx ;
  if( s < _sy ) s = _sy ;
  if( s < _sz ) s = _sz ;
  _sx /= s ;
  _sy /= s ;
  _sz /= s ;
  
  return true ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// data data access
bool data_file::need_refine( const Cube &c )
//-----------------------------------------------------------------------------
{
    int l = c.lv() ;
  if( l < 2 ) return true ;
  if( l >= _max_level ) return false ;
  
  // ajustar se o dado nao for cubico!
  real sz = c.sz() ;
  uint i  = floor( (c.cx()-sz) * _max_size ) ;
  uint iM = ceil ( (c.cx()+sz) * _max_size ) ;
  uint j  = floor( (c.cy()-sz) * _max_size ) ;
  uint jM = ceil ( (c.cy()+sz) * _max_size ) ;
  uint k  = floor( (c.cz()-sz) * _max_size ) ;
  uint kM = ceil ( (c.cz()+sz) * _max_size ) ;
  
  if( iM >= _size_x ) iM = _size_x-1 ;
  if( jM >= _size_y ) jM = _size_y-1 ;
  if( kM >= _size_z ) kM = _size_z-1 ;
  
  float iso_min = get_data(iM,jM,kM) ;
  float iso_max = iso_min ;
  for( ; i <= iM ; ++i )
  {
    for( ; j <= jM ; ++j )
    {
      for( ; k <= kM ; ++k )
      {
        float iso = get_data(i,j,k) ;
        if( iso_min > iso ) iso_min = iso ;
        if( iso_max < iso ) iso_max = iso ;
        
        if( iso_min*iso_max < 0 ) return true ;
      }
    }
  }
  
  return false ;
}
//_____________________________________________________________________________



