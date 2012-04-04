/**
 * @file    md_draw.cpp
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 *
 * @brief   MC_Draw Direct Draw Algorithm
 */
//________________________________________________


#if !defined(WIN32) || defined(__CYGWIN__)
#pragma implementation
#endif // WIN32

#include <math.h>
#include <time.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "mc_draw.h"
#include "LookUpTable.h"


//_____________________________________________________________________________
// print cube for debug
void MC_Draw::print_cube() { printf( "\t%f %f %f %f %f %f %f %f\n", _cube[0], _cube[1], _cube[2], _cube[3], _cube[4], _cube[5], _cube[6], _cube[7]) ; }
//_____________________________________________________________________________



//_____________________________________________________________________________
// main algorithm
bool MC_Draw::tesselate_cube( real iso )
//-----------------------------------------------------------------------------
{
  for( int p = 0 ; p < 8 ; ++p )
  {
    real &v = _cube[p] ;
    v -= iso ;
    if( fabs( v ) < R_EPSILON ) v = R_EPSILON ;
  }
  
  if( !compute_intersection_points() ) return false ;

  _lut_entry = 0 ;
  for( int p = 0 ; p < 8 ; ++p )
  {
    if( _cube[p] > 0 ) _lut_entry += 1 << p ;
  }

  return process_cube( ) ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
//_____________________________________________________________________________


//_____________________________________________________________________________
// Compute the intersection points
bool MC_Draw::compute_intersection_points()
//-----------------------------------------------------------------------------
{
  bool res = false ;
  res |= add_vertex( 0,1 ) ;
  res |= add_vertex( 1,2 ) ;
  res |= add_vertex( 2,3 ) ;
  res |= add_vertex( 3,0 ) ;

  res |= add_vertex( 4,5 ) ;
  res |= add_vertex( 5,6 ) ;
  res |= add_vertex( 6,7 ) ;
  res |= add_vertex( 7,4 ) ;

  res |= add_vertex( 0,4 ) ;
  res |= add_vertex( 1,5 ) ;
  res |= add_vertex( 2,6 ) ;
  res |= add_vertex( 3,7 ) ;
  
  return res ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Test a face
// if face>0 return true if the face contains a part of the surface
bool MC_Draw::test_face( schar face )
//-----------------------------------------------------------------------------
{
  real A,B,C,D ;

  switch( face )
  {
  case -1 : case 1 :  A = _cube[0] ;  B = _cube[4] ;  C = _cube[5] ;  D = _cube[1] ;  break ;
  case -2 : case 2 :  A = _cube[1] ;  B = _cube[5] ;  C = _cube[6] ;  D = _cube[2] ;  break ;
  case -3 : case 3 :  A = _cube[2] ;  B = _cube[6] ;  C = _cube[7] ;  D = _cube[3] ;  break ;
  case -4 : case 4 :  A = _cube[3] ;  B = _cube[7] ;  C = _cube[4] ;  D = _cube[0] ;  break ;
  case -5 : case 5 :  A = _cube[0] ;  B = _cube[3] ;  C = _cube[2] ;  D = _cube[1] ;  break ;
  case -6 : case 6 :  A = _cube[4] ;  B = _cube[7] ;  C = _cube[6] ;  D = _cube[5] ;  break ;
  default : printf( "Invalid face code %d\n", face ) ;  print_cube() ;  A = B = C = D = 0 ;
  };

  if( fabs( A*C - B*D ) < R_EPSILON )
    return face >= 0 ;
  return face * A * ( A*C - B*D ) >= 0  ;  // face and A invert signs
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Test the interior of a cube
// if s == 7, return true  if the interior is empty
// if s ==-7, return false if the interior is empty
bool MC_Draw::test_interior( schar s )
//-----------------------------------------------------------------------------
{
  real t, At=0, Bt=0, Ct=0, Dt=0, a, b ;
  char  test =  0 ;
  char  edge = -1 ; // reference edge of the triangulation

  switch( _case )
  {
  case  4 :
  case 10 :
    a = ( _cube[4] - _cube[0] ) * ( _cube[6] - _cube[2] ) - ( _cube[7] - _cube[3] ) * ( _cube[5] - _cube[1] ) ;
    b =  _cube[2] * ( _cube[4] - _cube[0] ) + _cube[0] * ( _cube[6] - _cube[2] )
             - _cube[1] * ( _cube[7] - _cube[3] ) - _cube[3] * ( _cube[5] - _cube[1] ) ;
    t = - b / (2*a) ;
    if( t<0 || t>1 ) return s>0 ;

    At = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
    Bt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
    Ct = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
    Dt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
    break ;

  case  6 :
  case  7 :
  case 12 :
  case 13 :
    switch( _case )
    {
    case  6 : edge = test6 [_config][2] ; break ;
    case  7 : edge = test7 [_config][4] ; break ;
    case 12 : edge = test12[_config][3] ; break ;
    case 13 : edge = tiling13_5_1[_config][_subconfig][0] ; break ;
    }
    switch( edge )
    {
    case  0 :
      t  = _cube[0] / ( _cube[0] - _cube[1] ) ;
      At = 0 ;
      Bt = _cube[3] + ( _cube[2] - _cube[3] ) * t ;
      Ct = _cube[7] + ( _cube[6] - _cube[7] ) * t ;
      Dt = _cube[4] + ( _cube[5] - _cube[4] ) * t ;
      break ;
    case  1 :
      t  = _cube[1] / ( _cube[1] - _cube[2] ) ;
      At = 0 ;
      Bt = _cube[0] + ( _cube[3] - _cube[0] ) * t ;
      Ct = _cube[4] + ( _cube[7] - _cube[4] ) * t ;
      Dt = _cube[5] + ( _cube[6] - _cube[5] ) * t ;
      break ;
    case  2 :
      t  = _cube[2] / ( _cube[2] - _cube[3] ) ;
      At = 0 ;
      Bt = _cube[1] + ( _cube[0] - _cube[1] ) * t ;
      Ct = _cube[5] + ( _cube[4] - _cube[5] ) * t ;
      Dt = _cube[6] + ( _cube[7] - _cube[6] ) * t ;
      break ;
    case  3 :
      t  = _cube[3] / ( _cube[3] - _cube[0] ) ;
      At = 0 ;
      Bt = _cube[2] + ( _cube[1] - _cube[2] ) * t ;
      Ct = _cube[6] + ( _cube[5] - _cube[6] ) * t ;
      Dt = _cube[7] + ( _cube[4] - _cube[7] ) * t ;
      break ;
    case  4 :
      t  = _cube[4] / ( _cube[4] - _cube[5] ) ;
      At = 0 ;
      Bt = _cube[7] + ( _cube[6] - _cube[7] ) * t ;
      Ct = _cube[3] + ( _cube[2] - _cube[3] ) * t ;
      Dt = _cube[0] + ( _cube[1] - _cube[0] ) * t ;
      break ;
    case  5 :
      t  = _cube[5] / ( _cube[5] - _cube[6] ) ;
      At = 0 ;
      Bt = _cube[4] + ( _cube[7] - _cube[4] ) * t ;
      Ct = _cube[0] + ( _cube[3] - _cube[0] ) * t ;
      Dt = _cube[1] + ( _cube[2] - _cube[1] ) * t ;
      break ;
    case  6 :
      t  = _cube[6] / ( _cube[6] - _cube[7] ) ;
      At = 0 ;
      Bt = _cube[5] + ( _cube[4] - _cube[5] ) * t ;
      Ct = _cube[1] + ( _cube[0] - _cube[1] ) * t ;
      Dt = _cube[2] + ( _cube[3] - _cube[2] ) * t ;
      break ;
    case  7 :
      t  = _cube[7] / ( _cube[7] - _cube[4] ) ;
      At = 0 ;
      Bt = _cube[6] + ( _cube[5] - _cube[6] ) * t ;
      Ct = _cube[2] + ( _cube[1] - _cube[2] ) * t ;
      Dt = _cube[3] + ( _cube[0] - _cube[3] ) * t ;
      break ;
    case  8 :
      t  = _cube[0] / ( _cube[0] - _cube[4] ) ;
      At = 0 ;
      Bt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
      Ct = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
      Dt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
      break ;
    case  9 :
      t  = _cube[1] / ( _cube[1] - _cube[5] ) ;
      At = 0 ;
      Bt = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
      Ct = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
      Dt = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
      break ;
    case 10 :
      t  = _cube[2] / ( _cube[2] - _cube[6] ) ;
      At = 0 ;
      Bt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
      Ct = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
      Dt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
      break ;
    case 11 :
      t  = _cube[3] / ( _cube[3] - _cube[7] ) ;
      At = 0 ;
      Bt = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
      Ct = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
      Dt = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
      break ;
    default : printf( "Invalid edge %d\n", edge ) ;  print_cube() ;  break ;
    }
    break ;

  default : printf( "Invalid ambiguous case %d\n", _case ) ;  print_cube() ;  break ;
  }

  if( At >= 0 ) test ++ ;
  if( Bt >= 0 ) test += 2 ;
  if( Ct >= 0 ) test += 4 ;
  if( Dt >= 0 ) test += 8 ;
  switch( test )
  {
  case  0 : return s>0 ;
  case  1 : return s>0 ;
  case  2 : return s>0 ;
  case  3 : return s>0 ;
  case  4 : return s>0 ;
  case  5 : if( At * Ct - Bt * Dt <  R_EPSILON ) return s>0 ; break ;
  case  6 : return s>0 ;
  case  7 : return s<0 ;
  case  8 : return s>0 ;
  case  9 : return s>0 ;
  case 10 : if( At * Ct - Bt * Dt >= R_EPSILON ) return s>0 ; break ;
  case 11 : return s<0 ;
  case 12 : return s>0 ;
  case 13 : return s<0 ;
  case 14 : return s<0 ;
  case 15 : return s<0 ;
  }

  return s<0 ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Process a unit cube
bool MC_Draw::process_cube( )
//-----------------------------------------------------------------------------
{
  if( _originalMC )
  {
    char nt = 0 ;
    while( casesClassic[_lut_entry][3*nt] != -1 ) nt++ ;
    draw_triangle( casesClassic[_lut_entry], nt ) ;
    return true ;
  }

  int   v12 = -1 ;
  _case   = cases[_lut_entry][0] ;
  _config = cases[_lut_entry][1] ;
  _subconfig = 0 ;

  switch( _case )
  {
  case  0 :
    break ;

  case  1 :
    draw_triangle( tiling1[_config], 1 ) ;
    break ;

  case  2 :
    draw_triangle( tiling2[_config], 2 ) ;
    break ;

  case  3 :
    if( test_face( test3[_config]) )
      draw_triangle( tiling3_2[_config], 4 ) ; // 3.2
    else
      draw_triangle( tiling3_1[_config], 2 ) ; // 3.1
    break ;

  case  4 :
    if( test_interior( test4[_config]) )
      draw_triangle( tiling4_1[_config], 2 ) ; // 4.1.1
    else
      draw_triangle( tiling4_2[_config], 6 ) ; // 4.1.2
    break ;

  case  5 :
    draw_triangle( tiling5[_config], 3 ) ;
    break ;

  case  6 :
    if( test_face( test6[_config][0]) )
      draw_triangle( tiling6_2[_config], 5 ) ; // 6.2
    else
    {
      if( test_interior( test6[_config][1]) )
        draw_triangle( tiling6_1_1[_config], 3 ) ; // 6.1.1
      else
	  {
        comp_c_vertex() ;
        draw_triangle( tiling6_1_2[_config], 9 , v12) ; // 6.1.2
      }
    }
    break ;

  case  7 :
    if( test_face( test7[_config][0] ) ) _subconfig +=  1 ;
    if( test_face( test7[_config][1] ) ) _subconfig +=  2 ;
    if( test_face( test7[_config][2] ) ) _subconfig +=  4 ;
    switch( _subconfig )
      {
      case 0 :
        draw_triangle( tiling7_1[_config], 3 ) ; break ;
      case 1 :
        draw_triangle( tiling7_2[_config][0], 5 ) ; break ;
      case 2 :
        draw_triangle( tiling7_2[_config][1], 5 ) ; break ;
      case 3 :
        comp_c_vertex() ;
        draw_triangle( tiling7_3[_config][0], 9, v12 ) ; break ;
      case 4 :
        draw_triangle( tiling7_2[_config][2], 5 ) ; break ;
      case 5 :
        comp_c_vertex() ;
        draw_triangle( tiling7_3[_config][1], 9, v12 ) ; break ;
      case 6 :
        comp_c_vertex() ;
        draw_triangle( tiling7_3[_config][2], 9, v12 ) ; break ;
      case 7 :
        if( test_interior( test7[_config][3]) )
          draw_triangle( tiling7_4_2[_config], 9 ) ;
        else
          draw_triangle( tiling7_4_1[_config], 5 ) ;
        break ;
      };
    break ;

  case  8 :
    draw_triangle( tiling8[_config], 2 ) ;
    break ;

  case  9 :
    draw_triangle( tiling9[_config], 4 ) ;
    break ;

  case 10 :
    if( test_face( test10[_config][0]) )
    {
      if( test_face( test10[_config][1]) )
        draw_triangle( tiling10_1_1_[_config], 4 ) ; // 10.1.1
      else
      {
        comp_c_vertex() ;
        draw_triangle( tiling10_2[_config], 8, v12 ) ; // 10.2
      }
    }
    else
    {
      if( test_face( test10[_config][1]) )
      {
        comp_c_vertex() ;
        draw_triangle( tiling10_2_[_config], 8, v12 ) ; // 10.2
      }
      else
      {
        if( test_interior( test10[_config][2]) )
          draw_triangle( tiling10_1_1[_config], 4 ) ; // 10.1.1
        else
          draw_triangle( tiling10_1_2[_config], 8 ) ; // 10.1.2
      }
    }
    break ;

  case 11 :
    draw_triangle( tiling11[_config], 4 ) ;
    break ;

  case 12 :
    if( test_face( test12[_config][0]) )
    {
      if( test_face( test12[_config][1]) )
        draw_triangle( tiling12_1_1_[_config], 4 ) ; // 12.1.1
      else
      {
        comp_c_vertex() ;
        draw_triangle( tiling12_2[_config], 8, v12 ) ; // 12.2
      }
    }
    else
    {
      if( test_face( test12[_config][1]) )
      {
        comp_c_vertex() ;
        draw_triangle( tiling12_2_[_config], 8, v12 ) ; // 12.2
      }
      else
      {
        if( test_interior( test12[_config][2]) )
          draw_triangle( tiling12_1_1[_config], 4 ) ; // 12.1.1
        else
          draw_triangle( tiling12_1_2[_config], 8 ) ; // 12.1.2
      }
    }
    break ;

  case 13 :
    if( test_face( test13[_config][0] ) ) _subconfig +=  1 ;
    if( test_face( test13[_config][1] ) ) _subconfig +=  2 ;
    if( test_face( test13[_config][2] ) ) _subconfig +=  4 ;
    if( test_face( test13[_config][3] ) ) _subconfig +=  8 ;
    if( test_face( test13[_config][4] ) ) _subconfig += 16 ;
    if( test_face( test13[_config][5] ) ) _subconfig += 32 ;
    switch( subconfig13[_subconfig] )
    {
      case 0 :/* 13.1 */
        draw_triangle( tiling13_1[_config], 4 ) ; break ;

      case 1 :/* 13.2 */
        draw_triangle( tiling13_2[_config][0], 6 ) ; break ;
      case 2 :/* 13.2 */
        draw_triangle( tiling13_2[_config][1], 6 ) ; break ;
      case 3 :/* 13.2 */
        draw_triangle( tiling13_2[_config][2], 6 ) ; break ;
      case 4 :/* 13.2 */
        draw_triangle( tiling13_2[_config][3], 6 ) ; break ;
      case 5 :/* 13.2 */
        draw_triangle( tiling13_2[_config][4], 6 ) ; break ;
      case 6 :/* 13.2 */
        draw_triangle( tiling13_2[_config][5], 6 ) ; break ;

      case 7 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][0], 10, v12 ) ; break ;
      case 8 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][1], 10, v12 ) ; break ;
      case 9 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][2], 10, v12 ) ; break ;
      case 10 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][3], 10, v12 ) ; break ;
      case 11 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][4], 10, v12 ) ; break ;
      case 12 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][5], 10, v12 ) ; break ;
      case 13 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][6], 10, v12 ) ; break ;
      case 14 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][7], 10, v12 ) ; break ;
      case 15 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][8], 10, v12 ) ; break ;
      case 16 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][9], 10, v12 ) ; break ;
      case 17 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][10], 10, v12 ) ; break ;
      case 18 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3[_config][11], 10, v12 ) ; break ;

      case 19 :/* 13.4 */
        comp_c_vertex() ;
        draw_triangle( tiling13_4[_config][0], 12, v12 ) ; break ;
      case 20 :/* 13.4 */
        comp_c_vertex() ;
        draw_triangle( tiling13_4[_config][1], 12, v12 ) ; break ;
      case 21 :/* 13.4 */
        comp_c_vertex() ;
        draw_triangle( tiling13_4[_config][2], 12, v12 ) ; break ;
      case 22 :/* 13.4 */
        comp_c_vertex() ;
        draw_triangle( tiling13_4[_config][3], 12, v12 ) ; break ;

      case 23 :/* 13.5 */
        _subconfig = 0 ;
        if( test_interior( test13[_config][6] ) )
          draw_triangle( tiling13_5_1[_config][0], 6 ) ;
        else
          draw_triangle( tiling13_5_2[_config][0], 10 ) ;
        break ;
      case 24 :/* 13.5 */
        _subconfig = 1 ;
        if( test_interior( test13[_config][6] ) )
          draw_triangle( tiling13_5_1[_config][1], 6 ) ;
        else
          draw_triangle( tiling13_5_2[_config][1], 10 ) ;
        break ;
      case 25 :/* 13.5 */
        _subconfig = 2 ;
        if( test_interior( test13[_config][6] ) )
          draw_triangle( tiling13_5_1[_config][2], 6 ) ;
        else
          draw_triangle( tiling13_5_2[_config][2], 10 ) ;
        break ;
      case 26 :/* 13.5 */
        _subconfig = 3 ;
        if( test_interior( test13[_config][6] ) )
          draw_triangle( tiling13_5_1[_config][3], 6 ) ;
        else
          draw_triangle( tiling13_5_2[_config][3], 10 ) ;
        break ;

      case 27 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][0], 10, v12 ) ; break ;
      case 28 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][1], 10, v12 ) ; break ;
      case 29 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][2], 10, v12 ) ; break ;
      case 30 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][3], 10, v12 ) ; break ;
      case 31 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][4], 10, v12 ) ; break ;
      case 32 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][5], 10, v12 ) ; break ;
      case 33 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][6], 10, v12 ) ; break ;
      case 34 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][7], 10, v12 ) ; break ;
      case 35 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][8], 10, v12 ) ; break ;
      case 36 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][9], 10, v12 ) ; break ;
      case 37 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][10], 10, v12 ) ; break ;
      case 38 :/* 13.3 */
        comp_c_vertex() ;
        draw_triangle( tiling13_3_[_config][11], 10, v12 ) ; break ;

      case 39 :/* 13.2 */
        draw_triangle( tiling13_2_[_config][0], 6 ) ; break ;
      case 40 :/* 13.2 */
        draw_triangle( tiling13_2_[_config][1], 6 ) ; break ;
      case 41 :/* 13.2 */
        draw_triangle( tiling13_2_[_config][2], 6 ) ; break ;
      case 42 :/* 13.2 */
        draw_triangle( tiling13_2_[_config][3], 6 ) ; break ;
      case 43 :/* 13.2 */
        draw_triangle( tiling13_2_[_config][4], 6 ) ; break ;
      case 44 :/* 13.2 */
        draw_triangle( tiling13_2_[_config][5], 6 ) ; break ;

      case 45 :/* 13.1 */
        draw_triangle( tiling13_1_[_config], 4 ) ; break ;

      default :
        printf("Marching Cubes: Impossible case 13?\n" ) ;  print_cube() ;
      }
      break ;

  case 14 :
    draw_triangle( tiling14[_config], 4 ) ;
    break ;
  };
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Adding triangles
void MC_Draw::draw_triangle( const char* trig_, char n, int v12 )
//-----------------------------------------------------------------------------
{
  for( int t = 0 ; t < 3*n ; t++ )
  {
    switch( trig_[t] )
    {
      case  0 : _verts[0][1].draw() ; break ;
      case  1 : _verts[1][2].draw() ; break ;
      case  2 : _verts[2][3].draw() ; break ;
      case  3 : _verts[3][0].draw() ; break ;
      case  4 : _verts[4][5].draw() ; break ;
      case  5 : _verts[5][6].draw() ; break ;
      case  6 : _verts[6][7].draw() ; break ;
      case  7 : _verts[7][4].draw() ; break ;
      case  8 : _verts[0][4].draw() ; break ;
      case  9 : _verts[1][5].draw() ; break ;
      case 10 : _verts[2][6].draw() ; break ;
      case 11 : _verts[3][7].draw() ; break ;
      case 12 : _c_vert     .draw() ; break ;
      default : break ;
    }
  }
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Adding vertices
bool MC_Draw::add_vertex( char i, char j )
//-----------------------------------------------------------------------------
{
  real ci = _cube[i] ;
  real cj = _cube[j] ;
  if( ci * cj >= 0 ) { _verts[i][j] = P_INV ; return false ; }
  
  _dat_access->interpolate( ci,cj, _space[i],_space[j], _verts[i][j] ) ;

  return true ;
}
//-----------------------------------------------------------------------------

bool MC_Draw::comp_c_vertex( )
//-----------------------------------------------------------------------------
{
  real u = 0 ;

  _c_vert.x() = _c_vert.y() = _c_vert.z() = 0 ;

  // Computes the average of the intersection points of the cube
  Point &v01 = _verts[0][1] ;  if( v01 != P_INV ) { ++u ; _c_vert += v01 ; }
  Point &v12 = _verts[1][2] ;  if( v12 != P_INV ) { ++u ; _c_vert += v12 ; }
  Point &v23 = _verts[2][3] ;  if( v23 != P_INV ) { ++u ; _c_vert += v23 ; }
  Point &v30 = _verts[3][0] ;  if( v30 != P_INV ) { ++u ; _c_vert += v30 ; }
  Point &v45 = _verts[4][5] ;  if( v45 != P_INV ) { ++u ; _c_vert += v45 ; }
  Point &v56 = _verts[5][6] ;  if( v56 != P_INV ) { ++u ; _c_vert += v56 ; }
  Point &v67 = _verts[6][7] ;  if( v67 != P_INV ) { ++u ; _c_vert += v67 ; }
  Point &v74 = _verts[7][4] ;  if( v74 != P_INV ) { ++u ; _c_vert += v74 ; }
  Point &v04 = _verts[0][4] ;  if( v04 != P_INV ) { ++u ; _c_vert += v04 ; }
  Point &v15 = _verts[1][5] ;  if( v15 != P_INV ) { ++u ; _c_vert += v15 ; }
  Point &v26 = _verts[2][6] ;  if( v26 != P_INV ) { ++u ; _c_vert += v26 ; }
  Point &v37 = _verts[3][7] ;  if( v37 != P_INV ) { ++u ; _c_vert += v37 ; }
  
  _c_vert  /= u ;

  return u > 0 ;
}
//_____________________________________________________________________________



