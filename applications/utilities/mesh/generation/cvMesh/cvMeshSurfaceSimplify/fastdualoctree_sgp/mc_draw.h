/**
 * @file    mc_draw.h
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 *
 * @brief   MarchingCubes Direct Draw
 */
//________________________________________________


#pragma once

#if !defined(WIN32) || defined(__CYGWIN__)
#pragma interface
#endif // WIN32

#include <map>
#include "point.h"
#include "data_access.h"


//_____________________________________________________________________________
/** Marching Cubes algorithm wrapper */
/** \class MarchingCubes Direct Draw
  * \brief Marching Cubes Direct Draw.
  */
class MC_Draw
//-----------------------------------------------------------------------------
{
// Constructors
public :
  /**
   * Main and default constructor
   * \brief constructor
   */
  MC_Draw () : _originalMC(true) {}
  /** Destructor */
  ~MC_Draw () {}

//-----------------------------------------------------------------------------
// Accessors
public :
  /**
   * selects wether the algorithm will use the enhanced topologically controlled lookup table or the original MarchingCubes
   * \param originalMC true for the original Marching Cubes
   */
  inline void set_method    ( const bool originalMC = false ) { _originalMC = originalMC ; }

//-----------------------------------------------------------------------------
// Algorithm
public :
  /** retrieves the isovalues at the cube vertices */
  real             *cube   () { return _cube ; }
  /** retrieves the geometry of the cube */
  Point            *space  () { return _space; }
  /** retrieves the data accessor */
  data_access     *&dat_access() { return _dat_access ; }

  /**
   * Main algorithm
   * \param iso isovalue
   */
  bool tesselate_cube( real iso ) ;

protected :
  /** tesselates one cube */
  bool process_cube ()             ;
  /** tests if the components of the tesselation of the cube should be connected by the interior of an ambiguous face */
  bool test_face    ( schar face ) ;
  /** tests if the components of the tesselation of the cube should be connected through the interior of the cube */
  bool test_interior( schar s )    ;


//-----------------------------------------------------------------------------
// Operations
protected :
  /** computes almost all the vertices of the mesh by interpolation along the cubes edges */
  bool compute_intersection_points() ;

  /** adds a vertex on edge cube[i] cube[j] */
  bool add_vertex( char i, char j ) ;

  /** compute a vertex inside the current cube */
  bool comp_c_vertex() ;

  /**
   * routine to add a triangle to the mesh
   * \param trig the code for the triangle as a sequence of edges index
   * \param n    the number of triangles to produce
   * \param v12  the index of the interior vertex to use, if necessary
   */
  void draw_triangle ( const char* trig, char n, int v12 = -1 ) ;

  /** prints cube for debug */
  void print_cube() ;

//-----------------------------------------------------------------------------
// Elements
protected :
  bool      _originalMC ;   /**< selects wether the algorithm will use the enhanced topologically controlled lookup table or the original MarchingCubes */

  real      _cube[8]     ;  /**< values of the implicit function on the active cube */
  Point     _space[8]    ;  /**< coordinates of the active cube */
  Point     _verts[8][8] ;  /**< coordinates of the vertices, per edge */
  Point     _c_vert      ;  /**< coordinates of the central vertex */
  uchar     _lut_entry   ;  /**< cube sign representation in [0..255] */
  uchar     _case        ;  /**< case of the active cube in [0..15] */
  uchar     _config      ;  /**< configuration of the active cube */
  uchar     _subconfig   ;  /**< subconfiguration of the active cube */
  data_access *_dat_access ;  /**< data accessor */
};
//_____________________________________________________________________________


