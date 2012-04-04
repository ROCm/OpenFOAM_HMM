/**
 * @file    MarchingCubes.h
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 *
 * @brief   MarchingCubes Algorithm
 */
//________________________________________________


#pragma once

#if !defined(WIN32) || defined(__CYGWIN__)
#pragma interface
#endif // WIN32

#include <map>
#include "morton.h"
#include "point.h"
#include "data_access.h"

//-----------------------------------------------------------------------------
// Triangle structure
/** \struct Triangle "MarchingCubes.h" MarchingCubes
 * Indices of the oriented triange vertices
 * \brief triangle structure
 * \param v1 First vertex index
 * \param v2 Second vertex index
 * \param v3 Third vertex index
 */
typedef struct
{
  int v1,v2,v3 ;  /**< Triangle vertices */
} Triangle ;
//_____________________________________________________________________________



//_____________________________________________________________________________
/** Marching Cubes algorithm wrapper */
/** \class MarchingCubes
  * \brief Marching Cubes algorithm.
  */
class MarchingCubes
//-----------------------------------------------------------------------------
{
// Constructors
public :
  /**
   * Main and default constructor
   * \brief constructor
   */
  MarchingCubes () ;
  /** Destructor */
  ~MarchingCubes() ;

//-----------------------------------------------------------------------------
// Accessors
public :
  /** accesses the number of vertices of the generated mesh */
  inline const int nverts() const { return _nverts ; }
  /** accesses the number of triangles of the generated mesh */
  inline const int ntrigs() const { return _ntrigs ; }
  /** accesses a specific vertex of the generated mesh */
  inline Point    * vert( const int i ) const { if( i < 0  || i >= _nverts ) return ( Point  *)NULL ; return _vertices  + i ; }
  /** accesses a specific triangle of the generated mesh */
  inline Triangle * trig( const int i ) const { if( i < 0  || i >= _ntrigs ) return (Triangle*)NULL ; return _triangles + i ; }

  /** accesses the vertex buffer of the generated mesh */
  inline Point    *vertices () { return _vertices  ; }
  /** accesses the triangle buffer of the generated mesh */
  inline Triangle *triangles() { return _triangles ; }

  /**
   * selects wether the algorithm will use the enhanced topologically controlled lookup table or the original MarchingCubes
   * \param originalMC true for the original Marching Cubes
   */
  inline void set_method    ( const bool originalMC = false ) { _originalMC = originalMC ; }

  // Data initialization
  /** inits all structures (must set sizes before call) : the temporary structures and the mesh buffers */
  void init_all   () ;
  /** clears temporary structures : the grid and the main */
  void clean_temps() ;
  /** clears all structures : the temporary structures and the mesh buffers */
  void clean_all  () ;


//-----------------------------------------------------------------------------
// Exportation
public :
  /**
   * OFF exportation of the generated mesh
   * \param fn  name of the IV file to create
   */
  void writeOFF( const char *fn ) ;


//-----------------------------------------------------------------------------
// Algorithm
public :
  /** retrieves the isovalues at the cube vertices */
  real             *cube   () { return _cube ; }
  /** retrieves the geometry of the cube */
  Point            *space  () { return _space; }
  /** retrieves the indexes of the cube */
  Key              *indexes() { return _indexes; }
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

  /**
   * routine to add a triangle to the mesh
   * \param trig the code for the triangle as a sequence of edges index
   * \param n    the number of triangles to produce
   * \param v12  the index of the interior vertex to use, if necessary
   */
  void add_triangle ( const char* trig, char n, int v12 = -1 ) ;

  /** tests and eventually doubles the vertex buffer capacity for a new vertex insertion */
  void test_vertex_addition() ;
  /** adds a vertex on edge cube[i] cube[j] */
  bool add_vertex( char i, char j ) ;
  /** gets the vertex of edge cube[i] cube[j] (no check) */
  int  get_vertex( char i, char j ) { std::map< std::pair< Key,Key >, int >::const_iterator it = _stored_vertices.find( std::make_pair( _indexes[i],_indexes[j] ) ) ; if( it == _stored_vertices.end() ) return -1 ; else return it->second ; }
  /** adds a vertex inside the current cube */
  int add_c_vertex() ;

  /** prints cube for debug */
  void print_cube() ;

public:  
  /** draw the surface using openGL */
  void draw_surf() { Triangle *ptr = _triangles ; for( int i = 0 ; i < _ntrigs ; ++i, ++ptr ) { _vertices[ptr->v1].draw() ; _vertices[ptr->v2].draw() ; _vertices[ptr->v3].draw() ; } }
  
//-----------------------------------------------------------------------------
// Elements
protected :
  bool      _originalMC ;   /**< selects wether the algorithm will use the enhanced topologically controlled lookup table or the original MarchingCubes */

  int       _nverts     ;  /**< number of allocated vertices  in the vertex   buffer */
  int       _ntrigs     ;  /**< number of allocated triangles in the triangle buffer */
  int       _Nverts     ;  /**< size of the vertex   buffer */
  int       _Ntrigs     ;  /**< size of the triangle buffer */
  Point    *_vertices   ;  /**< vertex   buffer */
  Triangle *_triangles  ;  /**< triangle buffer */

  std::map< std::pair< Key,Key >, int > _stored_vertices ;
  
  real      _cube[8]    ;  /**< values of the implicit function on the active cube */
  Point     _space[8]   ;  /**< coordinates of the active cube */
  Key       _indexes[8] ;  /**< indiexs of the active cube */
  uchar     _lut_entry  ;  /**< cube sign representation in [0..255] */
  uchar     _case       ;  /**< case of the active cube in [0..15] */
  uchar     _config     ;  /**< configuration of the active cube */
  uchar     _subconfig  ;  /**< subconfiguration of the active cube */
  data_access *_dat_access ;  /**< data accessor */
};
//_____________________________________________________________________________


