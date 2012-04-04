/**
 * \file    hash_octree.cpp
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 *  Octree structure with hashtable and hierarchical operations
 */
//_____________________________________________________________________________


#ifndef WIN32
#pragma implementation
#endif // WIN32


#include <stack>
#include "hash_octree.h"

/// invalid key
template<> Hash<real>::KeyData Hash<real>::KD_INV = { KEY_INV, R_INV } ;



//_____________________________________________________________________________
// Create the root
void HashOctree::init()
//-----------------------------------------------------------------------------
{
  _max_level = 0 ;
  _max_field = R_INV ;
  HashField::KeyData kd ;
  kd.key  =  1  ;
  kd.data = 0.0 ;
  _hash.insert( kd ) ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Check that each point scan or scribble is in the right node
bool HashOctree::check ()
//-----------------------------------------------------------------------------
{
  bool result = true ;

  for( cell_iterator it = cells_begin() ; it() ; ++it )
  {
    Key k = it.key() << 3 ;
    bool has_son = false ;
    bool all_sons = true  ;
    for( int i = 0 ; i < 8 ; ++i )
    {
      if( _hash[ k | i ].key != KEY_INV )
        has_son = true ;
      else
        all_sons = false ;
    }
    
    if( it.is_leaf() && has_son )
    {
      printf("(hash) Invalid leaf %d: value %f\n", (int)k, _hash[k].data ) ;
      result = false ;
    }

    if( has_son && !all_sons )
    {
      printf("(hash) Invalid subdivision %d: value %f\n", (int)k, _hash[k].data ) ;
      result = false ;
    }
  }

  return result ;
}
//_____________________________________________________________________________






//_____________________________________________________________________________
// Deletes the content of the octree
void HashOctree::clear_octree()
//-----------------------------------------------------------------------------
{
  _hash.reset() ;
  _max_level = 0 ;
  _max_field = R_INV ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// set the values of each leaf from the implicit function
bool HashOctree::set_impl( data_access *ref )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;
  
  _max_field = -FLT_MAX ;
  for( cell_iterator it = cells_begin() ; it() ; ++it )
  {
    if( it.is_leaf() )
    {
      real v = (*ref).value_at( it.top() ) ;
      *it = v ;
      v = fabs(v) ;
      if( _max_field < v ) _max_field = v ;
    }
//    else
//      *it = R_INV ;
  }
  
//  printf( "max field: %f\n", max_field() ) ;
  return true ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Refine the octree according to the data access
bool HashOctree::refine( data_access *ref /*= NULL*/ )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;

    clear() ;  init() ;
  
  // traversal stack: iterator has a pre-ordered stack, and cannot be recursed on
  geom_cell it = geom_root() ;
  std::stack< geom_cell > s ;
  s.push( it ) ;
  
  // although it may refine only the leaves
  // the leaf iterator may change during the refinement
  bool refined = false ;
  while( !s.empty() )
  {
    it = s.top() ;  s.pop() ;

    // test for subdivision
    if( !(*ref).need_refine( it ) ) continue ;
    refined = true ;
    
    //_____________________________________________________________________________
    // subdivides the cell in the tree
    _hash[it.key()].data = it._field = R_INV ;
    HashField::KeyData kd ;
    kd.data = 0.0 ;

    // recurse on the traversal
    Key sk = it.key() << 3 ;
    for( int i = 0 ; i < 8 ; ++i )
    {
      kd.key = sk | i ;
      if( !_hash.insert( kd ) )
      {
        printf( "hash table saturated!\n" ) ;
        exit(1) ;
        return false ;
      }
      
      s.push( geom_cell( kd.key, kd.data ) ) ;
    }
    
    Level lv = it.lv()+1 ;
    if( _max_level < lv ) _max_level = lv ;
  }
  
  return refined ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Refine the octree according to the data access
bool HashOctree::adapt( data_access *ref /*= NULL*/ )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;

#ifdef HASH_HAS_ERASE
    _max_level = 0 ;
  
  // traversal stack: iterator has a pre-ordered stack, and cannot be recursed on
  geom_cell it = geom_root() ;
  std::stack< geom_cell > s ;
  s.push( it ) ;
  
  // although it may refine only the leaves
  // the leaf iterator may change during the refinement
  geom_cell sons[8] ;
  bool refined = false ;
  while( !s.empty() )
  {
    it = s.top() ;  s.pop() ;
    
    // look for the leaves
    if( !it.is_leaf() )
    {
      if( !(*ref).need_refine( it ) )
      {
        // delete branch
        std::stack<Key> del_s ;
        del_s.push( it.key() ) ;
        while( !del_s.empty() )
        {
          Key n = del_s.top() << 3 ;
          del_s.pop() ;
          
          for( int i = 0 ; i < 8 ; ++i )
          {
            HashField::KeyData kd = _hash.erase( n|i ) ;
            if( kd == HashField::KD_INV ) continue ; // is_leaf(n)
            del_s.push( n|i ) ;
          }
        }
        
        refined = true ;
        _hash[it.key()].data = 0.0 ; // leaf now!
        if( _max_level < it.lv() ) _max_level = it.lv() ;
        continue ;
      }
      
      it.sons( sons, _hash ) ;
      for( int i = 0 ; i < 8 ; ++i )
        s.push( sons[i] ) ;
      
      if( _max_level < sons[0].lv() ) _max_level = sons[0].lv() ;
      continue ;
    }
    
    // test for subdivision
    if( !(*ref).need_refine( it ) ) continue ;
    refined = true ;
    
    //_____________________________________________________________________________
    // subdivides the cell in the tree
    _hash[it.key()].data = it._field = R_INV ;
    HashField::KeyData kd ;
    kd.data = 0.0 ;
    
    // recurse on the traversal
    Key sk = it.key() << 3 ;
    for( int i = 0 ; i < 8 ; ++i )
    {
      kd.key = sk | i ;
      if( !_hash.insert( kd ) )
      {
        printf( "hash table saturated!\n" ) ;
        exit(1) ;
        return false ;
      }
      
      s.push( geom_cell( kd.key, kd.data ) ) ;
    }
    
    Level lv = it.lv()+1 ;
    if( _max_level < lv ) _max_level = lv ;
  }
  
  return refined ;  

#else  // HASH_HAS_ERASE
  // no remotion in the hash table!
  return refine( ref ) ;
#endif // HASH_HAS_ERASE
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Draw the octree with wireframe
bool HashOctree::draw_wire()
//-----------------------------------------------------------------------------
{
  for( leaf_iterator it = leaves_begin() ; it() ; ++it )
    it.draw_wire() ;
  
  return true ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Draw the octree with dots
bool HashOctree::draw_centers()
//-----------------------------------------------------------------------------
{
  for( cell_iterator it = cells_begin() ; it() ; ++it )
  {
    ::glTexCoord1f( (GLfloat)fabs(*it / max_field()) ) ;
    it.top().draw() ;
  }
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Recursion structure
typedef struct
//-----------------------------------------------------------------------------
{
  enum { DUAL_VERTEX=0, DUAL_EDGE=1, DUAL_FACE=2, DUAL_CUBE=3 } type ;
  Key cells[8]  ;
  int direction ;
} hash_dual_iterator_struct ;
//-----------------------------------------------------------------------------
// edge pattern for cubes
// direction : x -> 0, y -> 1, z -> 2, other -> -1
static const int dir_edges[19][3] = {
  {0,1,0}, {0,2,1}, {0,4,2}, {1,3,1}, {1,5,2}, {2,3,0}, {2,6,2}, {3,7,2}, {4,5,0}, {4,6,1}, {5,7,1}, {6,7,0},
  {0,5,-1}, {0,6,-1}, {1,7,-1}, {2,7,-1}, {0,3,-1}, {4,7,-1}, {0,7,-1} } ;
//_____________________________________________________________________________



#define PRINT_DUAL_DEBUG 0

//_____________________________________________________________________________
// Debug print
void print_final_cube( Key *cells )
{
#if PRINT_DUAL_DEBUG
  printf( "cube %d %d %d %d %d %d %d %d \n",
         (int)(cells[0]),
         (int)(cells[1]),
         (int)(cells[2]),
         (int)(cells[3]),
         (int)(cells[4]),
         (int)(cells[5]),
         (int)(cells[6]),
         (int)(cells[7]) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//-----------------------------------------------------------------------------
void print_dual_cube( const char*ori, Key *cells )
{
#if PRINT_DUAL_DEBUG
  printf( "%s->cube %d %d %d %d %d %d %d %d \n", ori,
         (int)(cells[0]),
         (int)(cells[1]),
         (int)(cells[2]),
         (int)(cells[3]),
         (int)(cells[4]),
         (int)(cells[5]),
         (int)(cells[6]),
         (int)(cells[7]) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//-----------------------------------------------------------------------------
void print_dual_face( const char*ori, Key *cells )
{
#if PRINT_DUAL_DEBUG
  printf( "%s->face %d %d %d %d\n", ori,
         (int)(cells[0]),
         (int)(cells[1]),
         (int)(cells[2]),
         (int)(cells[3]) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//-----------------------------------------------------------------------------
void print_dual_edge( const char*ori, Key *cells )
{
#if PRINT_DUAL_DEBUG
  printf( "%s->edge %d %d\n", ori,
         (int)(cells[0]),
         (int)(cells[1]) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//-----------------------------------------------------------------------------
void print_dual_vertex( const char*ori, Key *cells )
{
#if PRINT_DUAL_DEBUG
  printf( "%s->vertex %d\n", ori,
         (int)(cells[0]) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// Rebuild the octree dual
bool HashOctree::dual_cubes_walk( hash_dual_walker &walker )
//-----------------------------------------------------------------------------
{
  bool result = true ;
  
  
  // face pattern for a cube
  // direction : yz -> 0, zx -> 1, xy -> 2
  static const int faces[6][5] = {
    {0,1,2,3, 2}, {4,5,6,7, 2}, {0,1,4,5, 1}, {2,3,6,7, 1}, {0,2,4,6, 0}, {1,3,5,7, 0} } ;
  
  // traversal stack: iterator has a pre-ordered stack, and cannot be recursed on
  hash_dual_iterator_struct dual_it ;
  dual_it.type     = hash_dual_iterator_struct::DUAL_VERTEX ;
  dual_it.cells[0] = (Key)1 ;
  std::stack< hash_dual_iterator_struct > s ;
  s.push( dual_it ) ;
  
  _dual_temp_memory = s.size() ;
  while( !s.empty() )
  {
    int s_size = s.size() ;
    if( s_size > _dual_temp_memory ) _dual_temp_memory = s_size ;
    
    dual_it = s.top() ;
    s.pop() ;
    
    //_____________________________________________________________________________
    //_____________________________________________________________________________
    // vertex case
    if( dual_it.type == hash_dual_iterator_struct::DUAL_VERTEX )
    {
      print_dual_vertex( "!", dual_it.cells ) ;

      Key v = dual_it.cells[0] ;
      if( is_leaf(v) ) continue ;
      v <<= 3 ;
      
      //_____________________________________________________________________________
      // vertex
      dual_it.type = hash_dual_iterator_struct::DUAL_VERTEX ;
      for( int i = 0 ; i < 8 ; ++i )
      {
        dual_it.cells[0] = v | i ;
        s.push( dual_it ) ;
        print_dual_vertex( "vertex", dual_it.cells ) ;
      }
      
      //_____________________________________________________________________________
      // edge
      dual_it.type = hash_dual_iterator_struct::DUAL_EDGE ;
      for( int i = 0 ; i < 12 ; ++i )
      {
        dual_it.cells[0] = v | dir_edges[i][0] ;
        dual_it.cells[1] = v | dir_edges[i][1] ;
        dual_it.direction = dir_edges[i][2] ;
        s.push( dual_it ) ;
        print_dual_edge( "vertex", dual_it.cells ) ;
      }
      
      //_____________________________________________________________________________
      // face
      dual_it.type = hash_dual_iterator_struct::DUAL_FACE ;
      for( int i = 0 ; i < 6 ; ++i )
      {
        dual_it.cells[0] = v | faces[i][0] ;
        dual_it.cells[1] = v | faces[i][1] ;
        dual_it.cells[2] = v | faces[i][2] ;
        dual_it.cells[3] = v | faces[i][3] ;
        dual_it.direction = faces[i][4] ;
        s.push( dual_it ) ;
        print_dual_face( "vertex", dual_it.cells ) ;
      }
      
      //_____________________________________________________________________________
      // cube
      dual_it.type = hash_dual_iterator_struct::DUAL_CUBE ;
      for( int i = 0 ; i < 8 ; ++i ) dual_it.cells[i] = v | i ;
      s.push( dual_it ) ;
      print_dual_cube( "vertex", dual_it.cells ) ;
      
      continue ;
    }
    
    
    //_____________________________________________________________________________
    //_____________________________________________________________________________
    // edge case
    else if( dual_it.type == hash_dual_iterator_struct::DUAL_EDGE )
    {
      print_dual_edge( "!", dual_it.cells ) ;

      Key all_sons[2][8] ;
      Key v0 = dual_it.cells[0] << 3 ;
      Key v1 = dual_it.cells[1] << 3 ;
      if( is_leaf( dual_it.cells[0] ) )
      {
        if( is_leaf( dual_it.cells[1] ) )
          continue ;
        
        // repeat father
        for( int i = 0 ; i < 8 ; ++i )
        {
          all_sons[0][i] = dual_it.cells[0] ;
          all_sons[1][i] = v1 | i ;
        }
      }
      else
      {
        if( is_leaf( dual_it.cells[1] ) )
        {
          // repeat father
          for( int i = 0 ; i < 8 ; ++i )
          {
            all_sons[0][i] = v0 | i ;
            all_sons[1][i] = dual_it.cells[1] ;
          }
        }
        else
        {
          for( int i = 0 ; i < 8 ; ++i )
          {
            all_sons[0][i] = v0 | i ;
            all_sons[1][i] = v1 | i ;
          }
        }
      }
      
      
      //_____________________________________________________________________________
      // local cube
      Key sons[8] ;
      int dir = dual_it.direction ;
      switch( dir )
      {
          // x-aligned edge
        case 0 :
          sons[0] = all_sons[0][1] ;
          sons[1] = all_sons[1][0] ;
          sons[2] = all_sons[0][3] ;
          sons[3] = all_sons[1][2] ;
          sons[4] = all_sons[0][5] ;
          sons[5] = all_sons[1][4] ;
          sons[6] = all_sons[0][7] ;
          sons[7] = all_sons[1][6] ;
          break ;
          
          // y-aligned edge
        case 1 :
          sons[0] = all_sons[0][2] ;
          sons[1] = all_sons[0][3] ;
          sons[2] = all_sons[1][0] ;
          sons[3] = all_sons[1][1] ;
          sons[4] = all_sons[0][6] ;
          sons[5] = all_sons[0][7] ;
          sons[6] = all_sons[1][4] ;
          sons[7] = all_sons[1][5] ;
          break ;
          
          // z-aligned edge
        case 2 :
          sons[0] = all_sons[0][4] ;
          sons[1] = all_sons[0][5] ;
          sons[2] = all_sons[0][6] ;
          sons[3] = all_sons[0][7] ;
          sons[4] = all_sons[1][0] ;
          sons[5] = all_sons[1][1] ;
          sons[6] = all_sons[1][2] ;
          sons[7] = all_sons[1][3] ;
          break ;
      }
      
      //_____________________________________________________________________________
      // edge
      dual_it.direction = dir ; // same direction
      dual_it.type = hash_dual_iterator_struct::DUAL_EDGE ;
      for( int i = 0 ; i < 12 ; ++i )
      {
        if( dir_edges[i][2] != dir ) continue ;
        dual_it.cells[0] = sons[ dir_edges[i][0] ] ;
        dual_it.cells[1] = sons[ dir_edges[i][1] ] ;
        s.push( dual_it ) ;
        print_dual_edge( "edge", dual_it.cells ) ;
      }
      
      //_____________________________________________________________________________
      // face
      dual_it.type = hash_dual_iterator_struct::DUAL_FACE ;
      for( int i = 0 ; i < 6 ; ++i )
      {
        if( faces[i][4] == dir ) continue ;
        dual_it.cells[0] = sons[ faces[i][0] ] ;
        dual_it.cells[1] = sons[ faces[i][1] ] ;
        dual_it.cells[2] = sons[ faces[i][2] ] ;
        dual_it.cells[3] = sons[ faces[i][3] ] ;
        dual_it.direction = faces[i][4] ;
        s.push( dual_it ) ;
        print_dual_face( "edge", dual_it.cells ) ;
      }
      
      //_____________________________________________________________________________
      // cube
      dual_it.type = hash_dual_iterator_struct::DUAL_CUBE ;
      for( int i = 0 ; i < 8 ; ++i ) dual_it.cells[i] = sons[i] ;
      s.push( dual_it ) ;
      print_dual_cube( "edge", dual_it.cells ) ;
      
      continue ;
    }
    
    
    //_____________________________________________________________________________
    //_____________________________________________________________________________
    // face case
    else if( dual_it.type == hash_dual_iterator_struct::DUAL_FACE )
    {
      print_dual_face( "!", dual_it.cells ) ;

      Key all_sons[4][8] ;
      bool leaf = true ;
      for( int i = 0 ; i < 4 ; ++i )
      {
        if( is_leaf( dual_it.cells[i] ) )
        {
          // repeat father
          for( int j = 0 ; j < 8 ; ++j )
            all_sons[i][j] = dual_it.cells[i] ;
        }
        else
        {
          leaf = false ;
          Key f = dual_it.cells[i] << 3 ;
          for( int j = 0 ; j < 8 ; ++j )
            all_sons[i][j] = f | j ;
        }
      }
      if( leaf ) continue ;
      
      
      //_____________________________________________________________________________
      // local cube
      Key sons[8] ;
      int dir = dual_it.direction ;
      switch( dir )
      {
          // yz-aligned face
        case 0 :
          sons[0] = all_sons[0][6] ;
          sons[1] = all_sons[0][7] ;
          sons[2] = all_sons[1][4] ;
          sons[3] = all_sons[1][5] ;
          sons[4] = all_sons[2][2] ;
          sons[5] = all_sons[2][3] ;
          sons[6] = all_sons[3][0] ;
          sons[7] = all_sons[3][1] ;
          break ;
          
          // zx-aligned edge
        case 1 :
          sons[0] = all_sons[0][5] ;
          sons[1] = all_sons[1][4] ;
          sons[2] = all_sons[0][7] ;
          sons[3] = all_sons[1][6] ;
          sons[4] = all_sons[2][1] ;
          sons[5] = all_sons[3][0] ;
          sons[6] = all_sons[2][3] ;
          sons[7] = all_sons[3][2] ;
          break ;
          
          // xy-aligned edge
        case 2 :
          sons[0] = all_sons[0][3] ;
          sons[1] = all_sons[1][2] ;
          sons[2] = all_sons[2][1] ;
          sons[3] = all_sons[3][0] ;
          sons[4] = all_sons[0][7] ;
          sons[5] = all_sons[1][6] ;
          sons[6] = all_sons[2][5] ;
          sons[7] = all_sons[3][4] ;
          break ;
      }
      
      //_____________________________________________________________________________
      // face
      // same direction
      dual_it.direction = dir ;
      dual_it.type = hash_dual_iterator_struct::DUAL_FACE ;
      for( int i = 0 ; i < 6 ; ++i )
      {
        if( faces[i][4] != dir ) continue ;
        dual_it.cells[0] = sons[ faces[i][0] ] ;
        dual_it.cells[1] = sons[ faces[i][1] ] ;
        dual_it.cells[2] = sons[ faces[i][2] ] ;
        dual_it.cells[3] = sons[ faces[i][3] ] ;
        //        dual_it.direction = faces[i][4] ;
        s.push( dual_it ) ;
        print_dual_face( "face", dual_it.cells ) ;
      }
      
      //_____________________________________________________________________________
      // cube
      dual_it.type = hash_dual_iterator_struct::DUAL_CUBE ;
      for( int i = 0 ; i < 8 ; ++i ) dual_it.cells[i] = sons[i] ;
      s.push( dual_it ) ;
      print_dual_cube( "face", dual_it.cells ) ;
      
      continue ;
    }
    
    
    //_____________________________________________________________________________
    //_____________________________________________________________________________
    // cube case
    else if( dual_it.type == hash_dual_iterator_struct::DUAL_CUBE )
    {
      print_dual_cube( "!", dual_it.cells ) ;

      bool leaf = true ;
      
      //_____________________________________________________________________________
      // cube
      for( int i = 0 ; i < 8 ; ++i )
      {
        // copy father
        if( is_leaf( dual_it.cells[i] ) )
          continue ;
        
        // else
        leaf = false ;
        dual_it.cells[i] = ( dual_it.cells[i]<<3 ) | 7-i ;
      }
      
      // finally add dual cube
      if( leaf )
      {
        result &= walker( *this, dual_it.cells ) ;
        print_final_cube( dual_it.cells ) ;
      }
      else
      {
        s.push( dual_it ) ;
        print_dual_cube( "cube", dual_it.cells ) ;
      }
      
      continue ;
    }
    
  }
  
  return result ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Isosurface walker
bool hash_isosurface_walker( HashOctree &fo, Key *keys )
//-----------------------------------------------------------------------------
{
  MarchingCubes &mc = fo.mc() ;
  
  /** isovalues at the cube vertices */
  real            *cube    = mc.cube() ;
  /** geometry of the cube */
  Point           *space   = mc.space() ;
  /** indexes of the cube */
  Key             *indexes = mc.indexes() ;
  
  for( int j = 0 ; j < 8 ; ++j )
  {
    int i = (j==2)?3:(j==3)?2:(j==7)?6:(j==6)?7:j ;
    const Key k = keys[j] ;
    HashOctree::geom_cell c = fo.geom_key(k) ;
    cube[i] = *c ;
    
    space[i] = (Point)c ;
    indexes[i] = k ;
  }
  
  mc.tesselate_cube( 0.0 ) ;
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Build the isosurface using dual marching cubes
bool HashOctree::build_isosurface( data_access *ref )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;
  _mc.dat_access() = ref ;

    printf( "Build isosurface... " ) ;
  bool result = false ;
  
  //_________________________________________________________________________
  // Dual Marching Cubes
  
  _mc.init_all() ; 
  
  result = dual_cubes_walk( hash_isosurface_walker ) ;
  
  _mc.clean_temps() ;
  
  printf( "generated %d vertices and %d faces!\n", _mc.nverts(), _mc.ntrigs() ) ;
  
  return result ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Isosurface drawing walker
bool hash_direct_isosurface_walker( HashOctree &fo, Key *keys )
//-----------------------------------------------------------------------------
{
  MC_Draw &mc = fo.mc_draw() ;
  
  /** isovalues at the cube vertices */
  real            *cube    = mc.cube() ;
  /** geometry of the cube */
  Point           *space   = mc.space() ;
  
  for( int j = 0 ; j < 8 ; ++j )
  {
    int i = (j==2)?3:(j==3)?2:(j==7)?6:(j==6)?7:j ;
    const HashOctree::geom_cell c = fo.geom_key( keys[j] ) ;
    cube[i] = *c ;
    space[i] = (Point)c ;
  }
  
  mc.tesselate_cube( 0.0 ) ;
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Build the isosurface using dual marching cubes
bool HashOctree::direct_draw_isosurface( data_access *ref )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;
  _mc_draw.dat_access() = ref ;
  
    bool result = false ;
  
  // Dual Marching Cubes
  ::glBegin( GL_TRIANGLES ) ;
  {
    result = dual_cubes_walk( hash_direct_isosurface_walker ) ;
  }
  ::glEnd() ; // GL_TRIANGLES
  
  return result ;
}
//_____________________________________________________________________________






//_____________________________________________________________________________
// Draw walker
bool hash_draw_walker( HashOctree &fo, Key *keys )
//-----------------------------------------------------------------------------
{
  Cube cells[8] ;
  for( int i = 0 ; i < 8 ; ++i )
    cells[i] = key2cube( keys[i] ) ;
    
  ((Point)cells[0]).draw() ;
  ((Point)cells[1]).draw() ;
  
  ((Point)cells[1]).draw() ;
  ((Point)cells[3]).draw() ;
  
  ((Point)cells[3]).draw() ;
  ((Point)cells[2]).draw() ;
  
  ((Point)cells[2]).draw() ;
  ((Point)cells[0]).draw() ;
  
  ((Point)cells[4]).draw() ;
  ((Point)cells[5]).draw() ;
  
  ((Point)cells[5]).draw() ;
  ((Point)cells[7]).draw() ;
  
  ((Point)cells[7]).draw() ;
  ((Point)cells[6]).draw() ;
  
  ((Point)cells[6]).draw() ;
  ((Point)cells[4]).draw() ;
  
  ((Point)cells[0]).draw() ;
  ((Point)cells[4]).draw() ;
  
  ((Point)cells[1]).draw() ;
  ((Point)cells[5]).draw() ;
  
  ((Point)cells[2]).draw() ;
  ((Point)cells[6]).draw() ;
  
  ((Point)cells[3]).draw() ;
  ((Point)cells[7]).draw() ;
  
  return true ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Draw the dual octree with wireframe
bool HashOctree::draw_dual()
//-----------------------------------------------------------------------------
{
  ::glBegin( GL_LINES ) ;
  bool result = dual_cubes_walk( hash_draw_walker ) ;
  ::glEnd() ;  // GL_LINES
  
  return result ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Do nothing walker, just for dual generation timing
bool hash_timing_walker( HashOctree &fo, Key *keys )
//-----------------------------------------------------------------------------
{
  return true ;
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// Do nothing, just for dual generation timing
bool HashOctree::dual_timing()
//-----------------------------------------------------------------------------
{
    return dual_cubes_walk( hash_timing_walker ) ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Prints statistics about the octree
void HashOctree::stats()
//-----------------------------------------------------------------------------
{
  uint level_dist[MAX_LEVEL+1] ;
  memset( level_dist, 0, (MAX_LEVEL+1)*sizeof(uint) ) ;
  
  int s = 0 ;
  for( cell_iterator it = cells_begin() ; it() ; ++it )
  {
    ++s ;
    if( it.is_leaf() ) level_dist[ it.lv() ] ++ ;
  }
  
  printf( " number of nodes:\t%d\n", s ) ;
#if USE_HASH_PTR
  printf( " total memory:\t%d\n", (int) ( s * sizeof(HashField::KeyBucket) ) ) ;
  _hash.stats() ;
#else  // USE_HASH_PTR 
  printf( " total memory:\t%d\n", (int) ( s * sizeof(HashField::KeyData  ) ) ) ;
#endif // USE_HASH_PTR
  
  printf( " leaves' levels:\t" ) ;
  for( Level l = 0 ; l <= MAX_LEVEL; ++l )
    printf( "\t%d", level_dist[l] ) ;
  printf( "\n" ) ;

  printf( " max size of the stack:\t%d\n", _dual_temp_memory ) ;
  printf( " total memory of stack:\t%d\n", (int) ( _dual_temp_memory * sizeof(hash_dual_iterator_struct) ) ) ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
//_____________________________________________________________________________
// search operations


//_____________________________________________________________________________
// Find cells of the octree at a given position
bool HashOctree::find_leaf( real x, real y, real z, geom_cell &n ) const
//-----------------------------------------------------------------------------
{
    
  n = geom_root() ;
  
  while( !n.is_leaf() )
  {
    ++n.lv() ;
    real sz = n.sz() ;
    
    // side test :
    //  0->back  lower left, 1->back  lower right, 2->back  upper left, 3->back  upper right
    //  4->front lower left, 5->front lower right, 6->front upper left, 7->front upper right
    int i = ( (n.cz() < z) << 2 ) | ( (n.cy() < y) << 1 ) | (n.cx() < x) ;
    n.cx() += (i&1) ? +sz : -sz ;
    n.cy() += (i&2) ? +sz : -sz ;
    n.cz() += (i&4) ? +sz : -sz ;
    
    geom_cell s ;
    n.son(i,s,_hash) ;
    n = s ;
  }
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Find cells of the octree inside a given box of center x,y and half side r
bool HashOctree::find_radius( real x, real y, real z, real r, List<geom_cell> &cells ) const
//-----------------------------------------------------------------------------
{
    
  cells.clear() ;
  
  geom_cell n = geom_root() ;
  
  std::stack< geom_cell > s ;
  s.push( n ) ;
  
  while( !s.empty() )
  {
    n = s.top() ;  s.pop() ;
//    printf( "%d\t%.4f\t%.4f\t%.4f\n", (int)n.lv(), n.x(), n.y(), n.z() ) ;
    if( n.is_leaf() )
    {
      cells.insert( n ) ;
      continue ;
    }
    
    geom_cell sons[8] ;
    n.sons( sons, _hash ) ;
    real sz_ = sons[0].sz() + r ;
    for( int i = 0 ; i < 8 ; ++i )
    {
      if( (fabs( sons[i].cx() - x ) < sz_ ) &&
          (fabs( sons[i].cy() - y ) < sz_ ) &&
          (fabs( sons[i].cz() - z ) < sz_ ) )
        s.push( sons[i] ) ;
    }
  }
//  printf( "\n" ) ;
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Find neighbor cells of the octree to a given cell
bool HashOctree::adjacent( const geom_cell &cell, List<geom_cell> &cells ) const
//-----------------------------------------------------------------------------
{
  return find_radius( cell.cx(), cell.cy(), cell.cz(), cell.sz() + (real) 0.5 / (1 << max_level() ) , cells ) ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
//_____________________________________________________________________________
// I/O

//_____________________________________________________________________________
// Draws the intersecting cells of a plane with the octree with color
void HashOctree::draw_plane ( real nx, real ny, real nz, real d )
//-----------------------------------------------------------------------------
{
  geom_cell n = geom_root() ;
  
  std::stack< geom_cell > s ;
  s.push( n ) ;
  
  //_____________________________________________________________________________
  // retrieve the index of the most positive cube corner
  bool fx = nx >= 0 ;
  bool fy = ny >= 0 ;
  bool fz = nz >= 0 ;
  
  //_____________________________________________________________________________
  
  ::glLineWidth( 1.0f ) ;
  ::glDisable( GL_LIGHTING ) ;
  
  real px,py,pz, mx,my,mz ;
  while( !s.empty() )
  {
    n = s.top() ;  s.pop() ;
    
    if( n.is_leaf() )
    {		
      
      ::glTexCoord1f( (GLfloat)fabs(*n/max_field()) ) ;
      ::glBegin( GL_LINES ) ;
      n.draw_wire() ;
      ::glEnd() ;  // GL_LINES		
      
      continue ;
    }
    
    //___________________________________________________________________________
    // recursion
    geom_cell sons[8] ;
    n.sons( sons, _hash ) ;
    real sz_ = sons[0].sz() ;
    for( int i = 0 ; i < 8 ; ++i )
    {
      // get the coordinates of the extremum vertices
      if( fx ) { px = sons[i].cx() + sz_ ;  mx = sons[i].cx() - sz_ ; }
      else     { px = sons[i].cx() - sz_ ;  mx = sons[i].cx() + sz_ ; }
      if( fy ) { py = sons[i].cy() + sz_ ;  my = sons[i].cy() - sz_ ; }
      else     { py = sons[i].cy() - sz_ ;  my = sons[i].cy() + sz_ ; }
      if( fz ) { pz = sons[i].cz() + sz_ ;  mz = sons[i].cz() - sz_ ; }
      else     { pz = sons[i].cz() - sz_ ;  mz = sons[i].cz() + sz_ ; }
      
      // plane/cube intersection test
      if( (nx * px + ny * py + nz * pz + d >= 0.0) && (nx * mx + ny * my + nz * mz + d <= 0.0)  )
        s.push( sons[i] ) ;
    }
  }
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Draws the intersection of a plane with the octree with color
void HashOctree::draw_slice ( real nx, real ny, real nz, real d, float alpha )
//-----------------------------------------------------------------------------
{
  geom_cell n = geom_root() ;
  
  std::stack< geom_cell > s ;
  s.push( n ) ;
  
  // retrieve the index of the most positive cube corner
  bool fx = nx >= 0 ;
  bool fy = ny >= 0 ;
  bool fz = nz >= 0 ;
  
  
  
  //_____________________________________________________________________________
  
  ::glDisable( GL_LIGHTING ) ;
  
  real px,py,pz, mx,my,mz ;
  while( !s.empty() )
  {
    n = s.top() ;  s.pop() ;
    if( n.is_leaf() )
    {
      // retrieves the signs of each cube corner
      real val[8] ;
      real v   = nx * n.cx() + ny * n.cy() + nz * n.cz() + d ;
      real sz_ = n.sz() ;
      int zeros = 0 ;
      for( int i = 0 ; i < 8 ; ++i )
      {
        val[i] = v + nx * ((i&1)?sz_:-sz_) + ny * ((i&2)?sz_:-sz_) + nz * ((i&4)?sz_:-sz_) ;
        if( fabs(val[i]) < R_EPSILON ) { val[i] = R_EPSILON ;  ++zeros ; }
      }
      
      //_________________________________________________________________________
      // degenerated intersection
      if( zeros > 0 )
      {
        int npos = 0, nneg = 0 ;
        for( int i = 0 ; i < 8 ; ++i )
        {
          if( fabs(val[i]) < R_EPSILON ) continue ;
          if( val[i] < 0 )
            ++nneg ;
          else
            ++npos ;
        }
        float def = (nneg>npos) ? (float)R_EPSILON : -(float)R_EPSILON ;
        
        for( int i = 0 ; i < 8 ; ++i )
        {
          if( fabs(val[i]) < R_EPSILON ) val[i] = def ;
        }
        //        continue ;
      }
      
      //_________________________________________________________________________
      // Continuation method
      int v0=-1, v1=-1 ;  // current edge vertices
      
      // get the first intersection
      static int edges[12][2] = { {0,1}, {0,2}, {0,4}, {1,3}, {1,5}, {2,3}, {2,6}, {3,7}, {4,5}, {4,6}, {5,7}, {6,7} } ;
      for( int e = 0 ; e < 12 ; ++e )
      {
        if( val[ edges[e][0] ] * val[ edges[e][1] ] < 0 )
        {
          v0 = edges[e][0] ;
          v1 = edges[e][1] ;
          break ;
        }
      }
      int  ov0 = v0, ov1 = v1 ;
      
      
      // equation of previous face adjacent to v0v1: v>>c & 1 == s (c-th coordinate equal to +/- 1)
      int pf = -1 ;  bool ps ;
      if       ( (v0 & 1) == (v1 & 1) )   { pf = 0 ; ps = ((v0 & 1) /*>> 0*/) == 1 ; }
      else   if( (v0 & 2) == (v1 & 2) )   { pf = 1 ; ps = ((v0 & 2)   >> 1  ) == 1 ; }
      else /*if( (v0 & 4) == (v1 & 4) )*/ { pf = 2 ; ps = ((v0 & 4)   >> 2  ) == 1 ; }
      
      //      printf( "cube %d %d %d %d %d %d %d %d, march : (%d,%d) ", val[0]>0, val[1]>0, val[2]>0, val[3]>0, val[4]>0, val[5]>0, val[6]>0, val[7]>0, v0, v1 ) ;
      
      //_________________________________________________________________________
      // marching
      
      GLfloat color = (GLfloat)fabs(*n / max_field()) ;
		  //cmap.gl_color( at(*n) ) ;	
      ::glBegin( GL_POLYGON ) ;
      {
        do
        {
          // equation of next face adjacent to v0v1
          int nf = -1 ;  bool ns ;
          if       ( pf != 0 && ((v0 & 1) == (v1 & 1)) )   { nf = 0 ; ns = ((v0 & 1) /*>> 0*/) == 1 ; }
          else   if( pf != 1 && ((v0 & 2) == (v1 & 2)) )   { nf = 1 ; ns = ((v0 & 2)   >> 1  ) == 1 ; }
          else /*if( pf != 2 && ((v0 & 4) == (v1 & 4)) )*/ { nf = 2 ; ns = ((v0 & 4)   >> 2  ) == 1 ; }
          
          /* equivalent construction
           // get the direction perpendicular to v0v1
           int of = (v0 ^ v1) & 7 ; // only bit of difference
           of = of==1 ? 0 : (of==2 ? 1 : 2) ;
           
           // compute the equation of the next face adjacent to v0v1
           int  nf = 3 - (pf+of) ;    // third option: {pf,nf,of}={0,1,2}
           bool ns = (v0 >> nf) & 1 ; // value of the edge on that direction
           */
          
          // compute the two next vertices (invert the pf-th bit
          int v2 = ps ? ( v0 & (~(1<<pf)) ) : ( v0 | (1<<pf) ) ;
          int v3 = ps ? ( v1 & (~(1<<pf)) ) : ( v1 | (1<<pf) ) ;
          
          // get the characteristic of the square
          int  c = (val[v3] > 0) | ((val[v2] > 0) << 1) | ((val[v1] > 0) << 2) | ((val[v0] > 0) << 3) ;
          if( c > 7 ) c = (~c)&7 ; // consider v0 as negative, and thus v1 as positive
          /*
           if( c < 4 )
           { // degenerated intesection
           printf( "degenerated intersection %d!\n", c ) ;
           }
           */
          
          //_______________________________________________________________________
          // method of continuity on the cube surface
          switch( c )
          {
            case  4 : v0 = v3 ;  break ;
            case  5 : v0 = v2 ;  v1 = v3 ;  break ;
            case  6 : /* impossible case */ break ;
            case  7 : v1 = v2 ;  break ;
          }
          
          
          
          
          real cx =  n.cx();
          real cy =  n.cy(); 
          real cz =  n.cz();
          real sz =  sz_;   
          real a0 = val[v0];  
          int i0 = v0; 
          real a1 = val[v1];  
          int i1  = v1;
          
          real w0 = a1 / (a1-a0) ;
          real w1 = a0 / (a0-a1) ;
          
          real x = cx + w0 * ((i0&1)?sz:-sz) + w1 * ((i1&1)?sz:-sz);
          real y = cy + w0 * ((i0&2)?sz:-sz) + w1 * ((i1&2)?sz:-sz);
          real z = cz + w0 * ((i0&4)?sz:-sz) + w1 * ((i1&4)?sz:-sz);
          
          ::glTexCoord1d( color ) ;
          ::glVertex3d(x, y, z) ;
          
          // next face
          pf = nf ;  ps = ns ;
          
          //          printf( "(%d,%d) ", v0, v1 ) ;
        } while( v0 != ov0 || v1 != ov1 ) ;
        //        printf( "\n" ) ;
      }
      ::glEnd() ; // GL_POLYGON
      
      continue ;
    }
    
    //___________________________________________________________________________
    // recursion
    geom_cell sons[8] ;
    n.sons( sons, _hash ) ;
    real sz_ = sons[0].sz() ;
    for( int i = 0 ; i < 8 ; ++i )
    {
      // get the coordinates of the extremum vertices
      if( fx ) { px = sons[i].cx() + sz_ ;  mx = sons[i].cx() - sz_ ; }
      else     { px = sons[i].cx() - sz_ ;  mx = sons[i].cx() + sz_ ; }
      if( fy ) { py = sons[i].cy() + sz_ ;  my = sons[i].cy() - sz_ ; }
      else     { py = sons[i].cy() - sz_ ;  my = sons[i].cy() + sz_ ; }
      if( fz ) { pz = sons[i].cz() + sz_ ;  mz = sons[i].cz() - sz_ ; }
      else     { pz = sons[i].cz() - sz_ ;  mz = sons[i].cz() + sz_ ; }
      
      // plane/cube intersection test
      if( (nx * px + ny * py + nz * pz + d >= 0.0) && (nx * mx + ny * my + nz * mz + d <= 0.0)  )
        s.push( sons[i] ) ;
    }
  }
}
//_____________________________________________________________________________





