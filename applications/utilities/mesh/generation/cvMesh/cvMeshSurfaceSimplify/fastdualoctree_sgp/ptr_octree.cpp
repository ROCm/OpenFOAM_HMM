/**
 * \file    ptr_octree.cpp
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 *  Octree structure with son-brother pointers
 */
//_____________________________________________________________________________


#ifndef WIN32
#pragma implementation
#endif // WIN32


#include "ptr_octree.h"


//_____________________________________________________________________________
// Create the root
void PtrOctree::init()
//-----------------------------------------------------------------------------
{
  _root.brother   = NULL ;
  _root.first_son = NULL ;
  _root.field     = 0.0  ;  
  
  _max_level = 0 ;
  _max_field = R_INV ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Check that each point scan or scribble is in the right node
bool PtrOctree::check ()
//-----------------------------------------------------------------------------
{
  bool result = true ;
  return result ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Delete the branch below a given node
void PtrOctree::clear_branch( cell *b )
//-----------------------------------------------------------------------------
{
  std::stack<cell*> s ;
  s.push( b ) ;
  while( !s.empty() )
  {
    cell *n = s.top() ;
    cell *son  = n->first_son ;
    cell *prev = n ;
    while( son != (cell*)NULL ) // look two levels to be able to update pointers
    {
      cell *next = son->brother ;
      if( !son->first_son ) // son is a leaf
      {
        delete son ;
        if( prev == n )
          n->first_son = next ;
        else
          prev->brother = next ;
      }
      else
      {
        prev = son ;
        s.push( son ) ;
      }
      son = next ;
    }
    if( !n->first_son ) // now it is a leaf
      s.pop() ;
  }
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// set the values of each leaf from the implicit function
bool PtrOctree::set_impl( data_access *ref )
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
    else
      *it = R_INV ;
  }
  
//  printf( "max field: %f\n", max_field() ) ;
  return true ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Refine the octree according to the data access
bool PtrOctree::refine( data_access *ref /*= NULL*/ )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;

    clear() ;  init() ;
  
  // traversal stack: iterator has a pre-ordered stack, and cannot be recursed on
  geom_cell it( &_root ) ;
  std::stack< geom_cell > s ;
  s.push( it ) ;
  
  // although it may refine only the leaves
  // the leaf iterator may change during the refinement
  bool refined = false ;
  while( !s.empty() )
  {
    it = s.top() ;  s.pop() ;
    geom_cell sons[8] ;
    
    // test for subdivision
    if( !(*ref).need_refine( it ) ) continue ;
    refined = true ;
    
    //_____________________________________________________________________________
    // subdivides the cell in the tree
    for( int i = 0 ; i < 8 ; ++i )
    {
      sons[i].cell     () = new cell ;
      sons[i].first_son() = (cell*)NULL ;
      sons[i].brother  () = (cell*)NULL ;
      sons[i].cell()->field  = 0.0 ;
      
      if( i == 0 ) it.first_son() = sons[i].cell() ;
      else sons[i-1].brother()    = sons[i].cell() ;
    }
    it.sons( sons ) ;
    //if( _max_level < it.lv() ) _max_level = it.lv() ;
    if( _max_level < sons[0].lv() ) _max_level = sons[0].lv() ;
    
    // recurse on the traversal
    for( int i = 0 ; i < 8 ; ++i )
      s.push( sons[i] ) ;
  }
  
  return refined ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Refine the octree according to the data access
bool PtrOctree::adapt( data_access *ref /*= NULL*/ )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;

    _max_level = 0 ;
  
  // traversal stack: iterator has a pre-ordered stack, and cannot be recursed on
  geom_cell it( &_root ) ;
  std::stack< geom_cell > s ;
  s.push( it ) ;
  
  // although it may refine only the leaves
  // the leaf iterator may change during the refinement
  bool refined = false ;
  while( !s.empty() )
  {
    it = s.top() ;  s.pop() ;
    geom_cell sons[8] ;

    // look for the leaves
    if( !it.is_leaf() )
    {
      if( !(*ref).need_refine( it ) )
      {
        clear_branch( it.cell() ) ;
        refined = true ;
        continue ;
      }

      it.sons( sons ) ;
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
    for( int i = 0 ; i < 8 ; ++i )
    {
      sons[i].cell     () = new cell ;
      sons[i].first_son() = (cell*)NULL ;
      sons[i].brother  () = (cell*)NULL ;
      sons[i].cell()->field  = 0.0 ;
      
      if( i == 0 ) it.first_son() = sons[i].cell() ;
      else sons[i-1].brother()    = sons[i].cell() ;
    }
    it.sons( sons ) ;
    //if( _max_level < it.lv() ) _max_level = it.lv() ;
    if( _max_level < sons[0].lv() ) _max_level = sons[0].lv() ;
    
    // recurse on the traversal
    for( int i = 0 ; i < 8 ; ++i )
      s.push( sons[i] ) ;
  }
  
  return refined ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Draw the octree with wireframe
bool PtrOctree::draw_wire()
//-----------------------------------------------------------------------------
{
  for( leaf_iterator it = leaves_begin() ; it() ; ++it )
    it.draw_wire() ;
  
  return true ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Draw the octree with dots
bool PtrOctree::draw_centers()
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
  PtrOctree::geom_cell  cells[8]  ;
  int                   direction ;
} ptr_dual_iterator_struct ;
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
void print_final_cube( PtrOctree::geom_cell *cells, float M )
{
#if PRINT_DUAL_DEBUG
  printf( "cube (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d) \n",
         (int)(cells[0].cx()*M), (int)(cells[0].cy()*M), (int)(cells[0].cz()*M),
         (int)(cells[1].cx()*M), (int)(cells[1].cy()*M), (int)(cells[1].cz()*M),
         (int)(cells[2].cx()*M), (int)(cells[2].cy()*M), (int)(cells[2].cz()*M),
         (int)(cells[3].cx()*M), (int)(cells[3].cy()*M), (int)(cells[3].cz()*M),
         (int)(cells[4].cx()*M), (int)(cells[4].cy()*M), (int)(cells[4].cz()*M),
         (int)(cells[5].cx()*M), (int)(cells[5].cy()*M), (int)(cells[5].cz()*M),
         (int)(cells[6].cx()*M), (int)(cells[6].cy()*M), (int)(cells[6].cz()*M),
         (int)(cells[7].cx()*M), (int)(cells[7].cy()*M), (int)(cells[7].cz()*M) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//-----------------------------------------------------------------------------
void print_dual_cube( const char*ori, PtrOctree::geom_cell *cells, float M )
{
#if PRINT_DUAL_DEBUG
  printf( "%s->cube (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d) \n", ori,
         (int)(cells[0].cx()*M), (int)(cells[0].cy()*M), (int)(cells[0].cz()*M),
         (int)(cells[1].cx()*M), (int)(cells[1].cy()*M), (int)(cells[1].cz()*M),
         (int)(cells[2].cx()*M), (int)(cells[2].cy()*M), (int)(cells[2].cz()*M),
         (int)(cells[3].cx()*M), (int)(cells[3].cy()*M), (int)(cells[3].cz()*M),
         (int)(cells[4].cx()*M), (int)(cells[4].cy()*M), (int)(cells[4].cz()*M),
         (int)(cells[5].cx()*M), (int)(cells[5].cy()*M), (int)(cells[5].cz()*M),
         (int)(cells[6].cx()*M), (int)(cells[6].cy()*M), (int)(cells[6].cz()*M),
         (int)(cells[7].cx()*M), (int)(cells[7].cy()*M), (int)(cells[7].cz()*M) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//-----------------------------------------------------------------------------
void print_dual_face( const char*ori, PtrOctree::geom_cell *cells, float M )
{
#if PRINT_DUAL_DEBUG
  printf( "%s->face (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)  (%d,%d,%d)\n", ori,
         (int)(cells[0].cx()*M), (int)(cells[0].cy()*M), (int)(cells[0].cz()*M),
         (int)(cells[1].cx()*M), (int)(cells[1].cy()*M), (int)(cells[1].cz()*M),
         (int)(cells[2].cx()*M), (int)(cells[2].cy()*M), (int)(cells[2].cz()*M),
         (int)(cells[3].cx()*M), (int)(cells[3].cy()*M), (int)(cells[3].cz()*M) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//-----------------------------------------------------------------------------
void print_dual_edge( const char*ori, PtrOctree::geom_cell *cells, float M )
{
#if PRINT_DUAL_DEBUG
  printf( "%s->edge (%d,%d,%d)  (%d,%d,%d)\n", ori,
         (int)(cells[0].cx()*M), (int)(cells[0].cy()*M), (int)(cells[0].cz()*M),
         (int)(cells[1].cx()*M), (int)(cells[1].cy()*M), (int)(cells[1].cz()*M) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//-----------------------------------------------------------------------------
void print_dual_vertex( const char*ori, PtrOctree::geom_cell *cells, float M )
{
#if PRINT_DUAL_DEBUG
  printf( "%s->vertex (%d,%d,%d)\n", ori,
         (int)(cells[0].cx()*M), (int)(cells[0].cy()*M), (int)(cells[0].cz()*M) ) ;
#endif  // PRINT_DUAL_DEBUG
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// Rebuild the octree dual
bool PtrOctree::dual_cubes_walk( ptr_dual_walker &walker )
//-----------------------------------------------------------------------------
{
  bool result = true ;
  
  float M = 1 << (max_level()+1) ;
  
  
  // face pattern for a cube
  // direction : yz -> 0, zx -> 1, xy -> 2
  static const int faces[6][5] = {
    {0,1,2,3, 2}, {4,5,6,7, 2}, {0,1,4,5, 1}, {2,3,6,7, 1}, {0,2,4,6, 0}, {1,3,5,7, 0} } ;
  
  // traversal stack: iterator has a pre-ordered stack, and cannot be recursed on
  ptr_dual_iterator_struct dual_it ;
  dual_it.type     = ptr_dual_iterator_struct::DUAL_VERTEX ;
  dual_it.cells[0] = geom_cell( &_root ) ;
  std::stack< ptr_dual_iterator_struct > s ;
  s.push( dual_it ) ;
  
  _dual_temp_memory = s.size() ;
  while( !s.empty() )
  {
    int s_size = s.size() ;
    if( s_size > _dual_temp_memory ) _dual_temp_memory = s_size ;

    dual_it = s.top() ;  s.pop() ;
    
    //_____________________________________________________________________________
    //_____________________________________________________________________________
    // vertex case
    if( dual_it.type == ptr_dual_iterator_struct::DUAL_VERTEX )
    {
      if( dual_it.cells[0].is_leaf() ) continue ;
      
      geom_cell sons[8] ;
      dual_it.cells[0].sons( sons ) ;
      
      //_____________________________________________________________________________
      // vertex
      dual_it.type = ptr_dual_iterator_struct::DUAL_VERTEX ;
      for( int i = 0 ; i < 8 ; ++i )
      {
        dual_it.cells[0] = sons[i] ;
        s.push( dual_it ) ;
        print_dual_vertex( "vertex", dual_it.cells, M ) ;
      }
      
      //_____________________________________________________________________________
      // edge
      dual_it.type = ptr_dual_iterator_struct::DUAL_EDGE ;
      for( int i = 0 ; i < 12 ; ++i )
      {
        dual_it.cells[0] = sons[ dir_edges[i][0] ] ;
        dual_it.cells[1] = sons[ dir_edges[i][1] ] ;
        dual_it.direction = dir_edges[i][2] ;
        s.push( dual_it ) ;
        print_dual_edge( "vertex", dual_it.cells, M ) ;
      }
      
      //_____________________________________________________________________________
      // face
      dual_it.type = ptr_dual_iterator_struct::DUAL_FACE ;
      for( int i = 0 ; i < 6 ; ++i )
      {
        dual_it.cells[0] = sons[ faces[i][0] ] ;
        dual_it.cells[1] = sons[ faces[i][1] ] ;
        dual_it.cells[2] = sons[ faces[i][2] ] ;
        dual_it.cells[3] = sons[ faces[i][3] ] ;
        dual_it.direction = faces[i][4] ;
        s.push( dual_it ) ;
        print_dual_face( "vertex", dual_it.cells, M ) ;
      }
      
      //_____________________________________________________________________________
      // cube
      dual_it.type = ptr_dual_iterator_struct::DUAL_CUBE ;
      for( int i = 0 ; i < 8 ; ++i ) dual_it.cells[i] = sons[i] ;
      s.push( dual_it ) ;
      print_dual_cube( "vertex", dual_it.cells, M ) ;
      
      continue ;
    }
    
    
    //_____________________________________________________________________________
    //_____________________________________________________________________________
    // edge case
    else if( dual_it.type == ptr_dual_iterator_struct::DUAL_EDGE )
    {
      geom_cell all_sons[2][8] ;
      if( dual_it.cells[0].is_leaf() )
      {
        if( dual_it.cells[1].is_leaf() )
          continue ;
        
        // repeat father
        for( int i = 0 ; i < 8 ; ++i )
          all_sons[0][i] = dual_it.cells[0] ;
        
        dual_it.cells[1].sons( all_sons[1] ) ;
      }
      else
      {
        dual_it.cells[0].sons( all_sons[0] ) ;
        
        if( dual_it.cells[1].is_leaf() )
        {
          // repeat father
          for( int i = 0 ; i < 8 ; ++i )
            all_sons[1][i] = dual_it.cells[1] ;
        }
        else
          dual_it.cells[1].sons( all_sons[1] ) ;
      }
      
      
      //_____________________________________________________________________________
      // local cube
      geom_cell sons[8] ;
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
      dual_it.type = ptr_dual_iterator_struct::DUAL_EDGE ;
      for( int i = 0 ; i < 12 ; ++i )
      {
        if( dir_edges[i][2] != dir ) continue ;
        dual_it.cells[0] = sons[ dir_edges[i][0] ] ;
        dual_it.cells[1] = sons[ dir_edges[i][1] ] ;
        s.push( dual_it ) ;
        print_dual_edge( "edge", dual_it.cells, M ) ;
      }
      
      //_____________________________________________________________________________
      // face
      dual_it.type = ptr_dual_iterator_struct::DUAL_FACE ;
      for( int i = 0 ; i < 6 ; ++i )
      {
        if( faces[i][4] == dir ) continue ;
        dual_it.cells[0] = sons[ faces[i][0] ] ;
        dual_it.cells[1] = sons[ faces[i][1] ] ;
        dual_it.cells[2] = sons[ faces[i][2] ] ;
        dual_it.cells[3] = sons[ faces[i][3] ] ;
        dual_it.direction = faces[i][4] ;
        s.push( dual_it ) ;
        print_dual_face( "edge", dual_it.cells, M ) ;
      }
      
      //_____________________________________________________________________________
      // cube
      dual_it.type = ptr_dual_iterator_struct::DUAL_CUBE ;
      for( int i = 0 ; i < 8 ; ++i ) dual_it.cells[i] = sons[i] ;
      s.push( dual_it ) ;
      print_dual_cube( "edge", dual_it.cells, M ) ;
      
      continue ;
    }
    
    
    //_____________________________________________________________________________
    //_____________________________________________________________________________
    // face case
    else if( dual_it.type == ptr_dual_iterator_struct::DUAL_FACE )
    {
      geom_cell all_sons[4][8] ;
      bool leaf = true ;
      for( int i = 0 ; i < 4 ; ++i )
      {
        if( dual_it.cells[i].is_leaf() )
        {
          // repeat father
          for( int j = 0 ; j < 8 ; ++j )
            all_sons[i][j] = dual_it.cells[i] ;
        }
        else
        {
          leaf = false ;
          dual_it.cells[i].sons( all_sons[i] ) ;
        }
      }
      if( leaf ) continue ;
      
      
      //_____________________________________________________________________________
      // local cube
      geom_cell sons[8] ;
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
      dual_it.type = ptr_dual_iterator_struct::DUAL_FACE ;
      for( int i = 0 ; i < 6 ; ++i )
      {
        if( faces[i][4] != dir ) continue ;
        dual_it.cells[0] = sons[ faces[i][0] ] ;
        dual_it.cells[1] = sons[ faces[i][1] ] ;
        dual_it.cells[2] = sons[ faces[i][2] ] ;
        dual_it.cells[3] = sons[ faces[i][3] ] ;
        //        dual_it.direction = faces[i][4] ;
        s.push( dual_it ) ;
        print_dual_face( "face", dual_it.cells, M ) ;
      }
      
      //_____________________________________________________________________________
      // cube
      dual_it.type = ptr_dual_iterator_struct::DUAL_CUBE ;
      for( int i = 0 ; i < 8 ; ++i ) dual_it.cells[i] = sons[i] ;
      s.push( dual_it ) ;
      print_dual_cube( "face", dual_it.cells, M ) ;
      
      continue ;
    }
    
    
    //_____________________________________________________________________________
    //_____________________________________________________________________________
    // cube case
    else if( dual_it.type == ptr_dual_iterator_struct::DUAL_CUBE )
    {
      bool leaf = true ;
      
      //_____________________________________________________________________________
      // cube
      for( int i = 0 ; i < 8 ; ++i )
      {
        // copy father
        if( dual_it.cells[i].is_leaf() )
          continue ;
        
        // else
        leaf = false ;
        geom_cell sons[8] ;
        dual_it.cells[i].sons( sons ) ;
        dual_it.cells[i] = sons[ 7-i ] ;
      }
      
      // finally add dual cube
      if( leaf )
      {
        result &= walker( *this, dual_it.cells ) ;
        print_final_cube( dual_it.cells, M ) ;
      }
      else
      {
        s.push( dual_it ) ;
        print_dual_cube( "cube", dual_it.cells, M ) ;
      }
      
      continue ;
    }
    
  }
  
  return result ;
}
//_____________________________________________________________________________





#include "morton.h"
//_____________________________________________________________________________
// Isosurface walker
bool ptr_isosurface_walker( PtrOctree &fo, PtrOctree::geom_cell *cells )
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
    const PtrOctree::geom_cell &c = cells[j] ;
    cube[i] = c.cell()->field ;
    
    space[i] = (Point&)c ;
    indexes[i] = cube2key((Cube&)c) ;
  }
  
  mc.tesselate_cube( 0.0 ) ;
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Build the isosurface using dual marching cubes
bool PtrOctree::build_isosurface( data_access *ref )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;
  _mc.dat_access() = ref ;
  
    printf( "Build isosurface... " ) ;
  bool result = false ;
  
  //_________________________________________________________________________
  // Dual Marching Cubes
  
  _mc.init_all() ;
  
  result = dual_cubes_walk( ptr_isosurface_walker ) ;
  
  _mc.clean_temps() ;
  
  printf( "generated %d vertices and %d faces!\n", _mc.nverts(), _mc.ntrigs() ) ;
  
  return result ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Isosurface drawing walker
bool ptr_direct_isosurface_walker( PtrOctree &fo, PtrOctree::geom_cell *cells )
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
    const PtrOctree::geom_cell &c = cells[j] ;
    cube[i] = c.cell()->field ;
    
    space[i] = (Point&)c ;
  }
  
  mc.tesselate_cube( 0.0 ) ;
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Build the isosurface using dual marching cubes
bool PtrOctree::direct_draw_isosurface( data_access *ref )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;
  _mc_draw.dat_access() = ref ;
  
    bool result = false ;
  
  // Dual Marching Cubes
  ::glBegin( GL_TRIANGLES ) ;
  {
    result = dual_cubes_walk( ptr_direct_isosurface_walker ) ;
  }
  ::glEnd() ; // GL_TRIANGLES
  
  return result ;
}
//_____________________________________________________________________________






//_____________________________________________________________________________
// Draw walker
bool ptr_draw_walker( PtrOctree &fo, PtrOctree::geom_cell *cells )
//-----------------------------------------------------------------------------
{
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
bool PtrOctree::draw_dual()
//-----------------------------------------------------------------------------
{
  ::glBegin( GL_LINES ) ;
  bool result = dual_cubes_walk( ptr_draw_walker ) ;
  ::glEnd() ;  // GL_LINES
  
  return result ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Do nothing walker, just for dual generation timing
bool ptr_timing_walker( PtrOctree &fo, PtrOctree::geom_cell *cells )
//-----------------------------------------------------------------------------
{
  return true ;
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// Do nothing, just for dual generation timing
bool PtrOctree::dual_timing()
//-----------------------------------------------------------------------------
{
    return dual_cubes_walk( ptr_timing_walker ) ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Prints statistics about the octree
void PtrOctree::stats()
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
  printf( " total memory:\t%d\n", (int) ( s * sizeof(cell) ) ) ;
  
  printf( " leaves' levels:\t" ) ;
  for( Level l = 0 ; l <= MAX_LEVEL; ++l )
    printf( "\t%d", level_dist[l] ) ;
  printf( "\n" ) ;
  
  printf( " max size of the stack:\t%d\n", _dual_temp_memory ) ;
  printf( " total memory of stack:\t%d\n", (int) ( _dual_temp_memory * sizeof(ptr_dual_iterator_struct) ) ) ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
//_____________________________________________________________________________
// search operations


//_____________________________________________________________________________
// Find cells of the octree at a given position
bool PtrOctree::find_leaf( real x, real y, real z, geom_cell &n )
//-----------------------------------------------------------------------------
{
  
  n = geom_cell( &_root ) ;
  
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
    
    n.cell() = n.first_son() ;
    while( --i > -1 ) n.cell() = n.brother() ;
  }
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Find cells of the octree inside a given box of center x,y and half side r
bool PtrOctree::find_radius( real x, real y, real z, real r, List<geom_cell> &cells )
//-----------------------------------------------------------------------------
{
    
  cells.clear() ;
  
  geom_cell n( &_root ) ;
  
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
    n.sons( sons ) ;
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
bool PtrOctree::adjacent( const geom_cell &cell, List<geom_cell> &cells )
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
void PtrOctree::draw_plane ( real nx, real ny, real nz, real d )
//-----------------------------------------------------------------------------
{
  geom_cell n( &_root ) ;
  
  std::stack< geom_cell > s ;
  s.push( n ) ;
  
  //_____________________________________________________________________________
  // retrieve the index of the most positive cube corner
  bool fx = nx >= 0 ;
  bool fy = ny >= 0 ;
  bool fz = nz >= 0 ;
  
  //_____________________________________________________________________________
  
  ::glLineWidth( (GLfloat)1.0 ) ;
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
    n.sons( sons ) ;
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
void PtrOctree::draw_slice ( real nx, real ny, real nz, real d, float alpha )
//-----------------------------------------------------------------------------
{
  geom_cell n( &_root ) ;
  
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
        float def = (nneg>npos) ? (GLfloat)R_EPSILON : (GLfloat)-R_EPSILON ;
        
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
    n.sons( sons ) ;
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

