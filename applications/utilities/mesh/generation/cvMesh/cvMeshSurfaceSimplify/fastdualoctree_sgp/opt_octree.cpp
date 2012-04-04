/**
 * \file    opt_octree.cpp
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 *  Octree structure with hashtable and optimized operations
 */
//_____________________________________________________________________________


#ifndef WIN32
#pragma implementation
#endif // WIN32


#include <stack>
#include "opt_octree.h"

/// invalid key
template<> Hash<Level>::KeyData Hash<Level>::KD_INV = { KEY_INV, L_INV } ;


//_____________________________________________________________________________
// Create the root
void OptOctree::init()
//-----------------------------------------------------------------------------
{
  _max_field = R_INV ;
  _opt_level = L_INV ;
  memset( _level_dist, 0, (MAX_LEVEL+1)*sizeof(uint) ) ;

  HashField::KeyData kd ;
  kd.key  =  1  ;
  kd.data = 0.0 ;
  _hash.insert( kd ) ;

  _level_dist[0] = 1 ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
/// Compute the optimal level
void OptOctree::compute_opt_level()
//-----------------------------------------------------------------------------
{
  Level left     = 0 ;
  Level right    = MAX_LEVEL ;
  uint  leftsum  =  0 ;
  uint  rightsum =  0 ;
  while( left < right )
  {
    uint leftsum_old = leftsum ;
    if( leftsum <= rightsum )
    {
      leftsum += _level_dist[ left ] ;
      ++left ;
    }
    if( leftsum_old >= rightsum )
    {
      rightsum += _level_dist[ right ] ;
      --right ;
    }
  }  
  _opt_level = right ;
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// Check that each point scan or scribble is in the right node
bool OptOctree::check ()
//-----------------------------------------------------------------------------
{
  bool result = true ;

  uint level_dist[MAX_LEVEL+1] ;
  memset( level_dist, 0, (MAX_LEVEL+1)*sizeof(uint) ) ;
  for( cell_iterator it = cells_begin() ; it() ; ++it )
  {
    if( is_leaf( it.key() ) )
      level_dist[ key_level(it.key()) ] ++ ;
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
      printf("(opt) Invalid leaf %d: value %f\n", (int)k, _hash[k].data ) ;
      result = false ;
    }

    if( has_son && !all_sons )
    {
      printf("(opt) Invalid subdivision %d: value %f\n", (int)k, _hash[k].data ) ;
      result = false ;
    }
  }

  for( Level lv = 0 ; lv <= MAX_LEVEL ; ++lv )
  {
    if( _level_dist[lv] != level_dist[lv] )
    {
      printf("(opt) Invalid level statistic %d: %d != %d\n", (int)lv, (int)_level_dist[lv], (int)level_dist[lv] ) ;
      result = false ;
    }
  }

  return result ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Deletes the content of the octree
void OptOctree::clear_octree()
//-----------------------------------------------------------------------------
{
  _hash .reset() ;
  _verts.reset() ;
  _max_field = R_INV ;
  _opt_level = L_INV ;
  memset( _level_dist, 0, (MAX_LEVEL+1)*sizeof(uint) ) ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// set the values of each leaf from the implicit function
bool OptOctree::set_impl( data_access *ref )
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
bool OptOctree::refine( data_access *ref /*= NULL*/ )
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
    Level lv = it.lv() ;
    _level_dist[ lv ] -= 1 ;
    _level_dist[lv+1] += 8 ;
  }
  
  compute_opt_level() ;

  return refined ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Refine the octree according to the data access
bool OptOctree::adapt( data_access *ref /*= NULL*/ )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;

#ifdef HASH_HAS_ERASE
    
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
//    printf( "opt key %d\n", (int) it.key() ) ;

    // look for the leaves
    if( !it.is_leaf() )
    {
      if( !(*ref).need_refine( it ) )
      {
        std::stack< std::pair<Key,Level> > del_s ; // avoid level recomputation
        del_s.push( std::make_pair( it.key(), it.lv() ) ) ;
        while( !del_s.empty() )
        {
          Key n    = del_s.top().first << 3 ;
          Level lv = del_s.top().second + 1 ;
          del_s.pop() ;
          
          for( int i = 0 ; i < 8 ; ++i )
          {
            HashField::KeyData kd = _hash.erase( n|i ) ;
            if( !is_inv( kd.data ) ) // is_leaf(n)
            {
              _level_dist[ lv ] -- ;
              continue ;
            }
            del_s.push( std::make_pair( n|i, lv ) ) ;
          }
        }

        refined = true ;
        _level_dist[it.lv()] += 1 ;
        _hash[it.key()].data = 0.0 ; // leaf now!
        continue ;
      }
      
      it.sons( sons, _hash ) ;
      for( int i = 0 ; i < 8 ; ++i )
        s.push( sons[i] ) ;
      
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
    
    Level lv = it.lv() ;
    _level_dist[ lv ] -= 1 ;
    _level_dist[lv+1] += 8 ;
  }
  
  compute_opt_level() ;
  return refined ;
  
#else  // HASH_HAS_ERASE
  // no remotion in the hash table!
  return refine( ref ) ;
#endif // HASH_HAS_ERASE
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Draw the octree with wireframe
bool OptOctree::draw_wire()
//-----------------------------------------------------------------------------
{
  for( leaf_iterator it = leaves_begin() ; it() ; ++it )
    it.draw_wire() ;
  
  return true ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Draw the octree with dots
bool OptOctree::draw_centers()
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
// Compute primal vertices
bool OptOctree::compute_primal_verts()
//-----------------------------------------------------------------------------
{
    
  // add all the inner vertices
  for( leaf_iterator l_it = leaves_begin() ; l_it() ; ++l_it )
  {
    Key verts_k[8] ;
    Level lv = vertices_keys( l_it.key(), verts_k ) ;
    
    HashVerts::KeyData kd ;
    kd.data = lv ;
    Key *v_ptr = verts_k ;
    for( int i = 0 ; i < 8 ; ++i, ++v_ptr )
    {
      kd.key = *v_ptr ;
      if( kd.key == KEY_INV ) continue ;
      
      HashVerts::KeyData &hk = _verts[ kd.key ] ;
      if( hk.key == KEY_INV )
        _verts.insert( kd ) ;
      else if( hk.data < kd.data )
        hk.data = kd.data ;
    }
  }
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Octree dual caller
bool OptOctree::dual_cubes_walk( opt_dual_walker &walker )
//-----------------------------------------------------------------------------
{
  bool result = compute_primal_verts() ;
  
  for( HashVerts::const_iterator v_it = _verts.cbegin() ; v_it() ; ++v_it )
  {
    const Level lv = *v_it ;
    Key dual_verts[8] ;
    dual_vertices_keys( v_it.key(), lv, dual_verts ) ;

    Key *dv_ptr = dual_verts ;
    for( int i = 0 ; i < 8 ; ++i, ++dv_ptr )
    {
      // find_at, only the level too high case
      Key &ki = *dv_ptr ;
      HashField::KeyData kd ;
      while( ki != 0 && !node_exists(ki,kd) )
        ki >>= 3;
    }
    result &= walker( *this, dual_verts ) ;
  }

  return result ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Isosurface walker
bool opt_isosurface_walker( OptOctree &fo, Key *keys )
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
    OptOctree::geom_cell c = fo.geom_key(k) ;
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
bool OptOctree::build_isosurface( data_access *ref )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;
  _mc.dat_access() = ref ;
  _verts.reset() ;
  
    printf( "Build isosurface... " ) ;
  bool result = false ;
  
  //_________________________________________________________________________
  // Dual Marching Cubes
  
  _mc.init_all() ; 
  
  result = dual_cubes_walk( opt_isosurface_walker ) ;
  
  _mc.clean_temps() ;
  
  printf( "generated %d vertices and %d faces!\n", _mc.nverts(), _mc.ntrigs() ) ;
  
  return result ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Isosurface drawing walker
bool opt_direct_isosurface_walker( OptOctree &fo, Key *keys )
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
    const OptOctree::geom_cell c = fo.geom_key( keys[j] ) ;
    cube[i] = *c ;
    
    space[i] = (Point)c ;
  }
  
  mc.tesselate_cube( 0.0 ) ;
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Build the isosurface using dual marching cubes
bool OptOctree::direct_draw_isosurface( data_access *ref )
//-----------------------------------------------------------------------------
{
  if( !ref ) return false ;
  _mc_draw.dat_access() = ref ;
  _verts.reset() ;
  
    bool result = false ;
  
  // Dual Marching Cubes
  ::glBegin( GL_TRIANGLES ) ;
  {
    result = dual_cubes_walk( opt_direct_isosurface_walker ) ;
  }
  ::glEnd() ; // GL_TRIANGLES
  
  return result ;
}
//_____________________________________________________________________________






//_____________________________________________________________________________
// Draw walker
bool opt_draw_walker( OptOctree &fo, Key *keys )
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
bool OptOctree::draw_dual()
//-----------------------------------------------------------------------------
{
  _verts.reset() ;

  ::glBegin( GL_LINES ) ;
  bool result = dual_cubes_walk( opt_draw_walker ) ;
  ::glEnd() ;  // GL_LINES
  
  return result ;
}
//_____________________________________________________________________________




//_____________________________________________________________________________
// Do nothing walker, just for dual generation timing
bool opt_timing_walker( OptOctree &fo, Key *keys )
//-----------------------------------------------------------------------------
{
  return true ;
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// Do nothing, just for dual generation timing
bool OptOctree::dual_timing()
//-----------------------------------------------------------------------------
{
  _verts.reset() ;

    return dual_cubes_walk( opt_timing_walker ) ;
}
//_____________________________________________________________________________





//_____________________________________________________________________________
// Prints statistics about the octree
void OptOctree::stats() const
//-----------------------------------------------------------------------------
{
  int s = _hash.size() ;
  printf( " number of nodes:\t%d\n", s ) ;
#if USE_HASH_PTR
  printf( " total memory:\t%d\n", (int) ( s * sizeof(HashField::KeyBucket) ) ) ;
  _hash.stats() ;
#else  // USE_HASH_PTR 
  printf( " total memory:\t%d\n", (int) ( s * sizeof(HashField::KeyData  ) ) ) ;
#endif // USE_HASH_PTR
  
  printf( " leaves' levels:\t" ) ;
  for( Level l = 0 ; l <= MAX_LEVEL; ++l )
    printf( "\t%d", _level_dist[l] ) ;
  printf( "\n" ) ;
  
  s = _verts.size() ;
  printf( " number of dual nodes:\t%d\n", s ) ;
#if USE_HASH_PTR
  printf( " total memory of dual:\t%d\n", (int) ( s * sizeof(HashVerts::KeyBucket) ) ) ;
  _verts.stats() ;
#else  // USE_HASH_PTR 
  printf( " total memory of dual:\t%d\n", (int) ( s * sizeof(HashVerts::KeyData  ) ) ) ;
#endif // USE_HASH_PTR
}
//_____________________________________________________________________________





//_____________________________________________________________________________
//_____________________________________________________________________________
// search operations


//_____________________________________________________________________________
// Find cells of the octree at a given position
bool OptOctree::find_leaf( Key k, Level o_lv, geom_cell &n ) const
//-----------------------------------------------------------------------------
{
    
  Level   shift = 3*(MAX_LEVEL-o_lv) ;
  Key     ki    = k >> shift ;
  HashField::KeyData kd ;
  
  // level too high
  if( !node_exists(ki,kd) )
  {
    do
    {
      ki >>= 3;
    }
    while( ki !=0 && !node_exists(ki,kd) ) ;
  }
  // level correct or too low
  else
  {
    HashField::KeyData kd_temp ;
    do
    {
      kd_temp = kd ;
      shift -= 3 ;
      ki = k >> shift ;
    }
    while( ki !=0 && node_exists(ki,kd) ) ;
    kd = kd_temp ;
  }
  n = geom_cell( kd ) ;

  return kd.key != KEY_INV ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Find cells of the octree inside a given box of center x,y and half side r
bool OptOctree::find_radius( real x, real y, real z, real r, List<geom_cell> &cells ) const
//-----------------------------------------------------------------------------
{
    
  cells.clear() ;

  Level lv_b = (Level)floor( log(r)/log(0.5) ) ;
  Key   k_b  = cube2key( Cube(x,y,z,lv_b) ) ;

  Key neighbors[27] ;
  neighbor_keys( k_b, neighbors ) ;

  std::stack<Key> s ;
  s.push( k_b ) ;
  
  Key *k_ptr = neighbors ;
  for( int dir = 0 ; dir < 27 ; ++dir, ++k_ptr )
    s.push( *k_ptr ) ;
  
  while( !s.empty() )
  {
    Key ki = s.top() ;
    s.pop() ;
    if( ki == KEY_INV ) continue ;
    
    HashField::KeyData kd ;
    if( node_exists(ki,kd) )
    {
      if( is_leaf( kd ) )
      {
        geom_cell gc(kd) ;
        real sz = gc.sz() + r ;

        // test directly with the key?
        if( fabs( x-gc.x() ) <= sz && fabs( y-gc.y() ) <= sz && fabs( z-gc.z() ) <= sz )
          cells.insert( gc ) ;
      }
      else 
      {
        Key ski = ki << 3 ;
        for( int i = 0 ; i < 8 ; ++i )
          // test directly with the key before insertion?
          s.push( ski | i ) ;
      }
    }        
    else
    {
      HashField::KeyData kd ;
      do
      {
        ki >>= 3 ;
      } while( ki !=0 && !node_exists(ki,kd) ) ;
      
      if( ki > 0 )
      {
        geom_cell gc(kd) ;
        real sz = gc.sz() + r ;

        // test directly with the key?
        if( fabs( x-gc.x() ) <= sz && fabs( y-gc.y() ) <= sz && fabs( z-gc.z() ) <= sz )
          cells.insert_unique( gc ) ;
      }
    }
  }
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// Find neighbor cells of the octree to a given cell
bool OptOctree::adjacent( Key k, List<geom_cell> &cells ) const
//-----------------------------------------------------------------------------
{
    
  cells.clear() ;

  Key neighbors[27] ;
  neighbor_keys( k, neighbors ) ;
  
  Key *k_ptr = neighbors ;
  for( int dir = 0 ; dir < 27 ; ++dir, ++k_ptr )
  {
    Key ki = *k_ptr ;
    if( ki == KEY_INV ) continue ;
    
    HashField::KeyData kd ;
    if( node_exists(ki,kd) )
    {
      if( is_leaf( kd ) ) cells.insert( geom_cell(kd) );        
      else 
      {
        std::stack<Key> subneighbors ;
        subneighbor_keys(ki, dir, subneighbors ) ;
        while( !subneighbors.empty() )
        {
          Key ski = subneighbors.top() ;
          subneighbors.pop() ;
          
          real sdi ;
          if( is_leaf( ski, sdi ) )
            cells.insert( geom_cell( ski, sdi ) ) ;
          else
            subneighbor_keys(ski, dir, subneighbors ) ;
        }
      }
    }        
    else
    {
      HashField::KeyData kd ;
      do
      {
        ki >>= 3 ;
      } while( ki !=0 && !node_exists(ki,kd) ) ;

      if( ki > 0 )
        cells.insert_unique( geom_cell(kd) ) ;
    }
  }
  
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
//_____________________________________________________________________________
// I/O

//_____________________________________________________________________________
// Draws the intersecting cells of a plane with the octree with color
void OptOctree::draw_plane ( real nx, real ny, real nz, real d )
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
void OptOctree::draw_slice ( real nx, real ny, real nz, real d, float alpha )
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

