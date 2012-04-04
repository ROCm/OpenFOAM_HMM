/**
 * \file    perf_main.cpp
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 *  Octree performances test
 */
//_____________________________________________________________________________


/*
with pointers
 level 10,12: hash 27, non hash 31
 level 8: hash 27, non hash 37
 
without pointers
 level 10: hash 21, non hash 31
 level 8: hash 21, non hash 37/40
*/


#define OCTREE_SWITCH 2

#ifdef OCTREE_SWITCH
# if   OCTREE_SWITCH==0
#  define OCTREE_PTR   1
# elif OCTREE_SWITCH==1
#  define OCTREE_HASH  1
# elif OCTREE_SWITCH==2
#  define OCTREE_OPT   1
# elif OCTREE_SWITCH==3
#  define OCTREE_LEAF  1
# elif OCTREE_SWITCH==4
#  define OCTREE_MEM   1
# endif // OCTREE_SWITCH
#endif // OCTREE_SWITCH


// switch values between pointer and hash octree
#if !OCTREE_PTR && !OCTREE_HASH && !OCTREE_OPT && !OCTREE_LEAF && !OCTREE_MEM
// # define OCTREE_PTR  1
// # define OCTREE_HASH 1
# define OCTREE_OPT  1
// # define OCTREE_LEAF 1
// # define OCTREE_MEM  1
#endif // !OCTREE_PTR && !OCTREE_HASH && !OCTREE_OPT && !OCTREE_LEAF && !OCTREE_MEM


#ifdef OCTREE_PTR
# undef OCTREE_HASH
# undef OCTREE_OPT  
# undef OCTREE_LEAF 
# undef OCTREE_MEM  
# define OCTREE_STRING "0 pointer"
# include "ptr_octree.h"
  /// octree type
  typedef PtrOctree Octree ;
#endif //OCTREE_PTR


#ifdef OCTREE_HASH
# undef OCTREE_PTR
# undef OCTREE_OPT  
# undef OCTREE_LEAF 
# undef OCTREE_MEM  
# define OCTREE_STRING "1 hash"
# include "hash_octree.h"
  /// octree type
  typedef HashOctree Octree ;
#endif //OCTREE_HASH


#ifdef OCTREE_OPT
# undef OCTREE_PTR
# undef OCTREE_HASH
# undef OCTREE_LEAF 
# undef OCTREE_MEM  
# define OCTREE_STRING "2 opt"
# include "opt_octree.h"
  /// octree type
  typedef OptOctree Octree ;
#endif //OCTREE_OPT


#ifdef OCTREE_LEAF
# undef OCTREE_PTR
# undef OCTREE_HASH
# undef OCTREE_OPT  
# undef OCTREE_MEM  
# define OCTREE_STRING "3 leaf"
# include "leaf_octree.h"
/// octree type
typedef LeafOctree Octree ;
#endif //OCTREE_LEAF


#ifdef OCTREE_MEM
# undef OCTREE_PTR
# undef OCTREE_HASH
# undef OCTREE_OPT  
# undef OCTREE_LEAF 
# define OCTREE_STRING "4 mem"
# include "mem_octree.h"
/// octree type
typedef MemOctree Octree ;
#endif //OCTREE_MEM


#ifndef NPOINTS
# define NPOINTS 100
#endif // NPOINTS

#ifndef RADIUS
# define RADIUS (0.003f)
#endif // RADIUS

#ifndef NDUALS
# define NDUALS 10
#endif // NDUALS

#ifndef MLEVEL
# define MLEVEL 5
#endif // MLEVEL

#ifndef IMPLFUN // set to 0 for random
# define IMPLFUN 1
#endif // IMPLFUN

#ifndef ISOVAL
# define ISOVAL 0.5
#endif // ISOVAL


#include "hash.h"
#include "implfuns.h"
#include <stdio.h>

//_____________________________________________________________________________
//
int main( int argc , char* argv[] )
//-----------------------------------------------------------------------------
{
#if OCTREE_LEAF && !USE_HASH_PTR
  return 0 ;
#endif // OCTREE_LEAF && !USE_HASH_PTR
  
  srand(77) ;
//  Chrono::ResetClocks() ;

  printf( "\n\n\n-----------------------------------------------------------------------------\n"
         OCTREE_STRING
         "\n\thash\t%s\tpointers\n"
         "\thash bits =\t%d\n"
         "\tnon hash bits=\t%d\n"
         "\tmax_level =\t%d\n"
         "\tdensity =\t%f\n"
         "\tnpoints =\t%d\n"
         "\tradius =\t%f\n"
         "\timplicit function =\t%s\n"
         "\tisovalue =\t%f\n"
         "\tnumber of duals=\t%d\n\n",
         (USE_HASH_PTR ? "with" : "without"), HASH_BITS, NON_HASH_BITS, MLEVEL, RAND_THRES, NPOINTS, RADIUS, fun_list[IMPLFUN], ISOVAL, NDUALS ) ;
  fflush(stdout);
  
  /// main data structure
  Octree octree ;
#if IMPLFUN==0
  data_rand ref( MLEVEL ) ;
  octree.refine( &ref ) ;
#else  // IMPLFUN
  data_func ref( MLEVEL, ISOVAL, fun_def[IMPLFUN] ) ;
  octree.refine( &ref ) ;
  octree.set_impl( &ref ) ;
#endif // IMPLFUN
  fflush(stdout);

  
/*  
  //-------------------------------------------------
  // search for random points
  Chrono::ResetClocks() ;
  printf( "\n\nRandom Points\n" ) ;
  for( int i = 0 ; i < NPOINTS ; ++i )
  {
    real x = (real)rand()/(real)RAND_MAX;
    real y = (real)rand()/(real)RAND_MAX;
    real z = (real)rand()/(real)RAND_MAX;
    Octree::geom_cell gc ;
    List<Octree::geom_cell> gc_list ;
    octree.find_leaf( x,y,z, gc ) ; 
    octree.find_radius( x,y,z, RADIUS, gc_list ) ; gc_list.clear() ;
    octree.adjacent( gc, gc_list ) ; gc_list.clear() ;
  }
  Chrono::PrintClocks() ;
  fflush(stdout);
  //-------------------------------------------------


  //-------------------------------------------------
  // search for all leaves once
  printf( "\nLeaves\n" ) ;
  Chrono::ResetClocks() ;
  for( Octree::leaf_iterator it = octree.leaves_begin() ; it() ; ++it )
  {
    Octree::geom_cell g_it = it.top() ;
    Octree::geom_cell gc ;
    List<Octree::geom_cell> gc_list ;
    octree.find_leaf( g_it.cx(), g_it.cy(), g_it.cz(), gc ) ; 
    octree.find_radius( g_it.cx(), g_it.cy(), g_it.cz(), RADIUS, gc_list ) ; gc_list.clear() ;
    octree.adjacent( gc, gc_list ) ; gc_list.clear() ;
  }
  Chrono::PrintClocks() ;
  fflush(stdout);
  //-------------------------------------------------
  
  
  //-------------------------------------------------
  // search for all points at the maximal level
  printf( "\n\nGrid\n" ) ;
  Level lv = octree.max_level() - 4 ;
  real step = 1.0 / (1 << lv) ;
  Chrono::ResetClocks() ;
  for( real x = 0.0 ; x <= 1.0 ; x += step )
  {
    for( real y = 0.0 ; y <= 1.0 ; y += step )
    {
      for( real z = 0.0 ; z <= 1.0 ; z += step )
      {
        Octree::geom_cell gc ;
        List<Octree::geom_cell> gc_list ;
        octree.find_leaf( x,y,z, gc ) ; 
        octree.find_radius( x,y,z, RADIUS, gc_list ) ; gc_list.clear() ;
        octree.adjacent( gc, gc_list ) ; gc_list.clear() ;
      }
    }
  }
  Chrono::PrintClocks() ;
  fflush(stdout);
  //-------------------------------------------------
*/

  //-------------------------------------------------
  // build many duals
  printf( "\n\nDuals\n" ) ;
//  Chrono::ResetClocks() ;
  for( int i = 0 ; i < NDUALS ; ++i )
  {
#if IMPLFUN==0
    octree.dual_timing() ;
#else  // IMPLFUN
    octree.mc().clean_all() ;
    octree.build_isosurface(&ref) ;
#endif // IMPLFUN
  }
//  Chrono::PrintClocks() ;
  octree.stats() ;
  fflush(stdout);
  //-------------------------------------------------

  
  return 0; 
}
//_____________________________________________________________________________

