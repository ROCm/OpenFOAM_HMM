/**
 * \file    morton.cpp
 * \author  Thomas Lewiner   <tomlew@puc-rio.br>
 * \author  Rener Pereira de Castro <rener@mat.puc-rio.br> 
 * \author  Matmidia Lab, Math Dept, PUC-Rio
 * \date    10/01/2010
 *
 * \brief  Morton codes manipulation
 *
 * Morton codes manipulation
 */
//_____________________________________________________________________________


#ifndef WIN32
#pragma implementation
#endif // WIN32

#include "morton.h"

//_____________________________________________________________________________
// static
#ifdef UNIX
real  R_INV = strtof("NAN(char-sequence)", NULL); //nanf(0) ;
#else
real  R_INV = nan(0) ;
#endif

Level L_INV = (Level)-1 ;
Point P_INV( R_INV, R_INV, R_INV ) ;
Cube  C_INV( P_INV, L_INV ) ;
Key   KEY_INV = 0 ;
//_____________________________________________________________________________


//_____________________________________________________________________________
// Static Members Aux


// 1 << ((dz==0)+(dy==0)+(dx==0))
static const char nsubneighbors[27] = { 1, // Direction (-1,-1,-1)
                                        2, // Direction (-1,-1, 0)
                                        1, // Direction (-1,-1, 1)
                                        2, // Direction (-1, 0,-1)
                                        4, // Direction (-1, 0, 0)
                                        2, // Direction (-1, 0, 1)
                                        1, // Direction (-1, 1,-1)
                                        2, // Direction (-1, 1, 0)
                                        1, // Direction (-1, 1, 1)
                                        2, // Direction ( 0,-1,-1)
                                        4, // Direction ( 0,-1, 0)
                                        2, // Direction ( 0,-1, 1)
                                        4, // Direction ( 0, 0,-1)
                                        0, // Direction ( 0, 0, 0)
                                        4, // Direction ( 0, 0, 1)
                                        2, // Direction ( 0, 1,-1)
                                        4, // Direction ( 0, 1, 0)
                                        2, // Direction ( 0, 1, 1)
                                        1, // Direction ( 1,-1,-1)
                                        2, // Direction ( 1,-1, 0)
                                        1, // Direction ( 1,-1, 1)
                                        2, // Direction ( 1, 0,-1)
                                        4, // Direction ( 1, 0, 0)
                                        2, // Direction ( 1, 0, 1)
                                        1, // Direction ( 1, 1,-1)
                                        2, // Direction ( 1, 1, 0)
                                        1, // Direction ( 1, 1, 1)                               
                                      };

static const char l3b[27][4] = { {7,-1,-1,-1},{6, 7,-1,-1},{6,-1,-1,-1},{5, 7,-1,-1},{4, 5, 6, 7},{4, 6,-1,-1},{ 5,-1,-1,-1},
                                 {4, 5,-1,-1},{4,-1,-1,-1},{3, 7,-1,-1},{2, 3, 6, 7},{2, 6,-1,-1},{1, 3, 5, 7},{-1,-1,-1,-1},
                                 {0, 2, 4, 6},{1, 5,-1,-1},{0, 1, 4, 5},{0, 4,-1,-1},{3,-1,-1,-1},{2, 3,-1,-1},{ 2,-1,-1,-1},
                                 {1, 3,-1,-1},{0, 1, 2, 3},{0, 2,-1,-1},{1,-1,-1,-1},{0, 1,-1,-1},{0,-1,-1,-1},};
//_____________________________________________________________________________
// Contraction and Dilatation Masks
static const uint three_2 = 9 ;
static const uint three_1 = 3 ;
static const uint three_0 = 1 ;
static const uint shift_2a = 18, shift_2b = 36 ;
static const uint shift_1a =  6, shift_1b = 12 ;
static const uint shift_0a =  2, shift_0b =  4 ;
static const Key dilate_mask_2 = (Key)0x7FC0000FF80001FFLL ;  // bits 1 from 0-9, 27-36 and 54-63
static const Key dilate_mask_1 = (Key)0x01C0E070381C0E07LL ;  // bits 1 from 0-3, 9-12, 18-21,...
static const Key dilate_mask_0 = (Key)0x9249249249249249LL ;  // bits 1 at 0,3,6,9,12,...
static const Key dilate_tz     = (Key)0x4924924924924924LL ;  // bits 1 at 2,5,8,11,14,...    rigth to left
static const Key dilate_ty     = (Key)0x2492492492492492LL ;  // bits 1 at 1,4,7,10,13,...
static const Key dilate_tx     = (Key)0x9249249249249249LL ;  // bits 1 at 0,3,6,9,12,...  
static const Key dilate_t1     = (Key)0xB6DB6DB6DB6DB6DBLL ;  // bits 0 at 2,5,8,11,14,...    = ~tz
static const Key dilate_t2     = (Key)0xDB6DB6DB6DB6DB6DLL ;  // bits 0 at 1,4,7,10,13,...    = ~ty
static const Key dilate_t3     = (Key)0x6DB6DB6DB6DB6DB6LL ;  // bits 0 at 0,3,6,9,12,...     = ~tx

//_____________________________________________________________________________
// Step Directions Find Neighbours - 26                              (dz,dy,dx)
static const Key dir_neigh[27] = { (Key)0xFFFFFFFFFFFFFFFFLL , // Direction (-1,-1,-1)
                                   (Key)0x6DB6DB6DB6DB6DB6LL , // Direction (-1,-1, 0)
                                   (Key)0x6DB6DB6DB6DB6DB7LL , // Direction (-1,-1, 1)
                                   (Key)0xDB6DB6DB6DB6DB6DLL , // Direction (-1, 0,-1)
                                   (Key)0x4924924924924924LL , // Direction (-1, 0, 0)
                                   (Key)0x4924924924924925LL , // Direction (-1, 0, 1)
                                   (Key)0xDB6DB6DB6DB6DB6FLL , // Direction (-1, 1,-1)
                                   (Key)0x4924924924924926LL , // Direction (-1, 1, 0)
                                   (Key)0x4924924924924927LL , // Direction (-1, 1, 1)
                                   (Key)0xB6DB6DB6DB6DB6DBLL , // Direction ( 0,-1,-1)
                                   (Key)0x2492492492492492LL , // Direction ( 0,-1, 0)
                                   (Key)0x2492492492492493LL , // Direction ( 0,-1, 1)
                                   (Key)0x9249249249249249LL , // Direction ( 0, 0,-1)
                                   (Key)0x0000000000000000LL , // Direction ( 0, 0, 0)
                                   (Key)0x0000000000000001LL , // Direction ( 0, 0, 1)
                                   (Key)0x924924924924924BLL , // Direction ( 0, 1,-1)
                                   (Key)0x0000000000000002LL , // Direction ( 0, 1, 0)
                                   (Key)0x0000000000000003LL , // Direction ( 0, 1, 1)
                                   (Key)0xB6DB6DB6DB6DB6DFLL , // Direction ( 1,-1,-1)
                                   (Key)0x2492492492492496LL , // Direction ( 1,-1, 0)
                                   (Key)0x2492492492492497LL , // Direction ( 1,-1, 1)
                                   (Key)0x924924924924924DLL , // Direction ( 1, 0,-1)
                                   (Key)0x0000000000000004LL , // Direction ( 1, 0, 0)
                                   (Key)0x0000000000000005LL , // Direction ( 1, 0, 1)
                                   (Key)0x924924924924924FLL , // Direction ( 1, 1,-1)
                                   (Key)0x0000000000000006LL , // Direction ( 1, 1, 0)
                                   (Key)0x0000000000000007LL , // Direction ( 1, 1, 1)                               
                                 };

//_____________________________________________________________________________





//_____________________________________________________________________________
// octree depth level associated to a morton key ( = (bit length - 1) / 3 )
Level key_level ( Key k )
//_____________________________________________________________________________
{
  static const real lg2_3 = 1.0 / (log(2)*3) ;
  if( k == KEY_INV ) return L_INV ;
  return (Level) floor( log(k) * lg2_3 ) ;
/*
  Level count = 0 ;
  do
  {
    k >>= 3 ;
    ++count;
  }
  while( k > 0 ) ;
	  
  return count-1 ;
*/
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// containment text implemented directly onf the key
bool key_contains( Key k_big, Key k_small )
//_____________________________________________________________________________
{
  while( (k_small!=KEY_INV) && ((k_small>>3) != k_big) )
    k_small >>= 3;
  return k_small != KEY_INV ;  
  // k_small & 7
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// convert an octree cube to a morton key using integer contraction
Key cube2key ( const Cube &c )
//_____________________________________________________________________________
{
  uint level = c.lv();
	Key  m     = 1 << level ;
	real sz    = (real) m ;
	
	Key bz = (Key) floor( c.z() * sz ) ;
	Key by = (Key) floor( c.y() * sz ) ;
	Key bx = (Key) floor( c.x() * sz ) ;
	/*
	 if( bz >= m ) bz = m-1 ;
	 if( by >= m ) by = m-1 ;
	 if( bx >= m ) bx = m-1 ;
	 */
	
	if( level > three_2 )
	{
		bz = ( bz | (bz << shift_2a) | (bz << shift_2b) ) & dilate_mask_2 ;
		by = ( by | (by << shift_2a) | (by << shift_2b) ) & dilate_mask_2 ;
		bx = ( bx | (bx << shift_2a) | (bx << shift_2b) ) & dilate_mask_2 ;
	}
	if( level > three_1 )
	{
		bz = ( bz | (bz << shift_1a) | (bz << shift_1b) ) & dilate_mask_1 ;
		by = ( by | (by << shift_1a) | (by << shift_1b) ) & dilate_mask_1 ;
		bx = ( bx | (bx << shift_1a) | (bx << shift_1b) ) & dilate_mask_1 ;
	}
	if( level > three_0 )
	{
		bz = ( bz | (bz << shift_0a) | (bz << shift_0b) ) & dilate_mask_0 ;
		by = ( by | (by << shift_0a) | (by << shift_0b) ) & dilate_mask_0 ;
		bx = ( bx | (bx << shift_0a) | (bx << shift_0b) ) & dilate_mask_0 ;
	}
  return (Key(1) << (3*level)) | (bz << 2) | (by << 1) | bx ;
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// convert a morton key to a key using integer dilation
const Cube key2cube ( Key  key )
//_____________________________________________________________________________
{
  if( key == KEY_INV )
    return C_INV ;
  
  Level level = key_level( key ) ;
  Key bz( (key >> 2) & dilate_mask_0 ) ;
  Key by( (key >> 1) & dilate_mask_0 ) ;
  Key bx( (   key  ) & dilate_mask_0 ) ;
  
  if( level > three_0 )
  {
    bz = ( bz | (bz >> shift_0a) | (bz >> shift_0b) ) & dilate_mask_1 ;
    by = ( by | (by >> shift_0a) | (by >> shift_0b) ) & dilate_mask_1 ;
    bx = ( bx | (bx >> shift_0a) | (bx >> shift_0b) ) & dilate_mask_1 ;
    
    if( level > three_1 )
    {
      bz = ( bz | (bz >> shift_1a) | (bz >> shift_1b) ) & dilate_mask_2 ;
      by = ( by | (by >> shift_1a) | (by >> shift_1b) ) & dilate_mask_2 ;
      bx = ( bx | (bx >> shift_1a) | (bx >> shift_1b) ) & dilate_mask_2 ;
      
      if( level > three_2 )
      {
        bz = bz | (bz >> shift_2a) | (bz >> shift_2b) ;
        by = by | (by >> shift_2a) | (by >> shift_2b) ;
        bx = bx | (bx >> shift_2a) | (bx >> shift_2b) ;
      }
    }
  }
  
  Key length_mask = (1 << level) -1 ;
  bz &= length_mask ;
  by &= length_mask ;
  bx &= length_mask ;
  
  Cube c ;
  c.lv() = level ;
  real sz  = c.sz() ;
  real sz2 = sz*2 ;
  c.z() = (real)(bz) * sz2 + sz ;
  c.y() = (real)(by) * sz2 + sz ;
  c.x() = (real)(bx) * sz2 + sz ;

  return c ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// addition in dilated integer
Key addition_locate_code( const Key& N , const Key& DN )
//-----------------------------------------------------------------------------
{
  return (
          ( ((N | dilate_t1) + (DN & dilate_tz)) & dilate_tz ) |
          ( ((N | dilate_t2) + (DN & dilate_ty)) & dilate_ty ) |
          ( ((N | dilate_t3) + (DN & dilate_tx)) & dilate_tx ) 
          );
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// addition in dilated integer
Key substraction_locate_code( const Key& N , const Key& DN )
//-----------------------------------------------------------------------------
{
  return (
          ( ((N & dilate_tz) - (DN & dilate_tz)) & dilate_tz ) |
          ( ((N & dilate_ty) - (DN & dilate_ty)) & dilate_ty ) |
          ( ((N & dilate_tx) - (DN & dilate_tx)) & dilate_tx ) 
         );
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// generate a list of 27 neighbor keys of same level, with their direction, indexed by 9*(dz+1) + 3*(dy+1) + (dx+1)
bool neighbor_keys( Key k, Key *neighbors )
//_____________________________________________________________________________
{
  Level nbits_neigh = 3*key_level(k) ;
  for( uint dir = 0 ; dir < 27 ; ++dir ) 
  {
    Key &nk = neighbors[dir] ;
    if( dir == 13 ) { nk = KEY_INV ; continue ; }
    nk = addition_locate_code( k , dir_neigh[dir] ) ; 
    if( (nk >> nbits_neigh) != 1 )
      nk = KEY_INV ;
  }
  return true ;
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// generate a list of the subneighbor keys from a neighbor direction
// returns the number of subneighbors
char subneighbor_keys( Key k, char dir, std::stack<Key> &subneighbors )
//_____________________________________________________________________________
{
  k <<= 3 ;
  char nsub = nsubneighbors[dir] ;
  for( char i = 0 ; i < nsub ; ++i )
  {
    subneighbors.push( k | l3b[dir][i] ) ;
  }
  return nsub ; 
}
//_____________________________________________________________________________



//_____________________________________________________________________________
void print_binary( const Key k )
//_____________________________________________________________________________
{
  printf( "0x" ) ;
  for( char i = 63 ; i >= 0 ; --i )
    printf( "%c", (k>>i)&1 ? '1' : '0' ) ;
}
//_____________________________________________________________________________


//_____________________________________________________________________________
// generate a list of 8 vertices keys
Level vertices_keys( Key k, Key *verts )
//_____________________________________________________________________________
{
  Level lv = key_level(k) ;
  Key   lv_k = 1 << (3*lv) ;
  Key   lv_K = lv_k << 1 ;

  for( int i = 0 ; i < 8 ; ++i )
  {
    Key vk = addition_locate_code( k,i ) ;

    if( vk >= lv_K || !((vk-lv_k)&dilate_tx) || !((vk-lv_k)&dilate_ty) || !((vk-lv_k)&dilate_tz) )
    {
      verts[i] = KEY_INV ;
      continue ;
    }
    
    verts[i] = vk << 3*(MAX_LEVEL - lv) ;
  }
  return lv ;

/*  
  Level lv = key_level(k) ;
  Key   lv_k = 1 << (3*lv) ;
  Key   lv_K = lv_k << 1 ;

  Cube c = key2cube(k) ;
  const real ix = c.x () ;
  const real iy = c.y () ;
  const real iz = c.z () ;
  const real is = c.sz() ;
  
  for( int i = 0 ; i < 8 ; ++i )
  {
    Key vk = addition_locate_code( k,i ) ;
    Cube v( (i&1) ? ix+is : ix-is,
            (i&2) ? iy+is : iy-is,
            (i&4) ? iz+is : iz-is,
            MAX_LEVEL ) ;
    bool invalid = (
                    v.x() <= 0.0 || v.x() >= 1 ||
                    v.y() <= 0.0 || v.y() >= 1 ||
                    v.z() <= 0.0 || v.z() >= 1
                   ) ;
    bool k_invalid = vk >= lv_K || !((vk-lv_k)&dilate_tx) || !((vk-lv_k)&dilate_ty) || !((vk-lv_k)&dilate_tz) ;
    if( invalid != k_invalid )
      printf( "problem\n" ) ;
    if( invalid )
    {
      verts[i] = KEY_INV ;
      continue ;
    }
    
    vk <<= 3*(MAX_LEVEL - lv) ;
    verts[i] = cube2key(v) ;
    if( vk != verts[i] )
      printf( "problem\n" ) ;
  }
  return c.lv() ;
*/
}
//_____________________________________________________________________________



//_____________________________________________________________________________
// generate a list of 8 dual vertices keys from a vertex key and the deepest neighbor level
bool dual_vertices_keys( Key k, Level lv, Key *dual_verts )
//_____________________________________________________________________________
{
  Key dk = k >> (3*(MAX_LEVEL - lv)) ;
  
  for( int i = 0 ; i < 8 ; ++i )
  {
    dual_verts[i] = substraction_locate_code( dk,i ) ;
  }
  
  return true ;

/*  
  Key dk = k >> (3*(MAX_LEVEL - lv)) ;
  
  Cube c = key2cube( k ) ;
  c.lv() = lv ;
  
  const real ix = c.x () ;
  const real iy = c.y () ;
  const real iz = c.z () ;
  const real is = c.sz() ;
  
//  print_binary( k ) ;
//  printf( " <- \n" ) ;
  for( int i = 0 ; i < 8 ; ++i )
  {
    Cube dv( (i&1) ? ix+is : ix-is,
             (i&2) ? iy+is : iy-is,
             (i&4) ? iz+is : iz-is,
             lv ) ;
    
    dual_verts[i] = cube2key(dv) ;
    if( substraction_locate_code( dk,7-i) != dual_verts[i] )
      printf( "problem\n" ) ;
//    print_binary( dual_verts[i] ) ;
//    printf( "\n" ) ;
//    print_binary( substraction_locate_code( dk,7-i) ) ;
//    printf( "\n" ) ;
  }
//  printf( "\n" ) ;
  
  return true ;
*/
}
//_____________________________________________________________________________


// testar que dois cubos se intersectam???


