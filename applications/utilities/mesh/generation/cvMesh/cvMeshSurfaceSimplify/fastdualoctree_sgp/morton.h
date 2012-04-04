/**
 * \file    morton.h
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



#pragma once

#ifndef WIN32
#pragma interface
#endif // WIN32


#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stack>
#include "cube.h"

typedef uint64_t Key;
extern Key KEY_INV ;
const uint MAX_LEVEL = (sizeof(Key)*8-1)/3 ;

//_____________________________________________________________________________

/// octree depth level associated to a morton key ( = (bit length - 1) / 3 )
Level key_level ( Key k ) ;

/// containment text implemented directly onf the key
bool key_contains( Key k_big, Key k_small ) ;

/// convert an octree cube to a morton key using integer contraction
Key cube2key ( const Cube &c ) ; 

/// convert a morton key to a key using integer dilation
const Cube key2cube ( Key  k ) ; 

/// generate a list of 27 neighbor keys of same level, with their direction, indexed by 9*(dx+1) + 3*(dy+1) + (dz+1)
bool neighbor_keys( Key k, Key *neighbors ) ;

/// generate a list of the subneighbor keys from a neighbor direction, and return the number of subneighbors
char subneighbor_keys( Key k, char dir, std::stack<Key> &subneighbors ) ;

/// generate a list of 8 vertices keys from a cell key
Level vertices_keys( Key k, Key *verts ) ;

/// generate a list of 8 dual vertices keys from a vertex key and the deepest neighbor level
bool dual_vertices_keys( Key k, Level lv, Key *dual_verts ) ;

//_____________________________________________________________________________



