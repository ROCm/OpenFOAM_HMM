/**
 * @file    hash.h
 * @author  Thomas Lewiner 
 * @version 
 * @date    
 *
 * @brief
 *
 */
//________________________________________________


#ifndef _HASH_H_
#define _HASH_H_

// may be overridden by the preprocessor
#ifndef USE_HASH_PTR
# define USE_HASH_PTR 1
#endif // USE_HASH_PTR

// may be overridden by the preprocessor
#ifndef HASH_BITS
# define HASH_BITS 22
#endif // HASH_BITS

#ifndef NON_HASH_BITS
# define NON_HASH_BITS 22
#endif // HASH_BITS




#if USE_HASH_PTR

// #define HASH_HAS_ERASE 1
# include "hash_ptr.h"

#else  // USE_HASH_PTR

# undef HASH_HAS_ERASE
# include "hash_noptr.h"

#endif // USE_HASH_PTR


#endif // _HASH_H_
