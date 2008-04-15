#ifndef __VECTOR_H
#define __VECTOR_H

/*@@
 *  Program: Star File Format Library  - $RCSfile: vector.h,v $
 *  Author:  Geoff Prewett
 *  Date:    August 12, 2003
 *
 *
 *  Star File Format Library - Copyright (C) 2003 by adapco, Ltd.
 *
 *  This program is the property of adapco, Ltd. and contains
 *  confidential and proprietary information.  The unauthorized use,
 *  distribution, or duplication of this program is prohibited.
 *  All rights reserved.
 *
 *  $Id: vector.h,v 1.3 2004/10/25 19:11:55 wayne Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

struct _Vector {
    int size;
    int alloc;
    int typeSize;
    int clear;		/* TRUE if need to clear memory when expanding */
    void *buffer;
    };

typedef struct _Vector* Vector;

/** Creates a growable array (or NULL if memory error).  If clear is TRUE,
    then unused memory will always be zero. */
extern Vector VCreate( int typeSize, int minSize, int clear );

/** Destroys a growable array (or does nothing if passed NULL) */
extern void VDestroy( Vector v );

/** Returns the number of elements currently used */
extern int VSize( Vector v );

/** Returns a pointer to the index i.  If i is larger than the array,
    it will be expanded to at least i.  If an expansion fails, returns NULL. */
extern void* VIndex( Vector v, int i );

#ifdef __cplusplus
}
#endif

#endif /* __VECTOR_H */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
