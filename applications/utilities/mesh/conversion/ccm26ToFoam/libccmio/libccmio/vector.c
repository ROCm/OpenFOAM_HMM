#ifndef __VECTOR_C
#define __VECTOR_C

/*@@
 *  Program: Star File Format Library  - $RCSfile: vector.c,v $
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
 *  $Id: vector.c,v 1.4 2005/01/11 21:51:20 wo Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MAKEDEPEND
#include <stdlib.h>
#include <string.h>
#endif

#include "vector.h"

Vector VCreate( int typeSize, int minSize, int clear )
{
    Vector v;

    if (typeSize <= 0)	/* Avoid bad mallocs and divide by zero in VSize */
	typeSize = 1;
    v = (Vector)malloc(sizeof(struct _Vector));
    if (v)
    {
	if (clear)
	    v->buffer = calloc(minSize, typeSize);
	else
	    v->buffer = malloc(minSize * typeSize);
	if (!v->buffer)
	{
	    free(v);
	    return(NULL);
	}
	v->typeSize = typeSize;
	v->alloc = minSize * typeSize;
	v->size = 0;
	v->clear = clear;
    }
    return(v);
}

void VDestroy( Vector v )
{
    if (!v)
	return;

    free(v->buffer);
    free(v);
}

int VSize( Vector v )
{
    if (!v || !v->buffer)
	return(0);
    return(v->size / v->typeSize);
}

void* VIndex( Vector v, int i )
{
    int offset, size;
    void *tmp;

    if (!v)
	return(NULL);

    offset = i * v->typeSize;
    if (offset + v->typeSize >= v->alloc)
    {
	size = v->alloc * 2;
	if (offset + v->typeSize > size)
	    size = offset + v->typeSize;
	if (v->clear)
	    tmp = calloc(size, 1);
	else
	    tmp = malloc(size);
	if (!tmp)
	    return(NULL);
	if (v->buffer)	/* In case somebody tried VCreate(n, 0); */
	{		/* in which case malloc(0) might return NULL */
	    memcpy(tmp, v->buffer, v->size);
	    free(v->buffer);
	}
	v->buffer = tmp;
	v->alloc = size;
    }
    
    if (offset >= v->size)
	v->size = offset + v->typeSize;
    return((char *)v->buffer + offset);
}

#ifdef __cplusplus
}
#endif

#endif /* __VECTOR_C */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
