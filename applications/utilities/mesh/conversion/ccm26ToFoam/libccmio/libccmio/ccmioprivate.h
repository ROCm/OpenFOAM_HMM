#ifndef CCMIO_PRIVATE_H
#define CCMIO_PRIVATE_H

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmioprivate.h,v $
 *  Author:  Geoff Prewett
 *  Date:    July 31, 2003
 *
 *
 *  Star File Format Library - Copyright (C) 2003 by adapco, Ltd.
 *
 *  This program is the property of adapco, Ltd. and contains
 *  confidential and proprietary information.  The unauthorized use,
 *  distribution, or duplication of this program is prohibited.
 *  All rights reserved.
 *
 *  $Id: ccmioprivate.h,v 1.7 2006/06/05 21:12:16 geoffp Exp $
 */

#ifndef MAKEDEPEND
#include <stdarg.h>
#endif

#include "libccmio/ccmio.h"

#ifdef __cplusplus
extern "C" {
#endif

#define StoreCStyleArrays	0	/* If set to 1, stores arrays in the
					   ADF file in C order (by reversing
					   the dimensions) */

#define kCCMIOMaxADFNodeSize	2147479552  /* 2 GB - 4096 (1 ADF block size)*/
/* Test the extended ADF feature with these sizes (keep the node size small!) */
/*#define kCCMIOMaxADFNodeSize	100000000 */
/*#define kCCMIOMaxADFNodeSize	7000000 */
/*#define kCCMIOMaxADFNodeSize	6000000 */
/*#define kCCMIOMaxADFNodeSize	128 */
/*#define kCCMIOMaxADFNodeSize	110 */
/*#define kCCMIOMaxADFNodeSize	36 */
/*#define kCCMIOMaxADFNodeSize	24 */  /* Use this for the ccmtest test also */
/*#define kCCMIOMaxADFNodeSize	100000 */
#define kCCMIOExtendedDataSize	"ExtendedSize"
#define kCCMIOExtendedDataName	"ExtendedData"
#define kCCMIOExtendedDataLabel	"extendedData"

typedef int ADFError;

/** Writes an invalid CCMIO node in 'node'. */
void MakeInvalidNode( CCMIONode *node );

/** Returns TRUE if 'node' is the root node, FALSE otherwise. */
int IsRootNode( CCMIONode node );

/** Returns the size of the specified datatype */
unsigned int CCMIOGetDataTypeSize( CCMIODataType type );

/** Returns an ADF typename string from the specified type */
char const* CCMIOGetDataTypeADFName( CCMIODataType type );

/** Returns an CCMIODataType from an ADF string. */
CCMIODataType CCMIOGetCCMIODataType( char const *dataStr );

/** Checks to see if the node has the same datatype and dimension
    size as specified */
CCMIOError IsSameFormat( CCMIONode node, CCMIODataType type, int dimSize );

/** Parses a va_list of dimension arguments */
int ParseArgs( va_list args, CCMIOIndex *out );

#if StoreCStyleArrays
/** Converts an array of Fortran dimensions to C dimensions (i.e. reverses
    them */
void FortranToCArray( int size, int *ary );
#endif /* StoreCStyleArrays */

/** Takes a coordinate and size of array and returns the offset
    of the index into the array (if StoreCStyleArrays is enabled then
    the array is assumed to be in C-ordering, otherwise it is in FORTRAN
    ordering). */
int CalcOffset( int n, int coord[], int dimWidth[] );

/** Like CCMIOGetDimensions, except only gives the dimensions of the node
    passed in;  does not take into account extended data */
CCMIOError GetADFNodeDimensions( CCMIOError *err, CCMIONode node,
				 int *nDims, int *dims );

/** Like CCMIOGetDataSize, except only returns the size (in bytes) of the node
    passed in;  does not take into account extended data */
CCMIOSize GetADFNodeDataSize( CCMIOError *err, CCMIONode node );

/** Returns 1 if adfErr is an ADF error, 0 otherwise. */
int IsADFError( int adfErr );

/** Returns the appropriate CCMIOError from the specified adfErr */
CCMIOError ADFToCCMIOError( int adfErr );

/* This is defined in ccmioread.c */
CCMIOError CCMIOExtendedADFIO( CCMIOError *err, CCMIONode node, 
			       CCMIOIOType ioType, CCMIODataType dataType,
			       int nDims, CCMIOIndex const *dims,
			       char *data, CCMIOIndex start, CCMIOIndex end );


#ifdef __cplusplus
}
#endif
#endif /* CCMIO_PRIVATE_H */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
