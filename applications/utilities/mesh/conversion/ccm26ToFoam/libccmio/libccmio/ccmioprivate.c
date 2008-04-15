#ifndef CCMIO_PRIVATE_C
#define CCMIO_PRIVATE_C

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmioprivate.c,v $
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
 *  $Id: ccmioprivate.c,v 1.7 2006/06/05 21:12:16 geoffp Exp $
 */

#ifndef MAKEDEPEND
#include <string.h>
#endif

#include "ccmiocore.h"
#include "ccmioprivate.h"
#include "libadf/ADF.h"

#if Debugging
#ifndef MAKEDEPEND
#include <stdio.h>
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int Sum( int n, int start, int vals[] );

/* -------------------------------------------------------------------------- */
struct _CCMIOTypeInfo {
    char name[8];
    unsigned int size;
};
struct _CCMIOTypeInfo gCCMIOTypeTable[(int)kCCMIOLastType] = {
  { "R4", 4 }, { "R8", 8 }, { "I4", 4 }, { "I8", 8 }, { "C1", 1 }, { "MT", 0 }};

unsigned int CCMIOGetDataTypeSize( CCMIODataType type )
{
    return(gCCMIOTypeTable[(int)type].size);
}

const char* CCMIOGetDataTypeADFName( CCMIODataType type )
{
    return(gCCMIOTypeTable[(int)type].name);
}

CCMIODataType CCMIOGetCCMIODataType( const char *dataStr )
{
    int i;
    char str[3];	/* We will not be using fancy ADF data strings... */

    memcpy(str, dataStr, 2);
    str[2] = '\0';

    for (i = 0;  i < (int)kCCMIOLastType;  ++i)
    {
	if (strcmp(str, gCCMIOTypeTable[i].name) == 0)
	    return((CCMIODataType)i);
    }
    return(kCCMIOBadType);
}

/* -------------------------------------------------------------------------- */
void MakeInvalidNode( CCMIONode *node )
{
    node->node = 0.0;
    node->parent = 0.0;
}

int IsRootNode( CCMIONode node )
{
    if (node.node == node.parent)
	return(TRUE);
    return(FALSE);
}

CCMIOError IsSameFormat( CCMIONode node, CCMIODataType type, int dimSize )
{
    char typeStr[ADF_DATA_TYPE_LENGTH + 1];
    int err, n;

    ADF_Get_Data_Type(node.node, typeStr, &err);
    if (IsADFError(err)) return(ADFToCCMIOError(err));
    if (strcmp(gCCMIOTypeTable[(int)type].name, typeStr) != 0)
	return(kCCMIOWrongDataTypeErr);

    ADF_Get_Number_of_Dimensions(node.node, &n, &err);
    if (IsADFError(err)) return(ADFToCCMIOError(err));
    if (n != dimSize && n > 0)
	/* n of zero means "no data" to ADF, but to us it means
	   one dimension with zero size;  dimSize is never < 1. */
	return(kCCMIOWrongDataTypeErr);
   return(kCCMIONoErr);
}

int IsADFError( int adfErr )
{
    if (adfErr == -1)
	return(FALSE);
    return(TRUE);
}

int ParseArgs( va_list args, CCMIOIndex *out )
{
    int n = 0, val;
    while((val = va_arg(args, CCMIOIndex)) != kCCMIOEndArgs &&
	  n < kCCMIOMaxDimension)
	out[n++] = val;
    if (n == 0)
	return(n);

#if StoreCStyleArrays
    FortranToCArray(n, out);
#endif
#if Debugging
    if (n >= kCCMIOMaxDimension)
	fprintf(stderr, "ParseArgs:  probably did not end a call to a var-arg CCMIO function with kCCMIOEndArgs!\n");
#endif

    return(n);
}

#if StoreCStyleArrays
void FortranToCArray( int size, int *ary )
{
    int i, tmp;

    for (i = 0;  i < size / 2;  ++i)
    {
	tmp = ary[i];
	ary[i] = ary[size - i - 1];
	ary[size - i - 1] = tmp;
    }
}
#endif /* StoreCStyleArrays */

int CalcOffset( int n, int coord[], int dimWidth[] )
{
    int i, offset;

#if StoreCStyleArrays
    offset = coord[n - 1];
    for (i = n - 2;  i >= 0;  --i)
	offset += coord[i] * Sum(n - i, i + 1, dimWidth);
    return(offset);
#else
    offset = coord[0];
    for (i = 1;  i < n;  ++i)
	offset += coord[i] * Sum(i, i - 1, dimWidth);
    return(offset);
#endif /* StoreCStyleArrays */
}

CCMIOError GetADFNodeDimensions( CCMIOError *err, CCMIONode node,
				 int *nDims, int *dims )
{
    int i, adfDims[ADF_MAX_DIMENSIONS];
    ADFError adfErr;

    CHECK_ERROR(err);

    ADF_Get_Number_of_Dimensions(node.node, nDims, &adfErr);
    if (IsADFError(adfErr))
	return(*err = ADFToCCMIOError(adfErr));

    if (dims)
    {
	ADF_Get_Dimension_Values(node.node, adfDims, &adfErr);
	if (adfErr == 32 /* No data */ || adfErr == 27 /* Dimension is zero */)
	{
	    *nDims = 0;
	    dims[0] = 0;
	}
	else if (IsADFError(adfErr))
	    return(*err = ADFToCCMIOError(adfErr));
	
	for (i = 0;  i < *nDims;  ++i)
	    dims[i] = (unsigned long)adfDims[i];
    }

    return(*err);
}

CCMIOSize GetADFNodeDataSize( CCMIOError *err, CCMIONode node )
{
    int i, n, dims[ADF_MAX_DIMENSIONS];
    CCMIOSize bytes;
    CCMIODataType type;

    CHECK_ERROR(err);
    
    bytes = 0;
    if (CCMIOGetDataType(err, node, &type))
	return(*err);
    bytes = CCMIOGetDataTypeSize(type);

    GetADFNodeDimensions(err, node, &n, dims);
    if (*err != kCCMIONoErr || n == 0)
    {
	bytes = 0;
	return(*err = kCCMIONoErr);
    }
    for (i = 0;  i < n;  i++)
	bytes *= dims[i];
    return(bytes);
}

CCMIOError ADFToCCMIOError( int adfErr )
{
    static int errorConversionTable[62][2] = {
	{ 0, kCCMIOInternalErr },	/* Not and ADF err;  placeholder only */
	{ 1, kCCMIOBadParameterErr },	/* Int value less than minimum */
	{ 2, kCCMIOBadParameterErr },	/* Int value greater than maximum */
	{ 3, kCCMIOBadParameterErr },	/* String too short */
	{ 4, kCCMIOBadParameterErr },	/* String too long */
	{ 5, kCCMIOCorruptFileErr },	/* String length corrupted */
	{ 6, kCCMIOInternalErr },	/* Too many ADF files open */
	{ 7, kCCMIOInternalErr },	/* ADF file status not recognized */
	{ 8, kCCMIOPermissionErr },	/* ADF file open err */
	{ 9, kCCMIOBadParameterErr }, 	/* ADF file not open */
	{ 10, kCCMIOBadParameterErr },	/* ADF file out of legal range */
	{ 11, kCCMIOInternalErr },	/* Block/offset out of legal range */
	{ 12, kCCMIOInternalErr },	/* NULL string pointer */
	{ 13, kCCMIOIOErr },		/* FSEEK error */
	{ 14, kCCMIOIOErr },		/* FWRITE error */
	{ 15, kCCMIOIOErr },		/* FREAD error */
	{ 16, kCCMIOInternalErr },	/* Memory corruption */
	{ 17, kCCMIOInternalErr },	/* Bad disk boundary */
	{ 18, kCCMIOInternalErr },	/* Requested new file already exists */
	{ 19, kCCMIOCorruptFileErr },	/* ADF file format not recognized */
	{ 20, kCCMIOCorruptFileErr }, 	/* Attempt to free RootNode info */
	{ 21, kCCMIOCorruptFileErr }, 	/* Attempt to free FreeChunkTable info*/
	{ 22, kCCMIOInternalErr },	/* Requested old file doesn't exist */
	{ 23, kCCMIOInternalErr },	/* Code unimplemented */
	{ 24, kCCMIOCorruptFileErr }, 	/* Bad subnode entries */
	{ 25, kCCMIONoMemoryErr },	/* Out of memory */
	{ 26, kCCMIODuplicateNodeErr }, /* Tried to create a duplicate child */
	{ 27, kCCMIOCorruptFileErr }, 	/* Node has no dimensions */
	{ 28, kCCMIOCorruptFileErr }, 	/* Node has illegal dimensions */
	{ 29, kCCMIONoNodeErr },	/* Child is not child of parent */
	{ 30, kCCMIOInternalErr },	/* Datatype is too long */
	{ 31, kCCMIOInternalErr },	/* Bad datatype */
	{ 32, kCCMIOInternalErr },	/* NULL pointer */
	{ 33, kCCMIONoDataErr },	/* Node has no data */
	{ 34, kCCMIOInternalErr },	/* Couldn't zero out memory */
	{ 35, kCCMIOInternalErr },	/* Not enough data for request */
	{ 36, kCCMIOInternalErr },	/* Bad end value */
	{ 37, kCCMIOInternalErr },	/* Bad stride value */
	{ 38, kCCMIOBadParameterErr }, 	/* Minimum val is greater than maximum*/
	{ 39, kCCMIOInternalErr },	/* Format of this machine is unknown */
	{ 40, kCCMIOInternalErr },	/* Conversion to/from unknown machine format */
	{ 41, kCCMIOInternalErr },	/* No conversion necessary */
	{ 42, kCCMIOInternalErr },	/* Data format not supported on this machine */
	{ 43, kCCMIOIOErr },		/* File close error */
	{ 44, kCCMIOInternalErr },	/* Numeric overflow/underflow in conversion */
	{ 45, kCCMIOInternalErr },	/* Bad start value */
	{ 46, kCCMIOBadParameterErr },	/* Value of zero not allowed */
	{ 47, kCCMIOBadParameterErr },	/* Bad dimension value */
	{ 48, kCCMIOInternalErr },	/* Error state must be zero or one */
	{ 49, kCCMIOInternalErr },	/* Dimension specs for memory and disk unequal*/
	{ 50, kCCMIOBadParameterErr },	 /* Link depth too large;  recursive link? */
	{ 51, kCCMIOInternalErr },	/* Node is not a link */
	{ 52, kCCMIOBadLinkErr },	/* Node linked to does not exist */
	{ 53, kCCMIOBadLinkErr },	/* ADF file linked to does not exist */
	{ 54, kCCMIOBadParameterErr }, 	/* Node ID of 0.0 is invalid */
	{ 55, kCCMIOCorruptFileErr }, 	/*Incomplete data when reading multiple blocks*/
	{ 56, kCCMIOBadParameterErr }, 	/* Node name contains invalid characters */
	{ 57, kCCMIOInternalErr },	/* ADF file version incompatible */
	{ 58, kCCMIOBadParameterErr }, 	/* Nodes are not from same file */
	{ 59, kCCMIOInternalErr },	/* Priority stack err */
	{ 60, kCCMIOInternalErr },	/* Machine and file formats incompatible */
	{ 61, kCCMIOIOErr }		/* Flush error */
    };

#if Debugging
    /* This didn't turn out to be as useful as I thought.  Mostly you just
       get a lot of ADF error 29, which is usually expected as it is often used
       to terminate loops */
/*
    if (IsADFError(adfErr))
	fprintf(stderr, "Converting ADF error %d to CCMIO error %d.\n", adfErr, 
		errorConversionTable[adfErr][1]);
*/
#endif
    if (adfErr == -1)
	return(kCCMIONoErr);
    return ((CCMIOError) errorConversionTable[adfErr][1]);
}

/* -------------------------------------------------------------------------- */
int Sum( int n, int start, int vals[] )
{
    int i, result = 0;

    for (i = start;  i < n;  ++i)
	result += vals[i];
    return(result);
}

#ifdef __cplusplus
}
#endif
#endif /* CCMIO_PRIVATE_C */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
