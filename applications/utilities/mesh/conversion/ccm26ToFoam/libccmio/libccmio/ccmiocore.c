#ifndef CCMIO_CORE_C
#define CCMIO_CORE_C

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmiocore.c,v $
 *  Author:  Geoff Prewett
 *  Date:    July 28, 2003
 *
 *
 *  Star File Format Library - Copyright (C) 2003 by adapco, Ltd.
 *
 *  This program is the property of adapco, Ltd. and contains
 *  confidential and proprietary information.  The unauthorized use,
 *  distribution, or duplication of this program is prohibited.
 *  All rights reserved.
 *
 *  $Id: ccmiocore.c,v 1.12 2006/06/06 15:50:44 geoffp Exp $
 */

#ifndef MAKEDEPEND
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#endif

#include "libadf/ADF.h"
#include "ccmio.h"
#include "ccmiocore.h"
#include "ccmioprivate.h"
#include "ccmioversion.h"
#include "ccmioutility.h"

const char kCCMIODefaultTitle[] = "";

#ifdef __cplusplus
extern "C" {
#endif

CCMIOError CCMIOGetDimensionsCore( CCMIOError *err, CCMIONode node, int *nDims,
				   CCMIOIndex *dims );
CCMIOError CCMIOSetDataTypeCore( CCMIOError *err, CCMIONode node,
				 CCMIODataType type, int nDims,
				 const CCMIOIndex *dims );

/** Opens the data file.  The file will be created if it does not already exist.
    \param filename	Name of the file to open.
    \param mode		Either kCCMIORead or kCCMIOWrite.
    \param root 	Returns the root node of the file. */
CCMIOError CCMIOOpen( const char *filename, CCMIOIOType mode, CCMIONode *root )
{
    int ver;
    ADFError adfErr;
    CCMIOError err = kCCMIONoErr;

    if (!root) return(kCCMIOBadParameterErr);

    ADF_Database_Open(filename, (mode == kCCMIORead) ? "READ_ONLY" : "UNKNOWN",
		      "NATIVE", &root->node, &adfErr);
    root->parent = root->node;
    if (IsADFError(adfErr))
	return(ADFToCCMIOError(adfErr));

    /* If there isn't a version, set one */
    if (CCMIOGetVersion(NULL, *root, &ver) != kCCMIONoErr)
    {
	if (mode == kCCMIORead)
	    return(err = kCCMIOCorruptFileErr);
	CCMIOSetVersion(&err, *root, kCCMIOVersion);
	CCMIOSetTitle(&err, *root, kCCMIODefaultTitle);
    }
    return(err);
}

CCMIOError CCMIOClose( CCMIONode root )
{
    ADFError adfErr;

    ADF_Database_Close(root.node, &adfErr);
    return(ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOGetNode( CCMIOError *err, CCMIONode parent, const char *path,
			 CCMIONode *node )
{
    char const *end;
    char *parentStr;
    int size;
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!node || !path) return(*err = kCCMIOBadParameterErr);

    if (*path == '/')	/* Always make paths relative to 'parent' */
	path++;

    if ((end = strrchr(path, '/')) == NULL)
	node->parent = parent.node;
    else
    {
	size = end - path;	/* The '/' will become the '\0' */
	parentStr = (char *)malloc((size + 1) * sizeof(char));
	memcpy(parentStr, path, size);
	parentStr[size] = '\0';
	ADF_Get_Node_ID(parent.node, parentStr, &node->parent, &adfErr);
	free(parentStr);
	if (IsADFError(adfErr))
	    return(*err = ADFToCCMIOError(adfErr));
    }
    ADF_Get_Node_ID(parent.node, path, &node->node, &adfErr);

    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOGetNumberOfChildren( CCMIOError *err, CCMIONode parent, int *n )
{
    ADFError adfErr;
    CHECK_ERROR_AND_CLEAR_PTR(err, n, 0);

    ADF_Number_of_Children(parent.node, n, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOGetNextChild( CCMIOError *err, CCMIONode parent, int *n,
			      CCMIONode *child )
{
    char name[ADF_NAME_LENGTH + 1];
    int nNames;
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!n || !child) return(*err = kCCMIOBadParameterErr);

    /* ADF uses FORTRAN indexing, so add 1 to *n to convert. */
    ADF_Children_Names(parent.node, *n + 1, 1, ADF_NAME_LENGTH, &nNames, name,
		       &adfErr);
    if (nNames == 0)
	return(*err = kCCMIONoNodeErr);
    else if (nNames != 1)
	return(*err = kCCMIOInternalErr);
    if (IsADFError(adfErr))
	return(*err = ADFToCCMIOError(adfErr));

    ADF_Get_Node_ID(parent.node, name, &child->node, &adfErr);
    child->parent = parent.node;
    *n += 1;

    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOGetName( CCMIOError *err, CCMIONode node, char *name )
{
    ADFError adfErr;

    CHECK_ERROR_AND_CLEAR_PTR(err, name, 0);
    if (!name) return(*err = kCCMIOBadParameterErr);

    *name = '\0';	/* in case of error... */
    ADF_Get_Name(node.node, name, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOSetName( CCMIOError *err, CCMIONode node, const char *name )
{
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!name) return(*err = kCCMIOBadParameterErr);

    ADF_Put_Name(node.parent, node.node, name, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOGetLabel( CCMIOError *err, CCMIONode node, char *label )
{
    ADFError adfErr;

    CHECK_ERROR_AND_CLEAR_PTR(err, label, 0);

    *label = '\0';	/* in case of error... */
    ADF_Get_Label(node.node, label, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOSetLabel( CCMIOError *err, CCMIONode node, const char *label )
{
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!label) return(*err = kCCMIOBadParameterErr);

    ADF_Set_Label(node.node, label, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOCreateNode( CCMIOError *err, CCMIONode parent, int openDup,
			    const char *name, const char *label,
			    CCMIONode *node )
{
    ADFError adfErr;
    CCMIONode tmpNode;

    CHECK_ERROR(err);
    if (!name) return(*err = kCCMIOBadParameterErr);
    if (!node) node = &tmpNode;

    node->parent = parent.node;
    ADF_Create(parent.node, name, &node->node, &adfErr);
    if (adfErr == 26 && openDup)  /* Duplicate node err */
    {
	if (!CCMIOGetNode(err, parent, name, node))
	    return(*err);
    }
    else if (IsADFError(adfErr))
	return(*err = ADFToCCMIOError(adfErr));

    if (label)
	ADF_Set_Label(node->node, label, &adfErr);

    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOCreateLink( CCMIOError *err, CCMIONode parent, const char *name,
			    const char *filename, const char *destName,
			    CCMIONode *node )
{
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!name || !filename || !destName || !node)
	return(*err = kCCMIOBadParameterErr);

    ADF_Link(parent.node, name, filename, destName, &node->node, &adfErr);
    node->parent = parent.node;
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIODeleteNode( CCMIOError *err, CCMIONode node )
{
    ADFError adfErr;
    CHECK_ERROR(err);

    ADF_Delete(node.parent, node.node, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIODeleteAllChildren( CCMIOError *err, CCMIONode node )
{
    int i, n;
    CCMIONode child;
    CHECK_ERROR(err);
    
    CCMIOGetNumberOfChildren(err, node, &n);
    for (i = 0;  n > 0 && !(*err);  )
    {
	CCMIOGetNextChild(err, node, &i, &child);
	CCMIODeleteNode(err, child);
	i = 0;	/* CCMIOGetNextChild increments this, but we want to delete
		   the first child each time */
	CCMIOGetNumberOfChildren(err, node, &n);
    }
    return(*err);
}

CCMIOError CCMIOMoveNode( CCMIOError *err, CCMIONode node, CCMIONode newParent )
{
    ADFError adfErr;
    CHECK_ERROR(err);

    ADF_Move_Child(node.parent, node.node, newParent.node, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

/* Requires:	nDims is a valid pointer.
   		dims is an array of ADF_MAX_DIMENSIONS, or a NULL pointer
   Ensures:	dims is set to the dimensions of the node (if it is a regular
		node this will just be the dimensions of the node;  if it is
		an extended node, these will be the dimensions of the extended
		data.) */
CCMIOError CCMIOGetDimensionsCore( CCMIOError *err, CCMIONode node, int *nDims,
				   CCMIOIndex *dims )
{
    int i, adfDims[ADF_MAX_DIMENSIONS];
    ADFError adfErr;
    CCMIONode nodeEx;

    { /* CHECK_ERROR declares a variable, so it needs its own scope */
    CHECK_ERROR(err);
    if (!nDims) return(*err = kCCMIOBadParameterErr);
    *nDims = 0;

    if (CCMIOGetNode(NULL, node, kCCMIOExtendedDataSize, &nodeEx) ==kCCMIONoErr)
    { /* This node is an extendend node */
	/* Check that this node is really a 1D array like it should be */
	i = 0;
	ADF_Get_Number_of_Dimensions(nodeEx.node, &i, &adfErr);
	if (IsADFError(adfErr) || i != 1)
	    return(*err = kCCMIOCorruptFileErr);

	/* Now find out how big it is */
	ADF_Get_Dimension_Values(nodeEx.node, adfDims, &adfErr);
	if (IsADFError(adfErr))
	    return(*err = kCCMIOCorruptFileErr);
	*nDims = adfDims[0];

	/* Now read in the number of dimensions */
	if (dims)
	{
	    if (sizeof(CCMIOIndex) > sizeof(int))
		ADF_Read_All_Data(nodeEx.node, (char *)dims, &adfErr);
	    else
	    {
		unsigned long long dims64[ADF_MAX_DIMENSIONS];
		ADF_Read_All_Data(nodeEx.node, (char *)dims64, &adfErr);
		for (i = 0;  i < *nDims;  ++i)
		    dims[i] = dims64[i];
	    }
	    if (IsADFError(adfErr))
		return(*err = kCCMIOCorruptFileErr);
	}
    }
    else
    {   /* Just a regular node;  call the regular ADF functions */
	ADF_Get_Number_of_Dimensions(node.node, nDims, &adfErr);
	if (IsADFError(adfErr))
	    return(*err = ADFToCCMIOError(adfErr));

	if (dims)
	{
	    ADF_Get_Dimension_Values(node.node, adfDims, &adfErr);
	    if (adfErr == 32 /* No data */ ||
		adfErr == 27 /* Dimension is zero */)
	    {
		*nDims = 0;
		dims[0] = 0;
	    }
	    else if (IsADFError(adfErr))
		return(*err = ADFToCCMIOError(adfErr));

	    for (i = 0;  i < *nDims;  ++i)
		dims[i] = (CCMIOIndex)adfDims[i];
	}
    }

    return(*err);
    } /* end CHECK_ERROR scope */
}

CCMIOError CCMIOGetDimensions( CCMIOError *err, CCMIONode node, int *nDims,
			       CCMIOIndex **dims )
{
    if (dims)	*dims = NULL;
    { /* CHECK_ERROR declares a variable, so it needs its own scope */
    CHECK_ERROR(err);
    if (!nDims) return(*err = kCCMIOBadParameterErr);

    CCMIOGetDimensionsCore(err, node, nDims, NULL);
	
    if (dims)
    {
	*dims = (CCMIOIndex *)malloc((*nDims) * sizeof(CCMIOIndex));
	/* According to the online ADF documentation,
	   ADF_Get_Number_of_Dimensions (called inside CCMIOGetDimensionsCore)
	   returns err 27 if the number of dimensions is zero.  According to
	   the actual code, it doesn't.  This appears to not be a problem
	   anywhere except on AIX, where malloc(0) produces a NULL pointer.
	   Apparently all the other systems actually give us a pointer.
	   So, if *nDims is zero and we got a NULL pointer, allocate a
	   non-zero one. */
	if (*nDims == 0 && !(*dims))
	    *dims = (CCMIOIndex *)malloc(sizeof(CCMIOIndex));
	CCMIOGetDimensionsCore(err, node, nDims, *dims);
    }

    return(*err);
    } /* end CHECK_ERROR scope */
}

CCMIOError CCMIOGetDataType( CCMIOError *err, CCMIONode node,
			     CCMIODataType *type )
{
    char dataStr[ADF_DATA_TYPE_LENGTH + 1];
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!type) return(kCCMIOBadParameterErr);

    ADF_Get_Data_Type(node.node, dataStr, &adfErr);
    if (IsADFError(adfErr))
	return(*err = ADFToCCMIOError(adfErr));
    *type = CCMIOGetCCMIODataType(dataStr);

    return(*err);
}

CCMIOError CCMIOGetDataSize( CCMIOError *err, CCMIONode node,
			     CCMIOSize *bytes )
{
    int i, dims;
    CCMIOIndex n[ADF_MAX_DIMENSIONS];
    CCMIODataType type;

    CHECK_ERROR(err);
    if (!bytes) return(*err = kCCMIOBadParameterErr);
    
    *bytes = 0;
    if (CCMIOGetDataType(err, node, &type))
	return(*err);
    *bytes = CCMIOGetDataTypeSize(type);

    CCMIOGetDimensionsCore(err, node, &dims, n);
    if (*err != kCCMIONoErr || dims == 0)
    {
	*bytes = 0;
	return(*err = kCCMIONoErr);
    }
    for (i = 0;  i < dims;  i++)
	*bytes *= n[i];
    return(*err);
}

CCMIOError CCMIOSetDataType( CCMIOError *err, CCMIONode node,
			     CCMIODataType type, ... )
{
    va_list args;
    CHECK_ERROR(err);

    va_start(args, type);
    CCMIOvSetDataType(err, node, type, args);
    va_end(args);
    return(*err);
}

CCMIOError CCMIOSetDataTypeCore( CCMIOError *err, CCMIONode node,
				 CCMIODataType type, int nDims,
				 const CCMIOIndex *dims )
{
    int i, adfDims[ADF_MAX_DIMENSIONS];
    CCMIOIndex itemSize, itemsPerNode;
    CCMIOSize size;
    ADFError adfErr;
    CCMIONode nodeEx;

    CHECK_ERROR(err);

    /* ADF does not like zero-size dimensions, so empty strings need a bit
       of tweaking -- set the node to dim 0, i.e. no data. */
    if (type == kCCMIOString && nDims == 1 && dims[0] == 0)
	nDims = 0;

    size = CCMIOGetDataTypeSize(type);
    for (i = 0;  i < nDims;  ++i)
	size *= dims[i];
    if (nDims == 0 || dims[0] == 0)
    {
	itemSize = 0;
	itemsPerNode = 0;
    }
    else
    {
	itemSize = CCMIOGetDataTypeSize(type);
	if (nDims > 1)
	    itemSize *= dims[0];
	itemsPerNode = kCCMIOMaxADFNodeSize /itemSize;/*This will be truncated*/
	if (itemsPerNode == 0)
	    return(*err = kCCMIOArrayDimensionToLargeErr);
    }

    for (i = 0;  i < nDims;  ++i)
	adfDims[i] = (int)dims[i];

    /* Delete any extended data nodes.  They will be recreated if necessary when
       the data is written.  Note that resetting an existing node's data type
       can potentially be a bit slow. */
    i = 0;
    while (CCMIOGetNextChildWithLabel(NULL, node, kCCMIOExtendedDataLabel,
				      &i, &nodeEx) == kCCMIONoErr)
	CCMIODeleteNode(err, nodeEx);

    /* Now we can set the size */
    if (size <= kCCMIOMaxADFNodeSize)
    {
	ADF_Put_Dimension_Information(node.node, CCMIOGetDataTypeADFName(type),
				      nDims, adfDims, &adfErr);
	*err = ADFToCCMIOError(adfErr);
    }
    else
    {
	int j;
	/* This can't be a CCMIOIndex because it *must* be 64-bits (since
	   that's the data type of the node, and CCMIOIndex might change) */
	unsigned long long dimArray64[ADF_MAX_DIMENSIONS];

	/* Setup the main node */
	adfDims[nDims - 1] = itemsPerNode;
	ADF_Put_Dimension_Information(node.node, CCMIOGetDataTypeADFName(type),
				      nDims, adfDims, &adfErr);
	

	/* Create the extended size node */
	CCMIOCreateNode(err, node, TRUE, kCCMIOExtendedDataSize,
			kCCMIOExtendedDataLabel, &nodeEx);
	if (*err != kCCMIONoErr)
	    return(*err);
	ADF_Put_Dimension_Information(nodeEx.node,
				      CCMIOGetDataTypeADFName(kCCMIOInt64),
				      1, &nDims, &adfErr);
	if (IsADFError(adfErr))
	    return(*err = ADFToCCMIOError(adfErr));
	for (j = 0;  j < ADF_MAX_DIMENSIONS && j < nDims;  ++j)
	    dimArray64[j] = dims[j];
	ADF_Write_All_Data(nodeEx.node, (char *)dimArray64, &adfErr);
	if (IsADFError(adfErr))
	    return(*err = ADFToCCMIOError(adfErr));
    }
    
    return(*err);
}

CCMIOError CCMIOvSetDataType( CCMIOError *err, CCMIONode node,
			      CCMIODataType type, va_list args )
{
    int n = 0;
    CCMIOIndex dims[ADF_MAX_DIMENSIONS + 1];

    CHECK_ERROR(err);
    n = ParseArgs(args, dims);
    if (n == 0)
	return(*err = kCCMIOBadParameterErr);

    /* ADF does not like zero-size dimensions, so empty strings need a bit
       of tweaking -- set the node to dim 0, i.e. no data. */
    if (type == kCCMIOString && n == 1 && dims[0] == 0)
	n = 0;

    CCMIOSetDataTypeCore(err, node, type, n, dims);

    return(*err);
}

CCMIOError CCMIOSetDataTypev( CCMIOError *err, CCMIONode node,
			      CCMIODataType type, int nDims,
			      const CCMIOIndex *dims )
{
    int i;
    CCMIOIndex dimsCopy[ADF_MAX_DIMENSIONS+1];
    CHECK_ERROR(err);
    if (nDims > ADF_MAX_DIMENSIONS)  return(*err = kCCMIOBadParameterErr);

    for (i = 0;  i < nDims;  ++i)
	dimsCopy[i] = dims[i];

#if StoreCStyleArrays
    FortranToCArray(n, dimsCopy);
#endif
    /* ADF does not like zero-size dimensions, so empty strings need a bit
       of tweaking -- set the node to dim 0, i.e. no data. */
    if (type == kCCMIOString && nDims == 1 && dimsCopy[0] == 0)
	nDims = 0;

    CCMIOSetDataTypeCore(err, node, type, nDims, dims);
    return(*err);
}
/*
CCMIOError CCMIOReadDataPoint( CCMIOError *err, CCMIONode node, void *data, ...)
{
    int i, n = 0;
    int start[ADF_MAX_DIMENSIONS + 1], end[ADF_MAX_DIMENSIONS + 1];
    int stride[ADF_MAX_DIMENSIONS + 1];
    int nDims = 1, dim = 1, memStart = 1, memEnd = 1, memStride = 1;
    va_list argp;
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!data) return(*err = kCCMIOBadParameterErr);

    va_start(argp, data);
    n = ParseArgs(argp, start);
    va_end(argp);

    if (n == 0)
	return(*err = kCCMIOBadParameterErr);

    for (i = 0;  i < n;  ++i)
    {
	start[i]++;	// FORTRAN arrays begin at 1 not zero
	end[i] = start[i];
	stride[i] = 1;
    }
	     
    ADF_Read_Data(node.node, start, end, stride, nDims, &dim, &memStart,
		  &memEnd, &memStride, (char *) data, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOReadData( CCMIOError *err, CCMIONode node, void *data,
			  CCMIODataType expected, int dimsExpected )
{
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!data) return(*err = kCCMIOBadParameterErr);

    if (dimsExpected <= 0)	// We always expect at least one dimension,
	dimsExpected = 1;	// even if the dimension's size is zero
    *err = IsSameFormat(node, expected, dimsExpected);
    if (*err != kCCMIONoErr)
	return(*err);

    ADF_Read_All_Data(node.node, (char *) data, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOWriteDataPoint( CCMIOError *err, CCMIONode node,
				void *data, ... )
{
    int i, n = 0;
    int start[ADF_MAX_DIMENSIONS + 1], end[ADF_MAX_DIMENSIONS + 1];
    int stride[ADF_MAX_DIMENSIONS + 1];
    int nDims = 1, dim = 1, memStart = 1, memEnd = 1, memStride = 1;
    va_list argp;
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!data) return(*err = kCCMIOBadParameterErr);

    va_start(argp, data);
    n = ParseArgs(argp, start);
    va_end(argp);
    // n == 0 is not an error;  it just means that instead of leaving
    // the node with no data, you have set it to have a specific type
    // of data, of which you are giving it nothing.  (e.g. an empty
    // string, instead of no data).  ADF does not see the difference,
    // so return now.
    if (n == 0)
	return(kCCMIONoErr);

    for (i = 0;  i < n;  ++i)
    {
	start[i]++;	// FORTRAN arrays begin at 1 not zero 
	end[i] = start[i];
	stride[i] = 1;
    }

    ADF_Write_Data(node.node, start, end, stride, nDims, &dim, &memStart,
		   &memEnd, &memStride, (char *) data, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}

CCMIOError CCMIOWriteData( CCMIOError *err, CCMIONode node, const void *data )
{
    ADFError adfErr;

    CHECK_ERROR(err);
    if (!data) return(*err = kCCMIOBadParameterErr);

    ADF_Write_All_Data(node.node, (const char *) data, &adfErr);
    return(*err = ADFToCCMIOError(adfErr));
}
*/
int CCMIOAreNodesEqual( CCMIONode node1, CCMIONode node2 )
{
    return(node1.node == node2.node && node1.parent == node2.parent);
}

CCMIOError CCMIOGetRootNode( CCMIOError *err, CCMIONode node, CCMIONode *root )
{
    ADFError adfErr;
    CHECK_ERROR(err);
    if (!root) return(*err = kCCMIOBadParameterErr);

    ADF_Get_Root_ID(node.node, &root->node, &adfErr);
    if (IsADFError(adfErr))
	return(*err = kCCMIOBadParameterErr);
    root->parent = root->node;
    return(*err);
}

#ifdef __cplusplus
    }
#endif
#endif /* CCMIO_CORE_C */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
