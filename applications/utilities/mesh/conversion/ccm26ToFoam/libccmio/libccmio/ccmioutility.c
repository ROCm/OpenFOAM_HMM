#ifndef CCMIO_UTILITY_C
#define CCMIO_UTILITY_C

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmioutility.c,v $
 *  Author:  Geoff Prewett
 *  Date:    August 28, 2003
 *
 *
 *  Star File Format Library - Copyright (C) 2003 by adapco, Ltd.
 *
 *  This program is the property of adapco, Ltd. and contains
 *  confidential and proprietary information.  The unauthorized use,
 *  distribution, or duplication of this program is prohibited.
 *  All rights reserved.
 *
 *  $Id: ccmioutility.c,v 1.13 2006/06/06 12:09:54 jal Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "ccmioutility.h"
#include "ccmiocore.h"
#include "ccmioprivate.h"
#include "libadf/ADF.h"   /* Only for ADF_MAX_DIMENSIONS */

#ifndef MAKEDEPEND
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#ifdef _WIN32
    #include <direct.h>
    #include <io.h>
    #define getcwd _getcwd
    #define mktemp _mktemp
#else
    #include <unistd.h>
    #include <strings.h> /* for rindex() */
#endif
#endif

/* Visual C++ 6 and HP compilers can't tell that a static const int is the
   same as a #define */
#define kMaxPath	10000		/* This does not seem to be defined
					   in a reasonable location, so just
					   make it big enough for everyone. */

static CCMIOError CCMIOWriteNodeCore( CCMIOError *err, CCMIONode parent,
				      const char *name, CCMIODataType type,
				      void *value );
static CCMIOError CCMIOReadNodeCore( CCMIOError *err, CCMIONode parent,
				     const char *name, CCMIODataType type,
				     void *value );

CCMIOError CCMIOWriteNodei( CCMIOError *err, CCMIONode parent, const char *name,
			    int value )
{
    return(CCMIOWriteNodeCore(err, parent, name, kCCMIOInt32, (void *)&value));
}

CCMIOError CCMIOWriteNodef( CCMIOError *err, CCMIONode parent, const char *name,
			    float value )
{
    return(CCMIOWriteNodeCore(err, parent, name, kCCMIOFloat32, (void *)&value));
}

CCMIOError CCMIOWriteNoded( CCMIOError *err, CCMIONode parent, const char *name,
			    double value )
{
    return(CCMIOWriteNodeCore(err, parent, name, kCCMIOFloat64, (void *)&value));
}

CCMIOError CCMIOWriteNodeCore( CCMIOError *err, CCMIONode parent,
			       const char *name, CCMIODataType type,
			       void *value )
{
    CCMIOIndex dims = 1;
    CCMIONode child;
    CHECK_ERROR(err);

    CCMIOCreateNode(err, parent, TRUE, name, name, &child);
    CCMIOSetDataType(err, child, type, 1ul, kCCMIOEndArgs);
    CCMIOExtendedADFIO(err, child, kCCMIOWrite, type, 1, &dims,
		       (char *)value, kCCMIOStart, kCCMIOEnd);
    return(*err);
}

CCMIOError CCMIOWriteNodestr( CCMIOError *err, CCMIONode parent,
			      const char *name, const char *value )
{
    CCMIOIndex dims = strlen(value);
    CCMIONode child;
    CHECK_ERROR(err);

    CCMIOCreateNode(err, parent, TRUE, name, name, &child);
    CCMIOSetDataType(err, child, kCCMIOString, dims, kCCMIOEndArgs);
    CCMIOExtendedADFIO(err, child, kCCMIOWrite, kCCMIOString, 1, &dims,
		       (char *)value, kCCMIOStart, kCCMIOEnd);
    return(*err);
}

CCMIOError CCMIOReadNodei( CCMIOError *err, CCMIONode parent, const char *name,
			   int *value )
{
    if (value) *value = 0;
    return(CCMIOReadNodeCore(err, parent, name, kCCMIOInt32, (void *)value));
}

CCMIOError CCMIOReadNodef( CCMIOError *err, CCMIONode parent, const char *name,
			   float *value )
{
    if (value) *value = 0;
    return(CCMIOReadNodeCore(err, parent, name, kCCMIOFloat32, (void *)value));
}

CCMIOError CCMIOReadNoded( CCMIOError *err, CCMIONode parent, const char *name,
			   double *value )
{
    if (value) *value = 0;
    return(CCMIOReadNodeCore(err, parent, name, kCCMIOFloat64, (void *)value));
}

CCMIOError CCMIOReadNodeCore( CCMIOError *err, CCMIONode parent,
			      const char *name, CCMIODataType type, void *value)
{
    CCMIOIndex dims = 1;
    CCMIONode child;
    CHECK_ERROR(err);
    if (!value)  return(*err = kCCMIOBadParameterErr);

    CCMIOGetNode(err, parent, name, &child);
    CCMIOExtendedADFIO(err, child, kCCMIORead, type, 1, &dims,
		       (char *)value, kCCMIOStart, kCCMIOEnd);
    return(*err);
}

CCMIOError CCMIOReadNodestr( CCMIOError *err, CCMIONode parent,
			     const char *name, char **value )
{
    {
    CCMIONode child;
    CHECK_ERROR_AND_CLEAR_PTR(err, value, 0);

    CCMIOGetNode(err, parent, name, &child);
    if (*err == kCCMIONoErr)
    {
	CCMIOSize size;

	CCMIOGetDataSize(err, child, &size);
	*value = (char *)malloc((size + 1) * sizeof(char));
	CCMIOExtendedADFIO(err, child, kCCMIORead, kCCMIOString, 1, &size,
			   (char *)*value, kCCMIOStart, kCCMIOEnd);
	if (*err == kCCMIONoDataErr)
	{
	    (*value)[0] = '\0';
	    *err = kCCMIONoErr;
	}
	else
	    (*value)[size] = '\0';
    }
    return(*err);
    }
}

/* CCMIORead*() and CCMIOWrite*() are actually included in ccmio.c, in order to
   limit the amount of fake templating. */

CCMIOError CCMIOGetNextChildWithLabel( CCMIOError *err, CCMIONode parent,
				       const char *label, int *n,
				       CCMIONode *child )
{
    char nodeLabel[kCCMIOMaxStringLength + 1];
    CHECK_ERROR(err);
    if (!child || !n) return(*err = kCCMIOBadParameterErr);

    do {
	CCMIOGetNextChild(err, parent, n, child);
	if ((*err) == kCCMIONoErr)
	{
	    CCMIOGetLabel(err, *child, nodeLabel);
	    if (strcmp(label, nodeLabel) == 0)
		return(*err);
	}
    } while ((*err) == kCCMIONoErr);

    return(*err);
}

CCMIOError CCMIOCompress( CCMIOError *err, char *filename )
{
    char *dirPtr, *tmpFilename, *moved = NULL;
    char basename[kMaxPath + 1];
    unsigned int bytes = 0, i;
    CCMIONode origRoot, copyRoot;
    CHECK_ERROR(err);
    if (!filename) return(*err = kCCMIOBadParameterErr);

#if (defined(_WIN32) && !defined(__NUTC__))
    if ((dirPtr = strrchr(filename, '\\')) == NULL)
#else
    if ((dirPtr = rindex(filename, '/')) == NULL)
#endif
    {
	if (!getcwd(basename, kMaxPath))
	    return(*err = kCCMIOIOErr);
	bytes = strlen(basename);
    }
    else
    {
	bytes = dirPtr - filename;
	if (bytes > kMaxPath)
	    bytes = kMaxPath;
	strncpy(basename, filename, bytes);
	basename[bytes] = '\0';
    }

    basename[bytes++] = '/';
    for (i = 0;  i < 6;  ++i)
	basename[bytes++] = 'X';
    basename[bytes] = '\0';
    tmpFilename = strdup(basename);
    if (!mktemp(tmpFilename))
	return(*err = kCCMIOIOErr);

    *err = CCMIOOpen(filename, kCCMIORead, &origRoot);
    if ((*err) != kCCMIONoErr)
	goto error;
    *err = CCMIOOpen(tmpFilename, kCCMIOWrite, &copyRoot);
    if ((*err) != kCCMIONoErr)
    {
	CCMIOClose(origRoot);
	goto error;
    }

    CCMIOCopyNode(err, origRoot, copyRoot, TRUE);

    CCMIOClose(copyRoot);
    CCMIOClose(origRoot);

    /* Rename the original file (just in case the rename of the new file fails),
       rename the temporary file to the original name, and delete the moved
       file */
    if ((*err) == kCCMIONoErr)
    {
	moved = strdup(basename);
	if (!mktemp(moved))
	    goto error;
	errno = 0;
	rename(filename, moved);
	if (errno)
	    goto error;
	rename(tmpFilename, filename);
	if (errno)
	{
	    rename(moved, filename);
	    goto error;
	}
	remove(moved);
    }
    else
	goto error;

    free(moved);
    free(tmpFilename);
    return(*err);

 error:
    remove(tmpFilename);
    free(moved);
    free(tmpFilename);
    return(*err = kCCMIOIOErr);
}

CCMIOError CCMIOCopyNode( CCMIOError *err, CCMIONode origNode,
			  CCMIONode copyNode, int copyExists )
{
    char *bytes = NULL;
    char name[kCCMIOMaxStringLength + 1], label[kCCMIOMaxStringLength + 1];
    int adfErr, n = 0, nDims, dims[ADF_MAX_DIMENSIONS];
    unsigned int size;
    CCMIONode copyParent, origChild;
    CCMIODataType type;
    CHECK_ERROR(err);

    CCMIOGetName(err, origNode, name);
    CCMIOGetLabel(err, origNode, label);
    CCMIOGetDataType(err, origNode, &type);
    /* Can't call CCMIOGetDimensions or CCMIOGetDataSize here, because
       we are doing a low-level copy of just this node, not any of the
       extended data nodes.  We'll copy those when we copy the children. */
    GetADFNodeDimensions(err, origNode, &nDims, dims);
    size = GetADFNodeDataSize(err, origNode);

    if (size)
    {
      bytes = (char*)malloc(size);
	if (!bytes)
	{
	    free(dims);
	    return(*err = kCCMIONoMemoryErr);
	}
	ADF_Read_All_Data(origNode.node, bytes, &adfErr);
	if (IsADFError(adfErr))
	    return(*err = ADFToCCMIOError(adfErr));
    }

    if (copyExists)
    {
	/* Can't set the name of the root node, because it doesn't have a
	   parent */
	CCMIONode rootNode;
	CCMIOGetRootNode(err, copyNode, &rootNode);
	if (!CCMIOAreNodesEqual(copyNode, rootNode))
	    CCMIOSetName(err, copyNode, name);

	CCMIOSetLabel(err, copyNode, name);
    }
    else
    {
	copyParent = copyNode;
	/* Need to pass TRUE because CCMIOOpen always creates a "general"
	   node;  otherwise we will get duplicate node errors */
	CCMIOCreateNode(err, copyParent, TRUE, name, label, &copyNode);
    }
    /* Can't call CCMIOSetDataTypev() for the same reason as above */
    ADF_Put_Dimension_Information(copyNode.node, CCMIOGetDataTypeADFName(type),
				  nDims, dims, &adfErr);
    *err = ADFToCCMIOError(adfErr);
    if (size && *err == kCCMIONoErr)
    {
	ADF_Write_All_Data(copyNode.node, bytes, &adfErr);
	if (IsADFError(adfErr))
	    return(*err = ADFToCCMIOError(adfErr));
    }

    free(bytes);

    while (CCMIOGetNextChild(err, origNode, &n, &origChild) == kCCMIONoErr)
	CCMIOCopyNode(err, origChild, copyNode, FALSE);
    if ((*err) == kCCMIONoNodeErr)
	*err = kCCMIONoErr;

    return(*err);
}
#ifdef __cplusplus
}
#endif
#endif /* CCMIO_UTILITY_C */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
