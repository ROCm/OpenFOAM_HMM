#ifndef CCMIO_UTILITY_H
#define CCMIO_UTILITY_H

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmioutility.h,v $
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
 *  $Id: ccmioutility.h,v 1.9 2006/06/05 21:12:16 geoffp Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "libccmio/ccmiotypes.h"

/* \{
   \name Utility functions
   \brief Condenses common tasks
   \ingroup core */

/** Creates a child node of parent with given name and value */
extern CCMIOError CCMIOWriteNodei( CCMIOError *err, CCMIONode parent,
				   char const *name, int value );

/** Creates a child node of parent with given name and value */
extern CCMIOError CCMIOWriteNodef( CCMIOError *err, CCMIONode parent,
				   char const *name, float value );

/** Creates a child node of parent with given name and value */
extern CCMIOError CCMIOWriteNoded( CCMIOError *err, CCMIONode parent,
				   char const *name, double value );

/** Creates a child node of parent with given name and string value */
extern CCMIOError CCMIOWriteNodestr( CCMIOError *err, CCMIONode parent,
				     char const *name, char const *value );

/** Reads the value of a child node;  returns \def kCCMIOWrongDataTypeErr if
    the type of the node is not \def kCCMIOInt32. */
extern CCMIOError CCMIOReadNodei( CCMIOError *err, CCMIONode parent,
				  char const *name, int *value );

/** Reads the value of a child node;  returns \def kCCMIOWrongDataTypeErr if
    the type of the node is not \def kCCMIOFloat32. */
extern CCMIOError CCMIOReadNodef( CCMIOError *err, CCMIONode parent,
				  char const *name, float *value );

/** Reads the value of a child node;  returns \def kCCMIOWrongDataTypeErr if
    the type of the node is not \def kCCMIOFloat64. */
extern CCMIOError CCMIOReadNoded( CCMIOError *err, CCMIONode parent,
				  char const *name, double *value );

/** Reads the value of a child node;  returns \def kCCMIOWrongDataTypeErr if
    the type of the node is not \def kCCMIOString.  The string returned is
    allocated by the library and must be freed by the application. */
extern CCMIOError CCMIOReadNodestr( CCMIOError *err, CCMIONode parent,
				    char const *name, char **value );

/* CCMIORead*() and CCMIOWrite*() are actually included in ccmio.c, in order
   to limit the amount of fake templating. */

/** Reads the entire contents of the node into the array provided.  The array
    must be the proper size (which can be determined with
    CCMIOGetDimensions()).  If actual node data is stored in a different format
    than requested, it will be converted.  For multidimensional arrays, the
    parameter isC specifies whether the data should be returned in Fortran or C
    order.  Since data is stored in Fortran order on disk, returning C order
    involves an extra copy. */
extern CCMIOError CCMIORead1i( CCMIOError *err, CCMIONode node, int *data,
			       CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIORead1f( CCMIOError *err, CCMIONode node, float *data,
			       CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIORead1d( CCMIOError *err, CCMIONode node, double *data,
			       CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIORead2i( CCMIOError *err, CCMIONode node, int *data,
			       CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIORead2f( CCMIOError *err, CCMIONode node, float *data,
			       CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIORead2d( CCMIOError *err, CCMIONode node, double *data,
			       CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIORead3i( CCMIOError *err, CCMIONode node, int *data,
			       CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIORead3f( CCMIOError *err, CCMIONode node, float *data,
			       CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIORead3d( CCMIOError *err, CCMIONode node, double *data,
			       CCMIOIndex start, CCMIOIndex end );

/** Writes the entire array to the node.  The node's size and data type will
    be automatically set to the proper values.  See CCMIORead*() for comments on
    the isC parameter. */
extern CCMIOError CCMIOWrite1i( CCMIOError *err, CCMIONode node, CCMIOIndex n,
				int const *data, CCMIOIndex start,
				CCMIOIndex end );
extern CCMIOError CCMIOWrite1f( CCMIOError *err, CCMIONode node, CCMIOIndex n,
				float const *data, CCMIOIndex start,
				CCMIOIndex end );
extern CCMIOError CCMIOWrite1d( CCMIOError *err, CCMIONode node, CCMIOIndex n,
				double const *data, CCMIOIndex start,
				CCMIOIndex end );
extern CCMIOError CCMIOWrite2i( CCMIOError *err, CCMIONode node,
				CCMIOIndex x, CCMIOIndex y,
				int const *data,
				CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIOWrite2f( CCMIOError *err, CCMIONode node,
				CCMIOIndex x, CCMIOIndex y,
				float const *data,
				CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIOWrite2d( CCMIOError *err, CCMIONode node,
				CCMIOIndex x, CCMIOIndex y,
				double const *data,
				CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIOWrite3i( CCMIOError *err, CCMIONode node,
				CCMIOIndex x, CCMIOIndex y,
				CCMIOIndex z, int const *data,
				CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIOWrite3f( CCMIOError *err, CCMIONode node,
				CCMIOIndex x, CCMIOIndex y,
				CCMIOIndex z, float const *data,
				CCMIOIndex start, CCMIOIndex end );
extern CCMIOError CCMIOWrite3d( CCMIOError *err, CCMIONode node,
				CCMIOIndex x, CCMIOIndex y,
				CCMIOIndex z, double const *data,
				CCMIOIndex start, CCMIOIndex end );

/** \internal */
CCMIOError CCMIOOldReadf( CCMIOError *err, CCMIONode node, int dimension,
			  int swapDims, float *data, CCMIOIndex start,
			  CCMIOIndex end );
/** \internal */
CCMIOError CCMIOOldReadd( CCMIOError *err, CCMIONode node, int dimension,
			  int swapDims, double *data, CCMIOIndex start,
			  CCMIOIndex end );
/** \internal */
CCMIOError CCMIOOldReadi( CCMIOError *err, CCMIONode node, int dimension,
			  int swapDims, int *data, CCMIOIndex start,
			  CCMIOIndex end );

/** Recursively copies a node.  origNode and copyNode need not be from
    the same file.
    \param origNode		The node to be copied.
    \param copyNode		If copyExists is TRUE, then copyNode is the
    				node to be copied over.  Otherwise, it
				is the parent, and a copy of origNode will be
				created as a child.
    \param copyExists		Determines the function of copyNode. */
extern CCMIOError CCMIOCopyNode( CCMIOError *err, CCMIONode origNode,
				 CCMIONode copyNode, int copyExists );
/** Same as CCMIOGetNextChild() (particularly with respect to the parameter n)
    except that it only returns children with the specified label. */
extern CCMIOError CCMIOGetNextChildWithLabel( CCMIOError *err, CCMIONode parent,
					      char const *label, int *n,
					      CCMIONode *child );

/** Compresses the CCMIO file specified.  This is occasionally necessary becase
    ADF (the underlying storage format) leaks disk space;  ADF does not
    completely recover the space when a node is delete or is rewritten with
    less information.  The compression is performed by copying the data to
    a new ADF file, which is then renamed to the original name.
    This function will correctly compress CCMIO files, and will probably, but
    not necessarily, correctly compress ADF files.  Also note that this function
    requires temporary storage equal to the size of the original file in the
    same filesystem as the original file and should not be called frequently,
    as it may take a some time. */
CCMIOError CCMIOCompress( CCMIOError *err, char *filename );

/* \} */
#ifdef __cplusplus
}
#endif
#endif /* CCMIO_UTILITY_H */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
