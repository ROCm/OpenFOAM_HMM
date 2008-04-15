#ifndef CCMIO_CORE_H
#define CCMIO_CORE_H

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmiocore.h,v $
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
 *  $Id: ccmiocore.h,v 1.8 2006/06/06 12:09:54 jal Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

#define kCCMIOBadNode 0.0

#ifndef MAKEDEPEND
#include <stdarg.h>
#endif

#include "libccmio/ccmiotypes.h"

/* \{
   \name Core functions
   \brief The minimal set of functions.
   \ingroup core */

/** Opens the data file.  The file will be created if it does not already exist.
    \param filename	Name of the file to open.
    \param mode		Mode to open the file:  kCCMIORead or kCCMIOWrite.
    \param root 	Returns the root node of the file. */
extern CCMIOError CCMIOOpen( char const *filename, CCMIOIOType mode,
			     CCMIONode *root );

/** Closes the data file
    \param root		Root node of the file to be closed. */
extern CCMIOError CCMIOClose( CCMIONode root );

/** Finds a node, given a node path
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param parent	The node that the root of the path begins with (need
    			not be the actual root node).
    \param path 	A unix-like path containing the names of nodes.
    			The beginning / is optional.
			(e.g. "/solver/brep/topology/endpoints").
    \param node		Returns the node (invalid if it does not exist). */
extern CCMIOError CCMIOGetNode( CCMIOError *err, CCMIONode parent,
				char const *path, CCMIONode *node );

/** Returns the number of children in the parent node. */
extern CCMIOError CCMIOGetNumberOfChildren( CCMIOError *err, CCMIONode parent,
					    int *n );

/** Returns the next child node.  Note that there is no guarantee that
    the order of the children will always be the same.
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param parent	Parent node of the prospective children.
    \param n		Which child (i.e. return the nth child).  n will be
    			incremented upon successful completion of the function.
			If there is no nth child, a kCCMIONoNodeErr error will be
			returned.
    \param child	The nth child, or NULL if the nth child does not exist.
*/
extern CCMIOError CCMIOGetNextChild( CCMIOError *err, CCMIONode parent, int *n,
				     CCMIONode *child );

/** Returns the name of the node
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param node		The node whose name is to be returned.
    \param name		Returns the name of the node.  Must be at least
    			kCCMIOMaxStringLength + 1 bytes. */
extern CCMIOError CCMIOGetName( CCMIOError *err, CCMIONode node, char *name );

/** Sets the name of an existing node
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param node		The node whose name is to be set.
    \param name		The name of the node.  If the file type (e.g. ADF)
			limits lengths of names, the name will be truncated. */
extern CCMIOError CCMIOSetName( CCMIOError *err, CCMIONode node,
				char const *name );

/** Returns the label of the node
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param node		The node whose label is to be returned.
    \param label	Returns the label of the node.  This string must be
    			at least kCCMIOMaxStringLength + 1 bytes. */
extern CCMIOError CCMIOGetLabel( CCMIOError *err, CCMIONode node, char *label );

/** Sets the label of an existing node. */
extern CCMIOError CCMIOSetLabel( CCMIOError *err, CCMIONode node,
				 char const *label );

/** Creates a new node
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param parent	The parent node of the new node.
    \param openDup	If TRUE, a \ref kCCMIODuplicateNodeErr will not be returned
    			in an attempt to create a sibling with the same name;
			instead, CCMIOGetNode() will be called on the original node
			If FALSE, the error will be returned and the value of
			node is undefined.
    \param name		The name of the node.  If the file type (e.g. ADF)
    			limits lengths of names, the name will be truncated.
    \param label	The label of the name.  This may be trunctated like
    			the name.
    \param node		Returns the new node.  NULL may be passed in if no
			return is required. */
extern CCMIOError CCMIOCreateNode( CCMIOError *err, CCMIONode parent,
				   int openDup, char const *name,
				   char const *label, CCMIONode *node );

/** Creates a new link
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param parent	The parent node of the new link.
    \param name		The name of the node.  See CCMIOCreateNode() for more
    			details on this parameter.
    \param filename	Filename that the link points to.
    \param destName	Path to the node to be linked to.
    \param node		Returns the new node (invalid in case of an error).
			NULL may be passed in if no return is required. */
extern CCMIOError CCMIOCreateLink( CCMIOError *err, CCMIONode parent,
				   char const *name, char const *filename,
				   char const *destName, CCMIONode *node );

/** Deletes the node */
extern CCMIOError CCMIODeleteNode( CCMIOError *err, CCMIONode node );

/** Deletes all children of the node, but not the node itself */
extern CCMIOError CCMIODeleteAllChildren( CCMIOError *err, CCMIONode node );

/** Moves the node underneath newParent */
extern CCMIOError CCMIOMoveNode( CCMIOError *err, CCMIONode node,
				 CCMIONode newParent );

/** Returns the number of items in the node's data
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param node		The node to use
    \param nDims	Returns the number of dimensions.  0 means no data.
    \param dims		Returns the number of elements of each dimension.  Note
    			that this array is allocated by the library and must be
			freed by the user. */
extern CCMIOError CCMIOGetDimensions( CCMIOError *err, CCMIONode node,
				      int *nDims, CCMIOIndex **dims );

/** Returns the number of bytes of the data */
extern CCMIOError CCMIOGetDataSize( CCMIOError *err, CCMIONode node,
				    CCMIOSize *bytes );

/** Returns the data type for the node */
extern CCMIOError CCMIOGetDataType( CCMIOError *err, CCMIONode node,
				    CCMIODataType *type);

/** Sets the datatype of the node.
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param node		The node
    \param type		The datatype
    \param ...		The number and size of the dimensions (passed as
    			CCMIOIndex), terminated with a \ref kCCMIOEndArgs.
			For instance an array of "int ary[5][10];"
			would be "CCMIOReadDataPoint(node, data, 5ul, 10ul, kCCMIOEndArgs)".
			A single dimension of size 1 is equivalent to a scalar
			value. */
extern CCMIOError CCMIOSetDataType( CCMIOError *err, CCMIONode node,
				    CCMIODataType type, ... );
extern CCMIOError CCMIOvSetDataType( CCMIOError *err, CCMIONode node,
				     CCMIODataType type, va_list args );
extern CCMIOError CCMIOSetDataTypev( CCMIOError *err, CCMIONode node,
				     CCMIODataType type,
				     int nDims, CCMIOIndex const *dims );

/** Returns all the data in the node
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param node		The node to use
    \param data		Pointer to memory large enough to hold the entire
    			data;  note that this is \em not allocated by the
			library.  The data is returned as FORTRAN data.
    \param expected	The expected type.  A \ref kCCMIOWrongDataTypeErr error
    			is returned if the node's data is not of the expected
			type and no data is written.
    \param dimsExpected	The number of dimensions that are expected.  A
    			kCCMIOWrongDataType error is returned if the node's data
			has a different number of dimensions.  A scalar and
			a one-dimensional array both have dimensions of 1. 
extern CCMIOError CCMIOReadData( CCMIOError *err, CCMIONode node, void *data,
				 CCMIODataType expected, int dimsExpected );
*/
/** Returns one element of data.  No datatype checking is performed.
    \param err		Error input and output.  If the incoming error is
    			not \ref kCCMIONoErr, the function will do nothing.
			Is set to the function's return value on exit.
    \param node		The node whose data to return.
    \param data		Pointer to memory large enough to hold one element
    \param ...		Var-arg list that is the C-style index into the
    			array.  Scalar values can be read by using a single
			value of 0, although CCMIOReadData() is preferred.
			At least one parameter must be passed to
			the list.  Note that the list must be terminated
			with \ref kCCMIOEndArgs. 
extern CCMIOError CCMIOReadDataPoint( CCMIOError *err, CCMIONode node,
				      void *data, ...);
*/
/** Write one element of data.  No datatype checking is performed.
    See CCMIOReadDataPoint() for a discussion of the parameters.  The var-arg
    list must be terminated with \ref kCCMIOEndArgs.
    Note: for performance reasons, the datatype of 'data' is unknown to
      the function, so CCMIOWriteDataPoint() cannot be used to write strings, as it
      does not know when 'data' is a string.  For string data use CCMIOWriteData.
    Note: When writing scalar data with CCMIOWriteDataPoint(), a dimension of 0
      should be used:  "CCMIOWriteDataPoint(err, node, data, 0, kCCMIOEndArgs)",
      although CCMIOWriteData() is preferred.
extern CCMIOError CCMIOWriteDataPoint( CCMIOError *err, CCMIONode node,
				       void *data, ... );
*/
/** Writes the entire node's data.  Requires that the node's datatype and size
    were correctly set with CCMIOSetDataType().
CCMIOError CCMIOWriteData( CCMIOError *err, CCMIONode node, void const *data );
*/
/** Returns TRUE if node1 == node2, FALSE otherwise. */
int CCMIOAreNodesEqual( CCMIONode node1, CCMIONode node2 );

/** Returns the root node, or \ref kCCMIOBadParameterErr if 'node' is not valid.*/
CCMIOError CCMIOGetRootNode( CCMIOError *err, CCMIONode node, CCMIONode *root );

/*\}*/
#ifdef __cplusplus
    }
#endif
#endif /* CCMIO_CORE_H */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
