#ifndef CCMIO_H
#define CCMIO_H

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmio.h,v $
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
 *  $Id: ccmio.h,v 1.18 2006/06/05 21:12:16 geoffp Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "libccmio/ccmiotypes.h"

/* \{ \name Basic functions
      \ingroup intermediate */

/** Opens an CCMIO file.  Calling CCMIOOpenFile() multiple times on a file does
    not create multiple references to it;  thus if CCMIOOpenFile() is called
    five times on the same file, only one CCMIOCloseFile() is required. */
extern CCMIOError CCMIOOpenFile( CCMIOError *err, char const *file,
				 CCMIOIOType mode, CCMIOID *root );

/** Closes all references to the CCMIO file.  Closing an already closed file
    generates a kCCMIOBadParameterErr. */
extern CCMIOError CCMIOCloseFile( CCMIOError *err, CCMIOID root );

/** Returns the version number of the file. */
extern CCMIOError CCMIOGetVersion( CCMIOError *err, CCMIONode root,
				   int *version );

/** Sets the version number of the file. */
extern CCMIOError CCMIOSetVersion( CCMIOError *err, CCMIONode root,
				   int version );

/** Returns the title of the file.  Note that the returned string is allocated
    by the library and must be freed by the user.  It is null-terminated. */
extern CCMIOError CCMIOGetTitle( CCMIOError *err, CCMIONode root, char **title );

/** Sets the title of the file.  The string must be null-terminated. */
extern CCMIOError CCMIOSetTitle( CCMIOError *err, CCMIONode root,
				 char const *title );

/* \} */
/* \{ \name Functions for Optional Nodes
      \ingroup intermediate */

/** Creates a child node of parent with given name and scalar value.
    If the node already exists it will be overwritten.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param value	Scalar value to be written. */
extern CCMIOError CCMIOWriteOpti( CCMIOError *err, CCMIOID parent,
				  char const *name, int value );

/** Creates a child node of parent with given name and scalar value.
    If the node already exists it will be overwritten.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param value	Scalar value to be written. */
extern CCMIOError CCMIOWriteOptf( CCMIOError *err, CCMIOID parent,
				  char const *name, float value );

/** Creates a child node of parent with given name and scalar value.
    If the node already exists it will be overwritten.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param value	Scalar value to be written. */
extern CCMIOError CCMIOWriteOptd( CCMIOError *err, CCMIOID parent,
				  char const *name, double value );

/** Creates a child node of parent with given name and string value.
    If the node already exists it will be overwritten.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param value	Null-terminated string. */
extern CCMIOError CCMIOWriteOptstr( CCMIOError *err, CCMIOID parent,
				    char const *name, char const *value );

/** Reads the scalar value of a child node;  converts to an integer
    if the type of the node is not \def kCCMIOInt32.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param value	Pointer to the scalar data to be read in. */
extern CCMIOError CCMIOReadOpti( CCMIOError *err, CCMIOID parent,
				 char const *name, int *value );

/** Reads the scalar value of a child node;  converts to an float
    if the type of the node is not \def kCCMIOFloat32.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param value	Pointer to the scalar data to be read in. */
extern CCMIOError CCMIOReadOptf( CCMIOError *err, CCMIOID parent,
				 char const *name, float *value );

/** Reads the scalar value of a child node;  converts to an double
    if the type of the node is not \def kCCMIOFloat64.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param value	Pointer to the scalar data to be read in. */
extern CCMIOError CCMIOReadOptd( CCMIOError *err, CCMIOID parent,
				 char const *name, double *value );

/** Reads the string value of a child node;  returns \def kCCMIOWrongDataTypeErr
    if the type of the node is not \def kCCMIOString.  'value' must previously
    be allocated to an appropriate size.  (This can be determined by calling
    the function once and passing NULL for 'value'.  'size' does not include
    the C string terminator.).  Either 'size' or 'value' can be NULL. */
extern CCMIOError CCMIOReadOptstr( CCMIOError *err, CCMIOID parent,
				   char const *name, int *size, char *value );

/** Reads a one-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  See CCMIOReadOpt1f() for a complete
    description. */
extern CCMIOError CCMIOReadOpt1i( CCMIOError *err, CCMIOID parent,
				  char const *name, int *data,
				  CCMIOIndex start, CCMIOIndex end );
/** Reads a one-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  'data' must be the proper size
    (which can be determined with CCMIOGetDimensions()).  If actual node data is
    stored in a different format than requested (i.e. float or double), it will
    be converted.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param data		Data to be read into.  Must be non-NULL. */
extern CCMIOError CCMIOReadOpt1f( CCMIOError *err, CCMIOID parent,
				  char const *name, float *data,
				  CCMIOIndex start, CCMIOIndex end );
/** Reads a one-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  See CCMIOReadOpt1f() for a complete
    description. */
extern CCMIOError CCMIOReadOpt1d( CCMIOError *err, CCMIOID parent,
				  char const *name, double *data,
				  CCMIOIndex start, CCMIOIndex end );
/** Reads a two-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  See CCMIOReadOpt2f() for a complete
    description. */
extern CCMIOError CCMIOReadOpt2i( CCMIOError *err, CCMIOID parent,
				  char const *name, int *data,
				  CCMIOIndex start, CCMIOIndex end );
/** Reads a two-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  'data' must be the proper size
    (which can be determined with CCMIOGetDimensions()).  If actual node data is
    stored in a different format than requested (i.e. float or double), it will
    be converted.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param data		Data to be read into.  Must be non-NULL.
    \param start	The offset of the element pointed to by 'data'.  Must
			Be in units of of the first element;  so if
			  float data[10][3];
			then start can be in the range [0, 9].
    \param end		The element one beyond the end.  So to read from [2, 5],
			start = 2, end = 6, which reads from data[2][0]
			to data[5][3]. */
extern CCMIOError CCMIOReadOpt2f( CCMIOError *err, CCMIOID parent,
				  char const *name, float *data,
				  CCMIOIndex start, CCMIOIndex end );
/** Reads a two-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  See CCMIOReadOpt2f() for a complete
    description. */
extern CCMIOError CCMIOReadOpt2d( CCMIOError *err, CCMIOID parent,
				  char const *name, double *data,
				  CCMIOIndex start, CCMIOIndex end );
/** Reads a three-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  See CCMIOReadOpt3f() for a complete
    description. */
extern CCMIOError CCMIOReadOpt3i( CCMIOError *err, CCMIOID parent,
				  char const *name, int *data,
				  CCMIOIndex start, CCMIOIndex end );
/** Reads a three-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  'data' must be the proper size
    (which can be determined with CCMIOGetDimensions()).  If actual node data is
    stored in a different format than requested (i.e. float or double), it will
    be converted.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param data		Data to be read into.  Must be non-NULL.
    \param start	The offset of the element pointed to by 'data'.  Must
			Be in units of of the first element;  so if
			  float data[10][3][3];
			then start can be in the range [0, 9].
    \param end		The element one beyond the end.  So to read from [2, 5],
			start = 2, end = 6, which reads from data[2][0][0]
			to data[5][3][3]. */
extern CCMIOError CCMIOReadOpt3f( CCMIOError *err, CCMIOID parent,
				  char const *name, float *data,
				  CCMIOIndex start, CCMIOIndex end );
/** Reads a three-dimensional array from the child node of 'parent'
    with name 'name' into 'data'.  See CCMIOReadOpt3f() for a complete
    description. */
extern CCMIOError CCMIOReadOpt3d( CCMIOError *err, CCMIOID parent,
				  char const *name, double *data,
				  CCMIOIndex start, CCMIOIndex end );

/** Writes a one-dimensional array to a child node of 'parent' with name 'name'.
    See CCMIOWriteOpt1f() for a complete description. */
extern CCMIOError CCMIOWriteOpt1i( CCMIOError *err, CCMIOID const parent,
				   char const *name, CCMIOSize const n,
                                   int const *data,
				   CCMIOIndex start, CCMIOIndex end );
/** Writes a one-dimensional array to a child node of 'parent' with name 'name'.
    If the child node does not exist it will be created and in either case its
    size and data type will be automatically set to the proper values.
    \param parent	The parent entity.
    \param name		Name of the child node.
    \param n		The total number of elements to be written.
			If not writing the entire array this is the number
			of elements that will be ultmately be written.
    \param data		The data to be written.
    \param start	When writing partial arrays, where to start writing.
    \param end		When writing partial arrays, where to stop writing.
			(The write will be from data[start] to data[start-1])
   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to kCCMIOStart. */
extern CCMIOError CCMIOWriteOpt1f( CCMIOError *err, CCMIOID parent,
				   char const *name, CCMIOSize n,
				   float const *data,
				   CCMIOIndex start, CCMIOIndex end );
/** Writes a one-dimensional array to a child node of 'parent' with name 'name'.
    See CCMIOWriteOpt1f() for a complete description. */
extern CCMIOError CCMIOWriteOpt1d( CCMIOError *err, CCMIOID parent,
				   char const *name, CCMIOSize n,
				   double const *data,
				   CCMIOIndex start, CCMIOIndex end );
/** Writes a two-dimensional array to a child node of 'parent' with name 'name'.
    See CCMIOWriteOpt2f() for a complete description. */
extern CCMIOError CCMIOWriteOpt2i( CCMIOError *err, CCMIOID const parent,
				   char const *name, CCMIOIndex x, CCMIOIndex y,
				   int const *data,
				   CCMIOIndex start, CCMIOIndex end );
/** Writes a two-dimensional array to a child node of 'parent' with name 'name'.
    If the child node does not exist it will be created and in either case its
    size and data type will be automatically set to the proper values.
    'data' is expected to be an array of data[x][y] elements.  Note that
    although the CCM specification specifies dimension in FORTRAN ordering,
    the API expects C ordering.  Thus a C array data[x][y] must be called
    as CCMIOWriteOpt2f(..., x, y, ...), even though the dimensions will be
    reversed on disk.
    \param parent	The parent node.
    \param name		Name of the child entity.
    \param x		The size of the first C dimension.
    \param y		The size of the second C dimension.
    \param data		The data to be written.
    \param start	The offset of the element pointed to by 'data'.  Must
			Be in units of of the first element;  so if
			  float data[10][3];
			then start can be in the range [0, 9].
    \param end		The element one beyond the end. So to write from [2, 5],
			start = 2, end = 6, which writes from data[2][0]
			to data[5][3].
   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to kCCMIOStart. */
extern CCMIOError CCMIOWriteOpt2f( CCMIOError *err, CCMIOID parent,
				   char const *name, CCMIOIndex x, CCMIOIndex y,
				   float const *data,
				   CCMIOIndex start, CCMIOIndex end );
/** Writes a two-dimensional array to a child node of 'parent' with name 'name'.
    See CCMIOWriteOpt2f() for a complete description. */
extern CCMIOError CCMIOWriteOpt2d( CCMIOError *err, CCMIOID parent,
				   char const *name, CCMIOIndex x, CCMIOIndex y,
				   double const *data,
				   CCMIOIndex start, CCMIOIndex end );
/** Writes a three-dimensional array to a child node of 'parent' with name
    'name'.  See CCMIOWriteOpt3f() for a complete description. */
extern CCMIOError CCMIOWriteOpt3i( CCMIOError *err, CCMIOID parent,
				   char const *name,
				   CCMIOIndex x, CCMIOIndex y, CCMIOIndex z,
				   int const *data,
				   CCMIOIndex start, CCMIOIndex end );
/** Writes a three-dimensional array to a child node of 'parent' with name
    'name'. If the child node does not exist it will be created and in either
    case its size and data type will be automatically set to the proper values.
    'data' is expected to be an array of data[x][y] elements.  Note that
    although the CCM specification specifies dimension in FORTRAN ordering,
    the API expects C ordering.  Thus a C array data[x][y] must be called
    as CCMIOWriteOpt3f(..., x, y, z, ...), even though the dimensions will be
    reversed on disk.
    \param parent	The parent node.
    \param name		Name of the child entity.
    \param isC		TRUE if 'data' is in C array order,
			FALSE if in Fortran array order.
    \param x		The size of the first dimension.
    \param y		The size of the second dimension.
    \param z		The size of the third dimension.
    \param data		The data to be written.
    \param start	The offset of the element pointed to by 'data'.  Must
			be in units of of the first element;  so if
			  float data[10][3][3];
			then start can be in the range [0, 9].
    \param end		The element one beyond the end. So to write from [2, 5],
			start = 2, end = 6, which writes from data[2][0][0]
			to data[5][3][3].
   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to kCCMIOStart. */
extern CCMIOError CCMIOWriteOpt3f( CCMIOError *err, CCMIOID parent,
				   char const *name,
				   CCMIOIndex x, CCMIOIndex y, CCMIOIndex z,
				   float const *data, CCMIOIndex start,
				   CCMIOIndex end );
/** Writes a three-dimensional array to a child node of 'parent' with name
    'name'.  See CCMIOWriteOpt3f() for a complete description. */
extern CCMIOError CCMIOWriteOpt3d( CCMIOError *err, CCMIOID parent,
				   char const *name,
				   CCMIOIndex x, CCMIOIndex y, CCMIOIndex z,
				   double const *data, CCMIOIndex start,
				   CCMIOIndex end );

/** Returns information about the optional node.
    \param parant	The parent node.
    \param name		Name of the child entity.
    \param type		Returns the data type of the child.  Can be NULL if
			the information is not desired.
    \param x		Returns the size of the first dimension (1 if the
			node is a scalar).  Can be NULL.
    \param y		Returns the size of the third dimension (0 if the
			node is one dimensional).  Can be NULL.
    \param z		Returns the size of the first dimension (0 if the
			node is less than three dimensions).  Can be NULL. */
extern CCMIOError CCMIOGetOptInfo( CCMIOError *err, CCMIOID parent,
				   char const *name, CCMIODataType *type,
				   CCMIOIndex *x, CCMIOIndex *y,
				   CCMIOIndex *z );

/* \} */
/* \{ \name Intermediate API
      \ingroup intermediate */

/** Marks the entity as invalid.  Further CCMIO calls will fail, until
    a call (such as CCMIOGetEntity()) writes a valid entity into it. */
extern void CCMIOInvalidateEntity( CCMIOID *entity );

/** Returns TRUE if the entity is valid, FALSE otherwise.  Note that
    unitialized entities will likely return TRUE, even though they are not
    actually valid.  Use CCMIOInvalidateEntity() to initialize them.  */
extern int CCMIOIsValidEntity( CCMIOID entity );

/** Returns TRUE if the two entities are from the same file,
    FALSE otherwise. */
extern int CCMIOIsFromSameFile( CCMIOID entity1, CCMIOID entity2 );

/** Creates a new entity.
    \param parent	For top-level entities this may be any CCMIOID in the same
   			file, but for child entities this must be the actual
			parent.
    \param type		The type of entity to be created.
    \param description	An optional description of the entity.  Passing NULL
			will omit a description.
    \param id		The CCMIOID of the created entity. */
extern CCMIOError CCMIONewEntity( CCMIOError *err, CCMIOID parent,
				  CCMIOEntity type, char const *description,
				  CCMIOID *id );

/** Retrieves an entity.  Although it is mostly an internal function, it is
    necessary for retrieving the cell and internal face entities.
    \param parent	For top-level entities this may be any CCMIOID in the
			same file, but for child entities this must be the
			actual parent.
    \param type		The type of entity to be retrieved.
    \param idVal	For entites that may have siblings, this specifies
			the ID number of the entity (this is most useful
			internally).  For entities that may not have siblings
			(e.g. cells and internal faces) this value is ignored.
    \param id		The CCMIOID of the created entity. */
extern CCMIOError CCMIOGetEntity( CCMIOError *err, CCMIOID parent,
				  CCMIOEntity type, int idVal, CCMIOID *id );

/** Creates a new entity identified by an index number.
    Although almost all entities are internally identified by an ID number, this
    number is unimportant for most entities and is only used internally.  Some
    entities, however, are most easily identified with a number and for these
    entities (kCCMIOBoundaryFaces, kCCMIOCellType, and kCCMIOBoundaryRegion) the
    ID number is provided by the creator.  The entities, called indexed
    entities, can be created with this function.
    \param parent	The CCMIOID of the parent entity.
    \param idVal	The boundary ID number.
    \param which	May be kCCMIOBoundaryFaces, kCCMIOCellType, or
    			kCCMIOBoundaryRegion.
    \param description	An optional description (NULL for no description)
    \param id		The CCMIOID to the new entity. */
extern CCMIOError CCMIONewIndexedEntity( CCMIOError *err, CCMIOID parent,
				     CCMIOEntity which, int idVal,
				     char const *description, CCMIOID *id );

/** Returns the index (i.e. the id) of the given indexed entity.  (See
    CCMIONewIndexedEntity() for more details and a list of valid entities.) */
extern CCMIOError CCMIOGetEntityIndex( CCMIOError *err, CCMIOID id , int *n );

/** Creates a new state.  If a state with that name exists, its contents,
    including processor entities, will be cleared.
    \param root		The root node of the file (returned by CCMIOOpenFile()).
    \param name		The name of the new state (limited to
   			kCCMIOMaxStringLength characters, not including a
			terminating NULL).
    \param problemDescription	The CCMIOID of the problemDescription entity
			relevant to this state, or NULL if there is none.
    \param state	The CCMIOID of the new entity.
    \param description	An optional description string
    			(NULL for no description) */
extern CCMIOError CCMIONewState( CCMIOError *err, CCMIOID root,
				 char const *name, CCMIOID *problemDescription,
				 char const *description, CCMIOID *state );

/** Retrieves a state.
    \param root		The root node of the file (returned by CCMIOOpenFile()).
    \param name		The name of the state (limited to
   			\def kCCMIOMaxStringLength characters, not including a
			terminating NULL).
    \param problemDescription	Returns the CCMIOID of the problemDescription
   			entity relevant to this state, or an invalid entity
			(check with CCMIOIsValidEntity()) if there is none.
			If NULL is passed to this parameter, no
			problemDescription will be looked for.
    \param state	The CCMIOID of the entity. */
extern CCMIOError CCMIOGetState( CCMIOError *err, CCMIOID root,
				 char const *name, CCMIOID *problemDescription,
				 CCMIOID *state );

/** Writes the problem description node to the state.  Note that there is
    not an CCMIOReadState();  this functionality is incorporated in
    CCMIOGetState().
    \param description	Optional character string describing the state.
			Pass NULL for no description (note that this will
			remove any existing description). */
extern CCMIOError CCMIOWriteState( CCMIOError *err, CCMIOID state,
				   CCMIOID problemDescription,
				   char const *description );

/** Creates a new field entity.
    \param phase	The FieldPhase entity.
    \param name		The name of the new field (limited to
   			\def kCCMIOMaxStringLength characters, not including a
			terminating NULL).
    \param shortName	The short name that prostar will use.  Must be
			a unique string among all short names in the fieldSet.
			Only the first kCCMIOProstarShortNameLength characters
			will be written.
    \param dim		The number of dimensions (\def kCCMIOScalar,
    			\def kCCMIOVector, \def kCCMIOTensor).
    \param field	The CCMIOID of the new entity. */
extern CCMIOError CCMIONewField( CCMIOError *err, CCMIOID phase,
				 char const *name, char const *shortName,
				 CCMIODimensionality dim, CCMIOID *field );

/** Retrieves a field entity.
    \param phase	The FieldPhase entity.
    \param name		The name of the field (limited to
   			kCCMIOMaxStringLength characters, not including a
			terminating NULL).
    \param dim		Returns the number of dimensions
   			kCCMIOTensor).
    \param field	The CCMIOID of the field entity. */
extern CCMIOError CCMIOGetField( CCMIOError *err, CCMIOID phase,
				 char const *name, CCMIODimensionality *dim,
				 CCMIOID *field );

/** Reads information about a field.  Note that this does not read the \em
    data of the node;  use CCMIOReadFieldData(f|v|i)() for that.  This is for
    determining the name of the field when using CCMIONextEntity().
    \param field	The CCMIOID of the field.
    \param name		Returns the name of the field.  Must be allocated to at
    			least kCCMIOMaxStringLength + 1 characters.
    \param shortName	Returns the short name of the field.  Must be allocated
			to at least kCCMIOProstarShortName + 1 characters.
    \param dim		Returns the dimensions of the field.
    \param datatype	Returns the type of data.  Only valid for scalar data;
			Vector and tensor components are not required to be the
			same data type (\ref kCCMIOUnknownType will be returned
			for non-scalar data).  May be NULL if the
			information is not desired. */
extern CCMIOError CCMIOReadField( CCMIOError *err, CCMIOID field, char *name,
				  char *shortName, CCMIODimensionality *dim,
				  CCMIODataType *datatype );

/** Creates a new prostar set.  The sets can be written with CCMIOWriteOpti().
    If called on a pre-existing set, all contents of the original set will
    be removed.
    \param name		The name of the set.  This is limited to
			\ref kCCMIOMaxStringLength.
    \param longName	The name of the set;  length is not limited.
    			If no long name is desired, NULL may be passed here,
			and the short name will be used.
    \param setID	The CCMIOID of the new prostar set.  Can be NULL if
			not desired. */
CCMIOError CCMIONewProstarSet( CCMIOError *err, CCMIOID root, const char *name,
			       const char *longName, CCMIOID *setID );

/** Gets a prostar set.
    \param name		The name of the set.
    \param longNameSize	The size of the long name, not including the C string
			terminator.  Can be NULL if the size is not desired.
    \param longName	The long name.  Can be NULL if the information is not
			desired.
    \param setsAvailableFlag	Lists the available sets.  The sets can be read
			using CCMIOReadOpti().  Can be NULL if not desired.
    \parma setID	The CCMIOID of the prostar set.  Can be NULL if not
			desired. */
CCMIOError CCMIOGetProstarSet( CCMIOError *err, CCMIOID root, const char *name,
			       int *longNameSize, char *longName,
			       unsigned int *setsAvailableFlag,
			       CCMIOID *setID );

/** Deletes the entity and all of its children.  Note that deleting
    a processor without calling CCMIOClearProcessor() first will not delete
    entites that become unused as a result of the delete (e.g. the mesh and
    post data referred to by that processor) and may lead to wasted disk space.
    However, calling CCMIODeleteEntity on a state entity is safe. */
extern CCMIOError CCMIODeleteEntity( CCMIOError *err, CCMIOID id );

/** Returns the next child entity of the parent or \def kCCMIONoNodeError
    if there are no more.
    \param parent	The CCMIOID of the parent node.  When iterating over
    			top-level entites (like a Map), any CCMIOID from the
			same file may be used.
    \param type		The type of entity to iterator over.
    \param i		The counter value.  This should always be set to 0
    			(zero) to retrieve the first child.  It will be
			automatically incremented.
    \param next		Returns the CCMIOID of the next child entity. */
extern CCMIOError CCMIONextEntity( CCMIOError *err, CCMIOID parent,
				   CCMIOEntity type, int *i, CCMIOID *next );

/** Returns the number of elements and the maximum element ID in the specified
    entity.  Either 'n' or 'max' may be NULL if the information is not
    desired. */
extern CCMIOError CCMIOEntitySize( CCMIOError *err, CCMIOID id, CCMIOSize *n,
				   CCMIOIndex *max );

/** Gets the name of the entity.  Although this is not useful for most entities,
    it may be convenient for named entities (kCCMIOState, kCCMIOField).
    \param name		Returns the name of the entity.  Must be at least
			kCCMIOMaxStringLength + 1. */
extern CCMIOError CCMIOEntityName( CCMIOError *err, CCMIOID id, char *name );

/** Gets the label of the entity.
    \deprecated Use CCMIOEntityDescription() instead.
    \param size		Returns the length of the label (not including a
    			terminating NULL).  May be NULL if the information
			is not required.
    \param name		Returns the label.  Must be at least 'size' + 1 bytes.
			If there is no label, this will be the empty string.
			May be NULL if the information is not required. */
extern CCMIOError CCMIOEntityLabel( CCMIOError *err, CCMIOID id,
				    CCMIOSize *size, char *label );

/** Gets the description (if any) of the specified entity.
    \param size		Returns the size of the string.  May be NULL if the
			information is not desired.
    \param desc		Returns the description string.  This must be allocated
			to at least as large as the number of characters + 1
			and may be NULL if the information is not desired. */
extern CCMIOError CCMIOEntityDescription( CCMIOError *err, CCMIOID id,
					  CCMIOSize *size, char *desc );

/** Returns the CCMIO node corresponding to this entity. */
extern CCMIOError CCMIOGetEntityNode( CCMIOError *err, CCMIOID id,
				      CCMIONode *node );

/** Returns the data type of the main data of the entity.  Currently only
    valid for vertices. */
extern CCMIOError CCMIOEntityDataType( CCMIOError *err, CCMIOID id,
				       CCMIODataType *type );

/** Writes the map data.
    \param id		The CCMIOID of the map.
    \param n		The total number of elements to be written.
    \param max		The maximum element value.
    \param data		The map data.
    \param start        The offset of the starting data.
    \param end          The offset of the ending data.  Data from start to
    			end - 1 will be written. If end is kCCMIOEnd,
			data will be writen from 'start' to the end of the data.
  
   To write the whole array, set start to kCCMIOStart and end to kCCMIOEnd.
  
   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to kCCMIOStart. */
extern CCMIOError CCMIOWriteMap( CCMIOError *err, CCMIOID id,
                                 CCMIOSize n, CCMIOSize max,
                                 int const *data,
				 CCMIOIndex start, CCMIOIndex end );

/** Reads the map data.
    \param id		The CCMIOID of the map.
    \param data		The map data.  Must be pre-allocated to the correct
			number of elements (given by CCMIOEntitySize()). */
extern CCMIOError CCMIOReadMap( CCMIOError *err, CCMIOID id, int *data,
				CCMIOIndex start, CCMIOIndex end );

/** Reads the vertex data.
    \param id		The CCMIOID of the vertices.
    \param dims		Returns the dimensionality of the vertex (i.e. 2 or 3).
    			May be NULL if the information is not desired.
    \param scale	Returns the scaling factor.
    			May be NULL if the information is not desired.
    \param mapID	Returns the CCMIOID of the map corresponding to this
			entity.  May be NULL if the information is not desired.
    \param vertices	A two dimensional array of size [nVerts][dims].  Must be
			pre-allocated to the correct size.  May be NULL if
			the information is not desired.
    \param start	The offset (in units of vertices) of the starting
    			vertex.
    \param end		The offset (in units of vertices) of the ending vertex.
    			Data from start to end - 1 will be read. If end is
			kCCMIOEnd, data will be read from 'start' to the
			end of the data.*/
extern CCMIOError CCMIOReadVerticesf( CCMIOError *err, CCMIOID id,
				      int *dims, float *scale,
				      CCMIOID *mapID, float  *vertices,
				      CCMIOIndex start, CCMIOIndex end );
/** Reads the vertex data.
    \param id		The CCMIOID of the vertices.
    \param dims		Returns the dimensionality of the vertex (i.e. 2 or 3).
    			May be NULL if the information is not desired.
    \param scale	Returns the scaling factor.
    			May be NULL if the information is not desired.
    \param mapID	Returns the CCMIOID of the map corresponding to this
			entity.  May be NULL if the information is not desired.
    \param vertices	A two dimensional array of size [nVerts][dims].  Must be
			pre-allocated to the correct size.  May be NULL if
			the information is not desired.
    \param start	The offset (in units of vertices) of the starting
			vertex.
    \param end		The offset (in units of vertices) of the ending vertex.
    			Data from start to end - 1 will be read. If end is
			kCCMIOEnd, data will be read from 'start' to the
			end of the data.*/
extern CCMIOError CCMIOReadVerticesd( CCMIOError *err, CCMIOID id,
				      int *dims, float *scale,
				      CCMIOID *mapID, double *vertices,
				      CCMIOIndex start, CCMIOIndex end );
/** Writes the vertex data.
    \param id		The CCMIOID of the vertices.
    \param dims		Returns the dimensionality of the vertex (i.e. 2 or 3).
    \param scale	Returns the scaling factor.
    \param mapID	Returns the CCMIOID of the map corresponding to this
    			entity.
    \param vertices	A two dimensional array of size [nVerts][dims].
    \param start	The offset (in units of vertices) of the starting
  		        vertex.
    \param end		The offset (in units of vertices) of the ending
                        vertex.  Data from start to end - 1 will be read.  If
                        end is kCCMIOEnd, data will be written from 'start'
			to the end of the data.

   To write the whole array, set start to kCCMIOStart and end to kCCMIOEnd.

   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to zero. */
extern CCMIOError CCMIOWriteVerticesf( CCMIOError *err, CCMIOID id,
				       int dims, float scale,
				       CCMIOID mapID, float const *vertices,
				       CCMIOIndex start, CCMIOIndex end );
/** Writes the vertex data.
    \param id		The CCMIOID of the vertices.
    \param dims		Returns the dimensionality of the vertex (i.e. 2 or 3).
    \param scale	Returns the scaling factor.
    \param mapID	Returns the CCMIOID of the map corresponding to this
    			entity.
    \param vertices	A two dimensional array of size [nVerts][dims].
    \param start	The offset (in units of vertices) of the starting
  		        vertex.
    \param end		The offset (in units of vertices) of the ending
                        vertex.  Data from start to end - 1 will be read.  If
                        end is kCCMIOEnd, data will be written from 'start'
			to the end of the data.
  
   To write the whole array, set start to kCCMIOStart and end to kCCMIOEnd.
  
   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to zero. */
extern CCMIOError CCMIOWriteVerticesd( CCMIOError *err, CCMIOID id,
				       int dims, float scale,
				       CCMIOID mapID, double const *vertices,
				       CCMIOIndex start, CCMIOIndex end );

/** Reads the cell data.
    \param id		The CCMIOID of the cell entity.
    \param mapID	Returns the CCMIOID of the map corresponding to this
			entity.  May be NULL if the information is not desired.
    \param cellTypes	Returns an array of cell type information.  May be
			NULL if the information is not desired.
    \param start	The beginning cell.
    \param end		The ending cell.  Data from 'start' to 'end' -1 will
    			be read.  If end is zero, data will be read from 'start'
			to the end of the data. */
extern CCMIOError CCMIOReadCells( CCMIOError *err, CCMIOID id, CCMIOID *mapID,
				  int *cellTypes, CCMIOIndex start,
				  CCMIOIndex end );

/** Writes the cell data
    \param id		The CCMIOID of the cell entity.
    \param mapID	The CCMIOID of the map corresponding to this entity.
    \param cellTypes	Returns an array of cell type information.
    \param start	The beginning cell.
    \param end		The ending cell.  Data from 'start' to 'end' -1 will
    			be read.  If end is zero, data will be written from
			'start' to the end of the data.
   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to kCCMIOStart. */
extern CCMIOError CCMIOWriteCells( CCMIOError *err, CCMIOID id, CCMIOID mapID,
				   int *cellTypes, CCMIOIndex start,
                                   CCMIOIndex end );

/** Reads the face data.
    \param entity	The face entity.
    \param which	Either kCCMIOInternalFaces or kCCMIOBoundaryFaces.
    \param mapID	Returns the CCMIOID of the map associated with these
    			faces.  May be NULL if the information is not desired.
    \param streamSize	Returns the total number of elements in the
  		        vertexStream.  May be NULL if the information is not
  		        desired.
    \param vertexStream	Returns the vertices in the faces in the form<br>
  		        <tt>nVerts v1 v2 ... vn<br>
  		            nVerts v1 v2 ... vn<br>
  			    &nbsp;&nbsp; ...</tt><br>
  			Note that if the stream is being buffered (with
  			'start' and 'end'), only a piece of the stream will
  			be returned, and may be in the middle of a face.
    \param start	The beginning element.
    \param end		The ending element.  Data from 'start' to 'end' -1 will
    			be read.  If end is zero, data will be read from
                        'start' to the end of the data.
  
   To read the whole array, set start to kCCMIOStart and end to kCCMIOEnd. */
extern CCMIOError CCMIOReadFaces( CCMIOError *err, CCMIOID entity,
				  CCMIOEntity which, CCMIOID *mapID,
				  CCMIOSize *streamSize, int *vertexStream,
				  CCMIOIndex start, CCMIOIndex end );

/** Writes the face data.  See CCMIOReadFaces() for a description of the
    parameters.
    \param streamSize	The total number of integers in the vertex stream.
    \param vertexStream The stream of vertices that compose a face.  This
			should be in the form:<br>
 		        <tt>nVerts v1 v2 ... vn<br>
 		            nVerts v1 v2 ... vn<br>
 			    &nbsp;&nbsp; ...</tt><br>

    <b>Note</b>: When writing partial arrays, the first block of the array you
    write <i>must</i> have start equal to zero. */
extern CCMIOError CCMIOWriteFaces( CCMIOError *err, CCMIOID entity,
                                   CCMIOEntity which, CCMIOID mapID,
                                   CCMIOSize streamSize,
                                   int const *vertexStream,
                                   CCMIOIndex start, CCMIOIndex end );

/** Reads the faces' cell associations with the faces.
    \param entity       The face entity.
    \param which        Either kCCMIOInternalFaces or kCCMIOBoundaryFaces.
    \param mapID        Returns the CCMIOID of the map associated with these
                        faces.  May be NULL if the information is not desired.
    \param cells        If reading internal faces,  this is a two dimensional
                        array with size [nFaces][2].  If reading boundary faces
                        it is a one dimensional array of size nFaces.  Must be
                        pre-allocated to the correct size.  May be NULL if the
                        information is not desired.
    \param start        The beginning face.
    \param end          The ending face.  Data from 'start' to 'end' -1 will
                        be read.  If end is zero, data will be read from
                        'start' to the end of the data.
  
   To read the whole array, set start to kCCMIOStart and end to kCCMIOEnd.

   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to kCCMIOStart. */
CCMIOError CCMIOReadFaceCells( CCMIOError *err, CCMIOID entity,
			       CCMIOEntity which, int *cells,
			       CCMIOIndex start, CCMIOIndex end );

/** Writes the face data.  See CCMIOReadFaces() for a description of the
    parameters.
  
   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to zero. */
CCMIOError CCMIOWriteFaceCells( CCMIOError *err, CCMIOID entity,
				CCMIOEntity which, CCMIOID mapID,
				int const *cells, CCMIOIndex start,
				CCMIOIndex end );

/** Writes the processor information.
    \param processor		The CCMIOID of the processor entity.
    \param verticesFile		The name of the file that the vertices are
				stored in.  NULL indicates the current file.
    \param vertices		The CCMIOID of the vertices entity.  If NULL,
				vertices information will be unchanged.
    \param topologyFile		The name of the file that the mesh is
				stored in.  NULL indicates the current file.
    \param topology		The CCMIOID of the mesh entity.  If NULL,
				mesh information will be unchanged.
    \param initialFieldFile	The name of the file that the inital field is 
				stored in.  NULL indicates the current file.
    \param initialField		The CCMIOID of the initial field set entity.
    				If NULL, the initial field information will be
				unchanged.
    \param solutionFile		The name of the file that the solution data is 
				stored in.  NULL indicates the current file.
    \param vertices		The CCMIOID of the solution data field set entity.
				If NULL, vertices information will be unchanged.
*/
extern CCMIOError CCMIOWriteProcessor( CCMIOError *err, CCMIOID processor,
				   char const *verticesFile, CCMIOID *vertices,
				   char const *topologyFile, CCMIOID *topology,
				   char const *initialFieldFile, CCMIOID *initialField,
				   char const *solutionFile, CCMIOID *solution);

/** Reads the processor information.  The returned entities may be in a 
    a different file and this will need to be closed with CCMIOClose, but since
    the node might be the root of the current file, they should not be closed
    until the end.  'positions' and 'solution' may be NULL if the information
    is not desired. */
extern CCMIOError CCMIOReadProcessor( CCMIOError *err, CCMIOID processor,
				      CCMIOID *vertices, CCMIOID *topology,
				      CCMIOID *initialField, CCMIOID *solution);

/** Clears the relevant information from the processor.  If this information
    is not used elsewhere, it will be removed from the file.
    \param state		The state that is the parent of this processor.
    \param processor		The processor to clear.
    \param clearVertices	If TRUE, vertex information will be cleared.
    				If FALSE, it will be untouched.
    \param clearTopology	If TRUE, topology information will be cleared.
    				If FALSE, it will be untouched.
    \param clearInitialField	If TRUE, initial field information will be 
    				cleared.  If FALSE, it will be untouched.
    \param clearSolution	If TRUE, solution post data will be cleared.
    				If FALSE, it will be untouched.
    \param clearLagrangian	If TRUE, any Lagrangian data will be cleared.
				If FALSE, it will be untouched. */
extern CCMIOError CCMIOClearProcessor( CCMIOError *err, CCMIOID state,
				       CCMIOID processor,
				       int clearVertices, int clearTopology,
				       int clearInitialField, int clearSolution,
				       int clearLagrangian );

/** Writes Lagrangian data.
    \param positionsFile	The file that contains the positions node,
				NULL for this file.
    \param positions		The node that contains the positions.  This node
				must be a kCCMIOVertices node.  (Vertices and
				positions ultimately store the same data)  This
				node may be NULL if no position data is to be
				written.
    \param solutionFile		The file containing the solution node, or
				NULL for this file.
    \param solution		The solution node (field set entity).  May be
				NULL if not solution data is to be written. */
extern CCMIOError CCMIOWriteLagrangianData( CCMIOError *err, CCMIOID lagrangian,
					    char const *positionsFile,
					    CCMIOID *positions,
					    char const *solutionFile,
					    CCMIOID *solution );

/** Reads the Langrangian information.  The returned entities may be in a 
    a different file and this will need to be closed with CCMIOClose, but since
    the node might be the root of the current file, they should not be closed
    until the end.  'positions' and 'solution' may be NULL if the information
    is not desired. */
extern CCMIOError CCMIOReadLagrangianData( CCMIOError *err, CCMIOID lagrangian,
					   CCMIOID *positions,
					   CCMIOID *solution );

/** Writes one component of a vector or tensor data field.
    \param fieldID		The postdata field.
    \param component		The component (\ref kCCMIOVectorX to
    				\ref kCCMIOVectorZ and \ref kCCMIOTensorXX
				to \ref kCCMIOTensorZZ).  All components
				should be specified.  If a component
				is not specified, the application is free
				to use whatever default value it feels is
				correct.
    \param componentField	The component field ID. */
extern CCMIOError CCMIOWriteMultiDimensionalFieldData( CCMIOError *err,
						       CCMIOID fieldID,
						       CCMIOComponent component,
						       CCMIOID componentField );

/** Gets the field ID of one component of a vector or tensor data field.
    Although applications <i>should</i> specify all components, if there
    is a \ref kCCMIONoNodeErr reading the entity, the application
    should assume that this component was not written and should use
    a proper default value (presumably 0.0).  If the file is from an older
    version of the library that stores multiple dimensions in one node,
    a \ref kCCMIOVersionErr is returned.  In this case CCMIOReadFieldData*()
    will read a multiple dimensional array, instead of just one component.
    
    \param fieldID		The postdata field.
    \param component		The component (\ref kCCMIOVectorX to
    				\ref kCCMIOVectorZ and \ref kCCMIOTensorXX
				to \ref kCCMIOTensorZZ).
    \param componentField	The component field ID. */
extern CCMIOError CCMIOReadMultiDimensionalFieldData( CCMIOError *err,
						      CCMIOID fieldID,
						      CCMIOComponent component,
						      CCMIOID *componentField );

/** Writes scalar data for a field.  Vector or tensor data should be written
    as individual components and each component written with
    CCMIOWriteMultiDimensionalFieldData().  Note that if
    CCMIOWriteMultiDimensionalFieldData() returns a kCCMIOVersionErr,
    this function will read an old-style field, which will return all components
    of the array at once, so the 'data' array will need to be 3 or 9
    (depending on whether vector or tensor data is being read) times longer.
    \param fieldData	The CCMIOID of the field entity.
    \param mapID	The CCMIOID of the map corresponding to this field.
    \param loc		What type of data this node is (cell, vertex, face)
    \param data		The data.
    \param start	The offset of the element pointed to by 'data'.  Must
			Be in units of of the first element;  so if
			  float data[10][3][3];
			then start can be in the range [0, 9].
    \param end		The element one beyond the end. So to write from [2, 5],
			start = 2, end = 6, which writes from data[2][0][0]
			to data[5][3][3].

   <b>Note</b>: When writing partial arrays, the first block of the array you
   write <i>must</i> have start equal to kCCMIOStart. */
extern CCMIOError CCMIOWriteFieldDataf( CCMIOError *err, CCMIOID fieldData,
					CCMIOID mapID, CCMIODataLocation loc,
					float *data, CCMIOIndex start,
					CCMIOIndex end  );
/** Writes the data for a field.  See CCMIOWriteFieldDataf() for a description
    of the paramaters. */
extern CCMIOError CCMIOWriteFieldDatad( CCMIOError *err, CCMIOID fieldData,
					CCMIOID mapID, CCMIODataLocation loc,
					double *data, CCMIOIndex start,
					CCMIOIndex end  );

/** Writes the data for a field.  See CCMIOWriteFieldDataf() for a description
    of the paramaters. */
extern CCMIOError CCMIOWriteFieldDatai( CCMIOError *err, CCMIOID fieldData,
					CCMIOID mapID, CCMIODataLocation loc,
					int *data, CCMIOIndex start,
					CCMIOIndex end  );

/** Writes data for a field that is constant.  Note that it is possible for
    one type of data (e.g. \ref kCCMIOCellValues) to be constant and another
    type to be specified.  (Whether this is meaningful is another issue.)
    \param fieldData	The CCMIOID of the field entity.
    \param mapID	The CCMIOID of the map corresponding to this field.
    \param loc		What type of data this node is (cell, vertex, face)
    \param value	The constant value. */
extern CCMIOError CCMIOWriteConstantFieldDataf( CCMIOError *err,
						CCMIOID fieldData,
						CCMIOID mapID,
						CCMIODataLocation loc,
						float value );

/** Writes data for a field that is constant.  See
    CCMIOWriteConstantFieldDataf() for a description of the parameters. */
extern CCMIOError CCMIOWriteConstantFieldDatad( CCMIOError *err,
						CCMIOID fieldData,
						CCMIOID mapID,
						CCMIODataLocation loc,
						double value );

/** Writes data for a field that is constant.  See
    CCMIOWriteConstantFieldDataf() for a description of the parameters. */
extern CCMIOError CCMIOWriteConstantFieldDatai( CCMIOError *err,
						CCMIOID fieldData,
						CCMIOID mapID,
						CCMIODataLocation loc,
						int value );

/** Reads scalar data from a field.  The data is returned in single-precision
    so if the original data is double-precision it will be converted.  If
    the data is a constant value, the entire array will be filled with that
    value.<br><br>
    <b>Note</b>: If the file pre-dates vector and tensor data being stored as
    separate components (that is, if CCMIOReadMultiDimensionalArray()
    returns a \ref kCCMIOVersion error), CCMIOReadFieldDataf() may be used
    to read the data.  If this is the case, CCMIOReadFieldDataf() will
    automatically return the multidimensional data into 'data', and 'data'
    must be big enough to accomodate it.
    \param fieldData	The CCMIOID of the field entity.
    \param mapID	Returns the CCMIOID of the map corresponding to this
    			field.  May be NULL if not desired.
    \param loc		Returns the type of data this node is (cell, vertex,
			face).  May be NULL if not desired.
    \param data		The data, must be preallocated to [n], [n][3],
			or [n][3][3] for kCCMIOScalar, kCCMIOVector, and
			kCCMIOTensor data, respectively.
			May be NULL if not desired.
    \param start	The offset of the element pointed to by 'data'.  Must
			Be in units of of the first element;  so if
			  float data[10][3][3];
			then start can be in the range [0, 9].
    \param end		The element one beyond the end.  So to read from [2, 5],
			start = 2, end = 6, which reads from data[2][0][0]
			to data[5][3][3].  (Note that newer files will only
			have scalar data)*/
extern CCMIOError CCMIOReadFieldDataf( CCMIOError *err, CCMIOID fieldData,
				       CCMIOID *mapID, CCMIODataLocation *loc,
				       float *data, CCMIOIndex start,
				       CCMIOIndex end  );

/** Reads the data from a field.  The data is returned in double-precision
    so if the original data is single-precision it will be converted.
    See CCMIOReadFieldDataf() for a full description. */
extern CCMIOError CCMIOReadFieldDatad( CCMIOError *err, CCMIOID fieldData,
				       CCMIOID *mapID, CCMIODataLocation *loc,
				       double *data, CCMIOIndex start,
				       CCMIOIndex end  );

/** Reads the data from a field.  The data is returned in double-precision
    so if the original data is single-precision it will be converted.
    See CCMIOReadFieldDataf() for a full description. */
extern CCMIOError CCMIOReadFieldDatai( CCMIOError *err, CCMIOID fieldData,
				       CCMIOID *mapID, CCMIODataLocation *loc,
				       int *data, CCMIOIndex start,
				       CCMIOIndex end  );

/** Writes the solver restart entity.
    \param restartInfo	The restart entity.
    \param solverName	The name of the solver that is writing this node.
			only the first \ref kCCMIOMaxStringLength characters
			are meaningful.
    \param iteration	The iteration number.
    \param time		The time (in solver units, with 0.0 being the
			initial time).
    \param timeUnits	A string describing the units of time.  If NULL,
			the units default to "s".
    \param startAngle	The starting angle of the mesh.  Should be 0.0 if
			the mesh is not rotating. */
extern CCMIOError CCMIOWriteRestartInfo( CCMIOError *err, CCMIOID restartInfo,
					 char const *solverName, int iteration,
					 float time, char const *timeUnits,
					 float startAngle );

/** Reads the solver restart entity.  Any pointer may be NULL if the
    information is not desired. */
extern CCMIOError CCMIOReadRestartInfo( CCMIOError *err, CCMIOID restartInfo,
					char *solverName, int *iteration,
					float *time, char *timeUnits,
					float *startAngle );




  /* I changed this routine to include a parameter for the number of faces to
   * speed it up so that there isn't a database read for every buffer write.
   * WRO 22-Oct-2004 */

  /** \internal
      Writes the face data.  See CCMIOReadFaces() for a description of the
      parameters.
    
     <b>Note</b>: When writing partial arrays, the first block of the array you
     write <i>must</i> have start equal to zero. */
extern CCMIOError CCMIOV2WriteFaceCells(
                                  CCMIOError *err, CCMIOID entity,
 				  CCMIOEntity which,
                                  CCMIOSize nFace, int const *cells,
                                  CCMIOIndex start, CCMIOIndex end );

/* \} */

#ifdef __cplusplus
}
#endif
#endif /* CCMIO_H */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
