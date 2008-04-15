#ifndef CCMIO_TYPES_H
#define CCMIO_TYPES_H

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmiotypes.h,v $
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
 *  $Id: ccmiotypes.h,v 1.11 2006/06/05 21:12:16 geoffp Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef DEBUG
    #define Debugging	1
#else
    #define Debugging	0
#endif

#ifndef TRUE
    #define TRUE	1
    #define FALSE	0
#endif

/** Please see \ref errorpage for details on the use of CCMIOError. */
typedef enum {
    kCCMIONoErr = 0,
    kCCMIONoFileErr,
    kCCMIOPermissionErr,
    kCCMIOCorruptFileErr,
    kCCMIOBadLinkErr,
    kCCMIONoNodeErr,
    kCCMIODuplicateNodeErr,
    kCCMIOWrongDataTypeErr,
    kCCMIONoDataErr,
    kCCMIOWrongParentErr,
    kCCMIOBadParameterErr,
    kCCMIONoMemoryErr,
    kCCMIOIOErr,
    kCCMIOTooManyFacesErr,
    kCCMIOVersionErr,
    kCCMIOArrayDimensionToLargeErr,
    kCCMIOInternalErr
    } CCMIOError;

typedef enum { kCCMIOFloat32 = 0, kCCMIOFloat64, kCCMIOInt32,
	       kCCMIOInt64 /* Internal use only */,
               kCCMIOString, kCCMIOUnknownType, kCCMIOBadType,
               kCCMIOLastType } CCMIODataType;

struct _CCMIONode {
    double node;
    double parent;
};

#define kCCMIOEndArgs	-1

/*typedef unsigned long long CCMIOIndex;
  typedef unsigned long long CCMIOSize; */
typedef unsigned int CCMIOIndex;
typedef unsigned int CCMIOSize;

typedef enum { kCCMIORead, kCCMIOWrite } CCMIOBufferType;
typedef CCMIOBufferType CCMIOIOType;

typedef enum { kCCMIOReadData, kCCMIONewData, kCCMIOAddData } CCMIOOpenType;

typedef struct _CCMIONode CCMIONode;

#define IGNORE_ERROR(err)		\
	CCMIOError dummyErr = kCCMIONoErr;	\
	if (!err) err = &dummyErr;

#define CHECK_ERROR(err)		\
	CCMIOError dummyErr = kCCMIONoErr;	\
	if (!err) err = &dummyErr;	\
	if (*err != kCCMIONoErr) return(*err);

#define CHECK_ERROR_AND_CLEAR_PTR(err, ptr, null)	\
	CCMIOError dummyErr = kCCMIONoErr;		\
	if (!err) err = &dummyErr;			\
	if (!ptr) return(*err = kCCMIOBadParameterErr);	\
	else	*ptr = null;				\
	if (*err != kCCMIONoErr) return(*err);
	
#ifdef _WIN32
	#ifndef snprintf
                #if !defined(__NUTC__)
			#define snprintf _snprintf
                #endif
	#endif
#endif
#if defined(__hpux__) || defined(__hpux) || defined(hpux) || defined(__alpha) || defined(__alpha__)
	#define kHasSNPrintf	0
#else
	#define kHasSNPrintf	1
#endif

typedef enum { kCCMIODimNull = 0, kCCMIOScalar = 1, kCCMIOVector, kCCMIOTensor } CCMIODimensionality;
typedef enum { kCCMIOVertex = 0, kCCMIOCell, kCCMIOFace } CCMIODataLocation;

/* Defines for CCMIOGetProstarSet() */
#define kCCMIOCellSet		(1 << 0)
#define kCCMIOVertexSet 	(1 << 1)
#define kCCMIOBoundarySet	(1 << 2)
#define kCCMIOBlockSet		(1 << 3)
#define kCCMIOSplineSet		(1 << 4)
#define kCCMIOCoupleSet		(1 << 5)

/* Keep in sync with gEntityNames and gTypeNames in ccmio.c */
typedef enum { kCCMIONull = -1, kCCMIOMap = 0, kCCMIOVertices, kCCMIOTopology,
	       kCCMIOInternalFaces, kCCMIOBoundaryFaces, kCCMIOCells,
	       kCCMIOProblemDescription, kCCMIOFieldSet, kCCMIOField,
	       kCCMIOFieldData, kCCMIOState, kCCMIOProcessor, kCCMIOCellType,
	       kCCMIOBoundaryRegion, kCCMIOLagrangianData, kCCMIOInterfaces,
	       kCCMIOFieldPhase, kCCMIORestart, kCCMIORestartData,
	       kCCMIOReferenceData, kCCMIOModelConstants, kCCMIOProstarSet,
	       kCCMIOMaxEntity } CCMIOEntity;

typedef enum { kCCMIOVectorX = 0, kCCMIOVectorY, kCCMIOVectorZ,
	       kCCMIOTensorXX = 0, kCCMIOTensorXY, kCCMIOTensorXZ,
	       kCCMIOTensorYX, kCCMIOTensorYY, kCCMIOTensorYZ,
	       kCCMIOTensorZX, kCCMIOTensorZY, kCCMIOTensorZZ } CCMIOComponent;
	       
typedef struct {
    CCMIONode root;
    CCMIONode node;
    int id;
    CCMIOEntity type;
    int version;
} CCMIOID;

#define kCCMIOMaxDimension	4
#define kCCMIOMaxStringLength	32
#define kCCMIOProstarShortNameLength	8
#define kCCMIOStart		0ul
#define kCCMIOEnd		0ul

#ifdef __cplusplus
    }
#endif
#endif /* CCMIO_TYPES_H */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
