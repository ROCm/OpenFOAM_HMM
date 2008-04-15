#ifndef CCMIO_C
#define CCMIO_C

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmio.c,v $
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
 *  $Id: ccmio.c,v 1.29 2006/06/06 13:59:18 cmm Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MAKEDEPEND
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#endif

#include <libadf/ADF.h>
#include "ccmio.h"
#include "ccmiocore.h"
#include "ccmioprivate.h"
#include "ccmioutility.h"
#include "ccmioversion.h"

const char kGeneralNodeName[] = "General";
const char kVersionNodePath[] = "/General/Version";
const char kVersionNodeName[] = "Version";
const char kTitleNodePath[] = "/General/Title";
const char kTitleNodeName[] = "Title";
const char kEntityIDDelimiter = '-';
const char kLabelName[] = "Label";
const char kMapName[] = "MapId";
const char kMapMaxName[] = "MaximumMappedId";
const char kMapDataName[] = "IdMap";
const char kMeshesName[] = "Meshes";
const char kVertexDimensionName[] = "Dimension";
const char kVertexScaleFactorName[] = "ScaleFactor";
const char kVertexCoordinatesName[] = "Coordinates";
const char kNumFacesName[] = "NumFaces";
const char kFacesVertexDataName[] = "Vertices";
const char kFacesCellDataName[] = "Cells";
const char kCellDataName[] = "CellType";
const char kNumCellsName[] = "NumCells";
const char kNBoundariesName[] = "NumBoundaryTypes";
const char kVerticesFileName[] = "VerticesFile";
const char kTopologyFileName[] = "TopologyFile";
const char kInitialFieldFileName[] = "InitialFieldFile";
const char kSolutionFileName[] = "SolutionFileName";
const char kPositionsFileName[] = "PositionsFile";
const char kVerticesIDName[] = "VerticesId";
const char kTopologyIDName[] = "TopologyId";
const char kInitialFieldIDName[] = "InitialFieldId";
const char kSolutionIDName[] = "SolutionId";
const char kPositionsIDName[] = "PositionsId";
const char kProblemDescriptionName[] = "ProblemDescription";
const char kFieldTypeName[] = "FieldType";
const char kFieldShortName[] = "ShortName";
const char kFieldUnitsName[] = "Units";
const char kFieldDataTypeName[] = "DataType";
const char kFieldDataName[] = "Data";
const char kMapMarkName[] = "unused";
const char kTempName[] = "tmp";
const char kSolverName[] = "SolverName";
const char kIterationName[] = "Iteration";
const char kTimeName[] = "Time";
const char kTimeUnitsName[] = "TimeUnits";
const char kStartAngleName[] = "StartAngle";
const char kDefaultTimeUnits[] = "s";
const char kConstantDataName[] = "ConstantData";
const char kComponentName[] = "Component";
const char kProstarSetsName[] = "ProstarSets";
const char kCellSetName[] = "CellSet";
const char kVertexSetName[] = "VertexSet";
const char kBoundarySetName[] = "BoundarySet";
const char kSplineSetName[] = "SplineSet";
const char kBlockSetName[] = "BlockSet";
const char kCoupleSetName[] = "CoupleSet";


static char gEntityNames[][kCCMIOMaxStringLength + 1] =
    { "Map", "Vertices", "FaceBasedTopology", "InternalFaces",
      "BoundaryFaces", "Cells", "ProblemDescription",
      "FieldSet", "Field", "Data", "State", "Processor", "CellType",
      "BoundaryRegion", "LagrangianData", "Interfaces", "Phase", "RestartInfo",
      "RestartData", "ReferenceData", "ModelConstants", "ProstarSet", "" };

static char gTypeNames[][kCCMIOMaxStringLength + 1] = 
    { "map", "vertices", "topology", "internalFaces", "boundaryFaces",
      "cells", "problemDescription", "fieldSet", "field", "fieldData", "state",
      "processor", "cellType", "boundaryRegion", "lagrangianData", 
      "interfaces", "phase", "restartInfo", "restartData", "referenceData",
      "modelConstants", "prostarSet", "" };

static void ClearCCMIOID( CCMIOID *id );
static void MakeRootCCMIOID( CCMIONode root, CCMIOID *id  );
static int CheckEntityParent( CCMIOID parent, CCMIOEntity childType );
static CCMIOError GetEntityParent( CCMIOError *err, CCMIOID id, CCMIOEntity type,
				 CCMIOID *parent );
static int ParseEntityID( const char *name );
static void MakeEntityName( CCMIOEntity type, int id, char *name);
static void ClearStateProbDef( CCMIOError *err, CCMIOID state );
static CCMIOError MarkMaps( CCMIOError *err, CCMIOID anyID, CCMIONode top, int unused );
static CCMIOError RemoveUnusedMaps( CCMIOError *err, CCMIONode root );
static CCMIOError ReadProcessorItem( CCMIOError *err, CCMIOID processor, CCMIOID *id,
				   const char *fileNodeName,
				   const char *mapNodeName, CCMIOEntity type );
static int IsProcessorEntityUsed( CCMIOError *err, CCMIOID state, CCMIOID proc,
				  CCMIOID *lagrangian, const char *fileNodeName,
				  const char *idNodeName );

/* ------------ Entity Cache -------------- */
/* Creating entities automatically requires knowing that last entity
   id.  Walking through them all each time is too slow with a large
   number of nodes, particularly since ADF gets very slow after about
   300 child nodes.  This cache attempts to ameliorate the problem. */
struct MaxIndexCache
{
    CCMIONode root;
    unsigned int max[kCCMIOMaxEntity];
    CCMIONode parent[kCCMIOMaxEntity];
};
static struct MaxIndexCache gCache = { { 0, 0 }, { 0 }, { { 0, 0} } };

static int CacheNeedsClear( struct MaxIndexCache *cache, CCMIOID parent,
			    CCMIOEntity type )
{
    if (CCMIOAreNodesEqual(cache->root, parent.root))
	return(!CCMIOAreNodesEqual(cache->parent[type], parent.node));
    return(TRUE);
}

static void CacheClear( struct MaxIndexCache *cache, CCMIOID parent,
			CCMIOEntity type )
{
    int i;

    if (type == kCCMIONull)
    {
	cache->root = parent.root;
	for (i = 0;  i < kCCMIOMaxEntity;  ++i)
	{
	    cache->max[i] = 0;
	    cache->parent[i].node = 0.0;
	    cache->parent[i].parent = 0.0;
	}
    }
    else
    {
	cache->max[type] = 0;
	cache->parent[type] = parent.node;
    }
}

static int CacheGet( struct MaxIndexCache *cache, CCMIOEntity type )
{
    return(cache->max[(int)type]);
}

static int CacheIncrement( struct MaxIndexCache *cache, CCMIOEntity type )
{
    return(cache->max[(int)type]++);
}

static void CacheSet( struct MaxIndexCache *cache, CCMIOEntity type,
		      unsigned int maxID )
{
    cache->max[(int)type] = maxID;
}
/* ---------------------------------------- */

void ClearCCMIOID( CCMIOID *id )
{
    if (id)
    {
	id->root.parent = 0;
	id->root.node = 0;
	id->node.parent = 0;
	id->node.node = 0;
	id->id = 0;
	id->type = kCCMIOMaxEntity;
	id->version = 0;
    }
}

void MakeRootCCMIOID( CCMIONode root, CCMIOID *id )
{
    if (id)
    {
	id->root = root;
	id->node = root;
	id->id = 0;
	id->type = kCCMIONull;
	id->version = 0;
	CCMIOGetVersion(NULL, root, &id->version);
    }
}

CCMIOError CCMIOGetVersion( CCMIOError *err, CCMIONode root, int *version )
{
    CHECK_ERROR(err);
    if (!version) return(*err = kCCMIOBadParameterErr);

    CCMIOReadNodei(err, root, kVersionNodePath, version);
    if (*err)
	return(*err = kCCMIOCorruptFileErr);
    return(*err = kCCMIONoErr);
}

CCMIOError CCMIOSetVersion( CCMIOError *err, CCMIONode root, int version )
{
    CCMIONode general;

    CHECK_ERROR(err);

    if (CCMIOGetNode(NULL, root, kGeneralNodeName, &general))
	if (CCMIOCreateNode(err, root, FALSE, kGeneralNodeName, kGeneralNodeName,
			  &general))
	    return(*err = kCCMIOCorruptFileErr);
    CCMIOWriteNodei(err, general, kVersionNodeName, version);
    if (*err)
	return(*err = kCCMIOCorruptFileErr);
    return(*err = kCCMIONoErr);
}

CCMIOError CCMIOGetTitle( CCMIOError *err, CCMIONode root, char **title )
{
    CHECK_ERROR(err);
    if (!title) return(*err = kCCMIOBadParameterErr);

    CCMIOReadNodestr(err, root, kTitleNodePath, title);
    if (*err)
	return(*err = kCCMIOCorruptFileErr);
    return(*err = kCCMIONoErr);
}

CCMIOError CCMIOSetTitle( CCMIOError *err, CCMIONode root, const char *title )
{
    CCMIONode general;

    CHECK_ERROR(err);
    if (!title) return(*err = kCCMIOBadParameterErr);

    if (CCMIOGetNode(NULL, root, kGeneralNodeName, &general))
	if (CCMIOCreateNode(NULL, root, FALSE, kGeneralNodeName,
			    kGeneralNodeName, &general))
	    return(*err = kCCMIOCorruptFileErr);
    CCMIOWriteNodestr(err, general, kTitleNodeName, title);
    if (*err)
	return(*err = kCCMIOCorruptFileErr);
    return(*err = kCCMIONoErr);
}

CCMIOError CCMIOOpenFile( CCMIOError *err, const char *file, CCMIOIOType mode,
			  CCMIOID *root )
{
    CCMIONode node;
    CHECK_ERROR(err);

    if (root)
	ClearCCMIOID(root);
    else
	return(*err = kCCMIOBadParameterErr);

    *err = CCMIOOpen(file, mode, &node);
    if (*err != kCCMIONoErr)
	return(*err);

    MakeRootCCMIOID(node, root);

    /* If you closed a file and opened a new one, ADF might give the new
       one the same root ID as the old one, so clear the cache, just in case.*/
    CacheClear(&gCache, *root, kCCMIONull);

    return(*err);
}

CCMIOError CCMIOCloseFile( CCMIOError *err, CCMIOID root )
{
    IGNORE_ERROR(err);

    CCMIOClose(root.root);
    return(*err);
}

CCMIOError CCMIOWriteOpti( CCMIOError *err, CCMIOID parent, const char *name,
		       int value )
{
    return(CCMIOWriteNodei(err, parent.node, name, value));
}

CCMIOError CCMIOWriteOptf( CCMIOError *err, CCMIOID parent, const char *name,
			   float value )
{
    return(CCMIOWriteNodef(err, parent.node, name, value));
}

CCMIOError CCMIOWriteOptd( CCMIOError *err, CCMIOID parent, const char *name,
			   double value )
{
    return(CCMIOWriteNoded(err, parent.node, name, value));
}

CCMIOError CCMIOWriteOptstr( CCMIOError *err, CCMIOID parent,
			     const char *name, const char *value )
{
    return(CCMIOWriteNodestr(err, parent.node, name, value));
}

CCMIOError CCMIOReadOptstr( CCMIOError *err, CCMIOID parent, const char *name,
			    int *size, char *value )
{
    if (size)	*size = 0;
    {
    char *tmp;
    CHECK_ERROR(err);

    CCMIOReadNodestr(err, parent.node, name, &tmp);
    if (tmp)
    {
	if (size)
	    *size = strlen(tmp);
	if (value)
	    strcpy(value, tmp);
    }
    else
    {
	if (size)
	    *size = 0;
	if (value)
	    *value = '\0';
    }
    free(tmp);

    return(*err);
    }
}

CCMIOError CCMIOGetOptInfo( CCMIOError *err, CCMIOID parent, const char *name,
			    CCMIODataType *type, CCMIOIndex *x,
			    CCMIOIndex *y, CCMIOIndex *z )
{
    if (type)	*type = kCCMIOBadType;
    if (x)	*x = 0;
    if (y)	*y = 0;
    if (z)	*z = 0;
    {
    char nodeName[kCCMIOMaxStringLength + 1];
    CCMIONode node;
    CHECK_ERROR(err);

    CCMIOGetName(err, parent.node, nodeName);
    CCMIOGetNode(err, parent.node, name, &node);
    if (type)
	CCMIOGetDataType(err, node, type);
    if (parent.version < 20300 && (x || y || z))
    {
	int nDims;
	CCMIOIndex *dims;

	CCMIOGetDimensions(err, node, &nDims, &dims);
	if (*err != kCCMIONoErr)
	    return(*err);
	if (nDims >= 1 && x)
	    *x = dims[0];
	if (nDims >= 2 && y)
	    *y = dims[1];
	if (nDims >= 2 && z)
	    *z = dims[2];
	free(dims);
    }
    else if (x || y || z)
    {
	int nDims;
	CCMIOIndex *dims;

	CCMIOGetDimensions(err, node, &nDims, &dims);
	if (*err != kCCMIONoErr)
	    return(*err);
	if (nDims >= 1 && x)
	    *x = dims[nDims - 1];
	if (nDims >= 2 && y)
	    *y = dims[nDims - 2];
	if (nDims >= 2 && z)
	    *z = dims[nDims - 3];
	free(dims);
    }
    }
    return(*err);
}

void CCMIOInvalidateEntity( CCMIOID *entity )
{
    ClearCCMIOID(entity);
}

int CCMIOIsValidEntity( CCMIOID entity )
{
    return(entity.type != kCCMIOMaxEntity);
}

int CCMIOIsFromSameFile( CCMIOID entity1, CCMIOID entity2 )
{
    return(CCMIOAreNodesEqual(entity1.root, entity2.root));
}

int CheckEntityParent (CCMIOID parent, CCMIOEntity childType)
{
  switch (childType)
  {
    case kCCMIOProcessor:
      return (parent.type != kCCMIOState) ? 0 : 1;

    case kCCMIOBoundaryFaces:
      return (parent.type != kCCMIOTopology) ? 0 : 1;

    case kCCMIOCellType:
      return (parent.type != kCCMIOProblemDescription) ? 0 : 1;

    case kCCMIOBoundaryRegion:
      return (parent.type != kCCMIOProblemDescription) ? 0 : 1;

    case kCCMIOInterfaces:
      return (parent.type != kCCMIOCells &&
              parent.type != kCCMIOTopology) ? 0 : 1;

    default:
      return 1;
  }
  return 1;
}

CCMIOError GetEntityParent( CCMIOError *err, CCMIOID id, CCMIOEntity type,
			    CCMIOID *parent )
{
    if (parent)		ClearCCMIOID(parent);
    {
    CHECK_ERROR(err);
    if (!parent)	return(*err = kCCMIOBadParameterErr);

    switch(type)
    {
	case kCCMIOVertices:
	case kCCMIOTopology:
	    CCMIOCreateNode(err, id.root, TRUE, kMeshesName, kMeshesName,
			  &parent->node);
	    break;
	case kCCMIOMap:
	case kCCMIOProblemDescription:
	case kCCMIOFieldSet:
	case kCCMIOField:
	case kCCMIOState:
	    {
	    char name[kCCMIOMaxStringLength + 1];

	    CCMIONode node = id.root;
	    if (type == kCCMIOField)
		node = id.node;
#if kHasSNPrintf
	    snprintf(name, kCCMIOMaxStringLength, "%ss",
		     gEntityNames[(int)type]);
#else
	    sprintf(name, "%ss", gEntityNames[(int)type]);
#endif
	    name[kCCMIOMaxStringLength] = '\0';
	    CCMIOCreateNode(err, node, TRUE, name, gTypeNames[(int)type],
			    &parent->node);
	    }
	    break;
	case kCCMIOLagrangianData:
	    CCMIOCreateNode(err, id.node, TRUE, gEntityNames[(int)type],
			    gTypeNames[(int)type], &parent->node);
	    break;
    	case kCCMIOProstarSet:
	    CCMIOCreateNode(err, id.root, TRUE, kProstarSetsName,
			    kProstarSetsName, &parent->node);
	    break;
	default:
	    if (!CheckEntityParent(id, type))
		return(*err = kCCMIOWrongParentErr);
	    parent->node = id.node;
	    break;
    }
    parent->root = id.root;
    parent->id = 0;
    parent->type = type;
    parent->version = id.version;
    return(*err);

    }
}

int ParseEntityID( const char *name )
{
    char const *pos = strrchr(name, kEntityIDDelimiter);

    if (pos == NULL)
	return(0);
    if (pos != name && (*(pos - 1) == '-'))	/* Negative ID number */
	--pos;
    return(atoi(++pos));
}

void MakeEntityName( CCMIOEntity type, int id, char *name )
{
    if (type == kCCMIOInternalFaces || type == kCCMIOCells ||
	type == kCCMIOInterfaces || type == kCCMIORestart ||
	type == kCCMIORestartData || type == kCCMIOReferenceData ||
	type == kCCMIOModelConstants)
#if kHasSNPrintf
	snprintf(name, kCCMIOMaxStringLength, gEntityNames[(int)type]);
#else
	sprintf(name, gEntityNames[(int)type]);
#endif
    else
#if kHasSNPrintf
	snprintf(name, kCCMIOMaxStringLength, "%s%c%d",
		 gEntityNames[(int)type], kEntityIDDelimiter, id);
#else
	sprintf(name, "%s%c%d",
		 gEntityNames[(int)type], kEntityIDDelimiter, id);
#endif
    name[kCCMIOMaxStringLength] = '\0';
}

CCMIOError GetEntityCore( CCMIOError *err, CCMIOID parent, int makeNew,
			  CCMIOEntity type, int idVal,
			  const char *description, CCMIOID *id )
{
    char name[kCCMIOMaxStringLength + 1];

    MakeEntityName(type, idVal, name);
    if (makeNew)
    {
	CCMIOCreateNode(err, parent.node, FALSE, name, gTypeNames[(int)type],
		      &id->node);
	if (description)
	    CCMIOWriteNodestr(err, id->node, kLabelName, description);
    }
    else
    {
	GetEntityParent(err, parent, type, &parent);
	CCMIOGetNode(err, parent.node, name, &id->node);
    }
    id->id = idVal;
    id->type = type;
    id->root = parent.root;
    id->version = parent.version;
    return(*err);
}

CCMIOError CCMIONewEntity( CCMIOError *err, CCMIOID parent, CCMIOEntity type,
			   const char *description, CCMIOID *id )
{
    /* Keep a cache of the last known maximum IDs which needs to be cleared
       if we start writing to a different file.  Only one file is cached because
       it is simpler and it is unlikely that entities will be alternatingly
       created. */
    if (CacheNeedsClear(&gCache, parent, type))
	CacheClear(&gCache, parent, type);

    if (id)	ClearCCMIOID(id);
    {
    char entityName[kCCMIOMaxStringLength + 1];
    int n, i = 0, maxID = 0;
    CCMIONode next;
    CCMIOID realParent;
    CHECK_ERROR(err);

    if (type == kCCMIOBoundaryFaces || type == kCCMIOState)
	return(*err = kCCMIOBadParameterErr);
    GetEntityParent(err, parent, type, &realParent);
    if (*err != kCCMIONoErr)
	return(*err);

    if (CacheGet(&gCache, type) != 0)
	maxID = CacheIncrement(&gCache, type);
    else
    {
	while (*err == kCCMIONoErr &&
	       CCMIOGetNextChildWithLabel(NULL,
					  realParent.node,gTypeNames[(int)type],
					  &i, &next) == kCCMIONoErr)
	{
	    CCMIOGetName(err, next, entityName);
	    n = ParseEntityID(entityName);
	    if (n > maxID)
		maxID = n;
	}
	CacheSet(&gCache, type, maxID + 1);
    }

    GetEntityCore(err, realParent, TRUE, type, maxID + 1, description, id);
    return(*err);
    }
}

CCMIOError CCMIONewIndexedEntity( CCMIOError *err, CCMIOID parent,
				  CCMIOEntity which, int idVal,
				  const char *description, CCMIOID *id )
{
    if (id)	ClearCCMIOID(id);
    {
    int nBoundaries;
    CHECK_ERROR(err);
    if (which != kCCMIOBoundaryFaces && which != kCCMIOCellType &&
	which != kCCMIOBoundaryRegion && which != kCCMIOFieldPhase)
	return(*err = kCCMIOWrongParentErr);

    GetEntityCore(err, parent, TRUE, which, idVal, description, id);
    if (*err != kCCMIONoErr)
	return(*err);
    if (which == kCCMIOBoundaryFaces)
    {
	CCMIOReadNodei(err, parent.node, kNBoundariesName, &nBoundaries);
	if (*err == kCCMIONoNodeErr)
	{
	    *err = kCCMIONoErr;
	    nBoundaries = 0;
	}
	CCMIOWriteNodei(err, parent.node, kNBoundariesName, nBoundaries + 1);
    }

    return(*err);
    }

}

CCMIOError CCMIOGetEntity( CCMIOError *err, CCMIOID parent, CCMIOEntity type,
			   int idVal, CCMIOID *id )
{
    if (id)	ClearCCMIOID(id);
    {
    CHECK_ERROR(err);
    
    if (!CheckEntityParent(parent, type))
	return(*err = kCCMIOWrongParentErr);
    if (type == kCCMIOState || type == kCCMIOField)
	return(*err = kCCMIOBadParameterErr);

    GetEntityCore(err, parent, FALSE, type, idVal, NULL, id);
    return(*err);
    }
}

void ClearStateProbDef( CCMIOError *err, CCMIOID state )
{
    int i = 0, idVal, idVal2;
    CCMIONode node;
    CCMIOID next, problem;

    if (CCMIOReadNodei(NULL, state.node, kProblemDescriptionName, &idVal)
								!= kCCMIONoErr)
	return;
    CCMIOGetEntity(err, state, kCCMIOProblemDescription, idVal, &problem);
    while (CCMIONextEntity(NULL, state, kCCMIOState, &i, &next) == kCCMIONoErr)
    {
	if (CCMIOAreNodesEqual(state.node, next.node))
	    continue;
	if (CCMIOReadNodei(NULL, next.node, kProblemDescriptionName, &idVal2)
	    				      == kCCMIONoErr && idVal == idVal2)
	    return;
    }
 
    CCMIODeleteEntity(err, problem);
    CCMIOGetNode(err, state.node, kProblemDescriptionName, &node);
    CCMIODeleteNode(err, node);
    return;
}

CCMIOError CCMIONewState( CCMIOError *err, CCMIOID root, const char *name,
			  CCMIOID *problemDescription,
			  const char *description, CCMIOID *state )
{
    CCMIOID tmpID;
    if (state)
	ClearCCMIOID(state);
    else
	state = &tmpID;
    {
    int i = 0;
    CCMIONode node;
    CCMIOID parent, proc;
    CHECK_ERROR(err);
    if (!name)	return(*err = kCCMIOBadParameterErr);

    GetEntityParent(err, root, kCCMIOState, &parent);
    CCMIOCreateNode(err, parent.node, TRUE, name, gTypeNames[(int)kCCMIOState],
		  &node);
    state->root = root.root;
    state->node = node;
    state->id = 0;
    state->type = kCCMIOState;

    /* In case this node already exists, remove everything */
    ClearStateProbDef(err, *state);
    while(CCMIONextEntity(NULL, *state, kCCMIOProcessor, &i, &proc)
								== kCCMIONoErr)
	CCMIOClearProcessor(NULL, *state, proc, TRUE, TRUE, TRUE, TRUE, TRUE);
    CCMIODeleteAllChildren(NULL, state->node);

    if (problemDescription)
	CCMIOWriteNodei(err, node, kProblemDescriptionName,
		      problemDescription->id);

    if (description)
	CCMIOWriteNodestr(err, state->node, kLabelName, description);

    return(*err);
    }
}

CCMIOError CCMIOGetState( CCMIOError *err, CCMIOID root, const char *name,
			  CCMIOID *problemDescription, CCMIOID *state)
{
    CCMIOID tmpState;
    if (problemDescription)	ClearCCMIOID(problemDescription);
    if (state)
	ClearCCMIOID(state);
    else
	state = &tmpState;
    {
    CCMIONode node;
    CCMIOID parent;
    CHECK_ERROR(err);
    if (!name)	return(*err = kCCMIOBadParameterErr);

    GetEntityParent(err, root, kCCMIOState, &parent);
    CCMIOGetNode(err, parent.node, name, &node);
    if (*err != kCCMIONoErr)
	return(*err);
    state->root = root.root;
    state->node = node;
    state->id = 0;
    state->type = kCCMIOState;
    state->version = root.version;
    if (problemDescription)
    {
	CCMIOReadNodei(err, state->node, kProblemDescriptionName,
		     &problemDescription->id);
	if (*err == kCCMIONoNodeErr)
	    *err = kCCMIONoErr;
	else
	    CCMIOGetEntity(err, *state, kCCMIOProblemDescription,
			 problemDescription->id, problemDescription);
    }
    }
    return(*err);
}

CCMIOError CCMIONewField( CCMIOError *err, CCMIOID phase, const char *name,
			  const char *shortName, CCMIODimensionality dim,
			  CCMIOID *field )
{
    if (field)	ClearCCMIOID(field);
    {
    CCMIOID parent;
    CHECK_ERROR(err);

    GetEntityParent(err, phase, kCCMIOField, &parent);
    CCMIOCreateNode(err, parent.node, TRUE, name, gTypeNames[(int)kCCMIOField],
		  &field->node);
    field->root = parent.root;
    field->id = 0;
    field->type = kCCMIOField;
    CCMIOWriteNodei(err, field->node, kFieldTypeName, dim);
    if (shortName)
    {
	char sname[kCCMIOProstarShortNameLength + 1] = "";
	snprintf(sname, kCCMIOProstarShortNameLength + 1, "%s", shortName);
	CCMIOWriteOptstr(err, *field, kFieldShortName, sname);
    }

    return(*err);
    }
}

CCMIOError CCMIOGetField( CCMIOError *err, CCMIOID phase, const char *name,
			  CCMIODimensionality *dim, CCMIOID *field )
{
    if (field)	ClearCCMIOID(field);
    {
    CCMIOID parent;
    CHECK_ERROR_AND_CLEAR_PTR(err, dim, kCCMIODimNull);
    if (phase.type != kCCMIOFieldPhase)
	return(*err = kCCMIOBadParameterErr);

    GetEntityParent(err, phase, kCCMIOField, &parent);
    CCMIOGetNode(err, parent.node, name, &field->node);
    field->root = parent.root;
    field->id = 0;
    field->type = kCCMIOField;
    field->version = phase.version;
    if (dim)
	CCMIOReadNodei(err, field->node, kFieldTypeName, (int *)dim);

    return(*err);
    }
}

CCMIOError CCMIOReadField( CCMIOError *err, CCMIOID field, char *name,
			   char *shortName, CCMIODimensionality *dim,
			   CCMIODataType *datatype )
{
    int dimensions;
    if (name)		name[0] = 0;
    if (shortName)	shortName[0] = 0;
    if (dim)		*dim = kCCMIODimNull;
    if (datatype)	*datatype = kCCMIOFloat32;
    {
    CHECK_ERROR(err);

    if (name)
	CCMIOGetName(err, field.node, name);
    if (shortName)
    {
	if (field.version >= 20200)
	    CCMIOReadOptstr(err, field, kFieldShortName, NULL, shortName);
    }
    CCMIOReadNodei(err, field.node, kFieldTypeName, &dimensions);
    if (dim)
	*dim = (CCMIODimensionality) dimensions;
    if (datatype && dimensions == 1)
    {
	int i = 0;
	CCMIONode node;
	CCMIOID fieldData;

	CCMIONextEntity(err, field, kCCMIOFieldData, &i, &fieldData);
	if (CCMIOGetNode(err, fieldData.node, kFieldDataName, &node)
							    == kCCMIONoNodeErr)
	{
	    /* Check if this is a constant data node, and if so,
	       get the data type of that node instead.*/
	    *err = kCCMIONoErr;
	    CCMIOGetNode(err, fieldData.node, kConstantDataName, &node);
	}
	CCMIOGetDataType(err, node, datatype);
    }
    else if (datatype)
	/* Vectors and tensors are just pointers to scalars, so there is
	   no surety that each component will have the same data type.*/
	*datatype = kCCMIOUnknownType;

    return(*err);
    }
    
}

CCMIOError CCMIONewProstarSet( CCMIOError *err, CCMIOID root, const char *name,
			       const char *longName, CCMIOID *setID )
{
    char adfName[kCCMIOMaxStringLength + 1];
    CCMIONode node;
    CCMIOID parent, tmpID;
    CHECK_ERROR(err);
    if (!name)	return(*err = kCCMIOBadParameterErr);
    if (!setID) setID = &tmpID;

    strncpy(adfName, name, kCCMIOMaxStringLength);
    adfName[kCCMIOMaxStringLength] = '\0';  /* Just in case */
    GetEntityParent(err, root, kCCMIOProstarSet, &parent);
    if (CCMIOGetNode(NULL, parent.node, adfName, &node) == kCCMIONoErr)
	CCMIODeleteNode(NULL, node);  /* Get rid of the sets */
    CCMIOCreateNode(err, parent.node, TRUE, adfName, adfName, &setID->node);
    if (longName)
	CCMIOWriteOptstr(err, *setID, kLabelName, longName);

    return(*err);
}

CCMIOError CCMIOGetProstarSet( CCMIOError *err, CCMIOID root, const char *name,
			       int *longNameSize, char *longName,
			       unsigned int *setsAvailableFlag,
			       CCMIOID *setID )
{
    char adfName[kCCMIOMaxStringLength + 1];
    CCMIONode node;
    CCMIOID parent, tmpID;
    CHECK_ERROR(err);
    if (!name)	return(*err = kCCMIOBadParameterErr);
    if (!setID) setID = &tmpID;

    strncpy(adfName, name, kCCMIOMaxStringLength);
    adfName[kCCMIOMaxStringLength] = '\0';  /* Just in case */
    GetEntityParent(err, root, kCCMIOProstarSet, &parent);
    MakeRootCCMIOID(root.root, setID);
    CCMIOGetNode(err, parent.node, name, &setID->node);

    if (longNameSize || longName)
    {
	CCMIOReadOptstr(err, *setID, kLabelName, longNameSize, longName);
	if (*err == kCCMIONoNodeErr)  /* If no label node, use the short name */
	{
	    if (longName)
		strncpy(longName, adfName, strlen(adfName) + 1);
	    if (longNameSize)
		*longNameSize = strlen(adfName);
	    *err = kCCMIONoErr;
	}
    }

    if (setsAvailableFlag && *err == kCCMIONoErr)
    {
	*setsAvailableFlag = 0;
	if (CCMIOGetNode(NULL, setID->node, kCellSetName, &node) == kCCMIONoErr)
	    *setsAvailableFlag |= kCCMIOCellSet;
	if (CCMIOGetNode(NULL, setID->node, kVertexSetName, &node) == kCCMIONoErr)
	    *setsAvailableFlag |= kCCMIOVertexSet;
	if (CCMIOGetNode(NULL, setID->node, kBoundarySetName, &node) == kCCMIONoErr)
	    *setsAvailableFlag |= kCCMIOBoundarySet;
	if (CCMIOGetNode(NULL, setID->node, kBlockSetName, &node) == kCCMIONoErr)
	    *setsAvailableFlag |= kCCMIOBlockSet;
	if (CCMIOGetNode(NULL, setID->node, kSplineSetName, &node) == kCCMIONoErr)
	    *setsAvailableFlag |= kCCMIOSplineSet;
	if (CCMIOGetNode(NULL, setID->node, kCoupleSetName, &node) == kCCMIONoErr)
	    *setsAvailableFlag |= kCCMIOCoupleSet;
    }

    return(*err);
}

CCMIOError CCMIOGetEntityIndex( CCMIOError *err, CCMIOID id, int *n )
{
    CHECK_ERROR_AND_CLEAR_PTR(err, n, 0);

    if (id.type != kCCMIOBoundaryFaces && id.type != kCCMIOCellType &&
	id.type != kCCMIOBoundaryRegion && id.type != kCCMIOFieldPhase &&
	id.type != kCCMIOProcessor)
	return(*err = kCCMIOBadParameterErr);

    *n = id.id;

    return(*err);
}

CCMIOError CCMIODeleteEntity( CCMIOError *err, CCMIOID id )
{
    int i = 0;
    CCMIOID child;
    CHECK_ERROR(err);

    if (id.type == kCCMIOState)
    {
	ClearStateProbDef(err, id);
	while(CCMIONextEntity(NULL, id, kCCMIOProcessor, &i, &child) == kCCMIONoErr)
	    CCMIOClearProcessor(NULL, id, child, TRUE, TRUE, TRUE, TRUE, TRUE);
    }

    CacheSet(&gCache, id.type, 0);
    CCMIODeleteNode(err, id.node);
    if (id.type == kCCMIOBoundaryFaces)
    {
	int nBoundaries;

	CCMIOReadNodei(err, id.node, kNBoundariesName, &nBoundaries);
	CCMIOWriteNodei(err, id.node, kNBoundariesName, nBoundaries - 1);
    }
    return(*err);
}

CCMIOError CCMIONextEntity( CCMIOError *err, CCMIOID parent, CCMIOEntity type,
			    int *i, CCMIOID *next )
{
    if (next)	ClearCCMIOID(next);
    {
    char name[kCCMIOMaxStringLength + 1];
    CHECK_ERROR(err);

    GetEntityParent(err, parent, type, &parent);

    if (CCMIOGetNextChildWithLabel(err, parent.node, gTypeNames[(int)type], i,
				 &next->node) == kCCMIONoErr)
    {
	next->root = parent.root;
	next->type = type;
	CCMIOGetName(err, next->node, name);
	next->id = ParseEntityID(name);
	next->version = parent.version;
    }
    return(*err);
    }
}

CCMIOError CCMIOEntitySize( CCMIOError *err, CCMIOID id, CCMIOSize *n,
			    CCMIOIndex *max )
{
    if (max)	*max = 0;
    if (n)	*n = 0;
    {
    int nDims;
    CCMIOIndex *dims;
    CCMIONode node;
    CHECK_ERROR(err);

    if (!(id.type == kCCMIOMap || id.type == kCCMIOVertices
	  || id.type == kCCMIOCells || id.type == kCCMIOInternalFaces
	  || id.type == kCCMIOBoundaryFaces || id.type == kCCMIOFieldData))
	return(*err);	/* No array data */

    if (id.type != kCCMIOMap)
    {
	int mapVal;
	CCMIOReadNodei(err, id.node, kMapName, &mapVal);
	CCMIOGetEntity(err, id, kCCMIOMap, mapVal, &id);
    }
    if (n)
    {
	CCMIOGetNode(err, id.node, kMapDataName, &node);
	CCMIOGetDimensions(err, node, &nDims, &dims);
	if (*err == kCCMIONoErr)
	    *n = dims[nDims - 1];
	free(dims);
    }
    if (max)
	CCMIOReadNodei(err, id.node, kMapMaxName, (int *) max);

    return(*err);
    }
}

CCMIOError CCMIOEntityName( CCMIOError *err, CCMIOID id, char *name )
{
    CHECK_ERROR(err);

    CCMIOGetName(err, id.node, name);
    return(*err);
}

CCMIOError CCMIOEntityLabel( CCMIOError *err, CCMIOID id, CCMIOSize *size,
			     char *label )
{
    if (size)	*size = 0;
    if (label)	*label = '\0';
    {
    CHECK_ERROR(err);

    CCMIOReadOptstr(err, id, kLabelName, (int *) size, label);
    if (*err == kCCMIONoNodeErr)
	*err = kCCMIONoErr;
    return(*err);
    }
}

CCMIOError CCMIOEntityDescription( CCMIOError *err, CCMIOID id,
				   CCMIOSize *size, char *desc )
{
    if (size)	*size = 0;
    if (desc)	*desc = '\0';
    {
    CHECK_ERROR(err);

    if (CCMIOReadOptstr(NULL, id, kLabelName, (int*)size, desc) != kCCMIONoErr)
    {
	if (size)
	    *size = 0;
	if (desc)
	    *desc = '\0';
    }
    return(*err);
    }
}

CCMIOError CCMIOGetEntityNode( CCMIOError *err, CCMIOID id, CCMIONode *node )
{
    CHECK_ERROR(err);

    *node = id.node;
    return(*err);
}

CCMIOError CCMIOEntityDataType( CCMIOError *err, CCMIOID id,
				CCMIODataType *type )
{
    CCMIONode node, dataNode;
    CHECK_ERROR_AND_CLEAR_PTR(err, type, kCCMIOFloat32);

    CCMIOGetEntityNode(err, id, &node);
    switch(id.type)
    {
	case kCCMIOVertices:
	    CCMIOGetNode(err, node, kVertexCoordinatesName, &dataNode);
	    CCMIOGetDataType(err, dataNode, type);
	    break;
	case kCCMIOFieldData:
	    CCMIOGetNode(err, node, kFieldDataName, &dataNode);
	    CCMIOGetDataType(err, dataNode, type);
	    break;
	default:
	    *type = kCCMIOInt32;
	    break;
    }
    return(*err);
}

CCMIOError CCMIOWriteState( CCMIOError *err, CCMIOID state,
			    CCMIOID problemDescription, const char *description)
{
    CHECK_ERROR(err);

    if (state.type != kCCMIOState ||
	problemDescription.type != kCCMIOProblemDescription)
	return(*err = kCCMIOBadParameterErr);

    /* Do not call ClearStateProbDef() here, otherwise it is impossible
       overwrite the file in such a way that the file size stays the same */

    CCMIOWriteNodei(err, state.node, kProblemDescriptionName,
		    problemDescription.id);

    if (description)
	CCMIOWriteNodestr(err, state.node, kLabelName, description);
    else
    {
	CCMIONode label;
	/* Node might not exist, but if it does, delete it */
	if (CCMIOGetNode(NULL, state.node, kLabelName, &label) == kCCMIONoErr)
	    CCMIODeleteNode(err, label);
    }

    return(*err);
}

CCMIOError CCMIOWriteMap( CCMIOError *err, CCMIOID id, CCMIOSize n,
			  CCMIOIndex max, int const *data,
                          CCMIOIndex start, CCMIOIndex end )
{
    static CCMIONode node;
    CHECK_ERROR(err);
    if (!data)	return(*err = kCCMIOBadParameterErr);

    if (start == kCCMIOStart)
    {
      CCMIOWriteNodei(err, id.node, kMapMaxName, max);
      CCMIOCreateNode(err, id.node, TRUE, kMapDataName, kMapDataName, &node);
    }
    /* The static declaration of node apparently is to avoid having to lookup
       the created node when writing maps in chunks, but if you interleave calls
       to WriteMap for different maps it will be invalid. For example if you 
       create one map for all cell ids and loop over regions writing chunks 
       of that map for the cells of each region, while also creating and 
       writing per-region maps in that loop, the global cell id map will be
       corrupted. Per-region cell maps are necessary to write solution data
       (eg, velocities only exist in fluid regions).
       This fixes that problem at the cost of re-introducing the node lookup.*/
    else if (node.parent != id.node.node)
    {
      CCMIOGetNode(err, id.node, kMapDataName, &node);
    }

    CCMIOWrite1i(err, node, n, data, start, end);

    return(*err);
}

CCMIOError CCMIOReadMap( CCMIOError *err, CCMIOID id, int *data,
			 CCMIOIndex start, CCMIOIndex end )
{
    CCMIONode node;
    CHECK_ERROR(err);
    if (!data) return(*err = kCCMIOBadParameterErr);

    CCMIOGetNode(err, id.node, kMapDataName, &node);
    CCMIORead1i(err, node, data, start, end);

    return(*err);
}

CCMIOError CCMIOReadVerticesCore( CCMIOError *err, CCMIOID id,
				  int *dims, float *scale, CCMIOID *mapID,
				  CCMIODataType type, void *vertices,
				  CCMIOIndex start, CCMIOIndex end  )
{
    if (dims)	*dims = 0;
    if (scale)	*scale = 0.0;
    if (mapID)	ClearCCMIOID(mapID);

    {
    int idVal;
    CCMIONode node;
    CHECK_ERROR(err);

    if (dims)
	CCMIOReadNodei(err, id.node, kVertexDimensionName, dims);
    if (scale)
	CCMIOReadNodef(err, id.node, kVertexScaleFactorName, scale);
    if (mapID)
    {
	CCMIOReadNodei(err, id.node, kMapName, &idVal);
	CCMIOGetEntity(err, id, kCCMIOMap, idVal, mapID);
    }
    if (vertices)
    {
	CCMIOGetNode(err, id.node, kVertexCoordinatesName, &node);
	if (type == kCCMIOFloat32)
	{
	    if (id.version < 20300)
		CCMIOOldReadf(err, node, 2, FALSE, (float *)vertices, start, end);
	    else
		CCMIORead2f(err, node, (float *)vertices, start, end);
	}
	else if (type == kCCMIOFloat64)
	{
	    if (id.version < 20300)
		CCMIOOldReadd(err, node, 2, FALSE, (double *)vertices, start, end);
	    else
		CCMIORead2d(err, node, (double *)vertices, start, end);
	}
	else
	    return(*err = kCCMIOBadParameterErr);
    }
    return(*err);
    }
}

CCMIOError CCMIOWriteVerticesCore( CCMIOError *err, CCMIOID id,
				   int dims, float scale, CCMIOID mapID,
				   CCMIODataType type, const void *vertices,
				   CCMIOIndex start, CCMIOIndex end )
{
  static CCMIOSize n;
  static CCMIONode node;
  CHECK_ERROR(err);

  if (start == kCCMIOStart)
  {
      CCMIOWriteNodei(err, id.node, kVertexDimensionName, dims);
      CCMIOWriteNodef(err, id.node, kVertexScaleFactorName, scale);
      CCMIOWriteNodei(err, id.node, kMapName, mapID.id);
      CCMIOCreateNode(err, id.node, TRUE, kVertexCoordinatesName,
		      kVertexCoordinatesName, &node);
      CCMIOEntitySize(err, mapID, &n, NULL);
  }
  /* See comment in WriteMap */
  else if (node.parent != id.node.node)
  {
      CCMIOGetNode(err, id.node, kVertexCoordinatesName, &node);
      CCMIOEntitySize(err, mapID, &n, NULL);
  }

  if (type == kCCMIOFloat32)
      CCMIOWrite2f(err, node, dims, n, (float *)vertices, start, end);
  else if (type == kCCMIOFloat64)
      CCMIOWrite2d(err, node, dims, n, (double *)vertices, start, end);
  else
      return(*err = kCCMIOBadParameterErr);
  return(*err);
}

CCMIOError CCMIOReadVerticesf( CCMIOError *err, CCMIOID id, int *dims,
			       float *scale, CCMIOID *mapID, float  *vertices,
			       CCMIOIndex start, CCMIOIndex end )
{
    return(CCMIOReadVerticesCore(err, id, dims, scale, mapID,
				 kCCMIOFloat32, vertices, start, end));
}

CCMIOError CCMIOReadVerticesd( CCMIOError *err, CCMIOID id, int *dims,
			       float *scale, CCMIOID *mapID, double *vertices,
			       CCMIOIndex start, CCMIOIndex end )
{
    return(CCMIOReadVerticesCore(err, id, dims, scale, mapID,
				 kCCMIOFloat64, vertices, start, end));
}
CCMIOError CCMIOWriteVerticesf( CCMIOError *err, CCMIOID id, int dims,
				float scale, CCMIOID mapID,
				const float *vertices,
				CCMIOIndex start, CCMIOIndex end )
{
    return(CCMIOWriteVerticesCore(err, id, dims, scale, mapID,
				  kCCMIOFloat32, vertices, start, end));
}

CCMIOError CCMIOWriteVerticesd( CCMIOError *err, CCMIOID id, int dims,
				float scale, CCMIOID mapID,
				const double *vertices,
				CCMIOIndex start, CCMIOIndex end )
{
    return(CCMIOWriteVerticesCore(err, id, dims, scale, mapID,
				  kCCMIOFloat64, vertices, start, end));
}

CCMIOError CCMIOReadCells( CCMIOError *err, CCMIOID id, CCMIOID *mapID,
			   int *cellTypes, CCMIOIndex start, CCMIOIndex end)
{
    if (mapID)	ClearCCMIOID(mapID);
    {
    CCMIONode node;
    CHECK_ERROR(err);

    if (mapID)
    {
	int idVal;

	CCMIOReadNodei(err, id.node, kMapName, &idVal);
	CCMIOGetEntity(err, id, kCCMIOMap, idVal, mapID);
    }
    if (cellTypes)
    {
	CCMIOGetNode(err, id.node, kCellDataName, &node);
	CCMIORead1i(err, node, cellTypes, start, end);
    }
    return(*err);
    }
}

CCMIOError CCMIOWriteCells( CCMIOError *err, CCMIOID id, CCMIOID mapID,
			    int *cellTypes, CCMIOIndex start,
			    CCMIOIndex end )
{
    static CCMIOSize n;
    static CCMIONode node;
    CHECK_ERROR(err);

    if (start == kCCMIOStart)
    {
	CCMIOEntitySize(err, mapID, &n, NULL);
	CCMIOWriteNodei(err, id.node, kNumCellsName, n);
	CCMIOWriteNodei(err, id.node, kMapName, mapID.id);
	CCMIOCreateNode(err, id.node, TRUE, kCellDataName, kCellDataName,&node);
    }
    /* See comment in WriteMap */
    else if (node.parent != id.node.node)
    {
	CCMIOGetNode(err, id.node, kCellDataName, &node);
	CCMIOEntitySize(err, mapID, &n, NULL);
    }

    CCMIOWrite1i(err, node, n, cellTypes, start, end);

    return(*err);
}

CCMIOError CCMIOReadFaces( CCMIOError *err, CCMIOID entity, CCMIOEntity which,
			   CCMIOID *mapID, CCMIOSize *streamSize,
			   int *vertexStream, CCMIOIndex start,
			    CCMIOIndex end )
{
    if (streamSize)	*streamSize = 0;
    if (mapID)		ClearCCMIOID(mapID);

    {
    CCMIONode node;
    CHECK_ERROR(err);
    if (which != kCCMIOInternalFaces && which != kCCMIOBoundaryFaces)
	return(*err = kCCMIOBadParameterErr);

    if (mapID)
    {
	int mapIDVal;
	CCMIOReadNodei(err, entity.node, kMapName, &mapIDVal);
	CCMIOGetEntity(err, entity, kCCMIOMap, mapIDVal, mapID);
    }
    if (streamSize || vertexStream)
    {
	CCMIOGetNode(err, entity.node, kFacesVertexDataName, &node);
	if (streamSize)
	{
	    int nDims;
	    CCMIOIndex *dims;
	    CCMIOGetDimensions(err, node, &nDims, &dims);
	    *streamSize = dims[0];
	    free(dims);
	}
	if (vertexStream)
	    CCMIORead1i(err, node, vertexStream, start, end);
    }

    return(*err);
    }
}

CCMIOError CCMIOWriteFaces( CCMIOError *err, CCMIOID entity,
                            CCMIOEntity which, CCMIOID mapID,
                            CCMIOSize streamSize, int const *vertexStream,
                            CCMIOIndex start, CCMIOIndex end )
{
    CCMIOSize nFaces;
    static CCMIONode node;
    CHECK_ERROR(err);
    if (!vertexStream || 
	(which != kCCMIOInternalFaces && which != kCCMIOBoundaryFaces))
	return(*err = kCCMIOBadParameterErr);

    if (start == kCCMIOStart)
    {
	CCMIOEntitySize(err, mapID, &nFaces, NULL);
	CCMIOWriteNodei(err, entity.node, kNumFacesName, nFaces);
	CCMIOWriteNodei(err, entity.node, kMapName, mapID.id);
	CCMIOCreateNode(err, entity.node, TRUE, kFacesVertexDataName,
			kFacesVertexDataName, &node);
    }
    /* See comment in WriteMap */
    else if (node.parent != entity.node.node)
    {
	CCMIOGetNode(err, entity.node, kFacesVertexDataName, &node);
    }

    CCMIOWrite1i(err, node, streamSize, vertexStream, start, end);

    return(*err);
}

CCMIOError CCMIOReadFaceCells( CCMIOError *err, CCMIOID entity,
			       CCMIOEntity which, int *cells,
			       CCMIOIndex start, CCMIOIndex end )
{
    CCMIONode node;
    CHECK_ERROR(err);
    if (!cells ||
	(which != kCCMIOInternalFaces && which != kCCMIOBoundaryFaces))
	return(*err = kCCMIOBadParameterErr);

    CCMIOGetNode(err, entity.node, kFacesCellDataName, &node);
    if (which == kCCMIOInternalFaces)
    {
	if (entity.version < 20300)
	    CCMIOOldReadi(err, node, 2, FALSE, cells, start, end);
	else
	    CCMIORead2i(err, node, cells, start, end);
    }
    else
	CCMIORead1i(err, node, cells, start, end);

    return(*err);
}

CCMIOError CCMIOWriteFaceCells( CCMIOError *err, CCMIOID entity,
				CCMIOEntity which, CCMIOID mapID,
				int const *cells, CCMIOIndex start,
				CCMIOIndex end )
{
    static CCMIOSize nFaces;
    static CCMIONode node;
    CHECK_ERROR(err);
    if (!cells ||
	(which != kCCMIOInternalFaces && which != kCCMIOBoundaryFaces))
	return(*err = kCCMIOBadParameterErr);

    if (start == kCCMIOStart)
    {
	CCMIOEntitySize(err, mapID, &nFaces, NULL);
	CCMIOCreateNode(err, entity.node, TRUE, kFacesCellDataName,
			kFacesCellDataName, &node);
    }
    /* See comment in WriteMap */
    else if (node.parent != entity.node.node)
    {
	CCMIOGetNode(err, entity.node, kFacesCellDataName, &node);
	CCMIOEntitySize(err, mapID, &nFaces, NULL);
    }

    if (which == kCCMIOInternalFaces)
	CCMIOWrite2i(err, node, 2, nFaces, cells, start, end);
    else
	CCMIOWrite1i(err, node, nFaces, cells, start, end);

    return(*err);
}

CCMIOError CCMIOWriteProcessor( CCMIOError *err, CCMIOID processor,
				const char *verticesFile, CCMIOID *vertices,
				const char *topologyFile, CCMIOID *topology,
				const char *initialFieldFile, CCMIOID *initialField,
				const char *solutionFile, CCMIOID *solution )
{
    CHECK_ERROR(err);

    if (!verticesFile)		verticesFile = "";
    if (!topologyFile)		topologyFile = "";
    if (!initialFieldFile)	initialFieldFile = "";
    if (!solutionFile)		solutionFile = "";

    if (vertices)
    {
	CCMIOWriteNodestr(err, processor.node, kVerticesFileName, verticesFile);
	CCMIOWriteNodei(err, processor.node, kVerticesIDName, vertices->id);
    }
    if (topology)
    {
	CCMIOWriteNodestr(err, processor.node, kTopologyFileName, topologyFile);
	CCMIOWriteNodei(err, processor.node, kTopologyIDName, topology->id);
    }
    if (initialField)
    {
	CCMIOWriteNodestr(err, processor.node, kInitialFieldFileName,
			initialFieldFile);
	CCMIOWriteNodei(err, processor.node, kInitialFieldIDName,
		      initialField->id);
    }
    if (solution)
    {
	CCMIOWriteNodestr(err, processor.node, kSolutionFileName, solutionFile);
	CCMIOWriteNodei(err, processor.node, kSolutionIDName, solution->id);
    }
    return(*err);
}

CCMIOError ReadProcessorItem( CCMIOError *err, CCMIOID processor, CCMIOID *id,
			      const char *fileNodeName, const char *mapNodeName,
			      CCMIOEntity type )
{
    char *name;
    int idVal;
    CCMIONode root;
    CCMIOID rootID;

    CCMIOReadNodestr(err, processor.node, fileNodeName, &name);
    if (*err == kCCMIONoNodeErr || (*err == kCCMIONoErr && strcmp(name, "") == 0))
    {
	root = processor.root;
	*err = kCCMIONoErr;
    }
    else if (*err != kCCMIONoErr)
	return(*err);
    else
	*err = CCMIOOpen(name, kCCMIORead, &root);
    free(name);
    CCMIOReadNodei(err, processor.node, mapNodeName, &idVal);
    MakeRootCCMIOID(root, &rootID);
    CCMIOGetEntity(err, rootID, type, idVal, id);

    return(*err);
}

CCMIOError CCMIOReadProcessor( CCMIOError *err, CCMIOID processor,
			   CCMIOID *vertices, CCMIOID *topology,
			   CCMIOID *initialField, CCMIOID *solution )
{
    if (vertices)	ClearCCMIOID(vertices);
    if (topology)	ClearCCMIOID(topology);
    if (initialField)	ClearCCMIOID(initialField);
    if (solution)	ClearCCMIOID(solution);

    {
    CHECK_ERROR(err);

    if (vertices)
	ReadProcessorItem(err, processor, vertices, kVerticesFileName,
			  kVerticesIDName, kCCMIOVertices);
    if (topology)
	ReadProcessorItem(err, processor, topology, kTopologyFileName,
			  kTopologyIDName, kCCMIOTopology);
    if (initialField)
	ReadProcessorItem(err, processor, initialField, kInitialFieldFileName,
			  kInitialFieldIDName, kCCMIOFieldSet);
    if (solution)
	ReadProcessorItem(err, processor, solution, kSolutionFileName,
			  kSolutionIDName, kCCMIOFieldSet);
    return(*err);
    }
}

CCMIOError CCMIOReadLagrangianData( CCMIOError *err, CCMIOID lagrangian,
				    CCMIOID *positions, CCMIOID *solution )
{
    if (positions)	ClearCCMIOID(positions);
    if (solution)	ClearCCMIOID(solution);

    {
    CHECK_ERROR(err);

    if (positions)
	ReadProcessorItem(err, lagrangian, positions, kPositionsFileName,
			  kPositionsIDName, kCCMIOVertices);
    if (solution)
	ReadProcessorItem(err, lagrangian, solution, kSolutionFileName,
			  kSolutionIDName, kCCMIOFieldSet);
    return(*err);
    }
}

CCMIOError CCMIOWriteLagrangianData( CCMIOError *err, CCMIOID lagrangian,
				  const char *positionsFile, CCMIOID *positions,
				  const char *solutionFile, CCMIOID *solution )
{
    CHECK_ERROR(err);

    if (!positionsFile)	positionsFile = "";
    if (!solutionFile)		solutionFile = "";

    if (positions)
    {
	CCMIOWriteNodestr(err, lagrangian.node, kPositionsFileName,
			  positionsFile);
	CCMIOWriteNodei(err, lagrangian.node, kPositionsIDName, positions->id);
    }
    if (solution)
    {
	CCMIOWriteNodestr(err, lagrangian.node, kSolutionFileName,
			  solutionFile);
	CCMIOWriteNodei(err, lagrangian.node, kSolutionIDName, solution->id);
    }
    return(*err);
}

int IsProcessorEntityUsed( CCMIOError *err, CCMIOID state, CCMIOID proc,
			   CCMIOID *lagrangian,
			   const char *fileNodeName, const char *idNodeName )
{
    char *name, *name2;
    int i = 0, j, k, idVal, idVal2;
    CCMIOID next, nextState, nextLagrangian;

    if (lagrangian)
    {
	CCMIOReadNodestr(err, lagrangian->node, fileNodeName, &name);
	CCMIOReadNodei(err, lagrangian->node, idNodeName, &idVal);
    }
    else
    {
	CCMIOReadNodestr(err, proc.node, fileNodeName, &name);
	CCMIOReadNodei(err, proc.node, idNodeName, &idVal);
    }
    while (CCMIONextEntity(NULL, state, kCCMIOState, &i, &nextState)
								 == kCCMIONoErr)
    {
	j = 0;
	while (CCMIONextEntity(NULL, nextState, kCCMIOProcessor, &j, &next)
								 == kCCMIONoErr)
	{
	    if (CCMIOAreNodesEqual(proc.node, next.node) &&
		CCMIOAreNodesEqual(proc.root, next.root))
		continue;
	    CCMIOReadNodestr(err, next.node, fileNodeName, &name2);
	    CCMIOReadNodei(err, next.node, idNodeName, &idVal2);
	    if (idVal == idVal2 &&
		((!name && !name2) ||
		 ((name && name2) && strcmp(name, name2) == 0)))
	    {
		free(name);
		free(name2);
		return(1);
	    }
	    free(name2);

	    /* Check in the Lagrangian nodes, too */
	    k = 0;
	    while (CCMIONextEntity(NULL, next, kCCMIOLagrangianData, &k,
				   &nextLagrangian) == kCCMIONoErr)
	    {
		CCMIOReadNodestr(err, nextLagrangian.node, fileNodeName,&name2);
		CCMIOReadNodei(err, nextLagrangian.node, idNodeName, &idVal2);
		if (idVal == idVal2 &&
		    ((!name && !name2) ||
		     ((name && name2) && strcmp(name, name2) ==0)))
		{
		    free(name);
		    free(name2);
		    return(1);
		}
		free(name2);
	    }
	}
    }
 
    free(name);
    return(0);
}

CCMIOError MarkMaps( CCMIOError *err, CCMIOID anyID, CCMIONode top, int unused )
{
    char name[kCCMIOMaxStringLength + 1];
    int i = 0, idVal;
    CCMIONode child, node;
    CCMIOID mapID;
    CHECK_ERROR(err);

    while (CCMIOGetNextChild(NULL, top, &i, &child) == kCCMIONoErr)
    {
	CCMIOGetName(err, child, name);
	if (strcmp(name, kMapName) == 0)
	{
	    CCMIOReadNodei(err, top, name, &idVal);
	    CCMIOGetEntity(err, anyID, kCCMIOMap, idVal, &mapID);
	    if (unused)
		CCMIOCreateNode(err, mapID.node, TRUE, kMapMarkName, kTempName,
			      &node);
	    else
		if (CCMIOGetNode(NULL, mapID.node, kMapMarkName, &node)
								   == kCCMIONoErr)
		    CCMIODeleteNode(err, node);
	}
	else
	    MarkMaps(err, anyID, child, unused);
    }
    return(*err);
}

CCMIOError RemoveUnusedMaps( CCMIOError *err, CCMIONode root )
{
    int i, n;
    CCMIONode node, markNode;
    CCMIOID rootID, id;
    CHECK_ERROR(err);

    /* The maps that are potentially unused have already been marked and their
       nodes deleted (requirement).  Now walk through what it left and unmark
       any maps that are used.  Note that the vertices and the topology are both
       under the same parent, so we can save a call to GetEntityParent. */
    MakeRootCCMIOID(root, &rootID);
    GetEntityParent(err, rootID, kCCMIOVertices, &id);
    MarkMaps(err, rootID, id.node, FALSE);
    GetEntityParent(err, rootID, kCCMIOFieldSet, &id);
    MarkMaps(err, rootID, id.node, FALSE);

    /* Now walk through the maps and delete everything that is marked.
       But go backwards, otherwise i will be wrong. */
    GetEntityParent(err, rootID, kCCMIOMap, &id);
    CCMIOGetNumberOfChildren(err, id.node, &n);
    i = n - 1;
    while (CCMIOGetNextChild(NULL, id.node, &i, &node) == kCCMIONoErr)
    {
	if (CCMIOGetNode(NULL, node, kMapMarkName, &markNode) == kCCMIONoErr)
	{
	    CCMIODeleteNode(err, node);
	    n--;
	}
	i -= 2;	/* CCMIOGetNextChild adds one to i, so need to subtract 2 */
    }
    return(*err);
}

CCMIOError CCMIOClearProcessor( CCMIOError *err, CCMIOID state,
				CCMIOID processor,
				int clearVertices, int clearTopology,
				int clearInitialField, int clearSolution,
				int clearLagrangian )
{
    CCMIOID id;
    CHECK_ERROR(err);

    /* Update the version.  This isn't entirely accurate, but if the processor
       is cleared, mostly likely the whole thing will be changed. */
    CCMIOSetVersion(err, state.root, kCCMIOVersion);
    state.version = kCCMIOVersion;

    if (clearVertices && !IsProcessorEntityUsed(err, state, processor, NULL,
					    kVerticesFileName, kVerticesIDName))
    {
	ReadProcessorItem(err, processor, &id, kVerticesFileName,
			  kVerticesIDName, kCCMIOVertices);
	if (CCMIOAreNodesEqual(state.root, id.root))
	{
	    MarkMaps(err, id, id.node, TRUE);
	    CCMIODeleteEntity(err, id);
	}
	else
	    CCMIOClose(id.root);
    }
    if (*err == kCCMIONoNodeErr)
	*err = kCCMIONoErr;
    if (clearTopology && !IsProcessorEntityUsed(err, state, processor, NULL,
					    kTopologyFileName, kTopologyIDName))
    {
	ReadProcessorItem(err, processor, &id, kTopologyFileName,
			  kTopologyIDName, kCCMIOTopology);
	if (CCMIOAreNodesEqual(state.root, id.root))
	{
	    MarkMaps(err, id, id.node, TRUE);
	    CCMIODeleteEntity(err, id);
	}
	else
	    CCMIOClose(id.root);
    }
    if (*err == kCCMIONoNodeErr)
	*err = kCCMIONoErr;
    if (clearInitialField && !IsProcessorEntityUsed(err, state, processor, NULL,
						    kInitialFieldFileName,
						    kInitialFieldIDName))
    {
	ReadProcessorItem(err, processor, &id, kInitialFieldFileName,
			  kInitialFieldIDName, kCCMIOFieldSet);
	if (CCMIOAreNodesEqual(state.root, id.root))
	{
	    MarkMaps(err, id, id.node, TRUE);
	    CCMIODeleteEntity(err, id);
	}
	else
	    CCMIOClose(id.root);
    }
    if (*err == kCCMIONoNodeErr)
	*err = kCCMIONoErr;
    if (clearSolution && !IsProcessorEntityUsed(err, state, processor, NULL,
					    kSolutionFileName, kSolutionIDName))
    {
	ReadProcessorItem(err, processor, &id, kSolutionFileName,
			  kSolutionIDName, kCCMIOFieldSet);
	if (CCMIOAreNodesEqual(state.root, id.root))
	{
	    MarkMaps(err, id, id.node, TRUE);
	    CCMIODeleteEntity(err, id);
	}
	else
	    CCMIOClose(id.root);
    }
    if (*err == kCCMIONoNodeErr)
	*err = kCCMIONoErr;
    if (clearLagrangian)
    {
	int i = 0;
	CCMIOID lagrangian;
	while (CCMIONextEntity(NULL, processor, kCCMIOLagrangianData, &i,
			       &lagrangian) == kCCMIONoErr)
	{
	    if (!IsProcessorEntityUsed(err, state, processor, &lagrangian,
				       kSolutionFileName, kSolutionIDName))
	    {
		ReadProcessorItem(err, lagrangian, &id, kSolutionFileName,
				  kSolutionIDName, kCCMIOFieldSet);
		if (CCMIOAreNodesEqual(state.root, id.root))
		{
		    MarkMaps(err, id, id.node, TRUE);
		    CCMIODeleteEntity(err, id);
		}
		else
		    CCMIOClose(id.root);
	    }
	    if (!IsProcessorEntityUsed(err, state, processor, &lagrangian,
				       kPositionsFileName, kPositionsIDName))
	    {
		ReadProcessorItem(err, lagrangian, &id, kPositionsFileName,
				  kPositionsIDName, kCCMIOVertices);
		if (CCMIOAreNodesEqual(state.root, id.root))
		{
		    MarkMaps(err, id, id.node, TRUE);
		    CCMIODeleteEntity(err, id);
		}
		else
		    CCMIOClose(id.root);
	    }
	}
    }
    if (*err == kCCMIONoNodeErr)
	*err = kCCMIONoErr;

    RemoveUnusedMaps(err, state.root);
    return(*err);
}

CCMIOError CCMIOWriteMultiDimensionalFieldData( CCMIOError *err,
						CCMIOID fieldID,
						CCMIOComponent component,
						CCMIOID componentField )
{
    char name[kCCMIOMaxStringLength+1], componentName[kCCMIOMaxStringLength+1];
    CHECK_ERROR(err);

    sprintf(name, "%s-%d", kComponentName, (int)component);
    CCMIOGetName(err, componentField.node, componentName);
    CCMIOWriteNodestr(err, fieldID.node, name, componentName);

    return(*err);
}

CCMIOError CCMIOReadMultiDimensionalFieldData( CCMIOError *err,
					       CCMIOID fieldID,
					       CCMIOComponent component,
					       CCMIOID *componentField )
{
    char name[kCCMIOMaxStringLength+1], *componentName;
    CHECK_ERROR(err);
    if (fieldID.version < 20400)
	return(*err = kCCMIOVersionErr);
    if (!componentField)
	return(*err = kCCMIOBadParameterErr);

    sprintf(name, "%s-%d", kComponentName, (int)component);
    CCMIOReadNodestr(err, fieldID.node, name, &componentName);
    ClearCCMIOID(componentField);
    if (componentName)
    {
	CCMIONode parent, node;

	parent.node = fieldID.node.parent;
	parent.parent = 0;
	if (CCMIOGetNode(err, parent, componentName, &node) == kCCMIONoErr)
	{
	    componentField->root = fieldID.root;
	    componentField->node = node;
	    componentField->id = 0;
	    componentField->type = kCCMIOField;
	    componentField->version = fieldID.version;
	}
    }
    free(componentName);

    return(*err);
}


CCMIOError CCMIOWriteFieldDataCore( CCMIOError *err, CCMIOID fieldData,
				    CCMIOID mapID, CCMIODataLocation loc,
				    CCMIODataType type, int isConst, void *data,
				    CCMIOIndex start, CCMIOIndex end )
{
    static CCMIOSize n;
    static CCMIONode node;
    CHECK_ERROR(err);

    if (isConst)
    {
	CCMIOWriteNodei(err, fieldData.node, kFieldDataTypeName, (int)loc);
	CCMIOWriteNodei(err, fieldData.node, kMapName, mapID.id);
	CCMIOEntitySize(err, mapID, &n, NULL);
	if (type == kCCMIOFloat32)
	    CCMIOWriteOptf(err, fieldData, kConstantDataName, *(float *)data);
	else if (type == kCCMIOFloat64)
	    CCMIOWriteOptd(err, fieldData, kConstantDataName, *(double *)data);
	else
	    CCMIOWriteOpti(err, fieldData, kConstantDataName, *(int *)data);
	node.node = 0.0;
	node.parent = 0.0;
    }
    else
    {
	if (start == kCCMIOStart)
	{
	    CCMIOWriteNodei(err, fieldData.node, kFieldDataTypeName, (int)loc);
	    CCMIOWriteNodei(err, fieldData.node, kMapName, mapID.id);
	    CCMIOEntitySize(err, mapID, &n, NULL);
	    CCMIOCreateNode(err, fieldData.node, TRUE, kFieldDataName,
			    kFieldDataName, &node);
	}
	/* See comment in WriteMap */
	else if (node.parent != fieldData.node.node)
	{
	    CCMIOGetNode(err, fieldData.node, kFieldDataName, &node);
	    CCMIOEntitySize(err, mapID, &n, NULL);
	}

	if (type == kCCMIOFloat32)
	    CCMIOWrite1f(err, node, n, (float *)data, start, end);
	else if (type == kCCMIOFloat64)
	    CCMIOWrite1d(err, node, n, (double *)data, start, end);
	else
	    CCMIOWrite1i(err, node, n, (int *)data, start, end);
    }

    return(*err);
}

void memfill(unsigned char *dest, unsigned char *src,
	     int nSrcBytes, int nDestCopies)
{
    int i, j;

    for (i = 0;  i < nDestCopies;  ++i)
    {
	for (j = 0;  j < nSrcBytes;  ++j)
	    *(dest + j) = *(src + j);
	dest += nSrcBytes;
    }
}

CCMIOError CCMIOReadFieldDataCore( CCMIOError *err, CCMIOID fieldData,
				   CCMIOID *mapID, CCMIODataLocation *loc,
				   CCMIODataType type, void *data,
				   CCMIOIndex start, CCMIOIndex end )
{
    int isOld = (fieldData.version < 20300);
    int isOldVector = (fieldData.version < 20400);
    CCMIONode node;
    CCMIOID map;
    CHECK_ERROR(err);

    if (loc)
	CCMIOReadNodei(err, fieldData.node, kFieldDataTypeName, (int *)loc);

    CCMIOReadNodei(err, fieldData.node, kMapName, &map.id);
    CCMIOGetEntity(err, fieldData, kCCMIOMap, map.id, &map);
    if (mapID)
	*mapID = map;

    if (data)
    {
	if (CCMIOGetNode(err, fieldData.node, kFieldDataName, &node)
							     == kCCMIONoNodeErr)
	{
	    int n;
	    CCMIOIndex nItems;

	    *err = kCCMIONoErr;
	    CCMIOEntitySize(err, map, &nItems, NULL);
	    if (end == 0 || end > nItems)
		end = nItems;
	    n = end - start;

	    type = kCCMIOFloat32;
	    CCMIOGetOptInfo(err, fieldData, kConstantDataName, &type,
			    NULL, NULL, NULL);
	    if (type == kCCMIOFloat32)
	    {
		float value;
		CCMIOReadOptf(err, fieldData, kConstantDataName, &value);
		memfill((unsigned char *) data, (unsigned char *) &value,
			CCMIOGetDataTypeSize(type), n);
	    }
	    else if (type == kCCMIOFloat64)
	    {
		double value;
		CCMIOReadOptd(err, fieldData, kConstantDataName, &value);
		memfill((unsigned char *) data, (unsigned char *) &value,
			CCMIOGetDataTypeSize(type), n);
	    }
	    else
	    {
		int value;
		CCMIOReadOpti(err, fieldData, kConstantDataName, &value);
		memfill((unsigned char *) data, (unsigned char *) &value,
			CCMIOGetDataTypeSize(type), n);
	    }
	}
	else if (*err == kCCMIONoNodeErr && fieldData.version >= 20400)
	    return(*err = kCCMIOVersionErr);
	else
	{
	    int nDims = 1;
	    CCMIOIndex *dims;
	    if (isOldVector)
	    {
		CCMIOGetDimensions(err, node, &nDims, &dims);
		free(dims);
	    }
	    switch (nDims)
	    {
		case 1:
		    {
		    if (type == kCCMIOFloat32)
			CCMIORead1f(err, node, (float *)data, start, end);
		    else if (type == kCCMIOFloat64)
			CCMIORead1d(err, node, (double *)data, start, end);
		    else
			CCMIORead1i(err, node, (int *)data, start, end);
		    }
		    break;
		case 2:
		    {
		    if (type == kCCMIOFloat32)
		    {
			if (isOld)
			    CCMIOOldReadf(err, node, 2, FALSE, (float *)data,
					  start, end);
			else
			    CCMIORead2f(err, node, (float *)data, start, end);
		    }
		    else if (type == kCCMIOFloat64)
		    {
			if (isOld)
			    CCMIOOldReadd(err, node, 2, FALSE, (double *)data,
					  start, end);
			else
			    CCMIORead2d(err, node, (double *)data, start, end);
		    }
		    else
		    {
			if (isOld)
			    CCMIOOldReadi(err, node, 2, FALSE, (int *)data,
					  start, end);
			else
			    CCMIORead2i(err, node, (int *)data, start, end);
		    }
		    }
		    break;
		case 3:
		    {
		    if (type == kCCMIOFloat32)
		    {
			if (isOld)
			    CCMIOOldReadf(err, node, 3, FALSE, (float *)data,
					  start, end);
			else
			    CCMIORead3f(err, node, (float *)data, start, end);
		    }
		    else if (type == kCCMIOFloat64)
		    {
			if (isOld)
			    CCMIOOldReadd(err, node, 3, FALSE, (double *)data,
					  start, end);
			else
			    CCMIORead3f(err, node, (float *)data, start, end);
		    }
		    else
		    {
			if (isOld)
			    CCMIOOldReadi(err, node, 3, FALSE, (int *)data,
					  start, end);
			else
			    CCMIORead3f(err, node, (float *)data, start, end);
		    }
		    break;
		    }
	    }
	}
    }
    return(*err);
}

CCMIOError CCMIOWriteFieldDataf( CCMIOError *err, CCMIOID fieldData,
				 CCMIOID mapID, CCMIODataLocation loc,
				 float *data, CCMIOIndex start,
				 CCMIOIndex end )
{
    return(CCMIOWriteFieldDataCore(err, fieldData, mapID, loc, kCCMIOFloat32,
				   FALSE, data, start, end));
}

CCMIOError CCMIOWriteFieldDatad( CCMIOError *err, CCMIOID fieldData,
				 CCMIOID mapID, CCMIODataLocation loc,
				 double *data, CCMIOIndex start,
				 CCMIOIndex end )
{
    return(CCMIOWriteFieldDataCore(err, fieldData, mapID, loc, kCCMIOFloat64,
				   FALSE, data, start, end));
}

CCMIOError CCMIOWriteFieldDatai( CCMIOError *err, CCMIOID fieldData,
				 CCMIOID mapID, CCMIODataLocation loc,
				 int *data, CCMIOIndex start,
				 CCMIOIndex end )
{
    return(CCMIOWriteFieldDataCore(err, fieldData, mapID, loc, kCCMIOInt32,
				   FALSE, data, start, end));
}

CCMIOError CCMIOWriteConstantFieldDataf( CCMIOError *err, CCMIOID fieldData,
					 CCMIOID mapID, CCMIODataLocation loc,
					 float value )
{
    return(CCMIOWriteFieldDataCore(err, fieldData, mapID, loc, kCCMIOFloat32,
				   TRUE, &value, kCCMIOStart, kCCMIOEnd));
}

CCMIOError CCMIOWriteConstantFieldDatad( CCMIOError *err, CCMIOID fieldData,
					 CCMIOID mapID, CCMIODataLocation loc,
					 double value )
{
    return(CCMIOWriteFieldDataCore(err, fieldData, mapID, loc, kCCMIOFloat64,
				   TRUE, &value, kCCMIOStart, kCCMIOEnd));
}

CCMIOError CCMIOWriteConstantFieldDatai( CCMIOError *err, CCMIOID fieldData,
					 CCMIOID mapID, CCMIODataLocation loc,
					 int value )
{
    return(CCMIOWriteFieldDataCore(err, fieldData, mapID, loc, kCCMIOInt32,
				   TRUE, &value, kCCMIOStart, kCCMIOEnd));
}

CCMIOError CCMIOReadFieldDataf( CCMIOError *err, CCMIOID fieldData,
				CCMIOID *mapID, CCMIODataLocation *loc,
				float *data, CCMIOIndex start,
				CCMIOIndex end )
{
    return(CCMIOReadFieldDataCore(err, fieldData, mapID, loc, kCCMIOFloat32,
				  data, start, end));
}

CCMIOError CCMIOReadFieldDatad( CCMIOError *err, CCMIOID fieldData,
				CCMIOID *mapID, CCMIODataLocation *loc,
				double *data, CCMIOIndex start,
				CCMIOIndex end )
{
    return(CCMIOReadFieldDataCore(err, fieldData, mapID, loc, kCCMIOFloat64,
				  data, start, end));
}

CCMIOError CCMIOReadFieldDatai( CCMIOError *err, CCMIOID fieldData,
				CCMIOID *mapID, CCMIODataLocation *loc,
				int *data, CCMIOIndex start,
				CCMIOIndex end )
{
    return(CCMIOReadFieldDataCore(err, fieldData, mapID, loc, kCCMIOInt32,
				  data, start, end));
}

CCMIOError CCMIOWriteRestartInfo( CCMIOError *err, CCMIOID restartInfo,
				  const char *solverName, int iteration,
				  float time, const char *timeUnits,
				  float startAngle )
{
    CHECK_ERROR(err);
    if (!solverName) return(*err = kCCMIOBadParameterErr);
    
    CCMIOWriteOptstr(err, restartInfo, kSolverName, solverName);
    CCMIOWriteOpti(err, restartInfo, kIterationName, iteration);
    CCMIOWriteOptf(err, restartInfo, kTimeName, time);
    if (timeUnits)
	CCMIOWriteOptstr(err, restartInfo, kTimeUnitsName, timeUnits);
    else
	CCMIOWriteOptstr(err, restartInfo, kTimeUnitsName, kDefaultTimeUnits);
    CCMIOWriteOptf(err, restartInfo, kStartAngleName, startAngle);

    return(*err);
}

CCMIOError CCMIOReadRestartInfo( CCMIOError *err, CCMIOID restartInfo,
				 char *solverName, int *iteration,
				 float *time, char *timeUnits,
				 float *startAngle )
{
    int size;
    char *tmpStr;
    CHECK_ERROR(err);

    if (solverName)
    {
	solverName[0] = '\0';
	CCMIOReadOptstr(NULL, restartInfo, kSolverName, &size, NULL);
	tmpStr = (char*) malloc(size + 1);
	CCMIOReadOptstr(NULL, restartInfo, kSolverName, NULL, tmpStr);
	strncpy(solverName, tmpStr, kCCMIOMaxStringLength);
	solverName[kCCMIOMaxStringLength] = '\0';
	free(tmpStr);
    }

    if (iteration && CCMIOReadOpti(err, restartInfo, kIterationName, iteration)
								 != kCCMIONoErr)
	*iteration = 0;

    if (time && CCMIOReadOptf(err, restartInfo, kTimeName, time) != kCCMIONoErr)
	*time = 0.0;

    if (startAngle &&
	CCMIOReadOptf(err, restartInfo, kStartAngleName, startAngle)
								 != kCCMIONoErr)
	*startAngle = 0.0;

    if (timeUnits)
    {
	timeUnits[0] = '\0';
	CCMIOReadOptstr(NULL, restartInfo, kTimeUnitsName, &size, NULL);
	tmpStr = (char*) malloc(size + 1);
	CCMIOReadOptstr(NULL, restartInfo, kTimeUnitsName, NULL, tmpStr);
	strncpy(timeUnits, tmpStr, kCCMIOMaxStringLength);
	timeUnits[kCCMIOMaxStringLength] = '\0';
	free(tmpStr);
    }

    return(*err);
}




/* I changed this routine to include a parameter for the number of faces to
 * speed it up so that there isn't a database read for every buffer write.
 * WRO 22-Oct-2004 */

CCMIOError CCMIOV2WriteFaceCells( CCMIOError *err, CCMIOID entity,
				  CCMIOEntity which,
                                  CCMIOSize nFaces, int const *cells,
                                  CCMIOIndex start, CCMIOIndex end )
{
    static CCMIONode node;
    CHECK_ERROR(err);
    if (!cells ||
	(which != kCCMIOInternalFaces && which != kCCMIOBoundaryFaces))
	return(*err = kCCMIOBadParameterErr);

    if (start == kCCMIOStart)
    {
      CCMIOCreateNode(err, entity.node, TRUE, kFacesCellDataName,
                      kFacesCellDataName, &node);
    }

    if (which == kCCMIOInternalFaces)
	CCMIOWrite2i(err, node, 2, nFaces, cells, start, end);
    else
	CCMIOWrite1i(err, node, nFaces, cells, start, end);

    return(*err);
}



/* Include the CCMIORead*() functions.  This uses #defines to accomplish
   template functions. */
#define TYPE		int
#define TYPE_ABBREV	i
#define CCMIOTYPE	kCCMIOInt32
#include "ccmioread.c"
#undef TYPE		/* Prevents redefinition warnings */
#undef TYPE_ABBREV
#undef CCMIOTYPE

#define TYPE		float
#define TYPE_ABBREV	f
#define CCMIOTYPE	kCCMIOFloat32
#include "ccmioread.c"
#undef TYPE		/* Prevents redefinition warnings */
#undef TYPE_ABBREV
#undef CCMIOTYPE

#define TYPE		double
#define TYPE_ABBREV	d
#define CCMIOTYPE	kCCMIOFloat64
#include "ccmioread.c"

#ifdef __cplusplus
}
#endif
#endif /* CCMIO_C */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
