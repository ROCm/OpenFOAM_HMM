#ifndef CCMIOMESH_VERIFY_C
#define CCMIOMESH_VERIFY_C

/*@@
 *  Program: Star File Format Library  - $RCSfile: ccmioverifymesh.c,v $
 *  Author:  Geoff Prewett
 *  Date:    Sept 18, 2003
 *
 *
 *  Star File Format Library - Copyright (C) 2003 by adapco, Ltd.
 *
 *  This program is the property of adapco, Ltd. and contains
 *  confidential and proprietary information.  The unauthorized use,
 *  distribution, or duplication of this program is prohibited.
 *  All rights reserved.
 *
 *  $Id: ccmioverifymesh.c,v 1.4 2005/01/11 21:51:20 wo Exp $
 */

#ifndef MAKEDEPEND
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>	/* tolower */
#endif

#include "ccmio.h"
#include "ccmiocore.h"
#include "ccmioutility.h"
#include "ccmiomesh.h"
#include "ccmiopost.h"
#include "ccmioverifymesh.h"

typedef struct {
    int valid;
    int scratch;
    char meshID[32];
    union {
	CCMIOMeshVertexf vertex;
	CCMIOMeshCell cell;
	struct {
	    int isInternal;
	    CCMIOMeshFace face;
	    } face;
        } u;
} Item;

typedef struct {
    int size;
    Item *items;
} Data;

#define ERR	stdout
const float kBigCoord = 100000.0;
const float kSmallScale = 1e-4;
const float kBigScale = 1e4;
const char kMaxGlobalIDName[] = "maxGlobal%sID";
const char kNBoundsName[] = "numBoundTypes";
const float kMaxPostVal = 1e13;
const float kMinPostVal = -1e13;

static const char gEntityNames[kCCMIOMeshLast][32] = {
    "kCCMIOMeshVertexData", "kCCMIOMeshCellData",
    "kCCMIOMeshNextBoundaryFaceData", "kCCMIOMeshInternalFaceData" };

static const char gNames[kCCMIOMeshLast][32] =
    { "Vertex", "Cell", "Boundary face", "Internal face" };

static CCMIOError CountMeshEntities( CCMIOError *err, CCMIOProcessor proc,
				     int *nMeshes,
				     int *nVerts, int *maxVertID,
				     int *nCells, int *maxCellID,
				     int *nFaces, int *maxFaceID );
static CCMIOError LoadAndCheckEntity( CCMIOError *err, CCMIOMesh mesh,
				      CCMIOPost post, CCMIOMeshEntity which,
				      Data *globalData, const char *meshID );
static CCMIOError VerifyEntityPostData( CCMIOError *err, CCMIOPost post,
					CCMIOMeshEntity which, int boundaryID,
					int maxID, int *ids, int nExpected,
					const char *meshID );

CCMIOError CCMIOCheckProcessorMesh( CCMIOError *err, CCMIOProcessor proc )
{
    if (err && ((*err) != kCCMIONoErr))
	fprintf(ERR, "Error;  mesh checking not performed\n");

    {
    int i, j, nMeshes, nVerts, nCells, nFaces, nBounds = 0;
    int maxVertID, maxCellID, maxFaceID;
    int iDom = 0, iSDom = 0;
    char domName[kCCMIOMaxStringLength + 1], sdomName[kCCMIOMaxStringLength+1];
    Data vertices, cells, faces;
    CCMIODomain domain;
    CCMIOSubDomain sdom;
    CCMIOMesh mesh;
    CCMIOPost post;
    CHECK_ERROR(err);
    
    fprintf(ERR, "Verifying mesh...\n");

    /* Allocate all the arrays */
    CountMeshEntities(err, proc, &nMeshes, &nVerts, &maxVertID, &nCells,
		      &maxCellID, &nFaces, &maxFaceID);
    vertices.items = (Item *)calloc(sizeof(Item), maxVertID + 1);
    if (!vertices.items)
	return(*err = kCCMIONoMemoryErr);
    vertices.size = maxVertID;
    cells.items = (Item *)calloc(sizeof(Item), maxCellID + 1);
    if (!cells.items)
    {
	free(vertices.items);
	return(*err = kCCMIONoMemoryErr);
    }
    cells.size = maxCellID;
    faces.items = (Item *)calloc(sizeof(Item), maxFaceID + 1);
    if (!faces.items)
    {
	free(vertices.items);
	free(cells.items);
	return(*err = kCCMIONoMemoryErr);
    }
    faces.size = maxFaceID;

    /* Read in the data */
    while (CCMIOProcessorGetNextDomainName(NULL, proc, &iDom, domName)
	   							== kCCMIONoErr)
    {
	CCMIODomainOpen(err, proc, domName, FALSE, &domain);
	iSDom = 0;
	while (CCMIODomainGetNextSubDomainName(NULL, domain, &iSDom,
					     sdomName) == kCCMIONoErr)
	{
	    char *name = sdomName;
	    int dimension;
	    float scale;
/*	    
	    printf("Reading mesh from domain %d, subdomain %d\n", domN, sdomN);
*/
	    CCMIOSubDomainOpen(err, domain, sdomName, FALSE, &sdom);
	    CCMIOMeshOpen(err, sdom, kCCMIORead, &mesh);
	    CCMIOMeshEnable(err, mesh, kCCMIOMeshFloatVertex, TRUE);
	    CCMIOMeshGetDimension(err, mesh, &dimension);
	    if (dimension < 2 || dimension > 3)
		fprintf(ERR, "Mesh %s has invalid dimension (%d)\n", name,
			dimension);
	    CCMIOMeshGetScaleFactor(err, mesh, &scale);
	    if (scale < 0)
		fprintf(ERR, "Mesh %s has negative scale factor!\n", name);
	    else if (scale < kSmallScale || scale > kBigScale)
		fprintf(ERR, "Mesh %s has a suspicious scale factor (%f).\n",
			name, scale);

	    if (CCMIOPostOpen(NULL, sdom, kCCMIORead, &post) != kCCMIONoErr)
		post = NULL;
	    if ((*err) != kCCMIONoErr)
		fprintf(ERR, "Error %d after open and initial check of mesh.\n",
			*err);

	    LoadAndCheckEntity(err, mesh, post, kCCMIOMeshVertexData, &vertices,
			       name);
	    if ((*err) != kCCMIONoErr)
		fprintf(ERR, "Error %d after loading vertices of mesh.\n",*err);

	    LoadAndCheckEntity(err, mesh, post, kCCMIOMeshCellData, &cells, name);
	    if ((*err) != kCCMIONoErr)
		fprintf(ERR, "Error %d after loading cells of mesh.\n",*err);

	    LoadAndCheckEntity(err, mesh, post, kCCMIOMeshInternalFaceData,
			       &faces, name);
	    if ((*err) != kCCMIONoErr)
		fprintf(ERR, "Error %d after loading internal faces of mesh.\n",
			*err);

	    if ((*err) == kCCMIONoErr)
	    {
		int written;

		nBounds = 0;
		while (LoadAndCheckEntity(err, mesh, post,
					  kCCMIOMeshNextBoundaryFaceData,
					  &faces, name) == kCCMIONoErr)
		    nBounds++;
		*err = kCCMIONoErr;
		CCMIOReadNodei(err, mesh->node, kNBoundsName, &written);
		if (written != nBounds)
		    fprintf(ERR, "Mesh %s specifies %d boundaries but actually has %d.\n", name, written, nBounds);
		else if (nBounds == 0)
		    fprintf(ERR, "Mesh %s has no boundary faces.\n", name);
	    }
	    CCMIOMeshClose(err, mesh);
	    CCMIOSubDomainClose(err, sdom);
	    if ((*err) != kCCMIONoErr)
	    {
		fprintf(ERR, "Error reading mesh %s;  some checks not performed.\n", name);
		CCMIODomainClose(err, domain);
		goto error;
	    }
	}
	CCMIODomainClose(err, domain);
    }

    /* Check that all vertices are on some face */
    for (i = 0;  i <= vertices.size;  ++i)
	vertices.items[i].scratch = 0;
    for (i = 0;  i <= faces.size;  ++i)
    {
	if (!faces.items[i].valid)  continue;
	for (j = 0;  j < faces.items[i].u.face.face.nVerts;  ++j)
	{
	    int vertID = faces.items[i].u.face.face.vertices[j];

	    if (vertID > vertices.size)
		fprintf(ERR, "Face %d of mesh %s:  Vertex %d (ID = %d) is out of range\n", faces.items[i].u.face.face.id, faces.items[i].meshID, j, vertID);
	    else if (!vertices.items[vertID].valid)
		fprintf(ERR, "Face %d of mesh %s:  Vertex %d (ID = %d) does not exist\n", faces.items[i].u.face.face.id, faces.items[i].meshID, j, vertID);
	    else
		vertices.items[vertID].scratch++;
	}
    }
    for (i = 0;  i <= vertices.size;  ++i)
    {
	if (vertices.items[i].valid && vertices.items[i].scratch == 0)
	    fprintf(ERR, "Vertex ID %d (mesh %s) is unused.\n", i,
		    vertices.items[i].meshID);
    }

    /* Check that cell has at least four faces */
    for (i = 0;  i <= cells.size;  ++i)
	cells.items[i].scratch = 0;
    for (i = 0;  i <= faces.size;  ++i)
    {
	int cell1, cell2 = 0;
	if (!faces.items[i].valid)  continue;

	if (faces.items[i].u.face.isInternal)
	{
	    cell1 = faces.items[i].u.face.face.u.internal.cell[0];
	    cell2 = faces.items[i].u.face.face.u.internal.cell[1];
	}
	else
	    cell1 = faces.items[i].u.face.face.u.boundary.cell;
	
	if (cell1 < 0 || cell1 > cells.size)
	    fprintf(ERR, "Face %d of mesh %s:  Cell %d (ID = %d) is out of range\n", faces.items[i].u.face.face.id, faces.items[i].meshID, i, cell1);
	else if (!cells.items[cell1].valid)
	    fprintf(ERR, "Face %d of mesh %s:  Cell %d (ID = %d) does not exist\n", faces.items[i].u.face.face.id, faces.items[i].meshID, i, cell1);
	else
	    cells.items[cell1].scratch++;
	if (faces.items[i].u.face.isInternal)
	{
	    if (cell2 < 0 || cell2 > cells.size)
		fprintf(ERR, "Face %d of mesh %s:  Cell %d (ID = %d) is out of range\n", faces.items[i].u.face.face.id, faces.items[i].meshID, i, cell2);
	    else if (!cells.items[cell2].valid)
		fprintf(ERR, "Face %d of mesh %s:  Cell %d (ID = %d) does not exist\n", faces.items[i].u.face.face.id, faces.items[i].meshID, i, cell2);
	    else
		cells.items[cell2].scratch++;
	}
    }
    for (i = 0;  i <= cells.size;  ++i)
    {
	if (cells.items[i].valid && cells.items[i].scratch < 4)
	    fprintf(ERR, "Cell ID %d (mesh %s) has fewer than 4 faces (%d).\n",
		    i, cells.items[i].meshID, cells.items[i].scratch);
    }

 error:
    free(vertices.items);
    free(cells.items);
    free(faces.items);
    if ((*err) != kCCMIONoErr)
	fprintf(ERR, "Error %d occurred during verification.\n", *err);
    return(*err);
    } /* end of CHECK_ERROR scope */
}

CCMIOError LoadAndCheckEntity( CCMIOError *err, CCMIOMesh mesh, CCMIOPost post, 
			       CCMIOMeshEntity which, Data *globalData,
			       const char *meshID )
{
    int i, n, id = 0, max, maxID = 0, *ids = NULL, boundaryID = 0;
    CCMIONode node;
    CCMIODataSet entity;
    CCMIOMeshVertexf vertex;
    CCMIOMeshCell cell;
    CCMIOMeshFace face;

    CHECK_ERROR(err);
    CCMIOMeshReadEntityOpen(err, mesh, which, &n, &max, &entity);
    if ((*err) == kCCMIONoErr && post)
	ids = (int *)calloc(sizeof(int), max + 1);  /* ID list for post check */
    if ((*err) != kCCMIONoErr && which != kCCMIOMeshNextBoundaryFaceData)
    {
	fprintf(ERR, "Mesh %s has no %c%s data.\n", meshID, tolower(gNames[(int)which][0]), gNames[(int)which]+1);
	return(*err);
    }
    for (i = 0;  i < n && (*err) == kCCMIONoErr;  ++i)
    {
	switch(which)
	{
	    case kCCMIOMeshVertexData:
		CCMIOMeshReadNextVertexf(err, entity, &vertex);
		id = vertex.id;
		if (fabsf(vertex.coord[0]) > kBigCoord ||
		    fabsf(vertex.coord[1]) > kBigCoord ||
		    fabsf(vertex.coord[2]) > kBigCoord)
		    fprintf(ERR, "Vertex %d (index %d) in mesh %s is out of range:\n    (%f, %f, %f)\n", vertex.id, vertex.index, meshID, vertex.coord[0], vertex.coord[1], vertex.coord[2]);
		break;
	    case kCCMIOMeshCellData:
		CCMIOMeshReadNextCell(err, entity, &cell);
		id = cell.id;
		break;
	    case kCCMIOMeshInternalFaceData:
	    case kCCMIOMeshNextBoundaryFaceData:
		CCMIOMeshReadNextFace(err, entity, &face);
		id = face.id;
		if (which == kCCMIOMeshNextBoundaryFaceData)
		    boundaryID = face.u.boundary.boundary;
		if (face.nVerts < 3)
		    fprintf(ERR, "%s %d (index %d) in mesh %s has less than three vertices!\n", gNames[(int)which], face.id, face.index, meshID);
		break;
	}
	if (id > globalData->size)
	    fprintf(ERR, "%s %d in mesh %s has an ID greater than maximum ID of all the meshes (%d).\n", gNames[(int)which], id, meshID, globalData->size);
	else if (id == 0)
	    fprintf(ERR, "%s %d in mesh %s has ID of zero.\n", gNames[(int)which], id, meshID);
	else if (id < 0)
	    fprintf(ERR, "%s %d in mesh %s has negative ID.\n", gNames[(int)which], id, meshID);
	else
	{
	    /* The same vertex may be in multiple meshes (but check that the
	       coordinates are the same).  A boundary face may be in multiple
	       meshes if it originally was an internal face, but internal
	       faces and cells can only be in one mesh. */
	    if (globalData->items[id].valid == TRUE)
	    {
		if (which == kCCMIOMeshVertexData &&
		    (vertex.coord[0] != globalData->items[id].u.vertex.coord[0] ||
		     vertex.coord[1] != globalData->items[id].u.vertex.coord[1] ||
		     vertex.coord[2] != globalData->items[id].u.vertex.coord[2]))
		    fprintf(ERR, "Vertex %d (mesh %s) has a different value than defined in mesh %s.\n", id, meshID, globalData->items[id].meshID);
		else if (which == kCCMIOMeshCellData ||
			 which == kCCMIOMeshInternalFaceData)
		    fprintf(ERR, "%s %d in mesh %s is already defined by mesh %s\n", gNames[(int)which], id, meshID, globalData->items[id].meshID);
		else if (which == kCCMIOMeshNextBoundaryFaceData)
		{
		    /* If we are adding a boundary face that already exists,
		       it must be an internal face that was split into two
		       boundary faces.  Merge these back together (after
		       verifying that this is indeed what happened). */
		    int j, equal = 0, cell;
		    if (face.nVerts != globalData->items[id].u.face.face.nVerts)
		    {
			fprintf(ERR, "%s %d in mesh %s is already defined by mesh %s\n", gNames[(int)which], id, meshID, globalData->items[id].meshID);
			fprintf(ERR, "\t(Boundary faces have different number of vertices (%d and %d)).\n", face.nVerts, globalData->items[id].u.face.face.nVerts); 
			continue;
		    }
		    for (j = 0;  j < face.nVerts;  ++j)
		    {
			if (face.vertices[j] == globalData->items[id].u.face.face.vertices[j])
			    equal++;
		    }
		    if (equal != face.nVerts)
		    {
			fprintf(ERR, "%s %d in mesh %s is already defined by mesh %s\n", gNames[(int)which], id, meshID, globalData->items[id].meshID);
			fprintf(ERR, "\t(Boundary faces have different vertices.\n");
			continue;
		    }
		    /* Everything is identical, convert to internal face */
		    cell = globalData->items[id].u.face.face.u.boundary.cell;
		    globalData->items[id].u.face.isInternal = TRUE;
		    globalData->items[id].u.face.face.u.internal.cell[0] = cell;
		    globalData->items[id].u.face.face.u.internal.cell[1] = face.u.boundary.cell;
		    /* Continue, otherwise this will be overwritten */
		    continue;
		}
	    }
	    globalData->items[id].valid = TRUE;
	    strcpy(globalData->items[id].meshID, meshID);
	    if (which == kCCMIOMeshVertexData)
		memcpy(&globalData->items[id].u.vertex,&vertex,sizeof(vertex));
	    else if (which == kCCMIOMeshCellData)
		memcpy(&globalData->items[id].u.cell, &cell, sizeof(cell));
	    else
	    {
		globalData->items[id].u.face.isInternal = (which == kCCMIOMeshInternalFaceData);
		memcpy(&globalData->items[id].u.face.face, &face, sizeof(face));
	    }
	}
	if (id > maxID)
	    maxID = id;
	if (ids)
	    ids[id] = 1;
    }
    if (which == kCCMIOMeshVertexData || which == kCCMIOMeshCellData)
    {
	char name[kCCMIOMaxStringLength + 1];
	CCMIODataSetGetNode(err, entity, &node);
#if kHasSNPrintf
	snprintf(name, kCCMIOMaxStringLength + 1, kMaxGlobalIDName,
		 gNames[(int)which]);
#else
	sprintf(name, kMaxGlobalIDName, gNames[(int)which]);
#endif
	CCMIOReadNodei(err, node, name, &max);
	if (max != maxID && (*err) == kCCMIONoErr)
	    fprintf(ERR, "For %s in mesh %s:\n\t\tspecified maximum in file (%d) is not observed maximum (%d)\n", gEntityNames[(int)which], meshID, max, maxID);
    }
    CCMIOMeshEntityClose(err, entity);
    if ((*err) != kCCMIONoErr && (*err) != kCCMIONoNodeErr)
	fprintf(ERR, "Error reading entity %s.\n", gEntityNames[(int)which]);
    if (post && VerifyEntityPostData(NULL, post, which, boundaryID, maxID, ids,
				     n, meshID) != kCCMIONoErr)
	fprintf(ERR, "Error reading post data of mesh '%s'.\n", meshID);
    free(ids);
	       
    return(*err);
}

CCMIOError CountMeshEntities( CCMIOError *err, CCMIOProcessor proc,int *nMeshes,
			      int *nVerts, int *maxVertID,
			      int *nCells, int *maxCellID,
			      int *nFaces, int *maxFaceID )
{
    int n, maxID, iDom = 0, iSDom = 0;
    char domName[kCCMIOMaxStringLength + 1], sdomName[kCCMIOMaxStringLength + 1];
    CCMIODomain domain;
    CCMIOSubDomain sdom;
    CCMIOMesh mesh;

    if (nMeshes) *nMeshes = 0;
    if (nVerts) *nVerts = 0;
    if (maxVertID) *maxVertID = 0;
    if (nCells) *nCells = 0;
    if (maxCellID) *maxCellID = 0;
    if (nFaces) *nFaces = 0;
    if (maxFaceID) *maxFaceID = 0;

    {
    CHECK_ERROR(err);
    while (CCMIOProcessorGetNextDomainName(err, proc, &iDom, domName)
	   							== kCCMIONoErr)
    {
	CCMIODomainOpen(err, proc, domName, FALSE, &domain);
	iSDom = 0;
	while (CCMIODomainGetNextSubDomainName(err, domain, &iSDom,
						     sdomName) == kCCMIONoErr)
	{
	    CCMIOSubDomainOpen(err, domain, sdomName, FALSE, &sdom);
	    CCMIOMeshOpen(err, sdom, kCCMIORead, &mesh);
	    *nMeshes += 1;

	    CCMIOMeshGetEntitySize(err, mesh, kCCMIOMeshVertexData, &n, &maxID);
	    *nVerts += n;
	    if (maxID > *maxVertID)
		*maxVertID = maxID;

	    CCMIOMeshGetEntitySize(err, mesh, kCCMIOMeshCellData, &n, &maxID);
	    *nCells += n;
	    if (maxID > *maxCellID)
		*maxCellID = maxID;

	    CCMIOMeshGetEntitySize(err, mesh, kCCMIOMeshInternalFaceData,
				 &n, &maxID);
	    *nFaces += n;
	    if (maxID > *maxFaceID)
		*maxFaceID = maxID;

	    CCMIOMeshClose(err, mesh);
	    CCMIOSubDomainClose(err, sdom);
	}
	if ((*err) == kCCMIONoNodeErr)
	    *err = kCCMIONoErr;
	CCMIODomainClose(err, domain);
    }
    }	/* End of CHECK_ERROR scope */
    if ((*err) == kCCMIONoNodeErr)
	*err = kCCMIONoErr;
    
    return(*err);
}

CCMIOError VerifyEntityPostData( CCMIOError *err, CCMIOPost post,
				 CCMIOMeshEntity which,
				 int boundaryID, int maxID, int *ids,
				 int nExpected, const char *meshID )
{
    char name[kCCMIOMaxStringLength + 1];
    int i = 0, j, n = 0, count = 0, nDims, size, ndims, *dims;
    float data[9];
    CCMIONode node;
    CCMIOPostData pdata;
    static char ordinal[][4] =
	{ "st", "nd", "rd", "th", "th", "th", "th", "th", "th", "th" };

    CHECK_ERROR(err);

    while (CCMIOPostGetNextDataName(NULL, post, n++, name) == kCCMIONoErr)
    {
	CCMIOPostReadDataOpen(err, post, which, kCCMIOFloat32, name, boundaryID,
			    &nDims, &pdata);
	if ((*err) == kCCMIONoNodeErr)
	{
	    *err = kCCMIONoErr;
	    continue;
	}
	else if ((*err) != kCCMIONoErr)
	    return(*err);

	CCMIOBufferGetNode(err, pdata->buff, &node);
	CCMIOGetDimensions(err, node, &ndims, &dims);
	if ((*err) == kCCMIONoErr && dims[0] != nExpected)
	{
	    fprintf(ERR, "Post data field '%s' in subdomain '%s' was expected\n\tto have %d elements but has %d.\n", name, meshID, nExpected, dims[0]);
	    free(dims);
	    continue;
	}
	free(dims);

	if (nDims == 1)
	    size = 1;
	else if (nDims == 2)
	    size = 3;
	else if (nDims == 3)
	    size = 9;
	else
	{
	    fprintf(ERR, "Post data field '%s' in subdomain '%s' has unsupported dimensionality %d.\n", name, meshID, nDims);
	    continue;
	}
	count = 1;
	for (i = 1;  i <= maxID && (*err) == kCCMIONoErr;  ++i)
	{
	    if (ids[i])
	    {
		/* Read a data point */
		if (nDims == 1)
		    CCMIOPostDataReadScalarf(err, pdata, i, &data[0]);
		if (nDims == 2)
		{
		    for (j = 0;  j < size;  ++j)
			CCMIOPostDataReadVectorf(err, pdata, i, j, &data[j]);
		}
		if (nDims == 3)
		{
		    for (j = 0;  j < size;  ++j)
			CCMIOPostDataReadTensorf(err, pdata, i, j / 3, j % 3,
					       &data[j]);
		}

		/* Check if it appears to be valid */
		for (j = 0;  j < size;  ++j)
		{
		    if (data[j] > kMaxPostVal || data[j] < kMinPostVal)
			fprintf(ERR, "The %d%s element in post data field '%s' in subdomain '%s' is extreme;\n\tit may not have been written.\n", count, ordinal[(i % 10)], name, meshID);
		}
		count++;
	    }
	}

	CCMIOPostDataClose(err, pdata);
    }

    return(*err);
}

#endif /* CCMIOMESH_VERIFY_C */


/* Automatic setting of emacs local variables. */
/* Local Variables: */
/* mode: C++ */
/* tab-width: 8 */
/* End: */
