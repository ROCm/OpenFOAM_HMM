#include "ccmio.h"
#include "ccmioutility.h"
#include <vector>
#include <iostream>
#include <algorithm>  // for min()

// If you are unfamiliar with the Standard Template Library (STL), please read.
//   vector:	A dynamic array.  vector<int> is an array of ints.  Suppose
//		we have vector<int> iarr;
//		    iarr.push_back( ... ); adds an item to the end of the
//					   array, resizing if necessary.
//		    iarr.size();	   Returns the number of items.
//		    iarr[10];		   The 11th item in the array.
//					   No bounds checking is performed.
//		    iarr.resize(10);	   Makes the array 10 items long.

using namespace std;

struct Face {
    int nVerts;
    int vertices[4];
    int cells[2];
};
static int const nCells = 8;
int cells[] = { 1, 1, 1, 2, 1, 1, 1, 1 };
static int const nInternalFaces = 12;
Face internalFaces[] = { { 4, {2, 5, 14, 11}, {2, 1} },
			 { 4, {4, 5, 14, 13}, {1, 3} },
			 { 4, {5, 6, 15, 14}, {2, 4} },
			 { 4, {5, 8, 17, 14}, {4, 3} },
			 { 4, {10, 11, 14, 13}, {5, 1} },
			 { 4, {11, 12, 15, 14}, {6, 2} },
			 { 4, {11, 14, 23, 20}, {6, 5} },
			 { 4, {13, 14, 17, 16}, {7, 3} },
			 { 4, {13, 14, 23, 22}, {5, 7} },
			 { 4, {14, 15, 18, 17}, {8, 4} },
			 { 4, {14, 15, 24, 23}, {6, 8} },
			 { 4, {14, 17, 26, 23}, {8, 7} } };

static int const nBoundaryFaces = 24;
Face boundaryFaces[] = { { 4, {1, 2, 5, 4}, {1, 0} },
			 { 4, {1, 10, 11, 2}, {1, 0} },
			 { 4, {1, 4, 13, 10}, {1, 0} },
			 { 4, {2, 3, 6, 5}, {2, 0} },
			 { 4, {2, 11, 12, 3}, {2, 0} },
			 { 4, {3, 12, 15, 6}, {2, 0} },
			 { 4, {4, 5, 8, 7}, {3, 0} },
			 { 4, {4, 7, 16, 13}, {3, 0} },
			 { 4, {5, 6, 9, 8}, {4, 0} },
			 { 4, {6, 15, 18, 9}, {4, 0} },
			 { 4, {7, 8, 17, 16}, {3, 0} },
			 { 4, {8, 9, 18, 17}, {4, 0} },
			 { 4, {10, 19, 20, 11}, {5, 0} },
			 { 4, {10, 13, 22, 19}, {5, 0} },
			 { 4, {11, 20, 21, 12}, {6, 0} },
			 { 4, {12, 21, 24, 15}, {6, 0} },
			 { 4, {13, 16, 25, 22}, {7, 0} },
			 { 4, {15, 24, 27, 18}, {8, 0} },
			 { 4, {16, 17, 26, 25}, {7, 0} },
			 { 4, {17, 18, 27, 26}, {8, 0} },
			 { 4, {19, 22, 23, 20}, {5, 0} },
			 { 4, {20, 23, 24, 21}, {6, 0} },
			 { 4, {22, 25, 26, 23}, {7, 0} },
			 { 4, {23, 26, 27, 24}, {8, 0} } };
int nVertices = 27;
float vertices[] = { 0, 0, 0, 0, 0.5, 0, 0, 1, 0, 0.5, 0, 0, 0.5, 0.5, 0, 0.5, 1, 0, 1, 0, 0, 1, 0.5, 0, 1, 1, 0, 0, 0, 0.5, 0, 0.5, 0.5, 0, 1, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 1, 0, 0.5, 1, 0.5, 0.5, 1, 1, 0.5, 0, 0, 1, 0, 0.5, 1, 0, 1, 1, 0.5, 0, 1, 0.5, 0.5, 1, 0.5, 1, 1, 1, 0, 1, 1, 0.5, 1, 1, 1, 1 };
int mapData[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36 };

float postData[] = { 100.0, 100.0, 100.0, 200.0, 100.0, 100.0, 100.0, 100.0 };

static int nDroplets = 8;
float positions[] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,  0, 1, 1,  1, 0, 1,
		      1, 1, 0,  1, 1, 1 };
int types[] = { 1, 1, 2, 1, 2, 1, 2, 2 };

static char const kStateName[] = "default";
static char const kDataName[] = "Temperature";
static char const kShortDataName[] = "Temp";
static char const kUnitsName[] = "Units";
static int const kInc = 10;

bool HasSameDescription(CCMIOID entity, const char *desc);
CCMIOError GetOrCreateEntity(CCMIOError *err, CCMIOID parent, CCMIOEntity which,
			     const char *description, CCMIOID *id);
CCMIOError GetOrCreateIndexedEntity(CCMIOError *err, CCMIOID parent,
				    CCMIOEntity which, int idx,
				    const char *description, CCMIOID *id);
CCMIOError GetOrCreateField(CCMIOError *err, CCMIOID phase, const char *name,
			    const char *description,
			    CCMIODimensionality dim, CCMIOID *field);
CCMIOError WriteMonolithicFieldData(CCMIOError *err, CCMIOID field,
				    CCMIOID mapID, CCMIODataType type,
				    void *postData, void *constData);
CCMIOError WriteLagrangianPositions(CCMIOError *err, CCMIOID processorID,
				    CCMIOID *lagrangian,
				    CCMIOID *verticesID, CCMIOID *solutionID);

int main( int argc, char *argv[] )
{
    int i;
    CCMIOID stateID, processorID, rootID, verticesID, topologyID, solutionID;
    CCMIOID mapID, problem;
    CCMIOID id, cellMapID;
    CCMIOError error = kCCMIONoErr, *err = &error;

    if (argc < 2)
    {
	cout << "Usage:  " << argv[0] << " outputFile" << endl;
	return(0);
    }

    CCMIOOpenFile(err, argv[1], kCCMIOWrite, &rootID);

    // Create a new state (or re-use an existing one).
    if (CCMIOGetState(err, rootID, kStateName, &problem, &stateID) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewState(err, rootID, kStateName, NULL, NULL, &stateID);
    }

    // Attempt to find the first processor.  If we find it, re-use it, otherwise
    // create one.  Multiprocessor meshes will obviously need to do this
    // multiple times...
    // Note that we pass NULL for the error because it is not an error for
    // us if there is not a processor;  it just means that we need to create
    // one.  If we pass in 'err' then all the rest of the functions will
    // needlessly fail because we encountered an expected error.
    i = 0;
    while (CCMIONextEntity(NULL, stateID, kCCMIOProcessor, &i, &processorID)
	   							== kCCMIONoErr)
	if (HasSameDescription(processorID, "CPU 1"))
	    break;
    
    if (CCMIOIsValidEntity(processorID))  // Won't be if we didn't find one
	CCMIOReadProcessor(err, processorID, &verticesID, &topologyID, NULL,
			   &solutionID);
    else
    {
	CCMIONewEntity(err, stateID, kCCMIOProcessor, "CPU 1", &processorID);
	// Initialize the processor's sub-nodes to be invalid;  if they
	// are unitialized, checks below may fail.
	CCMIOInvalidateEntity(&topologyID);
	CCMIOInvalidateEntity(&verticesID);
	CCMIOInvalidateEntity(&solutionID);

	CCMIONewEntity(err, rootID, kCCMIOTopology, "Topology", &topologyID);
    }

    // Write the vertices.  First we will need to get (or create) the mapping
    // from the index in the data array to the vertex ID.  Then we can write
    // the actual data.
    if (CCMIOIsValidEntity(verticesID))
	CCMIOReadVerticesf(err, verticesID, NULL, NULL, &mapID, NULL, 0, 0);
    else
    {
	CCMIONewEntity(err, rootID, kCCMIOMap, "Vertex map", &mapID);
	CCMIONewEntity(err, rootID, kCCMIOVertices, "Vertices", &verticesID);
    }

    CCMIOWriteMap(err, mapID, nVertices, mapData[nVertices], mapData,
		  kCCMIOStart, kCCMIOEnd);
    // Write the vertices piecemeal by way of illustration
    for (i = 0;  i < nVertices;  i += kInc)
	CCMIOWriteVerticesf(err, verticesID, 3, 1.0, mapID, vertices + 3 * i,
			    i, i + kInc);

    // Write the cells
    if (CCMIOGetEntity(err, topologyID, kCCMIOCells, 0, &id) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewEntity(err, topologyID, kCCMIOCells, "Cells", &id);
	CCMIONewEntity(err, rootID, kCCMIOMap, "Cell map", &cellMapID);
    }
    else
	CCMIOReadCells(err, id, &cellMapID, NULL, 0, 0);

    CCMIOWriteMap(err, cellMapID, nCells, mapData[nVertices], mapData,
		  kCCMIOStart, kCCMIOEnd);
    for (i = 0;  i < nCells;  i += 3)
	CCMIOWriteCells(err, id, cellMapID, cells + i, i, i + 3);

    // Write the faces.  There are two kinds of faces, internal faces and
    // boundary faces.  The procedure is the same for each except that
    // internal faces' cell array is a two dimensional array because internal
    // faces have a cell on both sides.  Boundary faces, obviously, only
    // have one cell.
    vector<int> v, c;

    if (CCMIOGetEntity(err, topologyID, kCCMIOInternalFaces, 0, &id) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewEntity(err, rootID, kCCMIOMap, NULL, &mapID);
	CCMIONewEntity(err, topologyID, kCCMIOInternalFaces, "Internal faces",
		       &id);
    }
    else
	CCMIOReadFaces(err, id, kCCMIOInternalFaces, &mapID, NULL, NULL, 0, 0);

    CCMIOWriteMap(err, mapID, nInternalFaces,
		mapData[nBoundaryFaces + nInternalFaces],
		&mapData[nBoundaryFaces], kCCMIOStart, kCCMIOEnd);
    for (int f = 0;  f < nInternalFaces;  ++f)
    {
	v.push_back(internalFaces[f].nVerts);
	for (int i = 0;  i < internalFaces[f].nVerts;  ++i)
	    v.push_back(internalFaces[f].vertices[i]);
	c.push_back(internalFaces[f].cells[0]);
	c.push_back(internalFaces[f].cells[1]);  // Our other cell
    }
    CCMIOWriteFaces(err, id, kCCMIOInternalFaces, mapID, v.size(), &v[0],
		    kCCMIOStart, kCCMIOEnd);
    CCMIOWriteFaceCells(err, id, kCCMIOInternalFaces, mapID, &c[0], kCCMIOStart,
			kCCMIOEnd);
	
    // Write the boundary faces
    if (CCMIOGetEntity(err, topologyID, kCCMIOBoundaryFaces, 0, &id) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewEntity(err, rootID, kCCMIOMap, NULL, &mapID);
	CCMIONewIndexedEntity(err, topologyID, kCCMIOBoundaryFaces, 0,
			      "Boundary faces", &id);
    }
    else
	CCMIOReadFaces(err, id, kCCMIOBoundaryFaces, &mapID, NULL, NULL, 0, 0);

    v.clear();
    c.clear();
    CCMIOWriteMap(err, mapID, nBoundaryFaces, mapData[nBoundaryFaces], mapData,
		  kCCMIOStart, kCCMIOEnd);
    for (int f = 0;  f < nBoundaryFaces;  ++f)
    {
	v.push_back(boundaryFaces[f].nVerts);
	for (int i = 0;  i < boundaryFaces[f].nVerts;  ++i)
	    v.push_back(boundaryFaces[f].vertices[i]);
	c.push_back(boundaryFaces[f].cells[0]);
	// Boundary faces only have one cell, so cells[1] is unused.
    }
    CCMIOWriteFaces(err, id, kCCMIOBoundaryFaces, mapID, v.size(), &v[0],
		    kCCMIOStart, kCCMIOEnd);
    CCMIOWriteFaceCells(err, id, kCCMIOBoundaryFaces, mapID, &c[0], kCCMIOStart,
			kCCMIOEnd);

    // Write any prostar sets that we are using (this is optional)
    int set1[] = { 1, 2, 3, 4 }, set2[] = {2, 4, 6, 8};
    int vset[] = { 2, 3, 5, 7, 11, 13,15, 17, 19, 23};
    char shortName[] = "half";
    char longName[] = "Supercalifragilisticexpialidocious";
    char truncatedName[] = "SetWithAReallyLongNameThatWillGetTruncated";
    vector<string> setNames;
    CCMIOID setID;
    if (CCMIOGetProstarSet(err, id, shortName, NULL, NULL, NULL, &setID) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewProstarSet(err, id, shortName, shortName, &setID);
    }
    CCMIOWriteOpt1i(err, setID, "CellSet", 4, set1, kCCMIOStart, kCCMIOEnd);
    CCMIOWriteOpt1i(err, setID, "VertexSet", 10, vset, kCCMIOStart, kCCMIOEnd);
    setNames.push_back(shortName);
    // The name of the next set is more than 32 characters long.  
    if (CCMIOGetProstarSet(err, id, longName, NULL, NULL, NULL, &setID) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewProstarSet(err, id, longName, longName, &setID);
    }
    CCMIOWriteOpt1i(err, setID, "SplineSet", 4, set2, kCCMIOStart, kCCMIOEnd);
    setNames.push_back(longName);
    // But we don't need to specify a long name if we don't want to.
    if (CCMIOGetProstarSet(err, id, truncatedName, NULL, NULL, NULL, &setID) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewProstarSet(err, id, truncatedName, NULL, &setID);
    }
    CCMIOWriteOpt1i(err, setID, "BoundarySet", 4, set2, kCCMIOStart, kCCMIOEnd);
    setNames.push_back(truncatedName);

    // Now make a string for the monitoring set node later on
    vector<char> setStr;
    vector<string>::iterator sit;
    for (sit = setNames.begin();  sit != setNames.end();  ++sit)
    {
	int i = 0, size = sit->size();
	if (size == 0)  // A size of zero will be a premature end-of-string
	    continue;

	setStr.push_back(min(kCCMIOMaxStringLength, size));
	while (i < size && i < kCCMIOMaxStringLength)
	    setStr.push_back((*sit)[i++]);
    }
    setStr.push_back('\0');  // End the C string

    // Write out a dummy problem description.  If we happen to know that
    // there already is a problem description previously recorded that
    // is valid we could skip this step.
    CCMIOID constants;

    if (!CCMIOIsValidEntity(problem))
    {
	*err = kCCMIONoErr;
	CCMIONewEntity(err, rootID, kCCMIOProblemDescription,
		       "Dummy description", &problem);
    }
    GetOrCreateEntity(err, problem, kCCMIOModelConstants, "Constant values",
		      &constants);
    GetOrCreateIndexedEntity(err, problem, kCCMIOCellType, 1, "Dummy celltypes",
			     &id);
    CCMIOWriteOptstr(err, id, "MaterialType", "solid");
    GetOrCreateIndexedEntity(err, problem, kCCMIOCellType, 2, "Dummy celltypes",
			     &id);
    CCMIOWriteOptstr(err, id, "MaterialType", "solid");
    
    CCMIOWriteOptf(err, constants, "Gravity", 9.82);
    CCMIOWriteOptf(err, constants, "B.P. of water", 373);
    CCMIOWriteOptstr(err, problem, "MonitoringSets", &setStr[0]);

    // We have a problem description recorded but our state does not know
    // about it.  So tell the state that it has a problem description.
    CCMIOWriteState(err, stateID, problem, "Example state");

    // Write out some simple solution data
    float zero = 0.0, one = 1.0, half = 0.5;
    CCMIOID phase, field;
    CCMIODimensionality dim;
    if (!CCMIOIsValidEntity(solutionID))
	CCMIONewEntity(err, rootID, kCCMIOFieldSet, "Dummy post data",
		       &solutionID);
    GetOrCreateIndexedEntity(err, solutionID, kCCMIOFieldPhase, 0, NULL,&phase);
    if (CCMIOGetField(err, phase, kDataName, &dim, &field) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewField(err, phase, kDataName, kShortDataName, kCCMIOScalar, &field);
    }
    CCMIOWriteOptstr(err, field, kUnitsName, "°F");
    WriteMonolithicFieldData(err, field, cellMapID, kCCMIOFloat32,
			     (void *)postData, NULL);

    // Write out a piece of vector data (which also illustrates constant data)
    CCMIOID vectorField;
    GetOrCreateField(err, phase, "Velocity", "VELO", kCCMIOVector, &vectorField);

    GetOrCreateField(err, phase, "Velocity (U)", "U", kCCMIOScalar, &field);
    CCMIOWriteOptstr(err, field, kUnitsName, "m/s");
    WriteMonolithicFieldData(err, field, cellMapID, kCCMIOFloat32,
			     NULL, (void *)&zero);
    CCMIOWriteMultiDimensionalFieldData(err, vectorField, kCCMIOVectorX, field);

    GetOrCreateField(err, phase, "Velocity (V)", "V", kCCMIOScalar, &field);
    CCMIOWriteOptstr(err, field, kUnitsName, "m/s");
    WriteMonolithicFieldData(err, field, cellMapID, kCCMIOFloat32,
			     NULL, (void *)&one);
    CCMIOWriteMultiDimensionalFieldData(err, vectorField, kCCMIOVectorY, field);

    GetOrCreateField(err, phase, "Velocity (W)", "W", kCCMIOScalar, &field);
    CCMIOWriteOptstr(err, field, kUnitsName, "m/s");
    WriteMonolithicFieldData(err, field, cellMapID, kCCMIOFloat32,
			     NULL, (void *)&half);
    CCMIOWriteMultiDimensionalFieldData(err, vectorField, kCCMIOVectorZ, field);

    // Add in some sample restart info
    int iteration = 0;
    float time = 0.0, startAngle = 0.0;
    CCMIOID restart, restartData;
    GetOrCreateEntity(err, solutionID, kCCMIORestart, NULL, &restart);
    CCMIOWriteRestartInfo(err, restart, "writeexample", iteration, time, NULL,
			  startAngle);
    // The solver part can have anything the solver wants
    GetOrCreateEntity(err, restart, kCCMIORestartData, NULL, &restartData);
    CCMIOWriteOptf(err, restartData, "Convergence", 100);

    // Now we have the mesh (vertices and topology) and the post data written.
    // Since we now have their IDs, we can write out the processor information.
    CCMIOWriteProcessor(err, processorID, NULL, &verticesID, NULL, &topologyID,
		      NULL, NULL, NULL, &solutionID);

    // Write out some simple Lagrangian data
    int nLagrangians = 0;
    CCMIOID tmp, lagrangian;
    i = 0;
    while (CCMIONextEntity(NULL, processorID, kCCMIOLagrangianData, &i, &tmp) == kCCMIONoErr)
    {
	if (nLagrangians > 0)
	    printf("Too many Langrangians in this file; only expecting one.\n");
	else
	{
	    lagrangian = tmp;
	    WriteLagrangianPositions(err, processorID, &lagrangian,
				     &verticesID, &solutionID);
	}
	nLagrangians++;
    }
    if (nLagrangians == 0)
    {
	lagrangian = tmp;  // tmp was set to invalid by CCMIONextEntity
	WriteLagrangianPositions(err, processorID, &lagrangian,
				 &verticesID, &solutionID);
    }

    GetOrCreateIndexedEntity(err, solutionID, kCCMIOFieldPhase, 0, NULL,&phase);
    GetOrCreateField(err, phase, "Droplet type", "droplet", kCCMIOScalar,
		     &field);
    CCMIOWriteOptstr(err, field, kUnitsName, "none");
    WriteMonolithicFieldData(err, field, cellMapID, kCCMIOInt32,
			     (void *)types, NULL);

    CCMIOWriteLagrangianData(err, lagrangian, NULL, &verticesID,
			     NULL, &solutionID);

    // We're done!
    CCMIOCloseFile(err, rootID);

    if (*err != kCCMIONoErr)
    {
	cout << "Error " << *err << " writing mesh" << endl;
	return(0);
    }

    return(1);
}

bool HasSameDescription(CCMIOID entity, const char *desc)
{
    unsigned int size;
    if (CCMIOEntityDescription(NULL, entity, &size, NULL) != kCCMIONoErr ||
	size <= 0)
	return(false);

    bool eq = false;
    char *str = new char[size + 1];
    if (CCMIOEntityDescription(NULL, entity, &size, str) == kCCMIONoErr &&
	strcmp(desc, str) == 0)
	    eq = true;

    delete [] str;

    return(eq);
}

CCMIOError GetOrCreateEntity(CCMIOError *err, CCMIOID parent, CCMIOEntity which,
			     const char *description, CCMIOID *id)
{
    if (CCMIOGetEntity(err, parent, which, 0, id) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewEntity(err, parent, which, description, id);
    }
    return(*err);
}

CCMIOError GetOrCreateIndexedEntity(CCMIOError *err, CCMIOID parent,
				    CCMIOEntity which, int idx,
				    const char *description, CCMIOID *id)
{
    if (CCMIOGetEntity(err, parent, which, idx, id) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewIndexedEntity(err, parent, which, idx, description, id);
    }
    return(*err);
}

CCMIOError GetOrCreateField(CCMIOError *err, CCMIOID phase, const char *name,
			    const char *description,
			    CCMIODimensionality dim, CCMIOID *field)
{
    if (CCMIOGetField(err, phase, name, &dim, field) == kCCMIONoNodeErr)
    {
	*err = kCCMIONoErr;
	CCMIONewField(err, phase, name, description, kCCMIOScalar,
		      field);
    }
    return(*err);
}

CCMIOError WriteMonolithicFieldData(CCMIOError *err, CCMIOID field,
				    CCMIOID mapID, CCMIODataType type,
				    void *postData, void *constData)
{
    int n = 0, i = 0;
    CCMIOID id, tmp;

    // We need to be very careful here, because field data nodes are numbered,
    // but we don't know which number is which.  CCMIONextEntity() will traverse
    // them in the order that they were created.
    if (CCMIONextEntity(NULL, field, kCCMIOFieldData, &i, &id) == kCCMIONoNodeErr)
    {
	CCMIONewEntity(err, field, kCCMIOFieldData, NULL, &id);
	// Since we have a new node, the do-while loop will find it, so 
	// send it past it.
	while (CCMIONextEntity(NULL, field, kCCMIOFieldData, &i, &tmp) == kCCMIONoErr)
	    ;
    }
    do
    {
	if (n == 0)
	{
	    if (postData)
	    {
		if (type == kCCMIOFloat32)
		    CCMIOWriteFieldDataf(err, id, mapID, kCCMIOCell,
					 (float *)postData,
					 kCCMIOStart, kCCMIOEnd);
		else if (type == kCCMIOFloat64)
		    CCMIOWriteFieldDatad(err, id, mapID, kCCMIOCell,
					 (double *)postData,
					 kCCMIOStart, kCCMIOEnd);
		else if (type == kCCMIOInt32)
		    CCMIOWriteFieldDatai(err, id, mapID, kCCMIOCell,
					 (int *)postData,
					 kCCMIOStart, kCCMIOEnd);
		else
		    cout << "Bad type" << endl;
	    }
	    else
	    {
		if (type == kCCMIOFloat32)
		    CCMIOWriteConstantFieldDataf(err, id, mapID, kCCMIOCell,
						 *(float *)constData);
		else if (type == kCCMIOFloat64)
		    CCMIOWriteConstantFieldDatad(err, id, mapID, kCCMIOCell,
						 *(double *)constData);
		else if (type == kCCMIOInt32)
		    CCMIOWriteConstantFieldDatai(err, id, mapID, kCCMIOCell,
						 *(int *)constData);
		else
		    cout << "Bad type" << endl;
	    }
	}
	else
	    cout << "Too many field entities;  only expected one!" << endl;
	n++;
    } while(CCMIONextEntity(NULL, field, kCCMIOFieldData, &i, &id)
								== kCCMIONoErr);

    return(*err);
}


CCMIOError WriteLagrangianPositions(CCMIOError *err, CCMIOID processorID,
				    CCMIOID *lagrangian,
				    CCMIOID *verticesID, CCMIOID *solutionID)
{
    CCMIOID mapID;

    if (!CCMIOIsValidEntity(*lagrangian))
	GetOrCreateEntity(err, processorID, kCCMIOLagrangianData, NULL,
			  lagrangian);
    CCMIOReadLagrangianData(NULL, *lagrangian, verticesID, solutionID);
    if (!CCMIOIsValidEntity(*verticesID)) {
	CCMIONewEntity(err, *lagrangian, kCCMIOVertices, "Positions",
		       verticesID);
	CCMIONewEntity(err, *lagrangian, kCCMIOMap, "Positions map", &mapID);
    }
    else
	CCMIOReadVerticesf(err, *verticesID, NULL, NULL, &mapID, NULL, 0, 0);
    if (!CCMIOIsValidEntity(*solutionID))
	CCMIONewEntity(err, *lagrangian, kCCMIOFieldSet,
		       "Dummy Lagrangian post data", solutionID);

    CCMIOWriteMap(err, mapID, nDroplets, mapData[nDroplets], mapData,
		  kCCMIOStart, kCCMIOEnd);
    CCMIOWriteVerticesf(err, *verticesID, 3, 1.0, mapID, positions,
			kCCMIOStart, kCCMIOEnd);

    return(*err);
}
