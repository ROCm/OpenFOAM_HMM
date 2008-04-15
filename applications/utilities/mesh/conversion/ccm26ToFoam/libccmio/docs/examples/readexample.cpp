#include <stdlib.h>	// For exit()
#include <ccmio.h>
#include <vector>
#include <map>
#include <iostream>

typedef std::vector< unsigned int > UIntArray;
typedef std::map< int, std::vector< int > > FaceMap;

// If you are unfamiliar with the Standard Template Library (STL), please read.
//   vector:	A dynamic array.  vector<int> is an array of ints.  Suppose
//		we have vector<int> iarr;
//		    iarr.push_back( ... ); adds an item to the end of the
//					   array, resizing if necessary.
//		    iarr.size();	   Returns the number of items.
//		    iarr[10];		   The 11th item in the array.
//					   No bounds checking is performed.
//		    iarr.resize(10);	   Makes the array 10 items long.
//   map:	An associative container.   map< int, float > is an sorted
//		based on an integer key and contains a value float.
//		Suppose we have map<int, float> m;
//		    m[18] = 3.141592;      Associates Pi with the key 18.
//					   Technically m[18] returns a reference
//					   to the value associated with 18,
//					   which is then assigned to.
//					   (A reference is like a pointer that
//					   cannot change where it is pointed to
//					   and is treated syntactically like a
//					   normal variable).
//   iterators: The way to increment through the container.
//		    container.begin()	   Returns an iterator to the first item
//		    container.end()	   Returns an iterator just past the end
//		    *iter		   The item being iterated to.

enum DataType { kScalar, kVector, kVertex, kCell, kInternalFace, kBoundaryFace,
		kBoundaryData, kBoundaryFaceData, kCellType };

using namespace std;

static int const kNValues = 10;	// Number of values of each element to print
static char const kDefaultState[] = "default";
static char const kUnitsName[] = "Units";
static int const kVertOffset = 2;
static int const kCellInc = 4;

static void ReadMesh( CCMIOError &err, CCMIOID &vertices, CCMIOID &topology );
static void ReadVertices( CCMIOError &err, CCMIOID &vertices,
			  char const *counter, int offset = 0 );
static void ReadSolverInfo( CCMIOError &err, CCMIOID &solution);
static void ReadPost( CCMIOError &err, CCMIOID &solution );
static void ReadScalar( CCMIOError &err, CCMIOID field, vector<int> &mapData,
			vector<float> &data, bool readingVector = false );
static void ReadSets( CCMIOError &err, CCMIOID &problem );
static void CheckError( CCMIOError const &err, char const *str );
static void PrintData( int n, int id, DataType type, void *data,
		       void *data2 = NULL );

int main( int argc, char *argv[])
{
    bool hasSolution = true;
    int i = 0;
    CCMIOID root, state, processor, vertices, topology, solution, problem, next;
    CCMIOError err;

    // Open the file.  Because we did not initialize 'err' we need to pass
    // in NULL (which always means kCCMIONoErr) and then assign the return
    // value to 'err'.).
    err = CCMIOOpenFile(NULL, argv[1], kCCMIORead, &root);

    // We are going to assume that we have a state with a known name.
    // We could instead use CCMIONextEntity() to walk through all the states in
    // the file and present the list to the user for selection.
    CCMIOGetState(&err, root, kDefaultState, &problem, &state);
    if (err != kCCMIONoErr)
    {
	cout << "No state named '" << kDefaultState << "'" << endl;
	return(0);
    }
    
    unsigned int size;
    CCMIOEntityDescription(&err, state, &size, NULL);
    char *desc = new char[size + 1];
    CCMIOEntityDescription(&err, state, NULL, desc);
    cout << "Reading state '" << kDefaultState << "' (" << desc << ")" << endl;
    delete [] desc;

    // Find the first processor (i has previously been initialized to 0) and
    // read the mesh and solution information.
    CCMIONextEntity(&err, state, kCCMIOProcessor, &i, &processor);
    CCMIOReadProcessor(&err, processor, &vertices, &topology, NULL, &solution);

    if (err != kCCMIONoErr)
    {
	// Maybe no solution;  try again
	err = kCCMIONoErr;
	CCMIOReadProcessor(&err, processor, &vertices, &topology, NULL, NULL);
	if (err != kCCMIONoErr)
	{
	    cout << "Could not read the file." << endl;
	    return(0);
	}
	hasSolution = false;
    }

    // Print out the problem description (just cell types, and model constants)
    i = 0;
    if (CCMIOIsValidEntity(problem))   // if we have a problem description...
    {
	// ... walk through each cell type and print it...
	while (CCMIONextEntity(NULL, problem, kCCMIOCellType, &i, &next)
	       							 == kCCMIONoErr)
	{
	    char *name;
	    int size, cellType;

	    // ... if it has a material type.  (Note that we do not pass in
	    // an array to get the name because we do not know how long the
	    // string is yet.  Many parameters to CCMIO functions that return
	    // data can be NULL if that data is not needed.)
	    if (CCMIOReadOptstr(NULL, next, "MaterialType", &size, NULL)
								   == kCCMIONoErr)
	    {
		name = new char[size + 1];
		CCMIOReadOptstr(&err, next, "MaterialType", &size, name);
		CCMIOGetEntityIndex(&err, next, &cellType);
		PrintData(i, cellType, kCellType, name);
		delete [] name;
	    }
	}

	float value;
	CCMIOID constants;

	// Read some constants of particular interest, fall back to defaults
	// if they don't exist (or at least, that is what we would do if we
	// actually used them).  If you want to iterate over all the constants
	// instead, you will need to use CCMIOGetEntityNode() and use
	// CCMIOGetNextChild(), and CCMIOReadData() in the Core API.
	if (CCMIOGetEntity(NULL, problem, kCCMIOModelConstants, 0, &constants)
	    							 == kCCMIONoErr)
	{
	    if (CCMIOReadOptf(NULL, constants, "Gravity", &value) == kCCMIONoErr)
		cout << "Gravity is " << value << " m/s^2" << endl;
	    else
		cout << "Gravity is unspecified" << endl;

	    if (CCMIOReadOptf(NULL, constants, "Hubble constant", &value) == kCCMIONoErr)
		cout << "Hubble constant is " << value << endl;
	    else
		cout << "Hubble constant is unspecified" << endl;
	}

	// Read in the prostar sets.  (This is probably not necessary for
	// most applications, but included to show proper parsing.)
	ReadSets(err, problem);
    }
    ReadMesh(err, vertices, topology);
    if (hasSolution)
	{
	ReadSolverInfo(err, solution);
	ReadPost(err, solution);
	}
    else
	cout << "No post data." << endl;

    // Read Lagrangian data
    vector<CCMIOID> ids;
    CCMIOID lagrangian;
    i = 0;
    while (CCMIONextEntity(&err, processor, kCCMIOLagrangianData, &i,
			   &lagrangian) == kCCMIONoErr)
    {
	CCMIOID solutionID, positionsID;
	CCMIOReadLagrangianData(&err, lagrangian, &positionsID, &solutionID);

	CCMIOEntityDescription(&err, lagrangian, &size, NULL);
	char *label = new char[size + 1];
	CCMIOEntityDescription(&err, lagrangian, NULL, label);
	cout << "Lagrangian data '" << label << "'" << endl;
	delete [] label;

	ReadVertices(err, positionsID, "positions");
	ReadPost(err, solutionID);
	ids.push_back(positionsID);
	ids.push_back(solutionID);
    }

    for (vector<CCMIOID>::iterator it = ids.begin();  it != ids.end();  ++it)
	CCMIOCloseFile(&err, *it);

    CCMIOCloseFile(&err, vertices);
    CCMIOCloseFile(&err, topology);
    CCMIOCloseFile(&err, solution);
    CCMIOCloseFile(&err, root);
}

void ReadMesh( CCMIOError &err, CCMIOID &vertices, CCMIOID &topology )
{
    unsigned int i, j, nCells, nFaces, size;
    CCMIOID mapID, id;
    vector<int> mapData, faces, cells, faceCells;
    vector<float> verts;

    ReadVertices(err, vertices, "vertices", kVertOffset);

    // Read the cells.  Store cell IDs so that we know what cells are in
    // this post data.
    CCMIOGetEntity(&err, topology, kCCMIOCells, 0, &id);
    CCMIOEntitySize(&err, id, &nCells, NULL);
    cells.resize(nCells);
    mapData.resize(nCells);
    for (unsigned int k = 0;  k < nCells;  k += kCellInc)
    {
	CCMIOReadCells(&err, id, &mapID, &cells[k], k, k + kCellInc);
	CCMIOReadMap(&err, mapID, &mapData[k], k, k + kCellInc);
    }

    cout << "\t" << nCells << " cells" << endl;
    for (i = 0;  i < nCells && err == kCCMIONoErr;  ++i)
	PrintData(i, mapData[i], kCell, &cells[i]);
    CheckError(err, "Error reading cells");

    // Read the internal faces.
    CCMIOGetEntity(&err, topology, kCCMIOInternalFaces, 0, &id);
    CCMIOEntitySize(&err, id, &nFaces, NULL);
    mapData.resize(nFaces);
    faceCells.resize(2 * nFaces);
    CCMIOReadFaces(&err, id, kCCMIOInternalFaces, NULL, &size, NULL,
		   kCCMIOStart, kCCMIOEnd);
    faces.resize(size);
    CCMIOReadFaces(&err, id, kCCMIOInternalFaces, &mapID, NULL, &faces[0],
		   kCCMIOStart, kCCMIOEnd);
    CCMIOReadFaceCells(&err, id, kCCMIOInternalFaces, &faceCells[0],
		       kCCMIOStart, kCCMIOEnd);
    CCMIOReadMap(&err, mapID, &mapData[0], kCCMIOStart, kCCMIOEnd);
    
    cout << "\t" << nFaces << " faces" << endl;
    unsigned int pos = 0;
    i = 0;
    while (pos < faces.size())
    {
	PrintData(i, mapData[i], kInternalFace, &faces[pos], &faceCells[2 * i]);
	pos += faces[pos] + 1;
	i++;
    }
    CheckError(err, "Error reading internal faces");

    // Read the boundary faces.
    int index = 0;
    while (CCMIONextEntity(NULL, topology, kCCMIOBoundaryFaces, &index, &id)
								 == kCCMIONoErr)
    {
	int boundaryVal;
	vector<int> prostarIDs;

	CCMIOEntitySize(&err, id, &nFaces, NULL);
	mapData.resize(nFaces);
	faceCells.resize(nFaces);
	CCMIOReadFaces(&err, id, kCCMIOBoundaryFaces, NULL, &size, NULL,
		       kCCMIOStart, kCCMIOEnd);
	faces.resize(size);
	CCMIOReadFaces(&err, id, kCCMIOBoundaryFaces, &mapID, NULL, &faces[0],
		       kCCMIOStart, kCCMIOEnd);
	CCMIOReadFaceCells(&err, id, kCCMIOBoundaryFaces, &faceCells[0],
			   kCCMIOStart, kCCMIOEnd);
	CCMIOReadMap(&err, mapID, &mapData[0], kCCMIOStart, kCCMIOEnd);

	// If the optional ProstarFaceId is there, read it
	if (CCMIOReadOpt1i(NULL, id, "ProstarFaceId", NULL,
			   kCCMIOStart, kCCMIOEnd) == kCCMIONoErr)
	{
	    prostarIDs.resize(nFaces);
	    CCMIOReadOpt1i(NULL, id, "ProstarFaceId", &prostarIDs[0],
			   kCCMIOStart, kCCMIOEnd);
	}
    
	CCMIOGetEntityIndex(&err, id, &boundaryVal);
	cout << "\t" << nFaces << " boundary faces (boundary " << boundaryVal
	     << ")" << endl;
	unsigned int pos = 0;
	j = 0;
	while (pos < faces.size())
	{
	    PrintData(j, mapData[j], kBoundaryFace, &faces[pos],
		      &faceCells[j]);
	    j++;
	    pos += faces[pos] + 1;
	}
	CheckError(err, "Error reading boundary faces");
    }
}

void ReadVertices( CCMIOError &err, CCMIOID &vertices, char const *counter,
		   int offset )
{
    int dims = 1;
    unsigned int i, nVertices, size;
    float scale;
    int mapData[kNValues] = { 0 };
    float verts[3 * kNValues] = { 0 };
    CCMIOID mapID;

    // Read the vertices.  This involves reading both the vertex data and
    // the map, which maps the index into the data array with the ID number.
    // As we process the vertices we need to be sure to scale them by the
    // appropriate scaling factor.  The offset is just to show you can read
    // any chunk.  Normally this would be in a for loop.
    CCMIOEntitySize(&err, vertices, &nVertices, NULL);
    CCMIOReadVerticesf(&err, vertices, &dims, &scale, &mapID, verts,
		       offset, offset + kNValues);
    CCMIOReadMap(&err, mapID, mapData, offset, offset + kNValues);

    CCMIOEntityDescription(&err, vertices, &size, NULL);
    char *label = new char[size + 1];
    CCMIOEntityDescription(&err, vertices, NULL, label);
    cout << "label: '" << label << "'" << endl;
    delete [] label;

    cout << "\t" << nVertices << " " << counter << endl;
    if (nVertices > (unsigned int)kNValues)
	nVertices = kNValues;
    for (i = 0;  i < nVertices && err == kCCMIONoErr;  ++i)
    {
	verts[dims * i    ] *= scale;
	verts[dims * i + 1] *= scale;
	verts[dims * i + 2] *= scale;   // This example assumes three
					// dimensional vertices.
	PrintData(i, mapData[i], kVertex, &verts[dims * i]);
    }
    cout << "\t\t..." << endl;
    CheckError(err, "Error reading vertices");
}

void ReadSolverInfo( CCMIOError &err, CCMIOID &solution)
{
    char solver[kCCMIOMaxStringLength+1], timeUnits[kCCMIOMaxStringLength+1];
    int iterations;
    float time, angle;
    CCMIOID restart;

    if (CCMIOGetEntity(NULL, solution, kCCMIORestart, 0, &restart)!=kCCMIONoErr)
    {
	cout << "(No solver information)" << endl;
	return;
    }

    CCMIOReadRestartInfo(&err, restart, solver, &iterations, &time, timeUnits,
			 &angle);
    cout << "The following post data is the result of iteration "
	 << iterations << " (time: " << time << " " << timeUnits
	 << ")" << endl;
    cout << "of solver '" << solver << "':" << endl;

    // You might want to read more solver data here, in which case, do
    // CCMIOGetEntity(&err, restart, kCCMIORestartData, 0, &restartData);
    // and proceed to read the information.  Do likewise if you want to read
    // the reference data.
}

void ReadPost( CCMIOError &err, CCMIOID &solution )
{
    bool oldFile = false;
    char name[kCCMIOMaxStringLength + 1];
    int h = 0, i = 0, dimSizes[] = { 0, 1, 3, 9 };
    CCMIOID field, phase;
    CCMIODimensionality dims;
    vector<int> mapData;
    vector<float> data;

    oldFile = (CCMIONextEntity(NULL, solution, kCCMIOFieldPhase, &h, &phase)
	       							!= kCCMIONoErr);
    h = 0;
    while (oldFile ||
	   CCMIONextEntity(NULL, solution, kCCMIOFieldPhase, &h, &phase)
	   							== kCCMIONoErr)
    {
	if (oldFile)
	    phase = solution;
	else
	{
	    int phaseNum = 0;
	    CCMIOGetEntityIndex(NULL, phase, &phaseNum);
	    cout << "Phase " << phaseNum << endl;
	}
	// Walk through each field in this field set (i.e. the post data)
	// and print out each one.
	while (CCMIONextEntity(NULL, phase, kCCMIOField, &i, &field)
								== kCCMIONoErr)
	{
	    // Read the information about the field so that we know how to
	    // process this field.
	    char shortName[kCCMIOProstarShortNameLength+1], *units = NULL;
	    int usize;
	    CCMIODataType datatype;
	    CCMIOReadField(&err, field, name, shortName, &dims, &datatype);
	    if (CCMIOReadOptstr(NULL, field, kUnitsName, &usize, NULL)
								 == kCCMIONoErr)
	    {
		units = new char[usize + 1];
		CCMIOReadOptstr(&err, field, kUnitsName, NULL, units);
	    }
	    cout << "\tPost field '" << name << "'  (" << shortName << "):\t(";
	    if (dims == 1)    // datatype is not meaningful for vectors/tensors
	    {
		if (datatype == kCCMIOFloat32)
		    cout << "Float32";
		else if (datatype == kCCMIOFloat64)
		    cout << "Float64";
		else if (datatype == kCCMIOInt32)
		    cout << "Int32";
		else
		    cout << "Invalid datatype";
	    }
	    cout << ")\tUnits: " << ((units) ? units : "<none>") << endl;
	    delete [] units;
	    
	    switch (dims)
	    {
		case kCCMIOScalar:
		    {
		    ReadScalar(err, field, mapData, data);
		    for (unsigned int k = 0;  k < data.size();  ++k)
			PrintData(k, mapData[i], kScalar, &data[k]); 
		    }
		    break;
		case kCCMIOVector:
		    {
		    vector<float> u, v, w;
		    CCMIOID scalar;
		    CCMIOReadMultiDimensionalFieldData(&err, field,
						       kCCMIOVectorX, &scalar);
		    if (err == kCCMIOVersionErr)
		    {
			// If we are reading an older version of the file,
			// where vectors are stored as vectors, not components,
			// we need to call CCMIOReadFieldData*(), which is
			// all that ReadScalar() does.
			err = kCCMIONoErr;
			ReadScalar(err, field, mapData, data, true);
			unsigned int const size = data.size() / dimSizes[dims];
			for (unsigned int k = 0;  k < size;  ++k)
			    PrintData(k, mapData[i], kVector,
				  &data[dimSizes[dims] * k]);
		    }
		    else
		    {
			ReadScalar(err, scalar, mapData, u);
			CCMIOReadMultiDimensionalFieldData(&err, field,
							   kCCMIOVectorY,
							   &scalar);
			ReadScalar(err, scalar, mapData, v);
			CCMIOReadMultiDimensionalFieldData(&err, field,
							   kCCMIOVectorZ,
							   &scalar);
			ReadScalar(err, scalar, mapData, w);
			data.resize(3 * u.size());
			for (unsigned int k = 0;  k < u.size();  ++k)
			{
			    data[3 * k    ] = u[k];
			    data[3 * k + 1] = v[k];
			    data[3 * k + 2] = w[k];
			}
			for (unsigned int k = 0;  k < u.size();  ++k)
			    PrintData(k, mapData[i], kVector,
				      &data[dimSizes[dims] * k]);
		    }
		    }
		    break;
		case kCCMIOTensor:
		    cout << "Tensor data not supported.  Ignoring field "
			 << name << "." << endl;
		    continue;
	    }
	    CheckError(err, "Error reading post data");
	}
	oldFile = false;
    }
}

void ReadScalar( CCMIOError &err, CCMIOID field,
		 vector<int> &mapData, vector<float> &data,
		 bool readingVector /* = false */)
{
    CCMIOSize n;
    CCMIOIndex max;
    int j = 0;
    CCMIOID fieldData, mapID;
    CCMIODataLocation type;

    // Read each piece of field data
    while (CCMIONextEntity(NULL, field, kCCMIOFieldData, &j, &fieldData)
								 == kCCMIONoErr)
    {
	// Figure out how big this data is so we can read it. If we were
	// storing this information permanently we might use a sparse
	// array, in which case we would need to find the maximum ID and
	// make the array that size.
	CCMIOEntitySize(&err, fieldData, &n, &max);
	mapData.resize(n);
	CCMIOReadFieldDataf(&err, fieldData, &mapID, &type, NULL,
			    kCCMIOStart, kCCMIOEnd);
	CCMIOReadMap(&err, mapID, &mapData[0], kCCMIOStart, kCCMIOEnd);

	// We are only going to process cell data.  Vertex data would
	// be processed similarly. If your appliation has only one value
	// for boundary data, you would separate the face data into
	// each boundary and combine it together.  If your application
	// stores boundary data on each face you could read it in using
	// similar procedures for the cell and vertex data.  Note that
	// the file may not contain all types of data.
	if (type == kCCMIOCell)
	{
	    if (readingVector)
		data.resize(3 * n);
	    else
		data.resize(n);
	    // If we want double precision we should use
	    // CCMIOReadFieldDatad().
	    CCMIOReadFieldDataf(&err, fieldData, &mapID, NULL,
				&data[0], kCCMIOStart, kCCMIOEnd);
	}

	CheckError(err, "Error reading post data");
    }
}

void ReadSets( CCMIOError &err, CCMIOID &problem )
{
    int nChars;
	
    if (CCMIOReadOptstr(NULL, problem, "MonitoringSets", &nChars, NULL) == kCCMIONoErr)
    {
	char *setStr = new char[nChars + 1];
	char setNames[][16] = { "CellSet", "VertexSet", "BoundarySet",
				"BlockSet", "SplineSet", "CoupleSet" };
	int i = 0, j, len, size;
	unsigned int flags;
	string name;
	CCMIOID setID;
	    
	cout << "Monitoring prostar sets:" << endl;

	// See the CCM specification for official documentation of this node.
	// The format is "<size of string> string <size of string> string ..."
	// The size of the string is one (binary) character.
	CCMIOReadOptstr(&err, problem, "MonitoringSets", NULL, setStr);
	while (i < nChars)
	{
	    j = 0;
	    len = setStr[i++];	// First character is the length
	    name.clear();
	    for (j = 0;  j < len;  ++j)  // Then the string (no C-style \0)
		name.push_back(setStr[i++]);
	    CCMIOGetProstarSet(&err, problem, name.c_str(), &size, NULL,
			       NULL, NULL);
	    char *longName = new char[size + 1];
	    CCMIOGetProstarSet(&err, problem, name.c_str(), NULL, longName,
			       &flags, &setID);
	    cout << "\t" << name << " (" << longName << "):" << endl;
	    delete [] longName;

	    // Read each set that is available and print it out.
	    // There is no need to loop through the nodes if you know which
	    // set you want;  just check 'flags' to see if the node is there
	    // (or just access the node and check for an error).
	    int n = 0;
	    for (unsigned int s = 1;  s <= kCCMIOCoupleSet;  s = s << 1)
	    {
		if (flags & s)
		{
		    int *theSet;
		    unsigned int size;
		    CCMIOGetOptInfo(&err, setID, setNames[n], NULL,
				    &size, NULL, NULL);
		    theSet = new int[size];
		    CCMIOReadOpt1i(&err, setID, setNames[n], theSet,
				   kCCMIOStart, kCCMIOEnd);

		    // Print out the set
		    cout << "\t\t" << setNames[n] << ":  ";
		    for (unsigned int idx = 0;  idx < size;  ++idx)
			cout << theSet[idx] << "  ";
		    cout << endl;
		    delete [] theSet;
		}
		n++;
	    }
	}

	delete [] setStr;
    }
    else
	cout << "No prostar sets selected for monitoring" << endl;
    CheckError(err, "Error reading prostar sets");
}
//--------------------------- Helper functions --------------------------------
void CheckError( CCMIOError const &err, char const *str )
{
    if (err == kCCMIONoErr)
	return;

    cout << str << " (error " << err << ")" << endl;
    exit(1);
}

void PrintData( int n, int id, DataType type, void *data, void *data2 /*=NULL*/ )
{
    if (n > kNValues)
	return;
    if (n == kNValues)
    {
	cout << "\t\t..." << endl;
	return;
    }

    switch (type)
    {
	case kScalar:
	    cout << "\t\t" << *(float *)data << endl;
	    break;
	case kVector:
	    cout << "\t\t(" << ((float *)data)[0] << ", "
		 << ((float *)data)[1] << ", " << ((float *)data)[2] << ")"
		 << endl;
	    break;
	case kVertex:
	    cout << "\t\t" << "Vertex " << id << ":  ("
		 << ((float *)data)[0] << ", "
		 << ((float *)data)[1] << ", "
		 << ((float *)data)[2] << ")" << endl;
	    break;
	case kCell:
	    cout << "\t\t" << "Cell " << id << ":  type " << *(int *)data
		 << endl;
	    break;
	case kInternalFace:
	    {
		int nVerts = *(int *)data;
		cout << "\t\t" << "Face " << id << ": \tCells: "
		     << ((int *)data2)[0] << ", " << ((int *)data2)[1]
		     << "\tVertices: ";
		for (int i = 0;  i < nVerts;  ++i)
		    cout << *(((int *)data) + 1 + i) << " ";
		cout << endl;
	    }
	    break;
	case kBoundaryFace:
	    {
		int nVerts = *(int *)data;
		cout << "\t\t" << "Face " << id << ": \tCell: "
		     << *(int *)data2 << "\tVertices: ";
		for (int i = 0;  i < nVerts;  ++i)
		    cout << *(((int *)data) + 1 + i) << " ";
		cout << endl;
	    }
	    break;
	case kBoundaryData:
	    cout << "\t\tBoundary " << *(int *)data << "  priority "
		 << *(int *)data2 << endl;
	    break;
	case kBoundaryFaceData:
	    {
		cout << "\t\tBoundary data [ ";
		for (int i = 0;  i < *(int *)data2; ++i)
		    cout << ((float *)data)[i] << " ";
		cout << "]" << endl;
	    }
	    break;
	case kCellType:
	    cout << "Cell type " << id << ":  " << (char *)data << endl;
	    break;
    }
}
