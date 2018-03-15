/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "tecplotWriter.H"
#include "fvMesh.H"
#include "TECIO.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const int32_t Foam::tecplotWriter::tecConst_0 = 0;
const int32_t Foam::tecplotWriter::tecConst_1 = 1;
const int32_t Foam::tecplotWriter::tecConst_False = 0;
const int32_t Foam::tecplotWriter::tecConst_True  = 1;

const Foam::string Foam::tecplotWriter::XYZ = "X Y Z";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tecplotWriter::tecplotWriter(const Time& runTime)
:
    time_(runTime)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tecplotWriter::writeInit
(
    const word& name,
    const string& varNames,
    const fileName& fName,
    const dataFileType fileType
) const
{
    const int32_t FileType   = fileType;
    const int32_t FileFormat = 0; // 0 = binary (plt), 1 = subzone (.szplt)

    Pout<< nl << nl
        << "Name:" << name
        << " varNames:" << varNames
        << " to file:" << fName
        << " of type:" << int(fileType)
        << endl;

    if
    (
        tecini142
        (
            name.c_str(),       //< DataSet Title
            varNames.c_str(),   //< Variables List
            fName.c_str(),      //< FileName
            time_.path().c_str(), //< ScratchDir
            &FileFormat,        //< FileFormat
            &FileType,          //< FileType
            &tecConst_False,    //< Debug (0: no debug, 1: debug)
            &tecConst_False     //< VIsDouble (0: single, 1: double)
        )
    )
    {
        FatalErrorInFunction
            << "Error in tecini142."
            << exit(FatalError);
    }
}


void Foam::tecplotWriter::writePolyhedralZone
(
    const word& zoneName,
    const int32_t strandID,
    const fvMesh& mesh,
    const UList<int32_t>& varLocArray,
    const int32_t NumFaceNodes
) const
{
    const int32_t NumNodes = mesh.nPoints();    // Number of unique nodes
    const int32_t NumElems = mesh.nCells();     // Number of elements
    const int32_t NumFaces = mesh.nFaces();     // Number of unique faces
    const double  SolTime  = time_.value();     // Solution time

    const int32_t ParentZone = 0;   // Bool: 0 = no parent zone
    const int32_t ShrConn   = 0;
    const int32_t NumBConns = 0;    // No Boundary Connections
    const int32_t NumBItems = 0;    // No Boundary Items

    Pout<< "zoneName:" << zoneName
        //<< " varLocArray:" << varLocArray
        << " solTime:" << SolTime
        << " strand:"  << strandID
        << endl;

    const int32_t ZoneType  = ZONE_FEPOLYHEDRON;
    if
    (
        teczne142
        (
            zoneName.c_str(),   //< ZoneTitle
            &ZoneType,          //< ZoneType
            &NumNodes,          //< IMxOrNumPts
            &NumElems,          //< JMxOrNumElements
            &NumFaces,          //< KMxOrNumFaces
            &tecConst_0,        //< (unused set to zero) ICellMax
            &tecConst_0,        //< (unused set to zero) JCellMax
            &tecConst_0,        //< (unused set to zero) KCellMax
            &SolTime,           //< SolutionTime
            &strandID,          //< StrandID
            &ParentZone,        //< ParentZone
            &tecConst_True,     //< IsBlock
            &tecConst_0,        //< (unused) NumFaceConnections
            &tecConst_0,        //< (unused) FaceNeighborMode
            &NumFaceNodes,      //< TotalNumFaceNodes
            &NumBConns,         //< NumConnectedBoundaryFaces
            &NumBItems,         //< TotalNumBoundaryConnections
            nullptr,            //< PassiveVarList
            varLocArray.cdata(), //< ValueLocation
            nullptr,            //< ShareVarFromZone
            &ShrConn            //< ShareConnectivityFromZone
        )
    )
    {
        FatalErrorInFunction
            << "Error in teczne142 - writing polyhedron zones."
            << exit(FatalError);
    }
}


void Foam::tecplotWriter::writePolygonalZone
(
    const word& zoneName,
    const int32_t strandID,
    const indirectPrimitivePatch& pp,
    const UList<int32_t>& varLocArray
) const
{
    const int32_t NumNodes = pp.nPoints();      // Number of unique nodes
    const int32_t NumElems = pp.size();         // Number of elements
    const int32_t NumFaces = pp.nEdges();       // Number of unique faces
    const double  SolTime  = time_.value();     // Solution time

    const int32_t ParentZone = 0;   // Int: 0 = no parent zone
    const int32_t NumFaceNodes = 2*pp.nEdges();

    const int32_t ShrConn   = 0;
    const int32_t NumBConns = 0;    // No Boundary Connections
    const int32_t NumBItems = 0;    // No Boundary Items

    Pout<< "zoneName:" << zoneName
        << " strandID:" << strandID
        //<< " varLocArray:" << varLocArray
        << " solTime:" << SolTime
        << endl;

    const int32_t ZoneType = ZONE_FEPOLYGON;
    if
    (
        teczne142
        (
            zoneName.c_str(),   //< ZoneTitle
            &ZoneType,          //< ZoneType
            &NumNodes,          //< IMax or NumPts
            &NumElems,          //< JMax or NumElements
            &NumFaces,          //< KMax or NumFaces
            &tecConst_0,        //< (Unused set to zero) ICellMax
            &tecConst_0,        //< (Unused set to zero) JCellMax
            &tecConst_0,        //< (Unused set to zero) KCellMax
            &SolTime,           //< SolutionTime
            &strandID,          //< StrandID
            &ParentZone,        //< ParentZone
            &tecConst_True,     //< IsBlock
            &tecConst_0,        //< (Unused for polygon zone) NumFaceConnections
            &tecConst_0,        //< (Unused for polygon zone) FaceNeighborMode
            &NumFaceNodes,      //< TotalNumFaceNodes
            &NumBConns,         //< NumConnectedBoundaryFaces
            &NumBItems,         //< TotalNumBoundaryConnections
            nullptr,            //< PassiveVarList
            varLocArray.cdata(), //< ValueLocation
            nullptr,            //< ShareVarFromZone
            &ShrConn            //< ShareConnectivityFromZone
        )
    )
    {
        FatalErrorInFunction
            << "Error in teczne142 - writing polygon zones."
            << exit(FatalError);
    }
}


void Foam::tecplotWriter::writeOrderedZone
(
    const word& zoneName,
    const int32_t strandID,
    const label n,
    const UList<int32_t>& varLocArray
) const
{
    const int32_t IMax = n;     // Number in I direction
    const int32_t JMax = 1;     // Number in J direction
    const int32_t KMax = 1;     // Number in K direction
    const double  SolTime = time_.value();  // Solution time

    const int32_t ParentZone = 0;   // Bool: no parent zone
    const int32_t NFConns = 0;      // Unused for ordered zones
    const int32_t FNMode  = 0;      // Unused for ordered zones

    const int32_t ShrConn  = 0;
    const int32_t NumFaceNodes = 1;
    const int32_t NumBConns = 0;    // No Boundary Connections
    const int32_t NumBItems = 0;    // No Boundary Items

    Pout<< "zoneName:" << zoneName
        << " strandID:" << strandID
        //<< " varLocArray:" << varLocArray
        << " solTime:" << SolTime
        << endl;

    const int32_t ZoneType = ZONE_ORDERED;
    if
    (
        teczne142
        (
            zoneName.c_str(),   //< ZoneTitle
            &ZoneType,          //< ZoneType
            &IMax,              //< IMax or NumPts
            &JMax,              //< JMax or NumElements
            &KMax,              //< KMax or NumFaces
            &tecConst_0,        //< (Unused set to zero) ICellMax
            &tecConst_0,        //< (Unused set to zero) JCellMax
            &tecConst_0,        //< (Unused set to zero) KCellMax
            &SolTime,           //< SolutionTime
            &strandID,          //< StrandID
            &ParentZone,        //< ParentZone
            &tecConst_True,     //< IsBlock
            &NFConns,           //< NumFaceConnections
            &FNMode,            //< FaceNeighborMode
            &NumFaceNodes,      //< TotalNumFaceNodes
            &NumBConns,         //< NumConnectedBoundaryFaces
            &NumBItems,         //< TotalNumBoundaryConnections
            nullptr,            //< PassiveVarList
            varLocArray.cdata(), //< ValueLocation
            nullptr,            //< ShareVarFromZone
            &ShrConn            //< ShareConnectivityFromZone
        )
    )
    {
        FatalErrorInFunction
            << "Error in teczne142 - writing ordered zones."
            << exit(FatalError);
    }
}


void Foam::tecplotWriter::writeConnectivity(const fvMesh& mesh) const
{
    // first pass: get the sizes
    List<int32_t> FaceNodeCounts(mesh.nFaces());
    label nFaceNodes = 0;
    forAll(mesh.faces(), facei)
    {
        const face& f = mesh.faces()[facei];
        nFaceNodes += f.size();
        FaceNodeCounts[facei] = int32_t(f.size());
    }

    // second pass: get the nodes as a flat list
    List<int32_t> FaceNodes(nFaceNodes);
    label nodeI = 0;
    forAll(mesh.faces(), facei)
    {
        const face& f = mesh.faces()[facei];
        forAll(f, fp)
        {
            FaceNodes[nodeI++] = int32_t(f[fp]+1);
        }
    }


    List<int32_t> FaceLeftElems(mesh.nFaces());
    forAll(mesh.faceOwner(), facei)
    {
        FaceLeftElems[facei] = mesh.faceOwner()[facei]+1;
    }

    List<int32_t> FaceRightElems(mesh.nFaces());
    forAll(mesh.faceNeighbour(), facei)
    {
        FaceRightElems[facei] = mesh.faceNeighbour()[facei]+1;
    }
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        FaceRightElems[facei] = 0;
    }

    if
    (
        tecpoly142
        (
            FaceNodeCounts.cdata(), // The face node counts array
            FaceNodes.cdata(),      // The face nodes array
            FaceLeftElems.cdata(),  // The left elements array
            FaceRightElems.cdata(), // The right elements array
            nullptr,       // No face boundary connection counts
            nullptr,       // No face boundary connection elements
            nullptr        // No face boundary connection zones
        )
    )
    {
        FatalErrorInFunction
            << "Error in tecpoly142."
            << exit(FatalError);
    }
}

void Foam::tecplotWriter::writeConnectivity
(
    const indirectPrimitivePatch& pp
) const
{
    const int32_t NumFaces     = pp.nEdges();   // Number of unique faces
    const int32_t NumFaceNodes = 2*NumFaces;    // 2 nodes per edge

    // All faces (=edges) have 2 nodes
    List<int32_t> FaceNodeCounts(NumFaces);
    FaceNodeCounts = 2;

    List<int32_t> FaceNodes(NumFaceNodes);
    label nodeI = 0;
    forAll(pp.edges(), edgei)
    {
        edge e = pp.edges()[edgei];
        if (e[0] > e[1])
        {
            e.flip();
        }

        FaceNodes[nodeI++] = int32_t(e[0]+1);
        FaceNodes[nodeI++] = int32_t(e[1]+1);
    }

    /* Define the right and left elements of each face.
     *
     * The last step for writing out the polyhedral data is to
     * define the right and left neighboring elements for each
     * face.  The neighboring elements can be determined using the
     * right-hand rule.  For each face, place your right-hand along
     * the face which your fingers pointing the direction of
     * incrementing node numbers (i.e. from node 1 to node 2).
     * Your right thumb will point towards the right element; the
     * element on the other side of your hand is the left element.
     *
     * The number zero is used to indicate that there isn't an
     * element on that side of the face.
     *
     * Because of the way we numbered the nodes and faces, the
     * right element for every face is the element itself
     * (element 1) and the left element is "no-neighboring element"
     * (element 0).
     */

    List<int32_t> FaceLeftElems(NumFaces);
    List<int32_t> FaceRightElems(NumFaces);

    const labelListList& edgeFaces = pp.edgeFaces();
    forAll(edgeFaces, edgei)
    {
        const labelList& eFaces = edgeFaces[edgei];

        if (eFaces.size() == 1)
        {
            FaceLeftElems[edgei]  = 0;
            FaceRightElems[edgei] = eFaces[0]+1;
        }
        else if (eFaces.size() == 2)
        {
            edge e = pp.edges()[edgei];
            if (e[0] > e[1])
            {
                e.flip();
            }

            const face& f0 = pp.localFaces()[eFaces[0]];

            // The face that uses the vertices of e in increasing order
            // is the left face.

            const label fp = f0.find(e[0]);
            const bool f0IsLeft = (f0.nextLabel(fp) == e[1]);

            if (f0IsLeft)
            {
                FaceLeftElems[edgei]  = eFaces[0]+1;
                FaceRightElems[edgei] = eFaces[1]+1;
            }
            else
            {
                FaceLeftElems[edgei]  = eFaces[1]+1;
                FaceRightElems[edgei] = eFaces[0]+1;
            }
        }
        else
        {
            // non-manifold. Treat as if open.
            FaceLeftElems[edgei]  = 0;
            FaceRightElems[edgei] = eFaces[0]+1;
        }
    }

    // Write the face map (created above)
    if
    (
        tecpoly142
        (
            FaceNodeCounts.cdata(),     // Face node counts array
            FaceNodes.cdata(),          // Face nodes array
            FaceLeftElems.cdata(),      // Left elements array
            FaceRightElems.cdata(),     // Right elements array
            nullptr,                    // No boundary connection counts
            nullptr,                    // No boundary connection elements
            nullptr                     // No boundary connection zones
        )
    )
    {
        FatalErrorInFunction
            << "Error in tecpoly142."
            << exit(FatalError);
    }
}


void Foam::tecplotWriter::writeEnd() const
{
    Pout<< "writeEnd" << endl;

    if (tecend142())
    {
        FatalErrorInFunction
            << "Error in tecend142."
            << exit(FatalError);
    }
}


// ************************************************************************* //
