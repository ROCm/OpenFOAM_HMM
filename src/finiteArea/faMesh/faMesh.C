/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

#include "faMesh.H"
#include "faGlobalMeshData.H"
#include "Time.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "IndirectList.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "faMeshLduAddressing.H"
#include "wedgeFaPatch.H"
#include "faPatchData.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faMesh, 0);
}

Foam::word Foam::faMesh::meshSubDir = "faMesh";

const int Foam::faMesh::quadricsFit_ = 0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMesh::setPrimitiveMeshData()
{
    DebugInFunction << "Setting primitive data" << endl;

    const indirectPrimitivePatch& bp = patch();

    // Set faMesh edges
    edges_.setSize(bp.nEdges());

    label edgeI = -1;

    label nIntEdges = bp.nInternalEdges();

    for (label curEdge = 0; curEdge < nIntEdges; ++curEdge)
    {
        edges_[++edgeI] = bp.edges()[curEdge];
    }

    forAll(boundary(), patchI)
    {
        const labelList& curP = boundary()[patchI];

        forAll(curP, eI)
        {
            edges_[++edgeI] = bp.edges()[curP[eI]];
        }
    }

    nEdges_ = edges_.size();
    nInternalEdges_ = nIntEdges;


    // Set edge owner and neighbour
    edgeOwner_.setSize(nEdges());
    edgeNeighbour_.setSize(nInternalEdges());

    edgeI = -1;

    const labelListList& bpef = bp.edgeFaces();

    for (label curEdge = 0; curEdge < nIntEdges; ++curEdge)
    {
        edgeOwner_[++edgeI] = bpef[curEdge][0];

        edgeNeighbour_[edgeI] = bpef[curEdge][1];
    }

    forAll(boundary(), patchI)
    {
        const labelList& curP = boundary()[patchI];

        forAll(curP, eI)
        {
            edgeOwner_[++edgeI] = bpef[curP[eI]][0];
        }
    }

    // Set number of faces
    nFaces_ = bp.size();

    // Set number of points
    nPoints_ = bp.nPoints();
}


void Foam::faMesh::clearGeomNotAreas() const
{
    DebugInFunction << "Clearing geometry" << endl;

    deleteDemandDrivenData(SPtr_);
    deleteDemandDrivenData(patchPtr_);
    deleteDemandDrivenData(patchStartsPtr_);
    deleteDemandDrivenData(LePtr_);
    deleteDemandDrivenData(magLePtr_);
    deleteDemandDrivenData(centresPtr_);
    deleteDemandDrivenData(edgeCentresPtr_);
    deleteDemandDrivenData(faceAreaNormalsPtr_);
    deleteDemandDrivenData(edgeAreaNormalsPtr_);
    deleteDemandDrivenData(pointAreaNormalsPtr_);
    deleteDemandDrivenData(faceCurvaturesPtr_);
    deleteDemandDrivenData(edgeTransformTensorsPtr_);
}


void Foam::faMesh::clearGeom() const
{
    DebugInFunction << "Clearing geometry" << endl;

    clearGeomNotAreas();
    deleteDemandDrivenData(S0Ptr_);
    deleteDemandDrivenData(S00Ptr_);
    deleteDemandDrivenData(correctPatchPointNormalsPtr_);
}


void Foam::faMesh::clearAddressing() const
{
    DebugInFunction << "Clearing addressing" << endl;

    deleteDemandDrivenData(lduPtr_);
}


void Foam::faMesh::clearOut() const
{
    clearGeom();
    clearAddressing();
    deleteDemandDrivenData(globalMeshDataPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMesh::faMesh(const polyMesh& pMesh)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>(pMesh),
    edgeInterpolation(*this),
    faSchemes(mesh()),
    faSolution(mesh()),
    data(mesh()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            time().findInstance(meshDir(), "faceLabels"),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            time().findInstance(meshDir(), "faBoundary"),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    DebugInFunction << "Creating faMesh from IOobject" << endl;

    setPrimitiveMeshData();

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    if (isFile(pMesh.time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
}


Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const labelList& faceLabels
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>(pMesh),
    edgeInterpolation(*this),
    faSchemes(mesh()),
    faSolution(mesh()),
    data(mesh()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        faceLabels
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    DebugInFunction << "Creating faMesh from components" << endl;
}


Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const fileName& defFile
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>(pMesh),
    edgeInterpolation(*this),
    faSchemes(mesh()),
    faSolution(mesh()),
    data(mesh()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        List<label>(0)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    DebugInFunction << "Creating faMesh from definition file" << endl;

    // Reading faMeshDefinition dictionary
    IOdictionary faMeshDefinition
    (
        IOobject
        (
            defFile,
            mesh().time().constant(),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const wordList polyMeshPatches(faMeshDefinition.lookup("polyMeshPatches"));

    const dictionary& bndDict = faMeshDefinition.subDict("boundary");

    const wordList faPatchNames(bndDict.toc());

    List<faPatchData> faPatches(faPatchNames.size() + 1);

    const polyBoundaryMesh& pbm = pMesh.boundaryMesh();

    forAll(faPatchNames, patchI)
    {
        dictionary curPatchDict = bndDict.subDict(faPatchNames[patchI]);

        faPatches[patchI].name_ = faPatchNames[patchI];

        faPatches[patchI].type_ = word(curPatchDict.lookup("type"));

        faPatches[patchI].ownPolyPatchID_ =
            pbm.findPatchID(word(curPatchDict.lookup("ownerPolyPatch")));

        faPatches[patchI].ngbPolyPatchID_ =
            pbm.findPatchID(word(curPatchDict.lookup("neighbourPolyPatch")));
    }


    // Setting faceLabels list size
    label size = 0;

    labelList patchIDs(polyMeshPatches.size(), -1);

    forAll(polyMeshPatches, patchI)
    {
        patchIDs[patchI] = pbm.findPatchID(polyMeshPatches[patchI]);

        size += pbm[patchIDs[patchI]].size();
    }

    faceLabels_ = labelList(size, -1);


    // Filling of faceLabels list
    label faceI = -1;

    sort(patchIDs);

    forAll(polyMeshPatches, patchI)
    {
        label start = pbm[patchIDs[patchI]].start();
        label size  = pbm[patchIDs[patchI]].size();

        for (label i = 0; i < size; ++i)
        {
            faceLabels_[++faceI] = start + i;
        }
    }


    // Determination of faPatch ID for each boundary edge.
    // Result is in the bndEdgeFaPatchIDs list
    labelList faceCells(faceLabels_.size(), -1);

    forAll(faceCells, faceI)
    {
        label faceID = faceLabels_[faceI];

        faceCells[faceI] = mesh().faceOwner()[faceID];
    }

    labelList meshEdges =
        patch().meshEdges
        (
            mesh().edges(),
            mesh().cellEdges(),
            faceCells
        );

    const labelListList& edgeFaces = mesh().edgeFaces();

    const label nTotalEdges = patch().nEdges();
    const label nInternalEdges = patch().nInternalEdges();

    labelList bndEdgeFaPatchIDs(nTotalEdges - nInternalEdges, -1);

    for (label edgeI = nInternalEdges; edgeI < nTotalEdges; ++edgeI)
    {
        label curMeshEdge = meshEdges[edgeI];

        labelList curEdgePatchIDs(2, label(-1));

        label patchI = -1;

        forAll(edgeFaces[curMeshEdge], faceI)
        {
            label curFace = edgeFaces[curMeshEdge][faceI];

            label curPatchID = pbm.whichPatch(curFace);

            if (curPatchID != -1)
            {
                curEdgePatchIDs[++patchI] = curPatchID;
            }
        }

        for (label pI = 0; pI < faPatches.size() - 1; ++pI)
        {
            if
            (
                (
                    curEdgePatchIDs[0] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[1] == faPatches[pI].ngbPolyPatchID_
                )
             ||
                (
                    curEdgePatchIDs[1] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[0] == faPatches[pI].ngbPolyPatchID_
                )
            )
            {
                bndEdgeFaPatchIDs[edgeI - nInternalEdges] = pI;
                break;
            }
        }
    }


    // Set edgeLabels for each faPatch
    for (label pI = 0; pI < (faPatches.size() - 1); ++pI)
    {
        SLList<label> tmpList;

        forAll(bndEdgeFaPatchIDs, eI)
        {
            if (bndEdgeFaPatchIDs[eI] == pI)
            {
                tmpList.append(nInternalEdges + eI);
            }
        }

        faPatches[pI].edgeLabels_ = tmpList;
    }

    // Check for undefined edges
    SLList<label> tmpList;

    forAll(bndEdgeFaPatchIDs, eI)
    {
        if (bndEdgeFaPatchIDs[eI] == -1)
        {
            tmpList.append(nInternalEdges + eI);
        }
    }

    if (tmpList.size() > 0)
    {
        // Check for processor edges
        labelList allUndefEdges(tmpList);
        labelList ngbPolyPatch(allUndefEdges.size(), -1);
        forAll(ngbPolyPatch, edgeI)
        {
            label curEdge = allUndefEdges[edgeI];

            label curPMeshEdge = meshEdges[curEdge];

            forAll(edgeFaces[curPMeshEdge], faceI)
            {
                label curFace = edgeFaces[curPMeshEdge][faceI];

                if (findIndex(faceLabels_, curFace) == -1)
                {
                    label polyPatchID = pbm.whichPatch(curFace);

                    if (polyPatchID != -1)
                    {
                        ngbPolyPatch[edgeI] = polyPatchID;
                    }
                }
            }
        }

        // Count ngb processorPolyPatch-es
        labelHashSet processorPatchSet;
        forAll(ngbPolyPatch, edgeI)
        {
            if (ngbPolyPatch[edgeI] != -1)
            {
                if (isA<processorPolyPatch>(pbm[ngbPolyPatch[edgeI]]))
                {
                    processorPatchSet.insert(ngbPolyPatch[edgeI]);
                }
            }
        }
        labelList processorPatches(processorPatchSet.toc());
        faPatches.setSize(faPatches.size() + processorPatches.size());

        for (label i = 0; i < processorPatches.size(); ++i)
        {
            SLList<label> tmpLst;

            forAll(ngbPolyPatch, eI)
            {
                if (ngbPolyPatch[eI] == processorPatches[i])
                {
                    tmpLst.append(allUndefEdges[eI]);
                }
            }

            faPatches[faPatchNames.size() + i].edgeLabels_ = tmpLst;

            faPatches[faPatchNames.size() + i].name_ =
                pbm[processorPatches[i]].name();

            faPatches[faPatchNames.size() + i].type_ =
                processorFaPatch::typeName;

            faPatches[faPatchNames.size() + i].ngbPolyPatchID_ =
                processorPatches[i];
        }

        // Remaining undefined edges
        SLList<label> undefEdges;
        forAll(ngbPolyPatch, eI)
        {
            if (ngbPolyPatch[eI] == -1)
            {
                undefEdges.append(allUndefEdges[eI]);
            }
            else if (!isA<processorPolyPatch>(pbm[ngbPolyPatch[eI]]))
            {
                undefEdges.append(allUndefEdges[eI]);
            }
        }

        if (undefEdges.size() > 0)
        {
            label pI = faPatches.size()-1;

            faPatches[pI].name_ = "undefined";
            faPatches[pI].type_ = "patch";
            faPatches[pI].edgeLabels_ = undefEdges;
        }
        else
        {
            faPatches.setSize(faPatches.size() - 1);
        }
    }
    else
    {
        faPatches.setSize(faPatches.size() - 1);
    }


    // Reorder processorFaPatch using
    // ordering of ngb processorPolyPatch
    forAll(faPatches, patchI)
    {
        if (faPatches[patchI].type_ == processorFaPatch::typeName)
        {
            labelList ngbFaces(faPatches[patchI].edgeLabels_.size(), -1);

            forAll(ngbFaces, edgeI)
            {
                label curEdge = faPatches[patchI].edgeLabels_[edgeI];

                label curPMeshEdge = meshEdges[curEdge];

                forAll(edgeFaces[curPMeshEdge], faceI)
                {
                    label curFace = edgeFaces[curPMeshEdge][faceI];

                    label curPatchID = pbm.whichPatch(curFace);

                    if (curPatchID == faPatches[patchI].ngbPolyPatchID_)
                    {
                        ngbFaces[edgeI] = curFace;
                    }
                }
            }

            SortableList<label> sortedNgbFaces(ngbFaces);
            labelList reorderedEdgeLabels(ngbFaces.size(), -1);
            for (label i = 0; i < reorderedEdgeLabels.size(); ++i)
            {
                reorderedEdgeLabels[i] =
                    faPatches[patchI].edgeLabels_
                    [
                        sortedNgbFaces.indices()[i]
                    ];
            }

            faPatches[patchI].edgeLabels_ = reorderedEdgeLabels;
        }
    }


    // Add good patches to faMesh
    SLList<faPatch*> faPatchLst;

    for (label pI = 0; pI < faPatches.size(); ++pI)
    {
        faPatches[pI].dict_.add("type", faPatches[pI].type_);
        faPatches[pI].dict_.add("edgeLabels", faPatches[pI].edgeLabels_);
        faPatches[pI].dict_.add
        (
            "ngbPolyPatchIndex",
            faPatches[pI].ngbPolyPatchID_
        );

        if (faPatches[pI].type_ == processorFaPatch::typeName)
        {
            if (faPatches[pI].ngbPolyPatchID_ == -1)
            {
                FatalErrorInFunction
                    << "ngbPolyPatch is not defined for processorFaPatch: "
                    << faPatches[pI].name_
                    << abort(FatalError);
            }

            const processorPolyPatch& procPolyPatch =
                refCast<const processorPolyPatch>
                (
                    pbm[faPatches[pI].ngbPolyPatchID_]
                );

            faPatches[pI].dict_.add("myProcNo", procPolyPatch.myProcNo());
            faPatches[pI].dict_.add
            (
                "neighbProcNo",
                procPolyPatch.neighbProcNo()
            );
        }

        faPatchLst.append
        (
            faPatch::New
            (
                faPatches[pI].name_,
                faPatches[pI].dict_,
                pI,
                boundary()
            ).ptr()
        );
    }

    addFaPatches(List<faPatch*>(faPatchLst));

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    if (isFile(mesh().time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
}


Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const label polyPatchID
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>(pMesh),
    edgeInterpolation(*this),
    faSchemes(mesh()),
    faSolution(mesh()),
    data(mesh()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        labelList(pMesh.boundaryMesh()[polyPatchID].size(), -1)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    DebugInFunction << "Creating faMesh from polyPatch" << endl;

    const polyBoundaryMesh& pbm = pMesh.boundaryMesh();

    // Set faceLabels
    forAll(faceLabels_, faceI)
    {
        faceLabels_[faceI] = pbm[polyPatchID].start() + faceI;
    }

    // Add one faPatch
    const indirectPrimitivePatch& bp = patch();

    const label nTotalEdges = bp.nEdges();

    const label nInternalEdges = bp.nInternalEdges();

    labelList edgeLabels(nTotalEdges - nInternalEdges, -1);

    forAll(edgeLabels, edgeI)
    {
        edgeLabels[edgeI] = nInternalEdges + edgeI;
    }

    dictionary patchDict;

    patchDict.add("type", "patch");
    patchDict.add("edgeLabels", edgeLabels);
    patchDict.add("ngbPolyPatchIndex", -1);

    List<faPatch*> faPatchLst(1);

    faPatchLst[0] = faPatch::New("default", patchDict, 0, boundary()).ptr();

    addFaPatches(faPatchLst);

    setPrimitiveMeshData();

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faMesh::~faMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::faMesh::meshDir() const
{
    return mesh().dbDir()/meshSubDir;
}


const Foam::Time& Foam::faMesh::time() const
{
    return mesh().time();
}


const Foam::fileName& Foam::faMesh::pointsInstance() const
{
    return mesh().pointsInstance();
}


const Foam::fileName& Foam::faMesh::facesInstance() const
{
    return mesh().facesInstance();
}


const Foam::indirectPrimitivePatch& Foam::faMesh::patch() const
{
    if (!patchPtr_)
    {
        patchPtr_ = new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh().faces(),
                faceLabels_
            ),
            mesh().points()
        );
    }

    return *patchPtr_;
}


Foam::indirectPrimitivePatch& Foam::faMesh::patch()
{
    if (!patchPtr_)
    {
        patchPtr_ = new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh().faces(),
                faceLabels_
            ),
            mesh().points()
        );
    }

    return *patchPtr_;
}


const Foam::pointField& Foam::faMesh::points() const
{
    return patch().localPoints();
}


const Foam::edgeList& Foam::faMesh::edges() const
{
    return edges_;
}


const Foam::faceList& Foam::faMesh::faces() const
{
    return patch().localFaces();
}


void Foam::faMesh::addFaPatches(const List<faPatch*>& p)
{
    DebugInFunction << "Adding patches to faMesh" << endl;

    if (boundary().size() > 0)
    {
        FatalErrorInFunction
            << "boundary already exists"
            << abort(FatalError);
    }

    boundary_.setSize(p.size());

    forAll(p, patchI)
    {
        boundary_.set(patchI, p[patchI]);
    }

    setPrimitiveMeshData();

    boundary_.checkDefinition();
}


Foam::label Foam::faMesh::comm() const
{
    return comm_;
}


Foam::label& Foam::faMesh::comm()
{
    return comm_;
}


const Foam::objectRegistry& Foam::faMesh::thisDb() const
{
    return mesh().thisDb();
}


const Foam::faBoundaryMesh& Foam::faMesh::boundary() const
{
    return boundary_;
}


const Foam::labelList& Foam::faMesh::patchStarts() const
{
    if (!patchStartsPtr_)
    {
        calcPatchStarts();
    }

    return *patchStartsPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::Le() const
{
    if (!LePtr_)
    {
        calcLe();
    }

    return *LePtr_;
}


const Foam::edgeScalarField& Foam::faMesh::magLe() const
{
    if (!magLePtr_)
    {
        calcMagLe();
    }

    return *magLePtr_;
}


const Foam::areaVectorField& Foam::faMesh::areaCentres() const
{
    if (!centresPtr_)
    {
        calcAreaCentres();
    }

    return *centresPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::edgeCentres() const
{
    if (!edgeCentresPtr_)
    {
        calcEdgeCentres();
    }

    return *edgeCentresPtr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S() const
{
    if (!SPtr_)
    {
        calcS();
    }

    return *SPtr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S0() const
{
    if (!S0Ptr_)
    {
        FatalErrorInFunction
            << "S0 is not available"
            << abort(FatalError);
    }

    return *S0Ptr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S00() const
{
    if (!S00Ptr_)
    {
        S00Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S00",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            S0()
        );

        S0Ptr_->writeOpt() = IOobject::AUTO_WRITE;
    }

    return *S00Ptr_;
}


const Foam::areaVectorField& Foam::faMesh::faceAreaNormals() const
{
    if (!faceAreaNormalsPtr_)
    {
        calcFaceAreaNormals();
    }

    return *faceAreaNormalsPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::edgeAreaNormals() const
{
    if (!edgeAreaNormalsPtr_)
    {
        calcEdgeAreaNormals();
    }

    return *edgeAreaNormalsPtr_;
}


const Foam::vectorField& Foam::faMesh::pointAreaNormals() const
{
    if (!pointAreaNormalsPtr_)
    {
        calcPointAreaNormals();

        if (quadricsFit_ > 0)
        {
            calcPointAreaNormalsByQuadricsFit();
        }
    }

    return *pointAreaNormalsPtr_;
}


const Foam::areaScalarField& Foam::faMesh::faceCurvatures() const
{
    if (!faceCurvaturesPtr_)
    {
        calcFaceCurvatures();
    }

    return *faceCurvaturesPtr_;
}


const Foam::FieldField<Foam::Field, Foam::tensor>&
Foam::faMesh::edgeTransformTensors() const
{
    if (!edgeTransformTensorsPtr_)
    {
        calcEdgeTransformTensors();
    }

    return *edgeTransformTensorsPtr_;
}


const Foam::faGlobalMeshData& Foam::faMesh::globalData() const
{
    if (!globalMeshDataPtr_)
    {
        globalMeshDataPtr_ = new faGlobalMeshData(*this);
    }

    return *globalMeshDataPtr_;
}


const Foam::lduAddressing& Foam::faMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        calcLduAddressing();
    }

    return *lduPtr_;
}


bool Foam::faMesh::movePoints()
{
    // Grab point motion from polyMesh
    const vectorField& newPoints = mesh().points();

    // Grab old time areas if the time has been incremented
    if (curTimeIndex_ < time().timeIndex())
    {
        if (S00Ptr_ && S0Ptr_)
        {
            DebugInfo<< "Copy old-old S" << endl;
            *S00Ptr_ = *S0Ptr_;
        }

        if (S0Ptr_)
        {
            DebugInfo<< "Copy old S" << endl;
            *S0Ptr_ = S();
        }
        else
        {
            DebugInfo<< "Creating old cell volumes." << endl;

            S0Ptr_ = new DimensionedField<scalar, areaMesh>
            (
                IOobject
                (
                    "S0",
                    time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                S()
            );
        }

        curTimeIndex_ = time().timeIndex();
    }

    clearGeomNotAreas();

    // To satisfy the motion interface for MeshObject, const cast is needed
    if (patchPtr_)
    {
        patchPtr_->movePoints(newPoints);
    }

    // Move boundary points
    const_cast<faBoundaryMesh&>(boundary_).movePoints(newPoints);

    // Move interpolation
    const edgeInterpolation& cei = *this;
    const_cast<edgeInterpolation&>(cei).edgeInterpolation::movePoints();

    // Note: Fluxes were dummy?

    return true;
}


bool Foam::faMesh::correctPatchPointNormals(const label patchID) const
{
    if ((patchID < 0) || (patchID >= boundary().size()))
    {
        FatalErrorInFunction
            << "patchID is not in the valid range"
            << abort(FatalError);
    }

    if (correctPatchPointNormalsPtr_)
    {
        return (*correctPatchPointNormalsPtr_)[patchID];
    }

    return false;
}


Foam::boolList& Foam::faMesh::correctPatchPointNormals() const
{
    if (!correctPatchPointNormalsPtr_)
    {
        correctPatchPointNormalsPtr_ = new boolList(boundary().size(), false);
    }

    return *correctPatchPointNormalsPtr_;
}


bool Foam::faMesh::write(const bool valid) const
{
    faceLabels_.write();
    boundary_.write();

    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::faMesh::operator!=(const faMesh& m) const
{
    return &m != this;
}


bool Foam::faMesh::operator==(const faMesh& m) const
{
    return &m == this;
}


// ************************************************************************* //
