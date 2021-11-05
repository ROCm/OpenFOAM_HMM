/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
#include "faMeshBoundaryHalo.H"
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faMesh, 0);
}


const Foam::word Foam::faMesh::prefix("finite-area");

Foam::word Foam::faMesh::meshSubDir = "faMesh";

int Foam::faMesh::origPointAreaMethod_ = 0;  // Tuning

const int Foam::faMesh::quadricsFit_ = 0;  // Tuning


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Convert patch names to face labels. Preserve patch order
static labelList selectPatchFaces
(
    const polyBoundaryMesh& pbm,
    const wordRes& polyPatchNames
)
{
    const labelList patchIDs
    (
        pbm.patchSet
        (
            polyPatchNames,
            false,  // warnNotFound
            true    // useGroups
        ).sortedToc()
    );

    if (patchIDs.empty())
    {
        FatalErrorInFunction
            << "No matching patches: " << polyPatchNames << nl
            << exit(FatalError);
    }

    label nFaceLabels = 0;
    for (const label patchi : patchIDs)
    {
        nFaceLabels += pbm[patchi].size();
    }

    labelList faceLabels(nFaceLabels);

    nFaceLabels = 0;
    for (const label patchi : patchIDs)
    {
        for (const label facei : pbm[patchi].range())
        {
            faceLabels[nFaceLabels] = facei;
            ++nFaceLabels;
        }
    }

    return faceLabels;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMesh::checkBoundaryEdgeLabelRange
(
    const labelUList& edgeLabels
) const
{
    label nErrors = 0;

    for (const label edgei : edgeLabels)
    {
        if (edgei < nInternalEdges_ || edgei >= nEdges_)
        {
            if (!nErrors++)
            {
                FatalErrorInFunction
                    << "Boundary edge label out of range "
                    << nInternalEdges_ << ".." << (nEdges_-1) << nl
                    << "   ";
            }

            FatalError<< ' ' << edgei;
        }
    }

    if (nErrors)
    {
        FatalError << nl << exit(FatalError);
    }
}


void Foam::faMesh::initPatch() const
{
    patchPtr_.reset
    (
        new uindirectPrimitivePatch
        (
            UIndirectList<face>(mesh().faces(), faceLabels_),
            mesh().points()
        )
    );
    bndConnectPtr_.reset(nullptr);
    haloMapPtr_.reset(nullptr);
    haloFaceCentresPtr_.reset(nullptr);
    haloFaceNormalsPtr_.reset(nullptr);
}


void Foam::faMesh::setPrimitiveMeshData()
{
    DebugInFunction << "Setting primitive data" << endl;

    const uindirectPrimitivePatch& bp = patch();
    const labelListList& edgeFaces = bp.edgeFaces();

    // Dimensions

    nEdges_ = bp.nEdges();
    nInternalEdges_ = bp.nInternalEdges();
    nFaces_ = bp.size();
    nPoints_ = bp.nPoints();

    edges_.resize(nEdges_);
    edgeOwner_.resize(nEdges_);
    edgeNeighbour_.resize(nInternalEdges_);

    // Internal edges
    for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
    {
        edges_[edgei] = bp.edges()[edgei];

        edgeOwner_[edgei] = edgeFaces[edgei][0];

        edgeNeighbour_[edgei] = edgeFaces[edgei][1];
    }

    // Continue with boundary edges
    label edgei = nInternalEdges_;

    for (const faPatch& p : boundary())
    {
        for (const label patchEdgei : p.edgeLabels())
        {
            edges_[edgei] = bp.edges()[patchEdgei];

            edgeOwner_[edgei] = edgeFaces[patchEdgei][0];

            ++edgei;
        }
    }
}


void Foam::faMesh::clearHalo() const
{
    DebugInFunction << "Clearing halo information" << endl;

    haloMapPtr_.reset(nullptr);
    haloFaceCentresPtr_.reset(nullptr);
    haloFaceNormalsPtr_.reset(nullptr);
}


void Foam::faMesh::clearGeomNotAreas() const
{
    DebugInFunction << "Clearing geometry" << endl;

    clearHalo();
    patchPtr_.reset(nullptr);
    bndConnectPtr_.reset(nullptr);
    deleteDemandDrivenData(SPtr_);
    deleteDemandDrivenData(patchStartsPtr_);
    deleteDemandDrivenData(LePtr_);
    deleteDemandDrivenData(magLePtr_);
    deleteDemandDrivenData(centresPtr_);
    deleteDemandDrivenData(edgeCentresPtr_);
    deleteDemandDrivenData(faceAreaNormalsPtr_);
    deleteDemandDrivenData(edgeAreaNormalsPtr_);
    pointAreaNormalsPtr_.reset(nullptr);
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
    globalMeshDataPtr_.reset(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMesh::faMesh(const polyMesh& pMesh, const zero)
:
    faMesh(pMesh, labelList())
{}


Foam::faMesh::faMesh(const polyMesh& pMesh)
:
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
            faMesh::meshSubDir,
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
            faMesh::meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    bndConnectPtr_(nullptr),
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
    globalMeshDataPtr_(nullptr),

    haloMapPtr_(nullptr),
    haloFaceCentresPtr_(nullptr),
    haloFaceNormalsPtr_(nullptr)
{
    DebugInFunction << "Creating from IOobject" << endl;

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

    if (fileHandler().isFile(pMesh.time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                faMesh::meshSubDir,
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
    const UList<label>& faceLabels
)
:
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
            faMesh::meshSubDir,
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
            faMesh::meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        label(0)
    ),
    comm_(Pstream::worldComm),
    patchPtr_(nullptr),
    bndConnectPtr_(nullptr),
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
    globalMeshDataPtr_(nullptr),

    haloMapPtr_(nullptr),
    haloFaceCentresPtr_(nullptr),
    haloFaceNormalsPtr_(nullptr)
{}


Foam::faMesh::faMesh(const polyPatch& pp)
:
    faMesh
    (
        pp.boundaryMesh().mesh(),
        identity(pp.range())
    )
{
    DebugInFunction << "Creating from polyPatch:" << pp.name() << endl;

    // Add single faPatch "default", but with processor connections
    PtrList<faPatch> newPatches
    (
        createOnePatch("default")
    );

    addFaPatches(newPatches);

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
}


Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const dictionary& faMeshDefinition
)
:
    faMesh
    (
        pMesh,
        selectPatchFaces
        (
            pMesh.boundaryMesh(),
            faMeshDefinition.get<wordRes>("polyMeshPatches")
        )
    )
{
    DebugInFunction << "Creating from definition (dictionary)" << endl;

    PtrList<faPatch> newPatches
    (
        createPatchList
        (
            faMeshDefinition.subDict("boundary"),

            // Optional 'empty' patch
            faMeshDefinition.getOrDefault<word>("emptyPatch", word::null),

            // Optional specification for default patch
            faMeshDefinition.findDict("defaultPatch")
        )
    );


    addFaPatches(newPatches);

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    if (fileHandler().isFile(pMesh.time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                faMesh::meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faMesh::~faMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::faMesh::meshDir() const
{
    return mesh().dbDir()/faMesh::meshSubDir;
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


bool Foam::faMesh::hasDb() const
{
    return true;
}


const Foam::objectRegistry& Foam::faMesh::thisDb() const
{
    return mesh().thisDb();
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

        S0Ptr_->writeOpt(IOobject::AUTO_WRITE);
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
        pointAreaNormalsPtr_.reset(new vectorField(nPoints()));

        if (origPointAreaMethod_)
        {
            calcPointAreaNormals_orig(*pointAreaNormalsPtr_);
        }
        else
        {
            calcPointAreaNormals(*pointAreaNormalsPtr_);
        }

        if (quadricsFit_ > 0)
        {
            calcPointAreaNormalsByQuadricsFit(*pointAreaNormalsPtr_);
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
        globalMeshDataPtr_.reset(new faGlobalMeshData(*this));
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
