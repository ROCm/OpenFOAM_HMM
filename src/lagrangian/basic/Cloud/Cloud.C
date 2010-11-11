/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "Cloud.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "mapPolyMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
const Foam::scalar Foam::Cloud<ParticleType>::trackingCorrectionTol = 1e-5;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::calcCellWallFaces() const
{
    cellWallFacesPtr_.reset(new PackedBoolList(pMesh().nCells(), false));

    PackedBoolList& cellWallFaces = cellWallFacesPtr_();

    const polyBoundaryMesh& patches = pMesh().boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<wallPolyPatch>(patches[patchI]))
        {
            const polyPatch& patch = patches[patchI];

            const labelList& pFaceCells = patch.faceCells();

            forAll(pFaceCells, pFCI)
            {
                cellWallFaces[pFaceCells[pFCI]] = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    particleCount_(0),
    labels_(),
    cellTree_(),
    nTrackingRescues_(),
    cellWallFacesPtr_()
{
    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    IDLList<ParticleType>::operator=(particles);
}


template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh, cloudName),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    particleCount_(0),
    labels_(),
    cellTree_(),
    nTrackingRescues_(),
    cellWallFacesPtr_()
{
    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    IDLList<ParticleType>::operator=(particles);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::findCellFacePt
(
    const point& pt,
    label& cellI,
    label& tetFaceI,
    label& tetPtI
) const
{
    cellI = -1;
    tetFaceI = -1;
    tetPtI = -1;

    const indexedOctree<treeDataCell>& tree = cellTree();

    // Find nearest cell to the point

    pointIndexHit info = tree.findNearest(pt, sqr(GREAT));

    if (info.hit())
    {
        label nearestCellI = tree.shapes().cellLabels()[info.index()];

        // Check the nearest cell to see if the point is inside.
        findFacePt(nearestCellI, pt, tetFaceI, tetPtI);

        if (tetFaceI != -1)
        {
            // Point was in the nearest cell

            cellI = nearestCellI;

            return;
        }
        else
        {
            // Check the other possible cells that the point may be in

            labelList testCells = tree.findIndices(pt);

            forAll(testCells, pCI)
            {
                label testCellI = tree.shapes().cellLabels()[testCells[pCI]];

                if (testCellI == nearestCellI)
                {
                    // Don't retest the nearest cell

                    continue;
                }

                // Check the test cell to see if the point is inside.
                findFacePt(testCellI, pt, tetFaceI, tetPtI);

                if (tetFaceI != -1)
                {
                    // Point was in the test cell

                    cellI = testCellI;

                    return;
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::Cloud<ParticleType>::findCellFacePt"
            "("
                "const point& pt, "
                "label& cellI, "
                "label& tetFaceI, "
                "label& tetPtI"
            ") const"
        )   << "Did not find nearest cell in search tree."
            << abort(FatalError);
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::findFacePt
(
    label cellI,
    const point& pt,
    label& tetFaceI,
    label& tetPtI
) const
{
    tetFaceI = -1;
    tetPtI = -1;

    List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
    (
        polyMesh_,
        cellI
    );

    forAll(cellTets, tetI)
    {
        const tetIndices& cellTetIs = cellTets[tetI];

        if (inTet(pt, cellTetIs.tet(polyMesh_)))
        {
            tetFaceI = cellTetIs.face();
            tetPtI = cellTetIs.tetPt();

            return;
        }
    }
}


template<class ParticleType>
bool Foam::Cloud<ParticleType>::inTet
(
    const point& pt,
    const tetPointRef& tet
) const
{
    // For robustness, assuming that the point is in the tet unless
    // "definitively" shown otherwise by obtaining a positive dot
    // product greater than a tolerance of SMALL.

    // The tet is defined: tet(Cc, tetBasePt, pA, pB) where the normal
    // vectors and base points for the half-space planes are:
    // area[0] = tet.Sa();
    // area[1] = tet.Sb();
    // area[2] = tet.Sc();
    // area[3] = tet.Sd();
    // planeBase[0] = tetBasePt = tet.b()
    // planeBase[1] = ptA       = tet.c()
    // planeBase[2] = tetBasePt = tet.b()
    // planeBase[3] = tetBasePt = tet.b()

    vector n = vector::zero;

    {
        // 0, a
        const point& basePt = tet.b();

        n = tet.Sa();
        n /= (mag(n) + VSMALL);

        if (((pt - basePt) & n) > SMALL)
        {
            return false;
        }
    }

    {
        // 1, b
        const point& basePt = tet.c();

        n = tet.Sb();
        n /= (mag(n) + VSMALL);

        if (((pt - basePt) & n) > SMALL)
        {
            return false;
        }
    }

    {
        // 2, c
        const point& basePt = tet.b();

        n = tet.Sc();
        n /= (mag(n) + VSMALL);

        if (((pt - basePt) & n) > SMALL)
        {
            return false;
        }
    }

    {
        // 3, d
        const point& basePt = tet.b();

        n = tet.Sd();
        n /= (mag(n) + VSMALL);

        if (((pt - basePt) & n) > SMALL)
        {
            return false;
        }
    }

    return true;
}


template<class ParticleType>
const Foam::indexedOctree<Foam::treeDataCell>&
Foam::Cloud<ParticleType>::cellTree() const
{
    if (cellTree_.empty())
    {
        treeBoundBox overallBb(polyMesh_.points());

        Random rndGen(261782);

        overallBb = overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        cellTree_.reset
        (
            new indexedOctree<treeDataCell>
            (
                treeDataCell
                (
                    false,      // not cache bb
                    polyMesh_
                ),
                overallBb,
                8,              // maxLevel
                10,             // leafsize
                3.0             // duplicity
            )
        );
    }

    return cellTree_();
}


template<class ParticleType>
const Foam::PackedBoolList& Foam::Cloud<ParticleType>::cellHasWallFaces()
const
{
    if (!cellWallFacesPtr_.valid())
    {
        calcCellWallFaces();
    }

    return cellWallFacesPtr_();
}



template<class ParticleType>
Foam::label Foam::Cloud<ParticleType>::getNewParticleID() const
{
    label id = particleCount_++;

    if (id == labelMax)
    {
        WarningIn("Cloud<ParticleType>::getNewParticleID() const")
            << "Particle counter has overflowed. This might cause problems"
            << " when reconstructing particle tracks." << endl;
    }
    return id;
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::addParticle(ParticleType* pPtr)
{
    append(pPtr);
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteParticle(ParticleType& p)
{
    delete(this->remove(&p));
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::cloudReset(const Cloud<ParticleType>& c)
{
    // Reset particle cound and particles only
    // - not changing the cloud object registry or reference to the polyMesh
    particleCount_ = 0;
    IDLList<ParticleType>::operator=(c);
}


template<class ParticleType>
template<class TrackingData>
void Foam::Cloud<ParticleType>::move(TrackingData& td, const scalar trackTime)
{
    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();
    const globalMeshData& pData = polyMesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Indexing of patches into the procPatches list
    const labelList& procPatchIndices = pData.processorPatchIndices();

    // Indexing of equivalent patch on neighbour processor into the
    // procPatches list on the neighbour
    const labelList& procPatchNeighbours = pData.processorPatchNeighbours();

    // Which processors this processor is connected to
    const labelList& neighbourProcs = pData[Pstream::myProcNo()];

    // Indexing from the processor number into the neighbourProcs list
    labelList neighbourProcIndices(Pstream::nProcs(), -1);

    forAll(neighbourProcs, i)
    {
        neighbourProcIndices[neighbourProcs[i]] = i;
    }

    // Initialise the stepFraction moved for the particles
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        pIter().stepFraction() = 0;
    }

    // Reset nTrackingRescues
    nTrackingRescues_ = 0;

    // While there are particles to transfer
    while (true)
    {
        // List of lists of particles to be transfered for all of the
        // neighbour processors
        List<IDLList<ParticleType> > particleTransferLists
        (
            neighbourProcs.size()
        );

        // List of destination processorPatches indices for all of the
        // neighbour processors
        List<DynamicList<label> > patchIndexTransferLists
        (
            neighbourProcs.size()
        );

        // Loop over all particles
        forAllIter(typename Cloud<ParticleType>, *this, pIter)
        {
            ParticleType& p = pIter();

            // Move the particle
            bool keepParticle = p.move(td, trackTime);

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                // If we are running in parallel and the particle is on a
                // boundary face
                if (Pstream::parRun() && p.faceI_ >= pMesh().nInternalFaces())
                {
                    label patchI = pbm.whichPatch(p.faceI_);

                    // ... and the face is on a processor patch
                    // prepare it for transfer
                    if (procPatchIndices[patchI] != -1)
                    {
                        label n = neighbourProcIndices
                        [
                            refCast<const processorPolyPatch>
                            (
                                pbm[patchI]
                            ).neighbProcNo()
                        ];

                        p.prepareForParallelTransfer(patchI, td);

                        particleTransferLists[n].append(this->remove(&p));

                        patchIndexTransferLists[n].append
                        (
                            procPatchNeighbours[patchI]
                        );
                    }
                }
            }
            else
            {
                deleteParticle(p);
            }
        }

        if (!Pstream::parRun())
        {
            break;
        }

        // Allocate transfer buffers
        PstreamBuffers pBufs(Pstream::nonBlocking);

        // Stream into send buffers
        forAll(particleTransferLists, i)
        {
            if (particleTransferLists[i].size())
            {
                UOPstream particleStream
                (
                    neighbourProcs[i],
                    pBufs
                );

                particleStream
                    << patchIndexTransferLists[i]
                    << particleTransferLists[i];
            }
        }

        // Set up transfers when in non-blocking mode. Returns sizes (in bytes)
        // to be sent/received.
        labelListList allNTrans(Pstream::nProcs());

        pBufs.finishedSends(allNTrans);

        bool transfered = false;

        forAll(allNTrans, i)
        {
            forAll(allNTrans[i], j)
            {
                if (allNTrans[i][j])
                {
                    transfered = true;
                    break;
                }
            }
        }

        if (!transfered)
        {
            break;
        }

        // Retrieve from receive buffers
        forAll(neighbourProcs, i)
        {
            label neighbProci = neighbourProcs[i];

            label nRec = allNTrans[neighbProci][Pstream::myProcNo()];

            if (nRec)
            {
                UIPstream particleStream(neighbProci, pBufs);

                labelList receivePatchIndex(particleStream);

                IDLList<ParticleType> newParticles
                (
                    particleStream,
                    typename ParticleType::iNew(*this)
                );

                label pI = 0;

                forAllIter(typename Cloud<ParticleType>, newParticles, newpIter)
                {
                    ParticleType& newp = newpIter();

                    label patchI = procPatches[receivePatchIndex[pI++]];

                    newp.correctAfterParallelTransfer(patchI, td);

                    addParticle(newParticles.remove(&newp));
                }
            }
        }
    }

    reduce(nTrackingRescues_, sumOp<label>());

    if (nTrackingRescues_ > 0)
    {
        Info<< nTrackingRescues_ << " tracking rescue corrections" << endl;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::autoMap(const mapPolyMesh& mapper)
{
    if (cloud::debug)
    {
        Info<< "Cloud<ParticleType>::autoMap(const morphFieldMapper& map) "
               "for lagrangian cloud " << cloud::name() << endl;
    }

    const labelList& reverseCellMap = mapper.reverseCellMap();
    const labelList& reverseFaceMap = mapper.reverseFaceMap();

    // Reset stored data that relies on the mesh
    cellTree_.clear();
    cellWallFacesPtr_.clear();

    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        if (reverseCellMap[pIter().cellI_] >= 0)
        {
            pIter().cellI_ = reverseCellMap[pIter().cellI_];

            if (pIter().faceI_ >= 0 && reverseFaceMap[pIter().faceI_] >= 0)
            {
                pIter().faceI_ = reverseFaceMap[pIter().faceI_];
            }
            else
            {
                pIter().faceI_ = -1;
            }

            pIter().initCellFacePt();
        }
        else
        {
            label trackStartCell = mapper.mergedCell(pIter().cellI_);

            if (trackStartCell < 0)
            {
                trackStartCell = 0;
            }

            vector p = pIter().position();

            const_cast<vector&>(pIter().position()) =
                polyMesh_.cellCentres()[trackStartCell];

            pIter().stepFraction() = 0;

            pIter().initCellFacePt();

            pIter().track(p);
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/this->name() + "_positions.obj"
    );

    forAllConstIter(typename Cloud<ParticleType>, *this, pIter)
    {
        const ParticleType& p = pIter();
        pObj<< "v " << p.position().x() << " " << p.position().y() << " "
            << p.position().z() << nl;
    }

    pObj.flush();
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudIO.C"

// ************************************************************************* //
