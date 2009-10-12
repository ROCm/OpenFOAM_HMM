/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Cloud.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "mapPolyMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
const Foam::scalar Foam::Cloud<ParticleType>::minValidTrackFraction = 1e-6;

template<class ParticleType>
const Foam::scalar Foam::Cloud<ParticleType>::trackingRescueTolerance = 1e-4;

template<class ParticleType>
const Foam::scalar Foam::Cloud<ParticleType>::intersectionTolerance = 0;

template<class ParticleType>
const Foam::scalar Foam::Cloud<ParticleType>::planarCosAngle = (1 - 1e-6);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// Calculate using face properties only
template<class ParticleType>
bool Foam::Cloud<ParticleType>::isConcaveCell(const label cellI) const
{
    const cellList& cells = pMesh().cells();

    const labelList& fOwner = pMesh().faceOwner();
    const vectorField& fAreas = pMesh().faceAreas();
    const pointField& fCentres = pMesh().faceCentres();
    const pointField& cCentres = pMesh().cellCentres();
    const pointField& points = pMesh().points();

    const cell& cFaces = cells[cellI];

    bool concave = false;

    forAll(cFaces, i)
    {
        label f0I = cFaces[i];

        // Get f0 centre
        const point& fc0 = fCentres[f0I];

        const face& f0 = pMesh().faces()[f0I];

        vector fn0 = fAreas[f0I];
        fn0 /= (mag(fn0) + VSMALL);

        // Flip normal if required so that it is always pointing out of
        // the cell
        if (fOwner[f0I] != cellI)
        {
            fn0 *= -1;
        }

        // Check if the cell centre is visible from the plane of the face

        if (((fc0 - cCentres[cellI]) & fn0) < 0)
        {
            Pout<< "Face " << f0I
                << " is not visible from the centre of cell " << cellI
                << endl;
        }

        forAll(cFaces, j)
        {
            if (j != i)
            {
                label f1I = cFaces[j];

                const face& f1 = pMesh().faces()[f1I];

                // Is any vertex of f1 on wrong side of the plane of f0?

                forAll(f1, f1pI)
                {
                    label ptI = f1[f1pI];

                    // Skip points that are shared between f1 and f0
                    if (findIndex(f0, ptI) > -1)
                    {
                        continue;
                    }

                    const point& pt = points[ptI];

                    // If the cell is concave, the point on f1 will be
                    // on the outside of the plane of f0, defined by
                    // its centre and normal, and the angle between
                    // (fc0 -pt) and fn0 will be greater than 90
                    // degrees, so the dot product will be negative.

                    scalar d = ((fc0 - pt) & fn0);

                    if (d < 0)
                    {
                        // Concave face

                        concave = true;
                    }
                }

                // Check for co-planar faces, which are also treated
                // as concave, as they are not strictly convex.

                vector fn1 = fAreas[f1I];
                fn1 /= (mag(fn1) + VSMALL);

                // Flip normal if required so that it is always pointing out of
                // the cell
                if (fOwner[f1I] != cellI)
                {
                    fn1 *= -1;
                }

                if ((fn0 & fn1) > planarCosAngle)
                {
                    // Planar face

                    concave = true;
                }
            }
        }
    }

    return concave;
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::calcConcaveCells() const
{
    const cellList& cells = pMesh().cells();

    concaveCellPtr_.reset(new PackedBoolList(pMesh().nCells()));
    PackedBoolList& concaveCell = concaveCellPtr_();

    forAll(cells, cellI)
    {
        // if (isConcaveCell(cellI))
        // {
        //     concaveCell[cellI] = 1;
        // }

        concaveCell[cellI] = 1;
    }

    {
        // Write cells that are a problem to file

        DynamicList<label> tmpConcaveCells;

        forAll(cells, cellI)
        {
            if (concaveCell[cellI])
            {
                tmpConcaveCells.append(cellI);
            }
        }

        Pout<< "Cloud<ParticleType>::calcConcaveCells() :"
            << " overall cells in mesh : " << pMesh().nCells() << nl
            << "    of these concave : " << tmpConcaveCells.size() << endl;

        fileName fName = pMesh().time().path()/"concaveCells";

        Pout<< "    Writing " << fName.name() << endl;

        OFstream file(fName);

        file << tmpConcaveCells;

        file.flush();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const IDLList<ParticleType>& particles,
    const bool concaveCheck
)
:
    cloud(pMesh),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    concaveCheck_(concaveCheck),
    particleCount_(0)
{
    IDLList<ParticleType>::operator=(particles);
}


template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles,
    const bool concaveCheck
)
:
    cloud(pMesh, cloudName),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    concaveCheck_(concaveCheck),
    particleCount_(0)
{
    IDLList<ParticleType>::operator=(particles);
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::clearOut()
{
    concaveCellPtr_.clear();
    labels_.clearStorage();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
const Foam::PackedBoolList& Foam::Cloud<ParticleType>::concaveCell()
const
{
    if (!concaveCellPtr_.valid())
    {
        calcConcaveCells();
    }

    return concaveCellPtr_();
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


namespace Foam
{

class combineNsTransPs
{

public:

    void operator()(labelListList& x, const labelListList& y) const
    {
        forAll(y, i)
        {
            if (y[i].size())
            {
                x[i] = y[i];
            }
        }
    }
};

} // End namespace Foam


template<class ParticleType>
template<class TrackingData>
void Foam::Cloud<ParticleType>::move(TrackingData& td)
{
    const globalMeshData& pData = polyMesh_.globalData();
    const labelList& processorPatches = pData.processorPatches();
    const labelList& processorPatchIndices = pData.processorPatchIndices();
    const labelList& processorPatchNeighbours =
        pData.processorPatchNeighbours();

    // Initialise the setpFraction moved for the particles
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        pIter().stepFraction() = 0;
    }

    // Assume there will be particles to transfer
    bool transfered = true;

    // While there are particles to transfer
    while (transfered)
    {
        // List of lists of particles to be transfered for all the processor
        // patches
        List<IDLList<ParticleType> > transferList(processorPatches.size());

        // Loop over all particles
        forAllIter(typename Cloud<ParticleType>, *this, pIter)
        {
            ParticleType& p = pIter();

            // Move the particle
            bool keepParticle = p.move(td);

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                // If we are running in parallel and the particle is on a
                // boundary face
                if (Pstream::parRun() && p.facei_ >= pMesh().nInternalFaces())
                {
                    label patchi = pMesh().boundaryMesh().whichPatch(p.facei_);
                    label n = processorPatchIndices[patchi];

                    // ... and the face is on a processor patch
                    // prepare it for transfer
                    if (n != -1)
                    {
                        p.prepareForParallelTransfer(patchi, td);
                        transferList[n].append(this->remove(&p));
                    }
                }
            }
            else
            {
                deleteParticle(p);
            }
        }

        if (Pstream::parRun())
        {
            // List of the numbers of particles to be transfered across the
            // processor patches
            labelList nsTransPs(transferList.size());

            forAll(transferList, i)
            {
                nsTransPs[i] = transferList[i].size();
            }

            // List of the numbers of particles to be transfered across the
            // processor patches for all the processors
            labelListList allNTrans(Pstream::nProcs());
            allNTrans[Pstream::myProcNo()] = nsTransPs;
            combineReduce(allNTrans, combineNsTransPs());

            transfered = false;

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

            forAll(transferList, i)
            {
                if (transferList[i].size())
                {
                    OPstream particleStream
                    (
                        Pstream::blocking,
                        refCast<const processorPolyPatch>
                        (
                            pMesh().boundaryMesh()[processorPatches[i]]
                        ).neighbProcNo()
                    );

                    particleStream << transferList[i];
                }
            }

            forAll(processorPatches, i)
            {
                label patchi = processorPatches[i];

                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (pMesh().boundaryMesh()[patchi]);

                label neighbProci =
                    procPatch.neighbProcNo() - Pstream::masterNo();

                label neighbProcPatchi = processorPatchNeighbours[patchi];

                label nRecPs = allNTrans[neighbProci][neighbProcPatchi];

                if (nRecPs)
                {
                    IPstream particleStream
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    IDLList<ParticleType> newParticles
                    (
                        particleStream,
                        typename ParticleType::iNew(*this)
                    );

                    forAllIter
                    (
                        typename Cloud<ParticleType>,
                        newParticles,
                        newpIter
                    )
                    {
                        ParticleType& newp = newpIter();
                        newp.correctAfterParallelTransfer(patchi, td);
                        addParticle(newParticles.remove(&newp));
                    }
                }
            }
        }
        else
        {
            transfered = false;
        }
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

    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        if (reverseCellMap[pIter().celli_] >= 0)
        {
            pIter().celli_ = reverseCellMap[pIter().celli_];

            if (pIter().facei_ >= 0 && reverseFaceMap[pIter().facei_] >= 0)
            {
                pIter().facei_ = reverseFaceMap[pIter().facei_];
            }
            else
            {
                pIter().facei_ = -1;
            }
        }
        else
        {
            label trackStartCell = mapper.mergedCell(pIter().celli_);

            if (trackStartCell < 0)
            {
                trackStartCell = 0;
            }

            vector p = pIter().position();
            const_cast<vector&>(pIter().position()) =
                polyMesh_.cellCentres()[trackStartCell];
            pIter().stepFraction() = 0;
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
