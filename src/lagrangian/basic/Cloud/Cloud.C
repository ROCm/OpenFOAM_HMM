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
    particleCount_(0)
{
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
    particleCount_(0)
{
    IDLList<ParticleType>::operator=(particles);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
template<class TrackingData>
void Foam::Cloud<ParticleType>::move(TrackingData& td)
{
    const globalMeshData& pData = polyMesh_.globalData();
    const labelList& neighbourProcs = pData[Pstream::myProcNo()];
    const labelList& procPatches = pData.processorPatches();
    const labelList& procPatchIndices = pData.processorPatchIndices();
    const labelList& procPatchNeighbours = pData.processorPatchNeighbours();
    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();

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
            bool keepParticle = p.move(td);

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                // If we are running in parallel and the particle is on a
                // boundary face
                if (Pstream::parRun() && p.facei_ >= pMesh().nInternalFaces())
                {
                    label patchi = pbm.whichPatch(p.facei_);

                    // ... and the face is on a processor patch
                    // prepare it for transfer
                    if (procPatchIndices[patchi] != -1)
                    {
                        label n = neighbourProcIndices
                        [
                            refCast<const processorPolyPatch>
                            (
                                pbm[patchi]
                            ).neighbProcNo()
                        ];

                        p.prepareForParallelTransfer(patchi, td);

                        particleTransferLists[n].append(this->remove(&p));

                        patchIndexTransferLists[n].append
                        (
                            procPatchNeighbours[patchi]
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
                    << labelList(patchIndexTransferLists[i])
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

                    label patchi = procPatches[receivePatchIndex[pI++]];

                    newp.correctAfterParallelTransfer(patchi, td);

                    addParticle(newParticles.remove(&newp));
                }
            }
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
