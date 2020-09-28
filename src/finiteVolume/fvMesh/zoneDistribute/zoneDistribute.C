/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "zoneDistribute.H"
#include "dummyTransform.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "wedgePolyPatch.H"

#include "globalPoints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneDistribute, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::zoneDistribute::coupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<processorPolyPatch>(pp))
        {
            nCoupled += pp.size();
        }
    }
    labelList coupledFaces(nCoupled);
    nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<processorPolyPatch>(pp))
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                coupledFaces[nCoupled++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            coupledFaces
        ),
        mesh_.points()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneDistribute::zoneDistribute(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, zoneDistribute>(mesh),
    stencil_(mesh),
    coupledBoundaryPoints_(coupledFacesPatch()().meshPoints()),
    send_(Pstream::nProcs())
{
}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

Foam::zoneDistribute& Foam::zoneDistribute::New(const fvMesh& mesh)
{
    zoneDistribute* ptr = mesh.thisDb().getObjectPtr<zoneDistribute>
    (
        zoneDistribute::typeName
    );

    if (!ptr)
    {
        ptr = new zoneDistribute(mesh);

        regIOobject::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneDistribute::updateStencil(const boolList& zone)
{
    stencil_.updateStencil(zone);
}


void Foam::zoneDistribute::setUpCommforZone
(
    const boolList& zone,
    bool updateStencil
)
{
    if (updateStencil)
    {
        stencil_.updateStencil(zone);
    }

    const labelHashSet comms = stencil_.needsComm();

    List<labelHashSet> needed_(Pstream::nProcs());

    if (Pstream::parRun())
    {
        for (const label celli : comms)
        {
            if (zone[celli])
            {
                for (const label gblIdx : stencil_[celli])
                {
                    if (!globalNumbering().isLocal(gblIdx))
                    {
                        const label procID =
                            globalNumbering().whichProcID(gblIdx);
                        needed_[procID].insert(gblIdx);
                    }
                }
            }
        }

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Stream data into buffer
        for (const int domain : Pstream::allProcs())
        {
            if (domain != Pstream::myProcNo())
            {
                // Put data into send buffer
                UOPstream toDomain(domain, pBufs);

                toDomain << needed_[domain];
            }
        }

        // wait until everything is written.
        pBufs.finishedSends();

        for (const int domain : Pstream::allProcs())
        {
            send_[domain].clear();

            if (domain != Pstream::myProcNo())
            {
                // get data from send buffer
                UIPstream fromDomain(domain, pBufs);

                fromDomain >> send_[domain];
            }
        }
    }
}


void Foam::zoneDistribute::updateMesh(const mapPolyMesh& mpm)
{
    if (mesh_.topoChanging())
    {
        coupledBoundaryPoints_ = coupledFacesPatch()().meshPoints();
    }
}


// ************************************************************************* //
