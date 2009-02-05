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

#include "surfMesh.H"
#include "Time.H"
#include "cellIOList.H"
#include "SubList.H"
#include "OSspecific.H"
#include "MeshedSurface.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::surfMesh, 0);

Foam::word Foam::surfMesh::meshSubDir = "surfMesh";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMesh::oneZone()
{
    word zoneName;

    if (surfZones_.size())
    {
        zoneName = surfZones_[0].name();
    }
    if (zoneName.empty())
    {
        zoneName = "zone0";
    }

    // set single default zone
    surfZones_.setSize(1);
    surfZones_[0] = surfZone
    (
        zoneName,
        nFaces(),       // zone size
        0,              // zone start
        0               // zone index
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::surfMesh::surfMesh(const IOobject& io, const word& surfName)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    surfMeshAllocator
    (
        IOobject
        (
            "points",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        IOobject
        (
            "faces",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    MeshReference(storedFaces_, storedPoints_),
    surfZones_
    (
        IOobject
        (
            "surfZones",
            time().findInstance(meshDir(), "surfZones"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{}


Foam::surfMesh::surfMesh
(
    const IOobject& io,
    const Xfer<pointField>& pointLst,
    const Xfer<faceList>& faceLst,
    const word& surfName
)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    surfMeshAllocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointLst,
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        faceLst
    ),
    MeshReference(storedFaces_, storedPoints_),
    surfZones_
    (
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    )
{}


Foam::surfMesh::surfMesh
(
    const IOobject& io,
    const Xfer< MeshedSurface<face> >& surf,
    const word& surfName
)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    surfMeshAllocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointField(),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        faceList()
    ),
    MeshReference(storedFaces_, storedPoints_),
    surfZones_
    (
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        surfZoneList()
    )
{
    // We can also send Xfer<..>::null just to force initialization
    if (&surf)
    {
        transfer(surf());
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMesh::~surfMesh()
{
    //    clearOut();
    //    resetMotion();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::resetPrimitives
(
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<surfZoneList>& zones,
    const bool validate
)
{
    // Clear addressing.
    MeshReference::clearGeom();

    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (&points)
    {
        storedPoints_.transfer(points());
    }

    if (&faces)
    {
        storedFaces_.transfer(faces());
    }

    if (&zones)
    {
        surfZones_.transfer(zones());
    }

    if (validate)
    {
        checkZones();
    }
}


void Foam::surfMesh::transfer
(
    MeshedSurface<face>& surf
)
{
    // Clear addressing.
    MeshReference::clearGeom();

    storedPoints_.transfer(surf.storedPoints());
    storedFaces_.transfer(surf.storedFaces());
    surfZones_.transfer(surf.storedZones());
}


void Foam::surfMesh::rename(const word& newName)
{
    FatalErrorIn
    (
        "surfMesh::rename(const word&)\n"
    )
        << "rename does not work correctly\n";
}



Foam::fileName Foam::surfMesh::meshDir() const
{
    return dbDir()/meshSubDir;
}


const Foam::fileName& Foam::surfMesh::pointsInstance() const
{
    return storedPoints_.instance();
}


const Foam::fileName& Foam::surfMesh::facesInstance() const
{
    return storedFaces_.instance();
}


Foam::label Foam::surfMesh::nPoints() const
{
    return storedPoints_.size();
}

Foam::label Foam::surfMesh::nFaces() const
{
    return storedFaces_.size();
}

const Foam::pointField& Foam::surfMesh::points() const
{
    return storedPoints_;
}

const Foam::faceList& Foam::surfMesh::faces() const
{
    return storedFaces_;
}

void Foam::surfMesh::checkZones()
{
    // extra safety, ensure we have at some zones
    // and they cover all the faces - fix start silently
    if (surfZones_.size() <= 1)
    {
        oneZone();
    }
    else
    {
        label count = 0;
        forAll(surfZones_, zoneI)
        {
            surfZones_[zoneI].start() = count;
            count += surfZones_[zoneI].size();
        }

        if (count < nFaces())
        {
            WarningIn
            (
                "surfMesh::checkZones()\n"
            )
                << "more faces " << nFaces() << " than zones " << count
                << " ... extending final zone"
                << endl;

            surfZones_[surfZones_.size()-1].size() += count - nFaces();
        }
        else if (count > size())
        {
            FatalErrorIn
            (
                "surfMesh::checkZones()\n"
            )
                << "more zones " << count << " than faces " << nFaces()
                << exit(FatalError);
        }
    }
}


// Add boundary patches. Constructor helper
void Foam::surfMesh::addZones
(
    const surfZoneList& zones,
    const bool validate
)
{
    surfZones_.setSize(zones.size());

    forAll(surfZones_, zoneI)
    {
        surfZones_[zoneI] = surfZone(zones[zoneI], zoneI);
    }

    if (validate)
    {
        checkZones();
    }
}


// Remove all files and some subdirs (eg, sets)
void Foam::surfMesh::removeFiles(const fileName& instanceDir) const
{
    fileName meshFilesPath = db().path()/instanceDir/meshSubDir;

    rm(meshFilesPath/"points");
    rm(meshFilesPath/"faces");
    rm(meshFilesPath/"surfZones");
}

void Foam::surfMesh::removeFiles() const
{
    removeFiles(instance());
}

// ************************************************************************* //
