/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "surfMesh.H"
#include "meshedSurf.H"
#include "MeshedSurfaceProxy.H"

#include "Time.H"
#include "OSspecific.H"
#include "MeshedSurface.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMesh, 0);
}

Foam::word Foam::surfMesh::meshSubDir = "surfMesh";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void Foam::surfMesh::oneZone()
// {
//     word zoneName;
//
//     if (surfZones_.size())
//     {
//         zoneName = surfZones_[0].name();
//     }
//     if (zoneName.empty())
//     {
//         zoneName = "zone0";
//     }
//
//     // Set single default zone
//     surfZones_.resize(1);
//     surfZones_[0] = surfZone
//     (
//         zoneName,
//         nFaces(),       // zone size
//         0,              // zone start
//         0               // zone index
//     );
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMesh::surfMesh(const IOobject& io)
:
    surfMesh(io, io.name())
{}


Foam::surfMesh::surfMesh(const IOobject& io, const word& surfName)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
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
    MeshReference(this->storedIOFaces(), this->storedIOPoints()),

    surfZones_
    (
        IOobject
        (
            "surfZones",
            time().findInstance(meshDir(), "surfZones"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}


Foam::surfMesh::surfMesh
(
    const word& surfName,
    const objectRegistry& obr
)
:
    surfaceRegistry(obr, surfName),
    Allocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints()),

    surfZones_
    (
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    )
{}


Foam::surfMesh::surfMesh
(
    const IOobject& io,
    const MeshedSurface<face>& surf,
    const word& surfName
)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            io.writeOpt()
        ),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            io.writeOpt()
        )
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints()),

    surfZones_
    (
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            io.writeOpt()
        )
    )
{
    DebugInfo
        <<"IOobject: " << io.path() << nl
        <<"  name: " << io.name()
        <<"  instance: " << io.instance()
        <<"  local: " << io.local()
        <<"  dbDir: " << io.db().dbDir() << nl
        <<"creating surfMesh at instance " << instance() << endl;

    copySurface(surf);
}


Foam::surfMesh::surfMesh
(
    const IOobject& io,
    MeshedSurface<face>&& surf,
    const word& surfName
)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints()),

    surfZones_
    (
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    )
{
    DebugInfo
        <<"IOobject: " << io.path() << nl
        <<" name: " << io.name()
        <<" instance: " << io.instance()
        <<" local: " << io.local()
        <<" dbDir: " << io.db().dbDir() << nl
        <<"creating surfMesh at instance " << instance() << nl
        <<"timeName: " << instance() << endl;

    transfer(surf);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMesh::~surfMesh()
{
    clearOut(); // Clear addressing
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::updatePointsRef()
{
    // Assign the reference to the points (quite ugly)
    // points() are returned as Field but are actually stored as SubField
    reinterpret_cast<typename MeshReference::PointFieldType&>
    (
        const_cast<Field<point>&>(MeshReference::points())
    ).shallowCopy
    (
        this->storedPoints()
    );
}


void Foam::surfMesh::updateFacesRef()
{
    // Assign the reference to the faces (UList)
    static_cast<MeshReference&>(*this).shallowCopy
    (
        this->storedFaces()
    );
}


void Foam::surfMesh::updateRefs()
{
    this->updatePointsRef();
    this->updateFacesRef();
}


void Foam::surfMesh::copySurface
(
    const pointField& points,
    const faceList& faces,
    bool validate
)
{
    clearOut(); // Clear addressing

    if
    (
        this->nPoints() != points.size()
     || this->nFaces() != faces.size()
    )
    {
        // Geometry changed
        clearFields();
    }

    this->storedIOPoints() = points;
    this->storedIOFaces() = faces;
    surfZones_.clear();

    this->updateRefs();

    // No zones
}


void Foam::surfMesh::copySurface
(
    const meshedSurf& surf,
    bool validate
)
{
    clearOut(); // Clear addressing

    if
    (
        this->nPoints() != surf.points().size()
     || this->nFaces() != surf.faces().size()
    )
    {
        // Geometry changed
        clearFields();
    }

    this->storedIOPoints() = surf.points();
    this->storedIOFaces() = surf.faces();
    surfZones_.clear();

    this->updateRefs();

    // No zones
}


void Foam::surfMesh::copySurface
(
    const MeshedSurface<face>& surf,
    bool validate
)
{
    clearOut(); // Clear addressing

    if
    (
        this->nPoints() != surf.points().size()
     || this->nFaces() != surf.surfFaces().size()
    )
    {
        // Geometry changed
        clearFields();
    }

    this->storedIOPoints() = surf.points();
    this->storedIOFaces() = surf.surfFaces();
    surfZones_ = surf.surfZones();

    this->updateRefs();

    if (validate)
    {
        checkZones();
    }
}


void Foam::surfMesh::transfer
(
    MeshedSurface<face>& surf,
    bool validate
)
{
    clearOut(); // Clear addressing

    this->storedPoints().transfer(surf.storedPoints());
    this->storedFaces().transfer(surf.storedFaces());
    this->storedZones().transfer(surf.storedZones());

    this->updateRefs();

    if (validate)
    {
        checkZones();
    }
}


Foam::autoPtr<Foam::MeshedSurface<Foam::face>>
Foam::surfMesh::releaseGeom()
{
    // Start with an empty geometry
    auto aptr = autoPtr<MeshedSurface<face>>::New();

    // Transfer in content
    aptr->storedPoints().transfer(this->storedPoints());
    aptr->storedFaces().transfer(this->storedFaces());
    aptr->storedZones().transfer(this->storedZones());

    this->updateRefs(); // This may not be needed...
    clearOut(); // Clear addressing.
    clearFields();

    return aptr;
}


Foam::fileName Foam::surfMesh::meshDir() const
{
    return dbDir()/meshSubDir;
}


const Foam::fileName& Foam::surfMesh::pointsInstance() const
{
    return this->storedIOPoints().instance();
}


const Foam::fileName& Foam::surfMesh::facesInstance() const
{
    return this->storedIOFaces().instance();
}


Foam::label Foam::surfMesh::nPoints() const
{
    return this->points().size();
}


Foam::label Foam::surfMesh::nFaces() const
{
    return this->faces().size();
}


const Foam::pointField& Foam::surfMesh::points() const
{
    return this->storedIOPoints();
}


const Foam::faceList& Foam::surfMesh::faces() const
{
    return this->storedIOFaces();
}


void Foam::surfMesh::checkZones()
{
    // Extra safety, ensure we have at some zones
    // and they cover all the faces - fix start silently

    if (surfZones_.size() <= 1)
    {
        removeZones();
        return;
    }

    label count = 0;
    for (surfZone& zn : surfZones_)
    {
        zn.start() = count;
        count += zn.size();
    }

    if (count < nFaces())
    {
        WarningInFunction
            << "More faces " << nFaces() << " than zones " << count
            << " ... extending final zone"
            << endl;

        surfZones_.last().size() += count - nFaces();
    }
    else if (size() < count)
    {
        FatalErrorInFunction
            << "More zones " << count << " than faces " << nFaces()
            << exit(FatalError);
    }
}


// Add boundary patches. Constructor helper
void Foam::surfMesh::addZones
(
    const surfZoneList& zones,
    const bool validate
)
{
    removeZones();

    forAll(surfZones_, zonei)
    {
        surfZones_[zonei] = surfZone(zones[zonei], zonei);
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


void Foam::surfMesh::write
(
    const fileName& name,
    const dictionary& options
) const
{
    write(name, name.ext(), options);
}


void Foam::surfMesh::write
(
    const fileName& name,
    const word& ext,
    const dictionary& options
) const
{
    MeshedSurfaceProxy<face>
    (
        this->points(),
        this->faces(),
        this->surfZones()
    ).write(name, ext, options);
}


// ************************************************************************* //
