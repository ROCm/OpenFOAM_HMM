/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "polySurface.H"
#include "Time.H"
#include "ModifiableMeshedSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polySurface, 0);
}

const Foam::word Foam::polySurface::pointDataName("PointData");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polySurface::calculateZoneIds(const UList<surfZone>& zones)
{
    if (returnReduce(zones.empty(), andOp<bool>()))
    {
        zoneIds_.clear();
        return;
    }

    // Extra safety, ensure we have at some zones
    // and they cover all the faces - fix start silently

    zoneIds_.resize(size(), Zero);

    label off = 0;
    for (const surfZone& zn : zones)
    {
        const label sz = zn.size();

        SubList<label>(zoneIds_, sz, off) = zn.index();

        off += zn.size();
    }

    if (off < size())
    {
        WarningInFunction
            << "More faces " << size() << " than zones " << off << endl;

        SubList<label>(zoneIds_, size()-off, off) = zones.last().index();
    }
    else if (size() < off)
    {
        FatalErrorInFunction
            << "More zones " << off << " than faces " << size()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polySurface::polySurface(const IOobject& io, bool doCheckIn)
:
    objectRegistry
    (
        IOobject
        (
            io.name(),
            io.db().time().constant(),
            io.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    MeshReference(faceList(), pointField()),
    zoneIds_()
{
    // Created without a point field sub-registry

    if (doCheckIn)
    {
        this->store();
    }
}


Foam::polySurface::polySurface
(
    const word& surfName,
    const objectRegistry& obr,
    bool doCheckIn
)
:
    polySurface
    (
        IOobject
        (
            surfName,
            obr.time().constant(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        doCheckIn
    )
{}


Foam::polySurface::polySurface
(
    const IOobject& io,
    const MeshedSurface<face>& surf,
    bool doCheckIn
)
:
    polySurface(io, doCheckIn)
{
    copySurface(surf);
}


Foam::polySurface::polySurface
(
    const IOobject& io,
    MeshedSurface<face>&& surf,
    bool doCheckIn
)
:
    polySurface(io, doCheckIn)
{
    transfer(surf);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polySurface::~polySurface()
{
    ///  clearOut(); // Clear addressing
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::polySurface::nFaceData() const
{
    label count = objectRegistry::size();

    // Remove PointData sub-registry from being included in the count
    if (objectRegistry::foundObject<objectRegistry>(pointDataName))
    {
        --count;
    }

    return count;
}


Foam::label Foam::polySurface::nPointData() const
{
    const objectRegistry* subreg =
        objectRegistry::cfindObject<objectRegistry>(pointDataName);

    if (subreg)
    {
        return subreg->size();
    }

    return 0;
}


const Foam::objectRegistry& Foam::polySurface::faceData() const
{
    return static_cast<const objectRegistry&>(*this);
}


const Foam::objectRegistry& Foam::polySurface::pointData() const
{
    // Force create on access
    return objectRegistry::subRegistry(pointDataName, true);
}


Foam::polySurface::FieldAssociation
Foam::polySurface::queryFieldAssociation(const word& fieldName) const
{
    unsigned where(FieldAssociation::NO_DATA);

    // Face Data
    {
        const objectRegistry* regptr = this;

        if (regptr && regptr->found(fieldName))
        {
            where |= FieldAssociation::FACE_DATA;
        }
    }

    // Point Data
    {
        const objectRegistry* regptr =
            cfindObject<objectRegistry>(pointDataName);

        if (regptr && regptr->found(fieldName))
        {
            where |= FieldAssociation::POINT_DATA;
        }
    }

    return FieldAssociation(where);
}


const Foam::regIOobject* Foam::polySurface::findFieldObject
(
    const word& fieldName,
    enum FieldAssociation association
) const
{
    const unsigned where(association);


    const regIOobject* ioptr = nullptr;

    // Face Data
    if (where & FieldAssociation::FACE_DATA)
    {
        const objectRegistry* regptr = this;

        if (regptr && (ioptr = regptr->cfindObject<regIOobject>(fieldName)))
        {
            return ioptr;
        }
    }

    // Point Data
    if (where & FieldAssociation::POINT_DATA)
    {
        const objectRegistry* regptr =
            cfindObject<objectRegistry>(pointDataName);

        if (regptr && (ioptr = regptr->cfindObject<regIOobject>(fieldName)))
        {
            return ioptr;
        }
    }

    return ioptr;
}


void Foam::polySurface::copySurface
(
    const pointField& points,
    const faceList& faces,
    bool unused
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

    this->storedPoints() = points;
    this->storedFaces() = faces;

    zoneIds_.clear();

    // if (validate)
    // {
    //     checkZones();
    // }
}


void Foam::polySurface::copySurface
(
    const meshedSurf& surf,
    bool unused
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

    this->storedPoints() = surf.points();
    this->storedFaces() = surf.faces();

    zoneIds_ = surf.zoneIds();

    // if (validate)
    // {
    //     checkZones();
    // }
}


void Foam::polySurface::copySurface
(
    const MeshedSurface<face>& surf,
    bool unused
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

    this->storedPoints() = surf.points();
    this->storedFaces() = surf.surfFaces();

    calculateZoneIds(surf.surfZones());

    // if (validate)
    // {
    //     checkZones();
    // }
}


void Foam::polySurface::transfer
(
    pointField&& points,
    faceList&& faces,
    labelList&& zoneIds
)
{
    clearOut(); // Clear addressing
    clearFields();

    this->storedPoints().transfer(points);
    this->storedFaces().transfer(faces);
    zoneIds_.transfer(zoneIds);
}


void Foam::polySurface::transfer
(
    MeshedSurface<face>& surf,
    bool unused
)
{
    clearOut(); // Clear addressing
    clearFields();

    ModifiableMeshedSurface<face> tsurf(std::move(surf));

    this->storedPoints().transfer(tsurf.storedPoints());
    this->storedFaces().transfer(tsurf.storedFaces());

    calculateZoneIds(tsurf.surfZones());

    // if (validate)
    // {
    //     checkZones();
    // }
}


// void Foam::polySurface::checkZones()
// {
//     // Extra safety, ensure we have at some zones
//     // and they cover all the faces - fix start silently
//
//     if (surfZones_.size() <= 1)
//     {
//         removeZones();
//         return;
//     }
//
//     label count = 0;
//     for (surfZone& zn : surfZones_)
//     {
//         zn.start() = count;
//         count += zn.size();
//     }
//
//     if (count < nFaces())
//     {
//         WarningInFunction
//             << "More faces " << nFaces() << " than zones " << count
//             << " ... extending final zone"
//             << endl;
//
//         surfZones_.last().size() += count - nFaces();
//     }
//     else if (size() < count)
//     {
//         FatalErrorInFunction
//             << "More zones " << count << " than faces " << nFaces()
//             << exit(FatalError);
//     }
// }


// * * * * * * * * * * * * * * * Specializations * * * * * * * * * * * * * * //

namespace Foam
{

template<>
const regIOobject* polySurface::findFieldObject<polySurfaceGeoMesh>
(
    const word& fieldName
) const
{
    // Face Data (main registry)
    return cfindObject<regIOobject>(fieldName);
}


template<>
const regIOobject* polySurface::findFieldObject<polySurfacePointGeoMesh>
(
    const word& fieldName
) const
{
    // Point Data (sub-registry)

    const objectRegistry* subreg =
        objectRegistry::cfindObject<objectRegistry>(pointDataName);

    if (subreg)
    {
        return subreg->cfindObject<regIOobject>(fieldName);
    }

    return nullptr;
}



template<>
const objectRegistry* polySurface::whichRegistry<polySurfaceGeoMesh>
(
    const word& fieldName
) const
{
    // Face Data (main registry)
    const objectRegistry* subreg = this;

    if (subreg->found(fieldName))
    {
        return subreg;
    }

    return nullptr;
}


template<>
const objectRegistry* polySurface::whichRegistry<polySurfacePointGeoMesh>
(
    const word& fieldName
) const
{
    // Point Data (sub registry)

    const objectRegistry* subreg =
        objectRegistry::cfindObject<objectRegistry>(pointDataName);

    if (subreg && subreg->found(fieldName))
    {
        return subreg;
    }

    return nullptr;
}

} // End of namespace

// ************************************************************************* //
