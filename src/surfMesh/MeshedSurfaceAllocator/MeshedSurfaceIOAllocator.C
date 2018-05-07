/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "MeshedSurfaceIOAllocator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Detail::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const IOobject& ioFaces,
    const IOobject& ioZones
)
:
    points_(ioPoints),
    faces_(ioFaces),
    zones_(ioZones)
{}


Foam::Detail::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints, const pointField& points,
    const IOobject& ioFaces,  const faceList& faces,
    const IOobject& ioZones,  const surfZoneList& zones
)
:
    points_(ioPoints, points),
    faces_(ioFaces, faces),
    zones_(ioZones, zones)
{}


Foam::Detail::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints, pointField&& points,
    const IOobject& ioFaces,  faceList&& faces,
    const IOobject& ioZones,  surfZoneList&& zones
)
:
    points_(ioPoints, std::move(points)),
    faces_(ioFaces, std::move(faces)),
    zones_(ioZones, std::move(zones))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Detail::MeshedSurfaceIOAllocator::~MeshedSurfaceIOAllocator()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Detail::MeshedSurfaceIOAllocator::setInstance
(
    const fileName& inst
)
{
    points_.instance() = inst;
    faces_.instance()  = inst;
    zones_.instance()  = inst;
}


void Foam::Detail::MeshedSurfaceIOAllocator::setWriteOption
(
    IOobject::writeOption wOpt
)
{
    points_.writeOpt() = wOpt;
    faces_.writeOpt()  = wOpt;
    zones_.writeOpt()  = wOpt;
}


void Foam::Detail::MeshedSurfaceIOAllocator::clear()
{
    points_.clear();
    faces_.clear();
    zones_.clear();
}


bool Foam::Detail::MeshedSurfaceIOAllocator::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    return
    (
        points_.writeObject(fmt, ver, cmp, valid)
     && faces_.writeObject(fmt, ver, cmp, valid)
     && zones_.writeObject(fmt, ver, cmp, valid)
    );
}


// ************************************************************************* //
