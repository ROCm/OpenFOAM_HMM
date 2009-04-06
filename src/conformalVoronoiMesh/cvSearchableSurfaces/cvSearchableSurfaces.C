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

#include "cvSearchableSurfaces.H"
#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cvSearchableSurfaces::cvSearchableSurfaces
(
    const conformalVoronoiMesh& cvMesh,
    const dictionary& geometryDict
)
:
    searchableSurfaces
    (
        IOobject
        (
            "cvSearchableSurfacesDirectory",
            cvMesh.time().constant(),
            "triSurface",
            cvMesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict
    ),
    cvMesh_(cvMesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cvSearchableSurfaces::~cvSearchableSurfaces()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::cvSearchableSurfaces::inside(const point& pt) const
{

}


bool Foam::cvSearchableSurfaces::outside(const point& pt) const
{

}


bool Foam::cvSearchableSurfaces::wellInside
(
    const point& pt,
    const scalar dist2
) const
{

}


bool Foam::cvSearchableSurfaces::wellOutside
(
    const point& pt,
    const scalar dist2
) const
{

}


// ************************************************************************* //
