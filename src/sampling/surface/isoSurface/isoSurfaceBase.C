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

#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::isoSurfaceBase::algorithmType
>
Foam::isoSurfaceBase::algorithmNames
({
    { algorithmType::ALGO_CELL, "cell" },
    { algorithmType::ALGO_TOPO, "topo" },
    { algorithmType::ALGO_POINT, "point" },
});


const Foam::Enum
<
    Foam::isoSurfaceBase::filterType
>
Foam::isoSurfaceBase::filterNames
({
    { filterType::NONE, "none" },
    { filterType::CELL, "cell" },
    { filterType::DIAGCELL, "diagcell" },
    { filterType::PARTIAL, "partial" },
    { filterType::FULL, "full" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::isoSurfaceBase::filterType
Foam::isoSurfaceBase::getFilterType
(
    const dictionary& dict,
    const isoSurfaceBase::filterType deflt
)
{
    word filterName;

    if (!dict.readIfPresent("regularise", filterName, keyType::LITERAL))
    {
        return deflt;
    }

    // Try as bool
    Switch sw(filterName, true);

    if (sw.valid())
    {
        return
        (
            sw
          ? deflt
          : isoSurfaceBase::filterType::NONE
        );
    }

    // As enum
    if (!isoSurfaceBase::filterNames.found(filterName))
    {
        FatalIOErrorInFunction(dict)
            << filterName << " is not in enumeration: "
            << isoSurfaceBase::filterNames << nl
            << exit(FatalIOError);
    }

    return isoSurfaceBase::filterNames[filterName];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurfaceBase::isoSurfaceBase
(
    const scalar iso,
    const boundBox& bounds
)
:
    meshedSurface(),
    iso_(iso),
    bounds_(bounds)
{}


// ************************************************************************* //
