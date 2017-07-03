/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "searchableSurfaceModifier.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableSurfaceModifier, 0);
    defineRunTimeSelectionTable(searchableSurfaceModifier, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceModifier::searchableSurfaceModifier
(
    const searchableSurfaces& geometry,
    const dictionary& dict
)
:
    geometry_(geometry),
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceModifier::~searchableSurfaceModifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::searchableSurfaceModifier>
Foam::searchableSurfaceModifier::New
(
    const word& type,
    const searchableSurfaces& geometry,
    const dictionary& dict
)
{
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(type);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown searchableSurfaceModifier type "
            << type << nl << nl
            << "Valid searchableSurfaceModifier types :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<searchableSurfaceModifier>(cstrIter()(geometry, dict));
}


// ************************************************************************* //
