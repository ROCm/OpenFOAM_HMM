/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "coordinateRotation.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateRotation> Foam::coordinateRotation::New
(
    const dictionary& dict,
    const objectRegistry& obr
)
{
    const word rotType = dict.lookup("type");

    auto cstrIter = objectRegistryConstructorTablePtr_->cfind(rotType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown coordinateRotation type " << rotType << nl << nl
            << "Valid coordinateRotation types :" <<  nl
            << objectRegistryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateRotation>(cstrIter()(dict, obr));
}


Foam::autoPtr<Foam::coordinateRotation> Foam::coordinateRotation::New
(
    const dictionary& dict
)
{
    const word rotType = dict.lookup("type");

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(rotType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown coordinateRotation type " << rotType << nl << nl
            << "Valid coordinateRotation types :" <<  nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<coordinateRotation>(cstrIter()(dict));
}


// ************************************************************************* //
