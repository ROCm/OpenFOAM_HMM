/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "reactionRateFlameArea.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::reactionRateFlameArea> Foam::reactionRateFlameArea::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
{
    const word modelType
    (
        dict.lookup("reactionRateFlameArea")
    );

    Info<< "Selecting reaction rate flame area correlation "
        << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown reactionRateFlameArea type "
            << modelType << nl << nl
            << "Valid reaction rate flame area types :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word className = modelType.substr(0, modelType.find('<'));

    return autoPtr<reactionRateFlameArea>
        (cstrIter()(className, dict, mesh, combModel));
}


// ************************************************************************* //
