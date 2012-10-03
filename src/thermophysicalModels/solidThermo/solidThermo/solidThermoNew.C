/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const fvMesh& mesh
)
{
    return basicThermo::New<solidThermo>(mesh);
}


Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    return basicThermo::New<solidThermo>(mesh, dict);
}

/*
Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "solidThermo::New(const fvMesh&, const dictionary&): "
            << "constructing solidThermo"
            << endl;
    }

    const word thermoType = dict.lookup("thermoType");

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(thermoType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidThermo::New(const fvMesh&, const dictionary&)"
        )   << "Unknown solidThermo type " << thermoType
            << endl << endl
            << "Valid solidThermo types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidThermo>(cstrIter()(mesh, dict));
}
*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
