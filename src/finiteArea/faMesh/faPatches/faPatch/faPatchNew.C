/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

#include "faPatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::faPatch> Foam::faPatch::New
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
{
    DebugInFunction
        << "constructing faPatch" << endl;

    word patchType(dict.lookup("type"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(patchType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown faPatch type " << patchType << nl << nl
            << "Valid faPatch types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<faPatch>(cstrIter()(name, dict, index, bm));
}


// ************************************************************************* //
