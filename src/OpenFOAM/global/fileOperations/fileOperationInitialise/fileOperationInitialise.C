/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "fileOperationInitialise.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(fileOperationInitialise, 0);
    defineRunTimeSelectionTable(fileOperationInitialise, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::fileOperationInitialise::fileOperationInitialise
(
    int& argc,
    char**& argv
)
{}


Foam::autoPtr<Foam::fileOperations::fileOperationInitialise>
Foam::fileOperations::fileOperationInitialise::New
(
    const word& type,
    int& argc,
    char**& argv
)
{
    if (debug)
    {
        InfoInFunction << "Constructing fileOperationInitialise" << endl;
    }

    auto cstrIter = wordConstructorTablePtr_->cfind(type);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown fileOperationInitialise type "
            << type << nl << nl
            << "Valid fileOperationInitialise types are" << endl
            << wordConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<fileOperationInitialise>(cstrIter()(argc, argv));
}


// ************************************************************************* //
