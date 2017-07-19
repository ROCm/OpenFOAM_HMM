/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "patchDistMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchDistMethod, 0);
    defineRunTimeSelectionTable(patchDistMethod, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethod::patchDistMethod
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    mesh_(mesh),
    patchIDs_(patchIDs)
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::patchDistMethod> Foam::patchDistMethod::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
{
    const word methodType(dict.lookup("method"));

    Info<< "Selecting patchDistMethod " << methodType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(methodType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown patchDistMethod type "
            << methodType << endl << endl
            << "Valid patchDistMethod types : " << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, mesh, patchIDs);
}



Foam::autoPtr<Foam::patchDistMethod> Foam::patchDistMethod::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const word& defaultPatchDistMethod
)
{
    word methodType = defaultPatchDistMethod;
    dict.readIfPresent("method", methodType);

    Info<< "Selecting patchDistMethod " << methodType << endl;
    auto cstrIter = dictionaryConstructorTablePtr_->cfind(methodType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown patchDistMethodType type "
            << methodType << endl << endl
            << "Valid patchDistMethod types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, mesh, patchIDs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchDistMethod::~patchDistMethod()
{}


// ************************************************************************* //
