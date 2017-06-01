/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "surfaceFeaturesExtraction.H"
#include "dictionary.H"
#include "ListOps.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{
    defineTypeName(method);
    defineRunTimeSelectionTable
    (
        method,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFeaturesExtraction::method::method()
:
    includedAngle_(0),
    geometricTestOnly_(Switch::NO)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::surfaceFeaturesExtraction::method::~method()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceFeaturesExtraction::method>
Foam::surfaceFeaturesExtraction::method::New
(
    const dictionary& dict
)
{
    const word methodName = dict.lookup("extractionMethod");

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(methodName);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown extractionMethod " << methodName << nl << nl
            << "Valid extraction methods:" << nl
            << flatOutput(dictionaryConstructorTablePtr_->sortedToc())
            << exit(FatalIOError);
    }

    return autoPtr<method>(cstrIter.object()(dict));
}


// ************************************************************************* //
