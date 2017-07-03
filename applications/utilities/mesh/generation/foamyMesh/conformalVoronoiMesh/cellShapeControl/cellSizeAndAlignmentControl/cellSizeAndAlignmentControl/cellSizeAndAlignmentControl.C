/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "cellSizeAndAlignmentControl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellSizeAndAlignmentControl, 0);
defineRunTimeSelectionTable(cellSizeAndAlignmentControl, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSizeAndAlignmentControl::cellSizeAndAlignmentControl
(
    const Time& runTime,
    const word& name,
    const dictionary& controlFunctionDict,
    const conformationSurfaces& geometryToConformTo,
    const scalar& defaultCellSize
)
:
    runTime_(runTime),
    defaultCellSize_(defaultCellSize),
    forceInitialPointInsertion_
    (
        controlFunctionDict.lookupOrDefault<Switch>
        (
            "forceInitialPointInsertion",
            Switch::OFF
        )
    ),
    name_(name)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::cellSizeAndAlignmentControl>
Foam::cellSizeAndAlignmentControl::New
(
    const Time& runTime,
    const word& name,
    const dictionary& controlFunctionDict,
    const conformationSurfaces& geometryToConformTo,
    const scalar& defaultCellSize
)
{
    const word controlType(controlFunctionDict.lookup("type"));

    Info<< indent << "Selecting cellSizeAndAlignmentControl "
        << controlType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(controlType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown cellSizeAndAlignmentControl type "
            << controlType << nl << nl
            << "Valid cellSizeAndAlignmentControl types :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<cellSizeAndAlignmentControl>
    (
        cstrIter()
        (
            runTime,
            name,
            controlFunctionDict,
            geometryToConformTo,
            defaultCellSize
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSizeAndAlignmentControl::~cellSizeAndAlignmentControl()
{}


// ************************************************************************* //
