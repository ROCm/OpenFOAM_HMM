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

#include "solverFieldSelection.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::solverFieldSelection::solverFieldSelection
(
    const objectRegistry& obr
)
:
    volFieldSelection(obr)
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "Registry must be of type " << fvMesh::typeName
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::solverFieldSelection::~solverFieldSelection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::solverFieldSelection::updateSelection()
{
    wordHashSet oldSet;
    oldSet.swap(selection_);

    wordHashSet volFields;
    addRegisteredGeoFields<fvPatchField, volMesh>(volFields);

    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);

    const Foam::dictionary& solverDict = mesh.solverPerformanceDict();

    for (const word& fieldName : volFields)
    {
        if (solverDict.found(fieldName))
        {
            selection_.insert(fieldName);
        }
    }

    return selection_ != oldSet;
}


// ************************************************************************* //
