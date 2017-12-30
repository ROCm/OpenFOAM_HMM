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

#include "fieldSelection.H"
#include "objectRegistry.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldSelection::fieldSelection
(
    const objectRegistry& obr
)
:
    HashSet<wordRe>(),
    obr_(obr),
    selection_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldSelection::~fieldSelection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldSelection::read(const dictionary& dict)
{
    dict.lookup("fields") >> *this;

    return true;
}


bool Foam::functionObjects::fieldSelection::containsPattern() const
{
    for (const wordRe& fieldName : *this)
    {
        if (fieldName.isPattern())
        {
            return true;
        }
    }

    return false;
}


void Foam::functionObjects::fieldSelection::clearSelection()
{
    selection_.clear();
}


bool Foam::functionObjects::fieldSelection::updateSelection()
{
    return false;
}


// ************************************************************************* //
