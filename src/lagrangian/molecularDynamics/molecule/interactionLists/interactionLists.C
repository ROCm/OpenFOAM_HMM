/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "interactionLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const dataType Foam::interactionLists::staticData();


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interactionLists::interactionLists()
:
    baseClassName(),
    data_()
{}


Foam::interactionLists::interactionLists(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


Foam::interactionLists::interactionLists(const interactionLists&)
:
    baseClassName(),
    data_()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interactionLists> Foam::interactionLists::New()
{
    return autoPtr<interactionLists>(new interactionLists);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interactionLists::~interactionLists()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::interactionLists::operator=(const interactionLists& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::interactionLists::operator=(const Foam::interactionLists&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
