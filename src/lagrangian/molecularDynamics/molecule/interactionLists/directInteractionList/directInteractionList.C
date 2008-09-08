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

#include "directInteractionList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const dataType Foam::directInteractionList::staticData();


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directInteractionList::directInteractionList()
:
    baseClassName(),
    data_()
{}


Foam::directInteractionList::directInteractionList(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


Foam::directInteractionList::directInteractionList(const directInteractionList&)
:
    baseClassName(),
    data_()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::directInteractionList> Foam::directInteractionList::New()
{
    return autoPtr<directInteractionList>(new directInteractionList);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directInteractionList::~directInteractionList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::directInteractionList::operator=(const directInteractionList& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::directInteractionList::operator=(const Foam::directInteractionList&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
