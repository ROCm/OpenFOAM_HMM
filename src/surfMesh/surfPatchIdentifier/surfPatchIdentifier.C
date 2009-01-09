/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "surfPatchIdentifier.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfPatchIdentifier::surfPatchIdentifier()
:
    name_(word::null),
    boundaryIndex_(0),
    geometricType_(word::null)
{}


Foam::surfPatchIdentifier::surfPatchIdentifier
(
    const word& name,
    const label index,
    const word& geometricType
)
:
    name_(name),
    boundaryIndex_(index),
    geometricType_(geometricType)
{}


Foam::surfPatchIdentifier::surfPatchIdentifier
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    name_(name),
    boundaryIndex_(index)
{
    dict.readIfPresent("geometricType", geometricType_);
}


Foam::surfPatchIdentifier::surfPatchIdentifier
(
    const surfPatchIdentifier& p,
    const label index
)
:
    name_(p.name()),
    boundaryIndex_(index),
    geometricType_(p.geometricType())
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfPatchIdentifier::~surfPatchIdentifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::surfPatchIdentifier::write(Ostream& os) const
{
    if (geometricType_.size())
    {
        os.writeKeyword("geometricType") << geometricType_
            << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// bool Foam::surfPatchIdentifier::operator!=
// (
//     const surfPatchIdentifier& p
// ) const
// {
//     return !(*this == p);
// }
//
//
// bool Foam::surfPatchIdentifier::operator==
// (
//     const surfPatchIdentifier& p
// ) const
// {
//     return geometricType() == p.geometricType() && name() == p.name();
// }


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Foam::Istream& Foam::operator>>(Istream& is, surfPatchIdentifier& p)
// {
//     is >> p.name_ >> p.geometricType_;
//
//     return is;
// }


Foam::Ostream& Foam::operator<<(Ostream& os, const surfPatchIdentifier& p)
{
    p.write(os);
    os.check
    (
        "Ostream& operator<<(Ostream&, const surfPatchIdentifier&)"
    );
    return os;
}

// ************************************************************************* //
