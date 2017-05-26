/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "surfZoneIdentifier.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::surfZoneIdentifier::emptyType = "empty";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZoneIdentifier::surfZoneIdentifier()
:
    name_(),
    index_(0),
    geometricType_()
{}


Foam::surfZoneIdentifier::surfZoneIdentifier(const label index)
:
    name_(),
    index_(index),
    geometricType_()
{}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    const word& name,
    const label index,
    const word& geometricType
)
:
    name_(name),
    index_(index),
    geometricType_(geometricType)
{}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    name_(name),
    index_(index),
    geometricType_()
{
    dict.readIfPresent("geometricType", geometricType_);
}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
    const surfZoneIdentifier& p,
    const label index
)
:
    name_(p.name()),
    index_(index),
    geometricType_(p.geometricType())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfZoneIdentifier::~surfZoneIdentifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfZoneIdentifier::write(Ostream& os) const
{
    if (geometricType_.size())
    {
        os.writeEntry("geometricType", geometricType_);
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const surfZoneIdentifier& a, const surfZoneIdentifier& b)
{
    return
    (
        (a.index() == b.index())
     && (a.name()  == b.name())
     && (a.geometricType() == b.geometricType())
    );
}


bool Foam::operator!=(const surfZoneIdentifier& a, const surfZoneIdentifier& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, surfZoneIdentifier& obj)
{
    is >> obj.name_ >> obj.geometricType_;
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const surfZoneIdentifier& obj)
{
    // Newlines to separate, since that is what triSurface currently expects
    os  << nl << obj.name_ << nl << obj.geometricType_;
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
