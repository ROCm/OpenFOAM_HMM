/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "geometricSurfacePatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(geometricSurfacePatch, 0);
}

const Foam::word Foam::geometricSurfacePatch::emptyType = "empty";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricSurfacePatch::geometricSurfacePatch()
:
    geometricType_(emptyType),
    name_("patch"),
    index_(0)
{}


Foam::geometricSurfacePatch::geometricSurfacePatch(const label index)
:
    geometricType_(emptyType),
    name_("patch"),
    index_(index)
{}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    const word& name,
    const label index,
    const word& geometricType
)
:
    geometricType_(geometricType),
    name_(name),
    index_(index)
{
    if (geometricType_.empty())
    {
        geometricType_ = emptyType;
    }
}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    const word& geometricType,
    const word& name,
    const label index
)
:
    geometricType_(geometricType),
    name_(name),
    index_(index)
{
    if (geometricType_.empty())
    {
        geometricType_ = emptyType;
    }
}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    const surfZoneIdentifier& ident
)
:
    geometricType_(ident.geometricType()),
    name_(ident.name()),
    index_(ident.index())
{
    if (geometricType_.empty())
    {
        geometricType_ = emptyType;
    }
}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    Istream& is,
    const label index
)
:
    geometricType_(is),
    name_(is),
    index_(index)
{
    if (geometricType_.empty())
    {
        geometricType_ = emptyType;
    }
}


Foam::geometricSurfacePatch::geometricSurfacePatch
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    geometricType_(emptyType),
    name_(name),
    index_(index)
{
    dict.readIfPresent("geometricType", geometricType_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::geometricSurfacePatch::write(Ostream& os) const
{
    os  << nl << name_
        << nl << geometricType_;
}


void Foam::geometricSurfacePatch::writeDict(Ostream& os) const
{
    os.writeEntry("geometricType", geometricType_);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

bool Foam::operator==
(
    const geometricSurfacePatch& a,
    const geometricSurfacePatch& b
)
{
    return
    (
        (a.geometricType() == b.geometricType())
     && (a.name() == b.name())
    );
}


bool Foam::operator!=
(
    const geometricSurfacePatch& a,
    const geometricSurfacePatch& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, geometricSurfacePatch& p)
{
    is >> p.name_ >> p.geometricType_;
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const geometricSurfacePatch& p)
{
    p.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
