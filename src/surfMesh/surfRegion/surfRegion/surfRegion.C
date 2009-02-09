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

Description

\*---------------------------------------------------------------------------*/

#include "surfRegion.H"
#include "dictionary.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::surfRegion, 0);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfRegion::surfRegion()
:
    surfRegionIdentifier(),
    size_(0),
    start_(0)
{}



Foam::surfRegion::surfRegion
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& geometricType
)
:
    surfRegionIdentifier(name, index, geometricType),
    size_(size),
    start_(start)
{}


Foam::surfRegion::surfRegion(Istream& is, const label index)
:
    surfRegionIdentifier(),
    size_(0),
    start_(0)
{
    word name(is);
    dictionary dict(is);

    operator=(surfRegion(name, dict, index));
}


Foam::surfRegion::surfRegion
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    surfRegionIdentifier(name, dict, index),
    size_(readLabel(dict.lookup("nFaces"))),
    start_(readLabel(dict.lookup("startFace")))
{}


Foam::surfRegion::surfRegion(const surfRegion& reg)
:
    surfRegionIdentifier(reg, reg.index()),
    size_(reg.size()),
    start_(reg.start())
{}


Foam::surfRegion::surfRegion(const surfRegion& reg, const label index)
:
    surfRegionIdentifier(reg, index),
    size_(reg.size()),
    start_(reg.start())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfRegion::write(Ostream& os) const
{
    writeDict(os);
}


void Foam::surfRegion::writeDict(Ostream& os) const
{
    os  << indent << name() << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    surfRegionIdentifier::write(os);
    os.writeKeyword("nFaces") << size() << token::END_STATEMENT << nl;
    os.writeKeyword("startFace") << start() << token::END_STATEMENT << nl;

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::surfRegion::operator!=(const surfRegion& reg) const
{
    return !(*this == reg);
}


bool Foam::surfRegion::operator==(const surfRegion& reg) const
{
    return
    (
        (geometricType() == reg.geometricType())
     && (size() == reg.size())
     && (start() == reg.start())
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, surfRegion& reg)
{
    reg = surfRegion(is, 0);

    is.check("Istream& operator>>(Istream&, surfRegion&)");
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const surfRegion& reg)
{
    reg.write(os);
    os.check("Ostream& operator<<(Ostream&, const surfRegion&");
    return os;
}


// ************************************************************************* //
