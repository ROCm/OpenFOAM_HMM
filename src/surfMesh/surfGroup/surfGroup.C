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

#include "surfGroup.H"
#include "dictionary.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::surfGroup, 0);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfGroup::surfGroup()
:
    surfPatchIdentifier(),
    size_(0),
    start_(0)
{}



Foam::surfGroup::surfGroup
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& geometricType
)
:
    surfPatchIdentifier(name, index, geometricType),
    size_(size),
    start_(start)
{}


Foam::surfGroup::surfGroup(Istream& is, const label index)
:
    surfPatchIdentifier(),
    size_(0),
    start_(0)
{
    word name(is);
    dictionary dict(is);

    operator=(surfGroup(name, dict, index));
}


Foam::surfGroup::surfGroup
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    surfPatchIdentifier(name, dict, index),
    size_(readLabel(dict.lookup("nFaces"))),
    start_(readLabel(dict.lookup("startFace")))
{}


Foam::surfGroup::surfGroup(const surfGroup& p)
:
    surfPatchIdentifier(p, p.index()),
    size_(p.size()),
    start_(p.start())
{}


Foam::surfGroup::surfGroup(const surfGroup& p, const label index)
:
    surfPatchIdentifier(p, index),
    size_(p.size()),
    start_(p.start())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfGroup::write(Ostream& os) const
{
    writeDict(os);
}


void Foam::surfGroup::writeDict(Ostream& os) const
{
    os  << indent << name() << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    surfPatchIdentifier::write(os);
    os.writeKeyword("nFaces") << size() << token::END_STATEMENT << nl;
    os.writeKeyword("startFace") << start() << token::END_STATEMENT << nl;

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::surfGroup::operator!=(const surfGroup& p) const
{
    return !(*this == p);
}


bool Foam::surfGroup::operator==(const surfGroup& p) const
{
    return
    (
        (geometricType() == p.geometricType())
     && (size() == p.size())
     && (start() == p.start())
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, surfGroup& p)
{
    p = surfGroup(is, 0);

    is.check("Istream& operator>>(Istream&, surfGroup&)");
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const surfGroup& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const surfGroup& p");
    return os;
}


// ************************************************************************* //
