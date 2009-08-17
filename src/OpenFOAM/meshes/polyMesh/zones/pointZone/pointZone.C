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

#include "pointZone.H"
#include "addToRunTimeSelectionTable.H"
#include "pointZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointZone, 0);
    defineRunTimeSelectionTable(pointZone, dictionary);
    addToRunTimeSelectionTable(pointZone, pointZone, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointZone::pointZone
(
    const word& name,
    const labelList& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(name, addr, index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const word& name,
    const Xfer<labelList>& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(name, addr, index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const pointZoneMesh& zm
)
:
    zone("point", name, dict, index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const pointZone& pz,
    const labelList& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(pz, addr, index),
    zoneMesh_(zm)
{}


Foam::pointZone::pointZone
(
    const pointZone& pz,
    const Xfer<labelList>& addr,
    const label index,
    const pointZoneMesh& zm
)
:
    zone(pz, addr, index),
    zoneMesh_(zm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointZone::~pointZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointZoneMesh& Foam::pointZone::zoneMesh() const
{
    return zoneMesh_;
}


Foam::label Foam::pointZone::whichPoint(const label globalPointID) const
{
    return zone::localID(globalPointID);
}


bool Foam::pointZone::checkDefinition(const bool report) const
{
    return zone::checkDefinition(zoneMesh_.mesh().points().size(), report);
}


void Foam::pointZone::writeDict(Ostream& os) const
{
    os  << nl << name_ << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry("pointLabels", os);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pointZone::operator=(const pointZone& pz)
{
    clearAddressing();
    labelList::operator=(pz);
}


void Foam::pointZone::operator=(const labelList& addr)
{
    clearAddressing();
    labelList::operator=(addr);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pointZone& pz)
{
    pz.write(os);
    os.check("Ostream& operator<<(Ostream& os, const pointZone& pz");
    return os;
}


// ************************************************************************* //
