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
    A subset of mesh cells.

\*---------------------------------------------------------------------------*/

#include "cellZone.H"
#include "addToRunTimeSelectionTable.H"
#include "cellZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "IOstream.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellZone, 0);

    defineRunTimeSelectionTable(cellZone, dictionary);
    addToRunTimeSelectionTable(cellZone, cellZone, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellZone::cellZone
(
    const word& name,
    const labelList& addr,
    const label index,
    const cellZoneMesh& zm
)
:
    zone(name, addr, index),
    zoneMesh_(zm)
{}


Foam::cellZone::cellZone
(
    const word& name,
    const Xfer<labelList>& addr,
    const label index,
    const cellZoneMesh& zm
)
:
    zone(name, addr, index),
    zoneMesh_(zm)
{}


Foam::cellZone::cellZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const cellZoneMesh& zm
)
:
    zone("cell", name, dict, index),
    zoneMesh_(zm)
{}


Foam::cellZone::cellZone
(
    const cellZone& cz,
    const labelList& addr,
    const label index,
    const cellZoneMesh& zm
)
:
    zone(cz, addr, index),
    zoneMesh_(zm)
{}

Foam::cellZone::cellZone
(
    const cellZone& cz,
    const Xfer<labelList>& addr,
    const label index,
    const cellZoneMesh& zm
)
:
    zone(cz, addr, index),
    zoneMesh_(zm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellZone::~cellZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellZone::whichCell(const label globalCellID) const
{
    return zone::localID(globalCellID);
}


const Foam::cellZoneMesh& Foam::cellZone::zoneMesh() const
{
    return zoneMesh_;
}


bool Foam::cellZone::checkDefinition(const bool report) const
{
    return zone::checkDefinition(zoneMesh_.mesh().nCells(), report);
}


void Foam::cellZone::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry("cellLabels", os);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cellZone::operator=(const cellZone& cz)
{
    clearAddressing();
    labelList::operator=(cz);
}


void Foam::cellZone::operator=(const labelList& addr)
{
    clearAddressing();
    labelList::operator=(addr);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const cellZone& cz)
{
    cz.write(os);
    os.check("Ostream& operator<<(Ostream& os, const cellZone& cz");
    return os;
}


// ************************************************************************* //
