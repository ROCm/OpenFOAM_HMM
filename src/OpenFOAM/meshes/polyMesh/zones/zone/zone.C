/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "zone.H"
#include "IOstream.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zone, 0);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::Map<Foam::label>& Foam::zone::lookupMap() const
{
    if (!lookupMapPtr_)
    {
        calcLookupMap();
    }

    return *lookupMapPtr_;
}


void Foam::zone::calcLookupMap() const
{
    if (debug)
    {
        Info<< "void zone::calcLookupMap() const: "
            << "Calculating lookup map"
            << endl;
    }

    if (lookupMapPtr_)
    {
        FatalErrorIn("void zone::calcLookupMap() const")
            << "Lookup map already calculated" << nl
            << abort(FatalError);
    }

    const labelList& addr = *this;

    lookupMapPtr_ = new Map<label>(2*addr.size());
    Map<label>& lm = *lookupMapPtr_;

    forAll(addr, i)
    {
        lm.insert(addr[i], i);
    }

    if (debug)
    {
        Info<< "void zone::calcLookupMap() const: "
            << "Finished calculating lookup map"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zone::zone
(
    const word& name,
    const labelList& addr,
    const label index
)
:
    labelList(addr),
    name_(name),
    index_(index),
    lookupMapPtr_(NULL)
{}


Foam::zone::zone
(
    const word& name,
    const Xfer<labelList>& addr,
    const label index
)
:
    labelList(addr),
    name_(name),
    index_(index),
    lookupMapPtr_(NULL)
{}


Foam::zone::zone
(
    const word& zoneType,
    const word& name,
    const dictionary& dict,
    const label index
)
:
    labelList(dict.lookup(zoneType + "Labels")),
    name_(name),
    index_(index),
    lookupMapPtr_(NULL)
{}


Foam::zone::zone
(
    const zone& z,
    const labelList& addr,
    const label index
)
:
    labelList(addr),
    name_(z.name()),
    index_(index),
    lookupMapPtr_(NULL)
{}


Foam::zone::zone
(
    const zone& z,
    const Xfer<labelList>& addr,
    const label index
)
:
    labelList(addr),
    name_(z.name()),
    index_(index),
    lookupMapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zone::~zone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::zone::localID(const label globalCellID) const
{
    const Map<label>& lm = lookupMap();

    Map<label>::const_iterator lmIter = lm.find(globalCellID);

    if (lmIter == lm.end())
    {
        return -1;
    }
    else
    {
        return lmIter();
    }
}


void Foam::zone::clearAddressing()
{
    deleteDemandDrivenData(lookupMapPtr_);
}


void Foam::zone::write(Ostream& os) const
{
    os  << nl << name_
        << nl << static_cast<const labelList&>(*this);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const zone& z)
{
    z.write(os);
    os.check("Ostream& operator<<(Ostream& f, const zone& z");
    return os;
}


// ************************************************************************* //
