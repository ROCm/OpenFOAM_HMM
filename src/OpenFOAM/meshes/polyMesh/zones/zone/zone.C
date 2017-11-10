/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "zone.H"
#include "IOstream.H"
#include "demandDrivenData.H"
#include "HashSet.H"

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
        InfoInFunction << "Calculating lookup map" << endl;
    }

    if (lookupMapPtr_)
    {
        FatalErrorInFunction
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
        InfoInFunction << "Finished calculating lookup map" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zone::zone(const word& name, const label index)
:
    labelList(),
    name_(name),
    index_(index),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const word& name,
    const labelUList& addr,
    const label index
)
:
    labelList(addr),
    name_(name),
    index_(index),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const word& name,
    labelList&& addr,
    const label index
)
:
    labelList(std::move(addr)),
    name_(name),
    index_(index),
    lookupMapPtr_(nullptr)
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
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const word& name,
    const dictionary& dict,
    const word& labelsName,
    const label index
)
:
    labelList(dict.lookup(labelsName)),
    name_(name),
    index_(index),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const zone& origZone,
    const labelUList& addr,
    const label index
)
:
    labelList(addr),
    name_(origZone.name()),
    index_(index),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const zone& origZone,
    labelList&& addr,
    const label index
)
:
    labelList(std::move(addr)),
    name_(origZone.name()),
    index_(index),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const zone& origZone,
    const Xfer<labelList>& addr,
    const label index
)
:
    labelList(addr),
    name_(origZone.name()),
    index_(index),
    lookupMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zone::~zone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::zone::localID(const label globalCellID) const
{
    return lookupMap().lookup(globalCellID, -1);
}


void Foam::zone::clearAddressing()
{
    deleteDemandDrivenData(lookupMapPtr_);
}


bool Foam::zone::checkDefinition(const label maxSize, const bool report) const
{
    const labelList& addr = *this;

    bool hasError = false;

    // To check for duplicate entries
    labelHashSet elems(size());

    for (const label idx : addr)
    {
        if (idx < 0 || idx >= maxSize)
        {
            hasError = true;

            if (report)
            {
                SeriousErrorInFunction
                    << "Zone " << name_
                    << " contains invalid index label " << idx << nl
                    << "Valid index labels are 0.."
                    << maxSize-1 << endl;
            }
            else
            {
                // w/o report - can stop checking now
                break;
            }
        }
        else if (!elems.insert(idx))
        {
            if (report)
            {
                WarningInFunction
                    << "Zone " << name_
                    << " contains duplicate index label " << idx << endl;
            }
        }
    }

    return hasError;
}


void Foam::zone::write(Ostream& os) const
{
    os  << nl << name_
        << nl << static_cast<const labelList&>(*this);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const zone& zn)
{
    zn.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
