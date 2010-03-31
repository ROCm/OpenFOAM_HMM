/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
     \\/     M anipulation  |
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

\*----------------------------------------------------------------------------*/

#include "sendingReferralList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sendingReferralList::sendingReferralList()
:
    labelList(),
    destinationProc_(-1)
{}


Foam::sendingReferralList::sendingReferralList
(
    const label destinationProc,
    const labelList& cellsToSend
)
:
    labelList(cellsToSend),
    destinationProc_(destinationProc)
{}


Foam::sendingReferralList::sendingReferralList
(
    const sendingReferralList& rL
)
:
    labelList(rL),
    destinationProc_(rL.destinationProc())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sendingReferralList::~sendingReferralList()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::sendingReferralList::operator=(const sendingReferralList& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::sendingReferralList::" \
            "operator=(const Foam::sendingReferralList&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    labelList::operator=(rhs);

    destinationProc_ = rhs.destinationProc();
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    sendingReferralList& sRL
)
{
    is  >> sRL.destinationProc_ >> static_cast<labelList&>(sRL);

    is.check("Istream& operator<<(Istream& f, const sendingReferralList& sRL");

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sendingReferralList& rL
)
{
    os  << rL.destinationProc() << token::SPACE
        << static_cast< const labelList& >(rL);

    os.check("Ostream& operator<<(Ostream& f, const sendingReferralList& rL");

    return os;
}


// ************************************************************************* //
