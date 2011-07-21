/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "constPropSite.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constPropSite::constPropSite()
:
    siteReferencePosition_(vector::zero),
    siteMass_(0.0),
    siteCharge_(0.0),
    siteId_(0),
    name_(),
    pairPotentialSite_(false),
    electrostaticSite_(false)
{}


Foam::constPropSite::constPropSite
(
    const vector& siteReferencePosition,
    const scalar& siteMass,
    const scalar& siteCharge,
    const label& siteId,
    const word& name,
    const bool& pairPotentialSite,
    const bool& electrostaticSite
)
:
    siteReferencePosition_(siteReferencePosition),
    siteMass_(siteMass),
    siteCharge_(siteCharge),
    siteId_(siteId),
    name_(name),
    pairPotentialSite_(pairPotentialSite),
    electrostaticSite_(electrostaticSite)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constPropSite::~constPropSite()
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, constPropSite& cPS)
{
    is  >> cPS.siteReferencePosition_
        >> cPS.siteMass_
        >> cPS.siteCharge_
        >> cPS.siteId_
        >> cPS.name_
        >> cPS.pairPotentialSite_
        >> cPS.electrostaticSite_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::constPropSite&)"
    );

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const constPropSite& cPS)
{

    os  << token::SPACE << cPS.siteReferencePosition()
        << token::SPACE << cPS.siteMass()
        << token::SPACE << cPS.siteCharge()
        << token::SPACE << cPS.siteId()
        << token::SPACE << cPS.name()
        << token::SPACE << cPS.pairPotentialSite()
        << token::SPACE << cPS.electrostaticSite();

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::constPropSite&)"
    );

    return os;
}


// ************************************************************************* //
