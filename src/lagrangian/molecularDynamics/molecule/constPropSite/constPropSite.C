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

Class
    constPropSite

Description

\*----------------------------------------------------------------------------*/

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
