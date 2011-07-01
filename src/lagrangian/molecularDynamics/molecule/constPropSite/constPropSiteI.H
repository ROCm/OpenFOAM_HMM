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


Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::constPropSite::siteReferencePosition() const
{
    return siteReferencePosition_;
}


inline Foam::vector& Foam::constPropSite::siteReferencePosition()
{
    return siteReferencePosition_;
}


inline const Foam::scalar& Foam::constPropSite::siteMass() const
{
    return siteMass_;
}


inline Foam::scalar& Foam::constPropSite::siteMass()
{
    return siteMass_;
}


inline const Foam::scalar& Foam::constPropSite::siteCharge() const
{
    return siteCharge_;
}


inline Foam::scalar& Foam::constPropSite::siteCharge()
{
    return siteCharge_;
}


inline const Foam::label& Foam::constPropSite::siteId() const
{
    return siteId_;
}


inline Foam::label& Foam::constPropSite::siteId()
{
    return siteId_;
}


inline const Foam::word& Foam::constPropSite::name() const
{
    return name_;
}


inline Foam::word& Foam::constPropSite::name()
{
    return name_;
}


inline const bool& Foam::constPropSite::pairPotentialSite() const
{
    return pairPotentialSite_;
}


inline bool& Foam::constPropSite::pairPotentialSite()
{
    return pairPotentialSite_;
}


inline const bool& Foam::constPropSite::electrostaticSite() const
{
    return electrostaticSite_;
}


inline bool& Foam::constPropSite::electrostaticSite()
{
    return electrostaticSite_;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::constPropSite::operator==
(
    const constPropSite& rhs
) const
{
    return
       siteReferencePosition_ == rhs.siteReferencePosition_
    && siteMass_ == rhs.siteMass_
    && siteCharge_ == rhs.siteCharge_
    && siteId_ == rhs.siteId_
    && name_ == rhs.name_
    && pairPotentialSite_ == rhs.pairPotentialSite_
    && electrostaticSite_ == rhs.electrostaticSite_;
}


bool Foam::constPropSite::operator!=
(
    const constPropSite& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
