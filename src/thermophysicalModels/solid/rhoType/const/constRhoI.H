/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "constRho.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::constRho::constRho(const word& name, const constRho& pg)
:
    specie(name, pg),
    rho_(pg.rho_)
{}


inline Foam::constRho::constRho(const specie& t, const scalar rho)
:
    specie(t),
    rho_(rho)
{}


inline Foam::autoPtr<Foam::constRho> Foam::constRho::clone() const
{
    return autoPtr<constRho>(new constRho(*this));
}


inline Foam::autoPtr<Foam::constRho> Foam::constRho::New(const dictionary& is)
{
    return autoPtr<constRho>(new constRho(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::constRho::rho(scalar, scalar) const
{
    return rho_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::constRho& Foam::constRho::operator=(const constRho& ct)
{
    specie::operator=(ct);
    rho_ = ct.rho_;

    return *this;
}


inline void Foam::constRho::operator+=(const constRho& ct)
{
    scalar molr1 = this->nMoles();

    specie::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    rho_ = molr1*rho_ + molr2*ct.rho_;
}


inline void Foam::constRho::operator-=(const constRho& ct)
{
    scalar molr1 = this->nMoles();

    specie::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    rho_ = molr1*rho_ - molr2*ct.rho_;

}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


inline Foam::constRho Foam::operator*
(
    const scalar s,
    const constRho& ct
)
{
    return constRho
    (
        s*static_cast<const specie&>(ct),
        ct.rho_
    );
}

// ************************************************************************* //
