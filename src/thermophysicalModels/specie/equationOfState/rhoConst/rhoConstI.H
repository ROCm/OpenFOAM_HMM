/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "rhoConst.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::rhoConst::rhoConst
(
    const specie& sp,
    const scalar rho
)
:
    specie(sp),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::rhoConst::rhoConst
(
    const word& name,
    const rhoConst& ico
)
:
    specie(name, ico),
    rho_(ico.rho_)
{}


inline Foam::autoPtr<Foam::rhoConst>
Foam::rhoConst::clone() const
{
    return autoPtr<rhoConst>(new rhoConst(*this));
}


inline Foam::autoPtr<Foam::rhoConst>
Foam::rhoConst::New(Istream& is)
{
    return autoPtr<rhoConst>(new rhoConst(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::rhoConst::rho(scalar p, scalar T) const
{
    return rho_;
}


inline Foam::scalar Foam::rhoConst::psi(scalar, scalar T) const
{
    return 0.0;
}


inline Foam::scalar Foam::rhoConst::Z(scalar, scalar) const
{
    return 0.0;
}


inline Foam::scalar Foam::rhoConst::cpMcv(scalar, scalar) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::rhoConst::operator+=(const rhoConst& ico)
{
    scalar molr1 = this->nMoles();

    specie::operator+=(ico);

    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ + molr2*ico.rho_;
}


inline void Foam::rhoConst::operator-=(const rhoConst& ico)
{
    scalar molr1 = this->nMoles();

    specie::operator-=(ico);

    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho_ = molr1*rho_ - molr2*ico.rho_;
}


inline void Foam::rhoConst::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::rhoConst Foam::operator+
(
    const rhoConst& ico1,
    const rhoConst& ico2
)
{
    scalar nMoles = ico1.nMoles() + ico2.nMoles();
    scalar molr1 = ico1.nMoles()/nMoles;
    scalar molr2 = ico2.nMoles()/nMoles;

    return rhoConst
    (
        static_cast<const specie&>(ico1)
      + static_cast<const specie&>(ico2),
        molr1*ico1.rho_ + molr2*ico2.rho_
    );
}


inline Foam::rhoConst Foam::operator-
(
    const rhoConst& ico1,
    const rhoConst& ico2
)
{
    scalar nMoles = ico1.nMoles() + ico2.nMoles();
    scalar molr1 = ico1.nMoles()/nMoles;
    scalar molr2 = ico2.nMoles()/nMoles;

    return rhoConst
    (
        static_cast<const specie&>(ico1)
      - static_cast<const specie&>(ico2),
        molr1*ico1.rho_ - molr2*ico2.rho_
    );
}


inline Foam::rhoConst Foam::operator*
(
    const scalar s,
    const rhoConst& ico
)
{
    return rhoConst(s*static_cast<const specie&>(ico), ico.rho_);
}


inline Foam::rhoConst Foam::operator==
(
    const rhoConst& ico1,
    const rhoConst& ico2
)
{
    return ico2 - ico1;
}


// ************************************************************************* //
