/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "pairPotential.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pairPotential::pairPotential()
{}


// Construct from components
pairPotential::pairPotential
(
    const word& pairPotentialName,
    const word& pairPotentialType,
    const scalar sigma,
    const scalar epsilon,
    const scalar rCut
)
:
    pairPotentialName_(pairPotentialName),
    pairPotentialType_(pairPotentialType),
    sigma_(sigma),
    epsilon_(epsilon),
    rCut_(rCut),
    rCutSqr_(rCut*rCut)
{
    // rCutShiftedForcePotentialConstTerm
    rCutSFPotConst_ = (pow((rCut_/sigma_),-12) - pow((rCut_/sigma_),-6));

    // rCutShiftedForceConstTerm
    rCutSFConst_ =
        (pow((rCut_/sigma_),-14) - 0.5*pow((rCut_/sigma_),-8))*rCut_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pairPotential::~pairPotential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pairPotential::write(Ostream& os) const
{
    os << pairPotentialName_
       << nl << pairPotentialType_
       << nl << sigma_
       << nl << epsilon_
       << nl << rCut_
       << endl;
}


scalar pairPotential::force(const scalar rIJMag) const
{
    // (rIJ/sigma)^-2
    scalar rIJoSMagSqInv = (sigma_/rIJMag)*(sigma_/rIJMag);

    // (rIJ/sigma)^-6
    scalar rIJoSMagSqInvCu = rIJoSMagSqInv*rIJoSMagSqInv*rIJoSMagSqInv;

    scalar f = 48.0*epsilon_/(sigma_*sigma_)
               *(rIJoSMagSqInvCu*(rIJoSMagSqInvCu - 0.5)*rIJoSMagSqInv
             - rCutSFConst_/rIJMag);

    return f;
}


scalar pairPotential::energy(const scalar rIJMag) const
{
    // (rIJ/sigma)^-2
    scalar rIJoSMagSqInv = (sigma_/rIJMag)*(sigma_/rIJMag);

    // (rIJ/sigma)^-6
    scalar rIJoSMagSqInvCu = rIJoSMagSqInv*rIJoSMagSqInv*rIJoSMagSqInv;

    scalar e = 4*epsilon_
               *(
                    rIJoSMagSqInvCu*(rIJoSMagSqInvCu - 1.0)
                  - rCutSFPotConst_
                  + 12*(rIJMag - rCut_)/(sigma_*sigma_)*rCutSFConst_
                );

    return e;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const pairPotential& pot)
{
    pot.write(os);
    os.check("Ostream& operator<<(Ostream& f, const pairPotential& pot");
    return os;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
