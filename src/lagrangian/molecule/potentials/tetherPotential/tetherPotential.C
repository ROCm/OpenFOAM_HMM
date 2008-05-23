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

#include "tetherPotential.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tetherPotential::tetherPotential()
{}


tetherPotential::tetherPotential
(
    const word& tetherPotentialName,
    const word& tetherPotentialType,
    const scalar springConstant
)
:
    tetherPotentialName_(tetherPotentialName),
    tetherPotentialType_(tetherPotentialType),
    springConstant_(springConstant)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetherPotential::~tetherPotential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tetherPotential::write(Ostream& os) const
{
    os << tetherPotentialName_
       << nl << tetherPotentialType_
       << nl << springConstant_
       << endl;
}


scalar tetherPotential::force(const scalar rITMag) const
{
    return -springConstant_ * rITMag;
}


scalar tetherPotential::energy(const scalar rITMag) const
{
    return 0.5*springConstant_*rITMag*rITMag;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const tetherPotential& pot)
{
    pot.write(os);
    os.check("Ostream& operator<<(Ostream& f, const tetherPotential& pot");
    return os;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
