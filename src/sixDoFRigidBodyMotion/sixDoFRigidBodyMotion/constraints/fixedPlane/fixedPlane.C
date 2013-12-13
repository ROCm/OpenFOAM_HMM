/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fixedPlane.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionConstraints
{
    defineTypeNameAndDebug(fixedPlane, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionConstraint,
        fixedPlane,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::fixedPlane
(
    const word& name,
    const dictionary& sDoFRBMCDict
)
:
    sixDoFRigidBodyMotionConstraint(name, sDoFRBMCDict),
    normal_(vector::zero)
{
    read(sDoFRBMCDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::~fixedPlane()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::constrainTranslation
(
    pointConstraint& pc
) const
{
    pc.applyConstraint(normal_);
}


void Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::constrainRotation
(
    pointConstraint& pc
) const
{}


bool Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::read
(
    const dictionary& sDoFRBMCDict
)
{
    sixDoFRigidBodyMotionConstraint::read(sDoFRBMCDict);

    normal_ = sDoFRBMCCoeffs_.lookup("normal");

    return true;
}


void Foam::sixDoFRigidBodyMotionConstraints::fixedPlane::write
(
    Ostream& os
) const
{
    os.writeKeyword("normal")
        << normal_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
