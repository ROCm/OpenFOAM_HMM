/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "linearSphericalAngularSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(linearSphericalAngularSpring, 0);
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        linearSphericalAngularSpring,
        dictionary
    );
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearSphericalAngularSpring::
linearSphericalAngularSpring
(
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(sDoFRBMRDict),
    refDir_(),
    stiffness_(),
    damping_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearSphericalAngularSpring::
~linearSphericalAngularSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::sixDoFRigidBodyMotionRestraints::linearSphericalAngularSpring::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    vector curDir = motion.currentOrientation(refDir_);

    vector a = refDir_ ^ curDir;

    scalar magA = mag(a);

    scalar theta = 0;

    if (magA > VSMALL)
    {
        a /= magA;

        theta = acos(curDir & refDir_);
    }
    else
    {
        a = vector::zero;
    }

    restraintMoment = -stiffness_*theta*a - damping_*motion.omega();

    restraintForce = vector::zero;

    // Not needed to be altered as restraintForce is zero, but set to
    // centre of mass to be sure of no spurious moment
    restraintPosition = motion.centreOfMass();
}


bool Foam::sixDoFRigidBodyMotionRestraints::linearSphericalAngularSpring ::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.lookup("referenceDirection") >> refDir_;

    sDoFRBMRCoeffs_.lookup("stiffness") >> stiffness_;

    sDoFRBMRCoeffs_.lookup("damping") >> damping_;

    return true;
}

// ************************************************************************* //
