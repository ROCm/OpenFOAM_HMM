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

#include "tabulatedAxialAngularSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "transform.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(tabulatedAxialAngularSpring, 0);
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        tabulatedAxialAngularSpring,
        dictionary
    );
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::tabulatedAxialAngularSpring::
tabulatedAxialAngularSpring
(
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(sDoFRBMRDict),
    refQ_(),
    axis_(),
    stiffness_(),
    convertToDegrees_(),
    damping_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::tabulatedAxialAngularSpring::
~tabulatedAxialAngularSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::sixDoFRigidBodyMotionRestraints::tabulatedAxialAngularSpring::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    vector refDir = rotationTensor(vector(1, 0 ,0), axis_) & vector(0, 1, 0);

    vector oldDir = refQ_ & refDir;

    vector newDir = motion.currentOrientation(refDir);

    if (mag(oldDir & axis_) > 0.95 || mag(newDir & axis_) > 0.95)
    {
        // Directions getting close to the axis, change reference

        refDir = rotationTensor(vector(1, 0 ,0), axis_) & vector(0, 0, 1);

        vector oldDir = refQ_ & refDir;

        vector newDir = motion.currentOrientation(refDir);
    }

    // Removing any axis component from oldDir and newDir and normalising
    oldDir -= (axis_ & oldDir)*axis_;
    oldDir /= mag(oldDir);

    newDir -= (axis_ & newDir)*axis_;
    newDir /= mag(newDir);

    scalar theta = mag(acos(oldDir & newDir));

    // Temporary axis with sign information.
    vector a = (oldDir ^ newDir);

    // Remove any component that is not along axis that may creep in
    a = (a & axis_)*axis_;

    scalar magA = mag(a);

    if (magA > VSMALL)
    {
        a /= magA;
    }
    else
    {
        a = vector::zero;
    }

    scalar stiffness;

    if (convertToDegrees_)
    {
        stiffness = stiffness_(radToDeg(theta));
    }
    else
    {
        stiffness = stiffness_(theta);
    }

    // Damping of along axis angular velocity only
    restraintMoment = -stiffness*theta*a - damping_*(motion.omega() & a)*a;

    restraintForce = vector::zero;

    // Not needed to be altered as restraintForce is zero, but set to
    // centreOfMass to be sure of no spurious moment
    restraintPosition = motion.centreOfMass();
}


bool Foam::sixDoFRigidBodyMotionRestraints::tabulatedAxialAngularSpring::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    refQ_ = sDoFRBMRCoeffs_.lookupOrDefault<tensor>("referenceOrientation", I);

    if (mag(mag(refQ_) - sqrt(3.0)) > 1e-9)
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionRestraints::"
            "tabulatedAxialAngularSpring::read"
            "("
                "const dictionary& sDoFRBMRDict"
            ")"
        )
            << "referenceOrientation " << refQ_ << " is not a rotation tensor. "
            << "mag(referenceOrientation) - sqrt(3) = "
            << mag(refQ_) - sqrt(3.0) << nl
            << exit(FatalError);
    }

    axis_ = sDoFRBMRCoeffs_.lookup("axis");

    scalar magAxis(mag(axis_));

    if (magAxis > VSMALL)
    {
        axis_ /= magAxis;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionRestraints::"
            "tabulatedAxialAngularSpring::read"
            "("
                "const dictionary& sDoFRBMCDict"
            ")"
        )
            << "axis has zero length"
            << abort(FatalError);
    }

    stiffness_ = interpolationTable<scalar>(sDoFRBMRCoeffs_);

    word angleFormat = sDoFRBMRCoeffs_.lookup("angleFormat");

    if (angleFormat == "degrees" || angleFormat == "degree")
    {
        convertToDegrees_ = true;
    }
    else if (angleFormat == "radians" || angleFormat == "radian")
    {
        convertToDegrees_ = false;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::sixDoFRigidBodyMotionRestraints::"
            "tabulatedAxialAngularSpring::read"
            "("
                "const dictionary& sDoFRBMCDict"
            ")"
        )
            << "angleFormat must be degree, degrees, radian or radians"
            << abort(FatalError);
    }

    sDoFRBMRCoeffs_.lookup("damping") >> damping_;

    return true;
}

// ************************************************************************* //
