/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::applyRestraints()
{
    forAll(restraints_, rI)
    {
        restraints_[rI].restraintForce();
    }
}


void Foam::sixDoFRigidBodyMotion::applyConstraints()
{

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion()
:
    motionState_(),
    restraints_(),
    refCentreOfMass_(vector::zero),
    momentOfInertia_(diagTensor::one*VSMALL),
    mass_(VSMALL)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const point& centreOfMass,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau,
    scalar mass,
    const point& refCentreOfMass,
    const diagTensor& momentOfInertia
)
:
    motionState_
    (
        centreOfMass,
        Q,
        v,
        a,
        pi,
        tau
    ),
    restraints_(),
    refCentreOfMass_(refCentreOfMass),
    momentOfInertia_(momentOfInertia),
    mass_(mass)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion(const dictionary& dict)
:
    motionState_(dict),
    restraints_(),
    refCentreOfMass_(dict.lookupOrDefault("refCentreOfMass", centreOfMass())),
    momentOfInertia_(dict.lookup("momentOfInertia")),
    mass_(readScalar(dict.lookup("mass")))
{
    addRestraints(dict);

    addConstraints(dict);
}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const sixDoFRigidBodyMotion& sDoFRBM
)
:
    motionState_(sDoFRBM.motionState()),
    restraints_(sDoFRBM.restraints()),
    refCentreOfMass_(sDoFRBM.refCentreOfMass()),
    momentOfInertia_(sDoFRBM.momentOfInertia()),
    mass_(sDoFRBM.mass())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::~sixDoFRigidBodyMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::addRestraints
(
    const dictionary& dict
)
{
    if (dict.found("restraints"))
    {
        const dictionary& restraintDict = dict.subDict("restraints");

        label i = 0;

        restraints_.setSize(restraintDict.size());

        forAllConstIter(IDLList<entry>, restraintDict, iter)
        {
            if (iter().isDict())
            {
                Info<< "Adding restraint: " << iter().keyword() << endl;

                restraints_.set
                (
                    i,
                    sixDoFRigidBodyMotionRestraint::New(iter().dict())
                );
            }

            i++;
        }

        restraints_.setSize(i);
    }
}


void Foam::sixDoFRigidBodyMotion::addConstraints
(
    const dictionary& dict
)
{

}


void Foam::sixDoFRigidBodyMotion::updatePosition
(
    scalar deltaT
)
{
    // First leapfrog velocity adjust and motion part, required before
    // force calculation

    if (Pstream::master())
    {
        v() += 0.5*deltaT*a();

        pi() += 0.5*deltaT*tau();

        // Leapfrog move part
        centreOfMass() += deltaT*v();

        // Leapfrog orientation adjustment

        tensor R;

        R = rotationTensorX(0.5*deltaT*pi().x()/momentOfInertia_.xx());
        pi() = pi() & R;
        Q() = Q() & R;

        R = rotationTensorY(0.5*deltaT*pi().y()/momentOfInertia_.yy());
        pi() = pi() & R;
        Q() = Q() & R;

        R = rotationTensorZ(deltaT*pi().z()/momentOfInertia_.zz());
        pi() = pi() & R;
        Q() = Q() & R;

        R = rotationTensorY(0.5*deltaT*pi().y()/momentOfInertia_.yy());
        pi() = pi() & R;
        Q() = Q() & R;

        R = rotationTensorX(0.5*deltaT*pi().x()/momentOfInertia_.xx());
        pi() = pi() & R;
        Q() = Q() & R;

    }

    Pstream::scatter(motionState_);
}


void Foam::sixDoFRigidBodyMotion::updateForce
(
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT
)
{
    // Second leapfrog velocity adjust part, required after motion and
    // force calculation

    if (Pstream::master())
    {
        a() = fGlobal/mass_;

        tau() = (Q().T() & tauGlobal);

        applyRestraints();

        applyConstraints();

        v() += 0.5*deltaT*a();

        pi() += 0.5*deltaT*tau();
    }

    Pstream::scatter(motionState_);
}


void Foam::sixDoFRigidBodyMotion::updateForce
(
    const pointField& positions,
    const vectorField& forces,
    scalar deltaT
)
{
    vector a = vector::zero;

    vector tau = vector::zero;

    if (Pstream::master())
    {
        forAll(positions, i)
        {
            const vector& f = forces[i];

            a += f/mass_;

            tau += (positions[i] ^ (Q().T() & f));
        }
    }

    updateForce(a, tau, deltaT);
}


Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyMotion::generatePositions(const pointField& pts) const
{
    return (centreOfMass() + (Q() & (pts - refCentreOfMass_)));
}


// ************************************************************************* //
