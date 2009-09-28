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

#include "sixDofRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDofRigidBodyMotion::sixDofRigidBodyMotion()
:
    centreOfMass_(vector::zero),
    refCentreOfMass_(vector::zero),
    momentOfInertia_(diagTensor::one*VSMALL),
    mass_(VSMALL),
    Q_(I),
    v_(vector::zero),
    a_(vector::zero),
    pi_(vector::zero),
    tau_(vector::zero)
{}


Foam::sixDofRigidBodyMotion::sixDofRigidBodyMotion
(
    const point& centreOfMass,
    const point& refCentreOfMass,
    const diagTensor& momentOfInertia,
    scalar mass,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau
)
:
    centreOfMass_(centreOfMass),
    refCentreOfMass_(refCentreOfMass),
    momentOfInertia_(momentOfInertia),
    mass_(mass),
    Q_(Q),
    v_(v),
    a_(a),
    pi_(pi),
    tau_(tau)
{}


Foam::sixDofRigidBodyMotion::sixDofRigidBodyMotion(const dictionary& dict)
:
    centreOfMass_(dict.lookup("centreOfMass")),
    refCentreOfMass_(dict.lookupOrDefault("refCentreOfMass", centreOfMass_)),
    momentOfInertia_(dict.lookup("momentOfInertia")),
    mass_(readScalar(dict.lookup("mass"))),
    Q_(dict.lookupOrDefault("Q", tensor(I))),
    v_(dict.lookupOrDefault("v", vector::zero)),
    a_(dict.lookupOrDefault("a", vector::zero)),
    pi_(dict.lookupOrDefault("pi", vector::zero)),
    tau_(dict.lookupOrDefault("tau", vector::zero))
{}


Foam::sixDofRigidBodyMotion::sixDofRigidBodyMotion
(
    const sixDofRigidBodyMotion& sDofRBM
)
:
    centreOfMass_(sDofRBM.centreOfMass()),
    refCentreOfMass_(sDofRBM.refCentreOfMass()),
    momentOfInertia_(sDofRBM.momentOfInertia()),
    mass_(sDofRBM.mass()),
    Q_(sDofRBM.Q()),
    v_(sDofRBM.v()),
    a_(sDofRBM.a()),
    pi_(sDofRBM.pi()),
    tau_(sDofRBM.tau())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDofRigidBodyMotion::~sixDofRigidBodyMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDofRigidBodyMotion::updatePosition
(
    scalar deltaT
)
{
    // First leapfrog velocity adjust and motion part, required before
    // force calculation

    v_ += 0.5*deltaT*a_;

    pi_ += 0.5*deltaT*tau_;

    // Leapfrog move part
    centreOfMass_ += deltaT*v_;

    // Leapfrog orientation adjustment

    tensor R;

    R = rotationTensorX(0.5*deltaT*pi_.x()/momentOfInertia_.xx());
    pi_ = pi_ & R;
    Q_ = Q_ & R;

    R = rotationTensorY(0.5*deltaT*pi_.y()/momentOfInertia_.yy());
    pi_ = pi_ & R;
    Q_ = Q_ & R;

    R = rotationTensorZ(deltaT*pi_.z()/momentOfInertia_.zz());
    pi_ = pi_ & R;
    Q_ = Q_ & R;

    R = rotationTensorY(0.5*deltaT*pi_.y()/momentOfInertia_.yy());
    pi_ = pi_ & R;
    Q_ = Q_ & R;

    R = rotationTensorX(0.5*deltaT*pi_.x()/momentOfInertia_.xx());
    pi_ = pi_ & R;
    Q_ = Q_ & R;
}


void Foam::sixDofRigidBodyMotion::updateForce
(
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT
)
{
    // Second leapfrog velocity adjust part, required after motion and
    // force calculation part

    a_ = fGlobal/mass_;

    tau_ = (Q_.T() & tauGlobal);

    v_ += 0.5*deltaT*a_;

    pi_ += 0.5*deltaT*tau_;

}


void Foam::sixDofRigidBodyMotion::updateForce
(
    const pointField& positions,
    const vectorField& forces,
    scalar deltaT
)
{
    // Second leapfrog velocity adjust part, required after motion and
    // force calculation part

    a_ = vector::zero;

    tau_ = vector::zero;

    forAll(positions, i)
    {
        const vector& f = forces[i];

        a_ += f/mass_;

        tau_ += (positions[i] ^ (Q_.T() & f));
    }

    v_ += 0.5*deltaT*a_;

    pi_ += 0.5*deltaT*tau_;

}


Foam::tmp<Foam::pointField>
Foam::sixDofRigidBodyMotion::generatePositions(const pointField& pts) const
{
    return (centreOfMass_ + (Q_ & (pts - refCentreOfMass_)));
}


// ************************************************************************* //
