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
    if (Pstream::master())
    {
        forAll(restraints_, rI)
        {
            // restraint position
            point rP = vector::zero;

            // restraint force
            vector rF = vector::zero;

            // restraint moment
            vector rM = vector::zero;

            restraints_[rI].restrain(*this, rP, rF, rM);

            Info<< "Restraint " << rI << " force " << rF << endl;

            a() += rF/mass_;

            tau() += Q().T() & (rM + ((rP - centreOfMass()) ^ rF));
        }
    }
}


void Foam::sixDoFRigidBodyMotion::applyConstraints(scalar deltaT)
{
    if (Pstream::master())
    {
        label iter = 0;

        bool converged = true;

        // constraint force accumulator
        vector cFA = vector::zero;

        // constraint moment accumulator
        vector cMA = vector::zero;

        do
        {
            iter++;

            if (iter > maxConstraintIters_)
            {
                FatalErrorIn
                (
                    "Foam::sixDoFRigidBodyMotion::applyConstraints"
                    "(scalar deltaT)"
                )
                    << nl << "Maximum number of constraint iterations ("
                    << maxConstraintIters_ << ") exceeded." << nl
                    << exit(FatalError);

            }

            converged = true;

            forAll(constraints_, cI)
            {
                // Info<< "Constraint " << cI << endl;

                // constraint position
                point cP = vector::zero;

                // constraint force
                vector cF = vector::zero;

                // constraint moment
                vector cM = vector::zero;

                converged = converged && constraints_[cI].constrain
                (
                    *this,
                    cFA,
                    cMA,
                    deltaT,
                    cP,
                    cF,
                    cM
                );

                // Accumulate constraint force
                cFA += cF;

                // Accumulate constraint moment
                cMA += cM + ((cP - centreOfMass()) ^ cF);
            }

        } while(!converged);

        Info<< "Constraints converged in " << iter << " iterations" << nl
            << "Constraint force: " << cFA << nl
            << "Constraint moment: " << cMA
            << endl;

        // Add the constrain forces and moments to the motion state variables
        a() += cFA/mass_;

        // The moment of constraint forces has already been added
        // during accumulation
        tau() += Q().T() & cMA;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion()
:
    motionState_(),
    restraints_(),
    constraints_(),
    maxConstraintIters_(0),
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
    constraints_(),
    maxConstraintIters_(0),
    refCentreOfMass_(refCentreOfMass),
    momentOfInertia_(momentOfInertia),
    mass_(mass)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion(const dictionary& dict)
:
    motionState_(dict),
    restraints_(),
    constraints_(),
    maxConstraintIters_(0),
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
    constraints_(sDoFRBM.constraints()),
    maxConstraintIters_(sDoFRBM.maxConstraintIters()),
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
    if (dict.found("constraints"))
    {
        const dictionary& constraintDict = dict.subDict("constraints");

        label i = 0;

        constraints_.setSize(constraintDict.size());

        forAllConstIter(IDLList<entry>, constraintDict, iter)
        {
            if (iter().isDict())
            {
                Info<< "Adding constraint: " << iter().keyword() << endl;

                constraints_.set
                (
                    i++,
                    sixDoFRigidBodyMotionConstraint::New(iter().dict())
                );
            }
        }

        constraints_.setSize(i);

        if (constraints_.size())
        {
            maxConstraintIters_ = readLabel
            (
                constraintDict.lookup("maxIterations")
            );
        }
    }
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

        rotate(Q(), pi(), deltaT);
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

        applyConstraints(deltaT);

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

            tau += Q().T() & ((positions[i] - centreOfMass()) ^ f);
        }
    }

    updateForce(a, tau, deltaT);
}


Foam::point Foam::sixDoFRigidBodyMotion::predictedPosition
(
    const point& pt,
    const vector deltaForce,
    const vector deltaMoment,
    scalar deltaT
) const
{
    vector vTemp = v() + deltaT*(a() + deltaForce/mass_);

    vector piTemp = pi() + deltaT*(tau() + (Q().T() & deltaMoment));

    point centreOfMassTemp = centreOfMass() + deltaT*vTemp;

    tensor QTemp = Q();

    rotate(QTemp, piTemp, deltaT);

    return (centreOfMassTemp + (QTemp & (pt - refCentreOfMass_)));
}


// ************************************************************************* //
