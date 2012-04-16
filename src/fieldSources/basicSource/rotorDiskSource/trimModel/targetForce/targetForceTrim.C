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

#include "targetForceTrim.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(targetForceTrim, 0);

    addToRunTimeSelectionTable(trimModel, targetForceTrim, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vector Foam::targetForceTrim::calcForce
(
    const vectorField& U,
    const scalarField& thetag,
    vectorField& force    
) const
{
    rotor_.calculate(U, thetag, force, false, false);

    const labelList& cells = rotor_.cells();
    const vectorField& C = rotor_.mesh().C();

    const vector& origin = rotor_.coordSys().origin();
    const vector& rollAxis = rotor_.coordSys().e1();
    const vector& pitchAxis = rotor_.coordSys().e2();
    const vector& yawAxis = rotor_.coordSys().e3();

    vector f(vector::zero);
    forAll(cells, i)
    {
        label cellI = cells[i];

        vector moment = force[cellI]^(C[cellI] - origin);
        f[0] += force[cellI] & yawAxis;
        f[1] += moment & pitchAxis;
        f[2] += moment & rollAxis;
    }

    reduce(f, sumOp<vector>());
    
    return f;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::targetForceTrim::targetForceTrim
(
    const rotorDiskSource& rotor,
    const dictionary& dict
)
:
    trimModel(rotor, dict, typeName),
    calcFrequency_(-1),
    target_(vector::zero),
    theta_(vector::zero),
    nIter_(50),
    tol_(1e-8),
    relax_(1.0),
    dTheta_(degToRad(0.1))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::targetForceTrim::~targetForceTrim()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::targetForceTrim::read(const dictionary& dict)
{
    trimModel::read(dict);

    const dictionary& targetDict(coeffs_.subDict("target"));
    target_[0] = readScalar(targetDict.lookup("fThrust"));
    target_[1] = readScalar(targetDict.lookup("mPitch"));
    target_[2] = readScalar(targetDict.lookup("mRoll"));

    const dictionary& pitchAngleDict(coeffs_.subDict("pitchAngles"));
    theta_[0] = degToRad(readScalar(pitchAngleDict.lookup("theta0Ini")));
    theta_[1] = degToRad(readScalar(pitchAngleDict.lookup("theta1cIni")));
    theta_[2] = degToRad(readScalar(pitchAngleDict.lookup("theta1sIni")));

    coeffs_.lookup("calcFrequency") >> calcFrequency_;

    coeffs_.readIfPresent("nIter", nIter_);
    coeffs_.readIfPresent("tol", tol_);
    coeffs_.readIfPresent("relax", relax_);

    if (coeffs_.readIfPresent("dTheta", dTheta_))
    {
        dTheta_ = degToRad(dTheta_);
    }
}


Foam::tmp<Foam::scalarField> Foam::targetForceTrim::thetag() const
{
    const List<vector>& x = rotor_.x();

    tmp<scalarField> ttheta(new scalarField(x.size()));
    scalarField& t = ttheta();

    forAll(t, i)
    {
        scalar psi = x[i].y();
        if (psi < 0)
        {
            psi += mathematical::twoPi;
        }

        t[i] = theta_[0] + theta_[1]*cos(psi) + theta_[2]*sin(psi);
    }

    return ttheta;
}


void Foam::targetForceTrim::correct(const vectorField& U, vectorField& force)
{
    if (rotor_.mesh().time().timeIndex() % calcFrequency_ == 0)
    {
        // iterate to find new pitch angles to achieve target force
        scalar err = GREAT;
        label iter = 0;
        tensor J(tensor::zero);

        while ((err > tol_) && (iter < nIter_))
        {
            // cache initial theta vector
            vector theta0(theta_);

            // set initial values
            vector old = calcForce(U, thetag(), force);

            // construct Jacobian by perturbing the pitch angles
            // by +/-(dTheta_/2)
            for (label pitchI = 0; pitchI < 3; pitchI++)
            {
                theta_[pitchI] -= dTheta_/2.0;
                vector f0 = calcForce(U, thetag(), force);

                theta_[pitchI] += dTheta_;
                vector f1 = calcForce(U, thetag(), force);

                vector ddTheta = (f1 - f0)/dTheta_;

                J[pitchI + 0] = ddTheta[0];
                J[pitchI + 3] = ddTheta[1];
                J[pitchI + 6] = ddTheta[2];

                theta_ = theta0;
            }

            // calculate the change in pitch angle vector
            vector dt = inv(J) & (target_ - old);

            // update pitch angles
            vector thetaNew = theta_ + relax_*dt;

            // update error
            err = mag(thetaNew - theta_);

            // update for next iteration
            theta_ = thetaNew;
            iter++;
        }

        if (iter == nIter_)
        {
            WarningIn
            (
                "void Foam::targetForceTrim::correct"
                "("
                    "const vectorField&, "
                    "vectorField&"
                ")"
            )   << "Trim routine not converged in " << iter
                << " iterations, max residual = " << err << endl;
        }
        else
        {
            Info<< type() << ": converged in " << iter
                << " iterations" << endl;
        }

        Info<< "    new pitch angles:" << nl
            << "        theta0  = " << radToDeg(theta_[0]) << nl
            << "        theta1c = " << radToDeg(theta_[1]) << nl
            << "        theta1s = " << radToDeg(theta_[2]) << nl
            << endl;
    }
}


// ************************************************************************* //
