/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "ODESolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ODESolver, 0);
    defineRunTimeSelectionTable(ODESolver, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ODESolver::ODESolver(const ODESystem& ode, const dictionary& dict)
:
    n_(ode.nEqns()),
    absTol_(n_, dict.lookupOrDefault<scalar>("absTol", SMALL)),
    relTol_(n_, dict.lookupOrDefault<scalar>("relTol", 1e-6)),
    maxSteps_(10000)
{}


Foam::ODESolver::ODESolver
(
    const ODESystem& ode,
    const scalarField& absTol,
    const scalarField& relTol
)
:
    n_(ode.nEqns()),
    absTol_(absTol),
    relTol_(relTol),
    maxSteps_(10000)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::ODESolver::normalizeError
(
    const scalarField& y0,
    const scalarField& y,
    const scalarField& err
) const
{
    // Calculate the maximum error
    scalar maxErr = 0.0;
    forAll(err, i)
    {
        scalar tol = absTol_[i] + relTol_[i]*max(mag(y0[i]), mag(y[i]));
        maxErr = max(maxErr, mag(err[i])/tol);
    }

    return maxErr;
}


void Foam::ODESolver::solve
(
    const ODESystem& ode,
    const scalar xStart,
    const scalar xEnd,
    scalarField& y,
    scalar& dxEst
) const
{
    scalar x = xStart;
    bool truncated = false;

    for (label nStep=0; nStep<maxSteps_; nStep++)
    {
        // Store previous iteration dxEst
        scalar dxEst0 = dxEst;

        // Check if this is a truncated step and set dxEst to integrate to xEnd
        if ((x + dxEst - xEnd)*(x + dxEst - xStart) > 0)
        {
            truncated = true;
            dxEst = xEnd - x;
        }

        // Integrate as far as possible up to dxEst
        solve(ode, x, y, dxEst);

        // Check if reached xEnd
        if ((x - xEnd)*(xEnd - xStart) >= 0)
        {
            if (nStep > 0 && truncated)
            {
                dxEst = dxEst0;
            }

            return;
        }
    }

    FatalErrorIn
    (
        "ODESolver::solve"
        "(const ODESystem& ode, const scalar xStart, const scalar xEnd,"
        "scalarField& y, scalar& dxEst) const"
    )   << "Integration steps greater than maximum " << maxSteps_
        << exit(FatalError);
}


// ************************************************************************* //
