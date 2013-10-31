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

#include "adaptiveSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

const scalar
    adaptiveSolver::safeScale=0.9,
    adaptiveSolver::alphaInc=0.2,
    adaptiveSolver::alphaDec=0.25,
    adaptiveSolver::minScale=0.2,
    adaptiveSolver::maxScale=10;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveSolver::adaptiveSolver(const ODESystem& ode)
:
    yTemp_(ode.nEqns())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::adaptiveSolver::normalizeError
(
    const scalar eps,
    const scalarField& y0,
    const scalarField& y,
    const scalarField& err
) const
{
    scalar epsAbs = 1e-10;
    scalar epsRel = eps;

    // Calculate the maximum error
    scalar maxErr = 0.0;
    forAll(err, i)
    {
        scalar eps = epsAbs + epsRel*max(mag(y0[i]), mag(y[i]));
        maxErr = max(maxErr, mag(err[i])/eps);
    }

    return maxErr;
}


void Foam::adaptiveSolver::solve
(
    const ODESystem& odes,
    scalar& x,
    scalarField& y,
    scalarField& dydx,
    const scalar eps,
    const scalar dxTry,
    scalar& dxDid,
    scalar& dxNext
) const
{
    scalar dx = dxTry;
    scalar err = 0.0;

    // Loop over solver and adjust step-size as necessary
    // to achieve desired error
    do
    {
        // Solve step and provide error estimate
        err = solve(odes, x, y, dydx, dx, yTemp_, eps);

        // If error is large reduce dx
        if (err > 1)
        {
            scalar scale = max(safeScale*pow(err, -alphaDec), minScale);
            dx *= scale;

            if (dx < VSMALL)
            {
                FatalErrorIn("adaptiveSolver::solve")
                    << "stepsize underflow"
                    << exit(FatalError);
            }
        }
    } while (err > 1);

    // Update the state
    dxDid = dx;
    x += dx;
    y = yTemp_;

    // If the error is small increase the step-size
    dxNext = min(max(safeScale*pow(err, -alphaInc), minScale), maxScale)*dx;
}


// ************************************************************************* //
