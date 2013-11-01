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

#include "RKCK5.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RKCK5, 0);
    addToRunTimeSelectionTable(ODESolver, RKCK5, ODE);

const scalar
    RKCK5::a2 = 0.2, RKCK5::a3 = 0.3, RKCK5::a4 = 0.6, RKCK5::a5 = 1.0,
    RKCK5::a6 = 0.875,
    RKCK5::b21 = 0.2, RKCK5::b31 = 3.0/40.0, RKCK5::b32 = 9.0/40.0,
    RKCK5::b41 = 0.3, RKCK5::b42 = -0.9, RKCK5::b43 = 1.2,
    RKCK5::b51 = -11.0/54.0, RKCK5::b52 = 2.5, RKCK5::b53 = -70.0/27.0,
    RKCK5::b54 = 35.0/27.0,
    RKCK5::b61 = 1631.0/55296.0, RKCK5::b62 = 175.0/512.0,
    RKCK5::b63 = 575.0/13824.0, RKCK5::b64 = 44275.0/110592.0,
    RKCK5::b65 = 253.0/4096.0,
    RKCK5::c1 = 37.0/378.0, RKCK5::c3 = 250.0/621.0,
    RKCK5::c4 = 125.0/594.0, RKCK5::c6 = 512.0/1771.0,
    RKCK5::dc1 = RKCK5::c1 - 2825.0/27648.0,
    RKCK5::dc3 = RKCK5::c3 - 18575.0/48384.0,
    RKCK5::dc4 = RKCK5::c4 - 13525.0/55296.0, RKCK5::dc5 = -277.00/14336.0,
    RKCK5::dc6 = RKCK5::c6 - 0.25;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RKCK5::RKCK5(const ODESystem& ode)
:
    ODESolver(ode),
    adaptiveSolver(ode),
    yTemp_(n_),
    k2_(n_),
    k3_(n_),
    k4_(n_),
    k5_(n_),
    k6_(n_),
    err_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::RKCK5::solve
(
    const ODESystem& odes,
    const scalar x0,
    const scalarField& y0,
    const scalarField& dydx0,
    const scalar dx,
    scalarField& y,
    const scalar eps
) const
{
    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i] + b21*dx*dydx0[i];
    }

    odes.derivatives(x0 + a2*dx, yTemp_, k2_);

    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i] + dx*(b31*dydx0[i] + b32*k2_[i]);
    }

    odes.derivatives(x0 + a3*dx, yTemp_, k3_);

    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i] + dx*(b41*dydx0[i] + b42*k2_[i] + b43*k3_[i]);
    }

    odes.derivatives(x0 + a4*dx, yTemp_, k4_);

    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i]
          + dx*(b51*dydx0[i] + b52*k2_[i] + b53*k3_[i] + b54*k4_[i]);
    }

    odes.derivatives(x0 + a5*dx, yTemp_, k5_);

    forAll(yTemp_, i)
    {
        yTemp_[i] = y0[i]
          + dx
           *(b61*dydx0[i] + b62*k2_[i] + b63*k3_[i] + b64*k4_[i] + b65*k5_[i]);
    }

    odes.derivatives(x0 + a6*dx, yTemp_, k6_);

    forAll(y, i)
    {
        y[i] = y0[i]
          + dx*(c1*dydx0[i] + c3*k3_[i] + c4*k4_[i] + c6*k6_[i]);
    }

    forAll(err_, i)
    {
        err_[i] =
            dx
           *(dc1*dydx0[i] + dc3*k3_[i] + dc4*k4_[i] + dc5*k5_[i] + dc6*k6_[i]);
    }

    return normalizeError(eps, y0, y, err_);
}


void Foam::RKCK5::solve
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
    adaptiveSolver::solve(odes, x, y, dydx, eps, dxTry, dxDid, dxNext);
}


// ************************************************************************* //
