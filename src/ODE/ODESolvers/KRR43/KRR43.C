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

#include "KRR43.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(KRR43, 0);
    addToRunTimeSelectionTable(ODESolver, KRR43, dictionary);

const scalar
    KRR43::gamma = 1.0/2.0,
    KRR43::a21 = 2.0,
    KRR43::a31 = 48.0/25.0,
    KRR43::a32 = 6.0/25.0,
    KRR43::c21 = -8.0,
    KRR43::c31 = 372.0/25.0,
    KRR43::c32 = 12.0/5.0,
    KRR43::c41 = -112.0/125.0,
    KRR43::c42 = -54.0/125.0,
    KRR43::c43 = -2.0/5.0,
    KRR43::b1 = 19.0/9.0,
    KRR43::b2 = 1.0/2.0,
    KRR43::b3 = 25.0/108.0,
    KRR43::b4 = 125.0/108.0,
    KRR43::e1 = 17.0/54.0,
    KRR43::e2 = 7.0/36.0,
    KRR43::e3 = 0.0,
    KRR43::e4 = 125.0/108.0,
    KRR43::c1X = 1.0/2.0,
    KRR43::c2X = -3.0/2.0,
    KRR43::c3X = 121.0/50.0,
    KRR43::c4X = 29.0/250.0,
    KRR43::a2X = 1.0,
    KRR43::a3X = 3.0/5.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KRR43::KRR43(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    adaptiveSolver(ode, dict),
    k1_(n_),
    k2_(n_),
    k3_(n_),
    k4_(n_),
    err_(n_),
    dydx_(n_),
    dfdx_(n_),
    dfdy_(n_, n_),
    a_(n_, n_),
    pivotIndices_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::KRR43::solve
(
    const ODESystem& ode,
    const scalar x0,
    const scalarField& y0,
    const scalarField& dydx0,
    const scalar dx,
    scalarField& y
) const
{
    ode.jacobian(x0, y0, dfdx_, dfdy_);

    for (register label i=0; i<n_; i++)
    {
        for (register label j=0; j<n_; j++)
        {
            a_[i][j] = -dfdy_[i][j];
        }

        a_[i][i] += 1.0/(gamma*dx);
    }

    LUDecompose(a_, pivotIndices_);

    forAll(k1_, i)
    {
        k1_[i] = dydx0[i] + dx*c1X*dfdx_[i];
    }

    LUBacksubstitute(a_, pivotIndices_, k1_);

    forAll(y, i)
    {
        y[i] = y0[i] + a21*k1_[i];
    }

    ode.derivatives(x0 + a2X*dx, y, dydx_);

    forAll(k2_, i)
    {
        k2_[i] = dydx_[i] + dx*c2X*dfdx_[i] + c21*k1_[i]/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k2_);

    forAll(y, i)
    {
        y[i] = y0[i] + a31*k1_[i] + a32*k2_[i];
    }

    ode.derivatives(x0 + a3X*dx, y, dydx_);

    forAll(k3_, i)
    {
        k3_[i] = dydx_[i] + dx*c3X*dfdx_[i] + (c31*k1_[i] + c32*k2_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k3_);

    forAll(k4_, i)
    {
        k4_[i] = dydx_[i] + dx*c4X*dfdx_[i]
           + (c41*k1_[i] + c42*k2_[i] + c43*k3_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, k4_);

    forAll(y, i)
    {
        y[i] = y0[i] + b1*k1_[i] + b2*k2_[i] + b3*k3_[i] + b4*k4_[i];
        err_[i] = e1*k1_[i] + e2*k2_[i] + e3*k3_[i] + e4*k4_[i];
    }

    return normalizeError(y0, y, err_);
}


void Foam::KRR43::solve
(
    const ODESystem& odes,
    scalar& x,
    scalarField& y,
    scalar& dxTry
) const
{
    adaptiveSolver::solve(odes, x, y, dxTry);
}


// ************************************************************************* //
