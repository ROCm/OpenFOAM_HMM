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

#include "KRR4.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(KRR4, 0);
    addToRunTimeSelectionTable(ODESolver, KRR4, dictionary);

const scalar
    KRR4::gamma = 1.0/2.0,
    KRR4::a21 = 2.0, KRR4::a31 = 48.0/25.0, KRR4::a32 = 6.0/25.0,
    KRR4::c21 = -8.0, KRR4::c31 = 372.0/25.0, KRR4::c32 = 12.0/5.0,
    KRR4::c41 = -112.0/125.0, KRR4::c42 = -54.0/125.0, KRR4::c43 = -2.0/5.0,
    KRR4::b1 = 19.0/9.0, KRR4::b2 = 1.0/2.0, KRR4::b3 = 25.0/108.0,
    KRR4::b4 = 125.0/108.0,
    KRR4::e1 = 17.0/54.0, KRR4::e2 = 7.0/36.0, KRR4::e3 = 0.0,
    KRR4::e4 = 125.0/108.0,
    KRR4::c1X = 1.0/2.0, KRR4::c2X = -3.0/2.0, KRR4::c3X = 121.0/50.0,
    KRR4::c4X = 29.0/250.0,
    KRR4::a2X = 1.0, KRR4::a3X = 3.0/5.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KRR4::KRR4(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    adaptiveSolver(ode, dict),
    g1_(n_),
    g2_(n_),
    g3_(n_),
    g4_(n_),
    err_(n_),
    dydx_(n_),
    dfdx_(n_),
    dfdy_(n_, n_),
    a_(n_, n_),
    pivotIndices_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::KRR4::solve
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

    forAll(g1_, i)
    {
        g1_[i] = dydx0[i] + dx*c1X*dfdx_[i];
    }

    LUBacksubstitute(a_, pivotIndices_, g1_);

    forAll(y, i)
    {
        y[i] = y0[i] + a21*g1_[i];
    }

    ode.derivatives(x0 + a2X*dx, y, dydx_);

    forAll(g2_, i)
    {
        g2_[i] = dydx_[i] + dx*c2X*dfdx_[i] + c21*g1_[i]/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, g2_);

    forAll(y, i)
    {
        y[i] = y0[i] + a31*g1_[i] + a32*g2_[i];
    }

    ode.derivatives(x0 + a3X*dx, y, dydx_);

    forAll(g3_, i)
    {
        g3_[i] = dydx0[i] + dx*c3X*dfdx_[i] + (c31*g1_[i] + c32*g2_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, g3_);

    forAll(g4_, i)
    {
        g4_[i] = dydx_[i] + dx*c4X*dfdx_[i]
           + (c41*g1_[i] + c42*g2_[i] + c43*g3_[i])/dx;
    }

    LUBacksubstitute(a_, pivotIndices_, g4_);

    forAll(y, i)
    {
        y[i] = y0[i] + b1*g1_[i] + b2*g2_[i] + b3*g3_[i] + b4*g4_[i];
        err_[i] = e1*g1_[i] + e2*g2_[i] + e3*g3_[i] + e4*g4_[i];
    }

    return normalizeError(y0, y, err_);
}


void Foam::KRR4::solve
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
