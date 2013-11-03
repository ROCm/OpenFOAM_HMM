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

Class
    Foam::KRR43

Description
    4th Order Kaps Rentrop Rosenbrock stiff ODE solver

    References:
    \verbatim
        "Generalized Runge-Kutta methods of order four with stepsize control
         for stiff ordinary differential equations"
        Kaps, P.,
        Rentrop, P.,
        Numerische Thematic, vol. 33, 1979, pp. 55–68.
    \endverbatim

    Constants from:
    \verbatim
        "Implementation of Rosenbrock Methods"
        Shampine, L. F.,
        ACM Transactions on Mathematical Software, vol. 8, 1982, pp. 93–113.
    \endverbatim

SourceFiles
    KRR43.C

\*---------------------------------------------------------------------------*/

#ifndef KRR43_H
#define KRR43_H

#include "ODESolver.H"
#include "adaptiveSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class KRR43 Declaration
\*---------------------------------------------------------------------------*/

class KRR43
:
    public ODESolver,
    public adaptiveSolver
{
    // Private data

        mutable scalarField k1_;
        mutable scalarField k2_;
        mutable scalarField k3_;
        mutable scalarField k4_;
        mutable scalarField err_;
        mutable scalarField dydx_;
        mutable scalarField dfdx_;
        mutable scalarSquareMatrix dfdy_;
        mutable scalarSquareMatrix a_;
        mutable labelList pivotIndices_;

        static const scalar
            gamma,
            a21, a31, a32,
            c21, c31, c32, c41, c42, c43,
            b1, b2, b3, b4,
            e1, e2, e3, e4,
            c1X, c2X, c3X, c4X,
            a2X, a3X;


public:

    //- Runtime type information
    TypeName("KRR43");


    // Constructors

        //- Construct from ODE
        KRR43(const ODESystem& ode, const dictionary& dict);


    // Member Functions

        //- Solve a single step dx and return the error
        scalar solve
        (
            const ODESystem& ode,
            const scalar x0,
            const scalarField& y0,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y
        ) const;

        //- Solve the ODE system and the update the state
        void solve
        (
            const ODESystem& ode,
            scalar& x,
            scalarField& y,
            scalar& dxTry
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
