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
    Foam::RKCK5

Description
    5th Order Cash-Karp Runge-Kutta ODE solver

    References:
    \verbatim
        "A variable order Runge-Kutta method for initial value problems
         with rapidly varying right-hand sides"
        Cash, J.R.,
        Karp, A.H.
        ACM Transactions on Mathematical Software, vol. 16, 1990, pp. 201–222.
    \endverbatim

SourceFiles
    RKCK5.C

\*---------------------------------------------------------------------------*/

#ifndef RKCK5_H
#define RKCK5_H

#include "ODESolver.H"
#include "adaptiveSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class RKCK5 Declaration
\*---------------------------------------------------------------------------*/

class RKCK5
:
    public ODESolver,
    public adaptiveSolver
{
    // Private data

        //- RKCK Constants
        static const scalar
            a2, a3, a4, a5, a6,
            b21, b31, b32, b41, b42, b43,
            b51, b52, b53, b54, b61, b62, b63, b64, b65,
            c1, c3, c4, c6,
            dc1, dc3, dc4, dc5, dc6;

        // Temporary fields
        mutable scalarField yTemp_;
        mutable scalarField k2_;
        mutable scalarField k3_;
        mutable scalarField k4_;
        mutable scalarField k5_;
        mutable scalarField k6_;

        //- Error-estimate field
        mutable scalarField err_;


public:

    //- Runtime type information
    TypeName("RKCK5");


    // Constructors

        //- Construct from ODE
        RKCK5(const ODESystem& ode);


    // Member Functions

        //- Solve a single step dx and return the error
        scalar solve
        (
            const ODESystem& ode,
            const scalar x0,
            const scalarField& y0,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y,
            const scalar eps
        ) const;

        //- Solve the ODE system and the update the state
        void solve
        (
            const ODESystem& ode,
            scalar& x,
            scalarField& y,
            scalarField& dydx,
            const scalar eps,
            const scalar dxTry,
            scalar& dxDid,
            scalar& dxNext
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
