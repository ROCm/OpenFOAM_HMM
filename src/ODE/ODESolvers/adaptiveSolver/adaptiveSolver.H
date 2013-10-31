/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::adaptiveSolver

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
    adaptiveSolver.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveSolver_H
#define adaptiveSolver_H

#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class adaptiveSolver Declaration
\*---------------------------------------------------------------------------*/

class adaptiveSolver
{
    // Private data

        //- Step-size adjustment controls
        static const scalar safeScale, alphaInc, alphaDec, minScale, maxScale;

        mutable scalarField yTemp_;


public:

    // Constructors

        //- Construct from ODESystem
        adaptiveSolver(const ODESystem& ode);


    //- Destructor
    virtual ~adaptiveSolver()
    {}


    // Member Functions

        //- Return the nomalized scalar error
        scalar normalizeError
        (
            const scalar eps,
            const scalarField& y0,
            const scalarField& y,
            const scalarField& err
        ) const;

        //- Solve a single step dx and return the error
        virtual scalar solve
        (
            const ODESystem& ode,
            const scalar x0,
            const scalarField& y0,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y,
            const scalar eps
        ) const = 0;

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
