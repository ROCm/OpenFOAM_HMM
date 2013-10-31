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
    Foam::ODESolver

Description
    Abstract base-class for ODE system solvers

SourceFiles
    ODESolver.C

\*---------------------------------------------------------------------------*/

#ifndef ODESolver_H
#define ODESolver_H

#include "ODESystem.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ODESolver Declaration
\*---------------------------------------------------------------------------*/

class ODESolver
{

protected:

    // Private data

        label n_;
        mutable scalarField yScale_;
        mutable scalarField dydx_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ODESolver(const ODESolver&);

        //- Disallow default bitwise assignment
        void operator=(const ODESolver&);


public:

    //- Runtime type information
    TypeName("ODESolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ODESolver,
            ODE,
            (const ODESystem& ode),
            (ode)
        );


    // Constructors

        //- Construct for given ODE
        ODESolver(const ODESystem& ode);


    // Selectors

        //- Select null constructed
        static autoPtr<ODESolver> New
        (
            const word& ODESolverTypeName,
            const ODESystem& ode
        );


    //- Destructor
    virtual ~ODESolver()
    {}


    // Member Functions

        virtual void solve
        (
            const ODESystem& ode,
            scalar& x,
            scalarField& y,
            scalarField& dydx,
            const scalar eps,
            const scalarField& yScale,
            const scalar hTry,
            scalar& hDid,
            scalar& hNext
        ) const = 0;


        virtual void solve
        (
            const ODESystem& ode,
            const scalar xStart,
            const scalar xEnd,
            scalarField& y,
            const scalar eps,
            scalar& hEst
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
