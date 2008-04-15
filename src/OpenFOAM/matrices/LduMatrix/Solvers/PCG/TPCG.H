/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::TPCG

Description
    Preconditioned conjugate gradient solver for symmetric lduMatrices
    using a run-time selectable preconditiioner.

SourceFiles
    TPCG.C

\*---------------------------------------------------------------------------*/

#ifndef TPCG_H
#define TPCG_H

#include "LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TPCG Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class TPCG
:
    public LduMatrix<Type, DType, LUType>::solver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        TPCG(const TPCG&);

        //- Disallow default bitwise assignment
        void operator=(const TPCG&);


public:

    //- Runtime type information
    TypeName("PCG");


    // Constructors

        //- Construct from matrix components and solver data dictionary
        TPCG
        (
            const word& fieldName,
            const LduMatrix<Type, DType, LUType>& matrix,
            const dictionary& solverDict
        );


    // Destructor

        virtual ~TPCG()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual typename LduMatrix<Type, DType, LUType>::solverPerformance solve
        (
            Field<Type>& psi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TPCG.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
