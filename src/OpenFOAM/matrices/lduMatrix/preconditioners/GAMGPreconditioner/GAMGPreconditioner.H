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
    Foam::GAMGPreconditioner

Description
    Geometric agglomerated algebraic multigrid preconditioner.

See Also
    GAMGSolver for more details.

SourceFiles
    GAMGPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGPreconditioner_H
#define GAMGPreconditioner_H

#include "GAMGSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GAMGPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class GAMGPreconditioner
:
    public GAMGSolver,
    public lduMatrix::preconditioner
{
    // Private data

        //- Number of V-cycles to perform
        label nVcycles_;


    // Private member functions

        //- Read control parameters from the control dictionary
        void readControls();


public:

    //- Runtime type information
    TypeName("GAMG");


    // Constructors

        //- Construct for given solver and preconditioner data stream
        GAMGPreconditioner
        (
            const lduMatrix::solver& sol,
            Istream& preconditionerData
        );


    // Destructor

        virtual ~GAMGPreconditioner();


    // Member Functions

        //- Read and reset the preconditioner parameters from the given stream
        virtual void read(Istream& preconditionerData);

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
