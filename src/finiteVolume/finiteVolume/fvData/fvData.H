/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::fvData

Description
    Database for finite volume solution data, solver performance and
    other reduced data.  fvMesh is derived from fvData so that all fields have
    access to the fvData from the mesh reference they hold.

SourceFiles
    fvData.C

\*---------------------------------------------------------------------------*/

#ifndef fvData_H
#define fvData_H

#include "IOdictionary.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvData Declaration
\*---------------------------------------------------------------------------*/

class fvData
:
    public IOdictionary
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvData(const fvData&);

        //- Disallow default bitwise assignment
        void operator=(const fvData&);


public:

    //- Debug switch
    static int debug;


    // Constructors

        //- Construct for objectRegistry
        fvData(const objectRegistry& obr);


    // Member Functions

        // Access

            //- Return the dictionary of solver performance data
            //  which includes initial and final residuals for convergence
            //  checking
            const dictionary& solverPerformanceDict() const;

            //- Add/set the solverPerformance entry for the named field
            void setSolverPerformance
            (
                const word& name,
                const lduMatrix::solverPerformance&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
