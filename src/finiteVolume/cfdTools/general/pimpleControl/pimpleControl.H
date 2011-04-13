/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::pimpleControl

Description
    PIMPLE control class to supply convergence information/checks for
    the PIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef pimpleControl_H
#define pimpleControl_H

#include "fvMesh.H"
//#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pimpleControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleControl
{
    struct fieldData
    {
        wordRe name;
        scalar relTol;
        scalar absTol;
        scalar initialResidual;
    };


    // Private data

        //- Reference to the mesh database
        fvMesh& mesh_;

        //- Maximum number of PIMPLE correctors
        int nCorr_;

        //- Current PIMPLE corrector
        int corr_;

        //- List of residual data per field
        List<fieldData> residualControl_;


    // Private Member Functions

        //- Return index of field in residualControl_ if present
        label applyToField(const word& fieldName) const;

        //- Return true if all convergence checks are satified
        bool criteriaSatisfied();

        //- Disallow default bitwise copy construct
        pimpleControl(const pimpleControl&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("pimpleControl");


    // Constructors

        //- Construct from mesh
        pimpleControl(fvMesh& mesh);


    //- Destructor
    virtual ~pimpleControl();


    // Member Functions

        //- Read constrols from fvSolution dictionary
        void read();

        //- Loop start
        bool start();

        //- Loop loop
        bool loop();

        //- Helper function to identify final iteration
        bool finalIter() const;


    // Member Operators

        void operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
