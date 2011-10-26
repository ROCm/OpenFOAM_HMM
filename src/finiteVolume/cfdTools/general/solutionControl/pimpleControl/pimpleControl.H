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
    Foam::pimpleControl

Description
    PIMPLE control class to supply convergence information/checks for
    the PIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef pimpleControl_H
#define pimpleControl_H

#include "solutionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pimpleControl Declaration
\*---------------------------------------------------------------------------*/

class pimpleControl
:
    public solutionControl
{
protected:

    // Protected data

        // Solution controls

            //- Maximum number of PIMPLE correctors
            label nCorrPIMPLE_;

            //- Maximum number of PISO correctors
            label nCorrPISO_;

            //- Maximum number of non-orthogonal correctors
            label nCorrNonOrtho_;

            //- Current PIMPLE corrector
            label corrPIMPLE_;

            //- Current PISO corrector
            label corrPISO_;

            //- Current non-orthogonal corrector
            label corrNonOrtho_;

            //- Flag to indicate whether to only solve turbulence on final iter
            bool turbOnFinalIterOnly_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void read();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

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

        // Access

            //- Maximum number of PIMPLE correctors
            inline label nCorrPIMPLE() const;

            //- Maximum number of PISO correctors
            inline label nCorrPISO() const;

            //- Maximum number of non-orthogonal correctors
            inline label nCorrNonOrtho() const;

            //- Current PIMPLE corrector index
            inline label corrPIMPLE() const;

            //- Current PISO corrector index
            inline label corrPISO() const;

            //- Current non-orthogonal corrector index
            inline label corrNonOrtho() const;


        // Solution control

            //- PIMPLE loop
            inline bool loop();

            //- Corrector loop
            inline bool correct();

            //- Non-orthogonal corrector loop
            inline bool correctNonOrthogonal();

            //- Helper function to identify when to store the intial residuals
            inline bool storeInitialResiduals() const;

            //- Helper function to identify final PIMPLE (outer) iteration
            inline bool finalIter() const;

            //- Helper function to identify final non-orthogonal iteration
            inline bool finalNonOrthogonalIter() const;

            //- Helper function to identify final inner iteration
            inline bool finalInnerIter() const;

            //- Helper function to identify whether to solve for turbulence
            inline bool turbCorr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pimpleControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
