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
    Foam::solutionControl

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef solutionControl_H
#define solutionControl_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class solutionControl Declaration
\*---------------------------------------------------------------------------*/

class solutionControl
{
public:

    struct fieldData
    {
        wordRe name;
        scalar absTol;
        scalar relTol;
        scalar initialResidual;
    };


protected:

    // Protected data

        //- Reference to the mesh database
        fvMesh& mesh_;

        //- List of residual data per field
        List<fieldData> residualControl_;

        //- The dictionary name, e.g. SIMPLE, PIMPLE
        const word algorithmName_;


        // Solution controls

            //- Maximum number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Flag to indicate to solve for momentum
            bool momentumPredictor_;

            //- Flag to indictae to solve using transonic algorithm
            bool transonic_;


    // Protected Member Functions

        //- Read constrols from fvSolution dictionary
        virtual void read(const bool absTolOnly);

        //- Return index of field in residualControl_ if present
        virtual label applyToField(const word& fieldName) const;

        //- Return true if all convergence checks are satified
        virtual bool criteriaSatisfied() = 0;

        //- Disallow default bitwise copy construct
        solutionControl(const solutionControl&);

        //- Disallow default bitwise assignment
        void operator=(const solutionControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solutionControl");


    // Constructors

        //- Construct from mesh
        solutionControl(fvMesh& mesh, const word& algorithmName);


    //- Destructor
    virtual ~solutionControl();


    // Member Functions

        // Access

            //- Return the solution dictionary
            inline const dictionary& dict() const;


        // Solution control

            //- Maximum number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Flag to indicate to solve for momentum
            inline bool momentumPredictor() const;

            //- Flag to indictae to solve using transonic algorithm
            inline bool transonic() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
