/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::objectives::objectiveForceTarget

Description

SourceFiles
    objectiveForceTarget.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveForceTarget_H
#define objectiveForceTarget_H

#include "objectiveForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                    Class objectiveForceTarget Declaration
\*---------------------------------------------------------------------------*/

class objectiveForceTarget
:
    public objectiveForce
{
    // Private data

        scalar force_;

        scalar target_;


public:

    //- Runtime type information
    TypeName("forceTarget");


    // Constructors

        //- Construct from components
        objectiveForceTarget
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveForceTarget() = default;


    // Member Functions

        //- Return the objective function value
        scalar J();

        //- Write objective value
        void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
