/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR AUTHOR,AFFILIATION
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

Description
    Template for use with dynamic code generation of a
    points0MotionSolver.

SourceFiles
    codedPoints0MotionSolverTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef codedPoints0MotionSolverTemplate_H
#define codedPoints0MotionSolverTemplate_H

#include "points0MotionSolver.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     A templated codedPoints0MotionSolver
\*---------------------------------------------------------------------------*/

class ${typeName}Points0MotionSolver
:
    public points0MotionSolver
{
    // Private Member Functions

        //- No copy construct
        ${typeName}Points0MotionSolver
        (
            const ${typeName}Points0MotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const ${typeName}Points0MotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from mesh and dictionary
        ${typeName}Points0MotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary&
        );


    //- Destructor
    virtual ~${typeName}Points0MotionSolver();


    // Member Functions

        //- Provide current points for motion. Uses current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
