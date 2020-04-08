/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::functionObjects::momentumError

Group
    grpForcesFunctionObjects

Description
    Produces a balance terms result for the steady momentum equation

Usage
    Example of function object specification:
    \verbatim
    momErr
    {
        type        momentumError;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: momentumError| yes         |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        phi          | Flux field name         | no          | phi
    \endtable
See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    momentumError.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_momentumError_H
#define functionObjects_momentumError_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class momentumError Declaration
\*---------------------------------------------------------------------------*/

class momentumError
:
    public fvMeshFunctionObject
{

protected:

    // Protected data

        // Read from dictionary

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Flux
            word phiName_;


    // Protected Member Functions

         //- Return the effective viscous stress (laminar + turbulent).
        tmp<volVectorField> divDevRhoReff();

        //- No copy construct
        momentumError(const momentumError&) = delete;

        //- No copy assignment
        void operator=(const momentumError&) = delete;


public:

    //- Runtime type information
    TypeName("momentumError");


    // Constructors

        //- Construct from Time and dictionary
        momentumError
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~momentumError() = default;

    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();

        //- Calculate the momentum error
        void calcMomentError();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
