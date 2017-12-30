/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::functionObjects::setFlow

Group
    grpFieldFunctionObjects

Description
    Provides options to set the velocity and flux fields as a function of time

    Useful for testing advection behaviour of numerical schemes by e.g.
    imposing solid body rotation, vortex flows.  All types include a scaling
    Foam::Function1 type enabling the strength of the transformation to vary as
    a function of time.

Usage
    Example of function object specification:
    \verbatim
    setFlow1
    {
        type            setFlow;
        libs            ("libfieldFunctionObjects.so");
        ...
        mode            rotation;
        scale           1;
        reverseTime     1;
        omega           6.28318530718;
        origin          (0.5 0 0.5);
        refDir          (1 0 0);
        axis            (0 1 0);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description                | Required | Default value
        type         | Type name: setFlow         | yes      |
        U            | Name of velocity field     | no       | U
        rho          | Name of density field      | no       | none
        phi          | Name of flux field         | no       | phi
        mode         | operating mode - see below | yes      |
    \endtable

    Available \c mode types include:
    - function
    - rortation
    - vortex2D
    - vortex3D


See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    setFlow.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_setFlow_H
#define functionObjects_setFlow_H

#include "fvMeshFunctionObject.H"
#include "Function1.H"
#include "Enum.H"
#include "point.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class setFlow Declaration
\*---------------------------------------------------------------------------*/

class setFlow
:
    public fvMeshFunctionObject
{
    enum class modeType
    {
        FUNCTION,
        ROTATION,
        VORTEX2D,
        VORTEX3D
    };

   static const Enum<modeType> modeTypeNames;


   // Private Data

        //- Name of vecloity field, default = "U"
        word UName_;

        //- Name of density field, default = "none"
        word rhoName_;

        //- Name of flux field, default = "phi"
        word phiName_;

        //- Operating mode
        modeType mode_;

        //- Reverse time
        scalar reverseTime_;

        //- Scaling function
        autoPtr<Function1<scalar>> scalePtr_;


        // Rotation

            //- Origin
            point origin_;

            //- Rotation tensor for rotational mode
            tensor R_;

            //- Rotational speed function
            autoPtr<Function1<scalar>> omegaPtr_;


        // Function

            //- Velocity function
            autoPtr<Function1<vector>> velocityPtr_;


    // Private Member Functions

        //- Set the flux field
        void setPhi(const volVectorField& U);


public:

    //- Runtime type information
    TypeName("setFlow");


    // Constructors

        //- Construct from Time and dictionary
        setFlow
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~setFlow();


    virtual bool read(const dictionary& dict);

    virtual bool execute();

    virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
