/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::functionObjects::lambVector

Group
    grpFieldFunctionObjects

Description
    Calculates Lamb vector, i.e. the cross product of vorticity and velocity.
    Calculation of the divergence of the Lamb vector can be performed by using
    'div' functionObject on this 'lambVector' functionObject.

    The field is stored on the mesh database so that it can be retrieved
    and used for other applications.

Usage
    \verbatim
    lambVector1
    {
        type        lambVector;
        libs        ("libfieldFunctionObjects.so");
        field       UMean;
        ... Base options of fieldExpression ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: lambVector   | yes         |
        field        | Name of volVectorField  | no          | U
    \endtable

See also
    Foam::functionObjects::fieldExpression

SourceFiles
    lambVector.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_lambVector_H
#define functionObjects_lambVector_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class lambVector Declaration
\*---------------------------------------------------------------------------*/

class lambVector
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the lambVector field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("lambVector");


    // Constructors

        //- Construct from Time and dictionary
        lambVector
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~lambVector() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
