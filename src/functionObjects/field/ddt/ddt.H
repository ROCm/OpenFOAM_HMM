/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::functionObjects::ddt

Group
    grpFieldFunctionObjects

Description
    Calculates the Eulerian time derivative of a field.

    The operation can be applied to any volume field generating a field of the
    same type.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    ddt.C
    ddtTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ddt_H
#define functionObjects_ddt_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class ddt Declaration
\*---------------------------------------------------------------------------*/

class ddt
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate and store the result
        //  \return true if fieldName_ was known and had the correct type
        template<class Type>
        bool calcDdt();

        //- Calculate time derivative of the field
        //  \return true on success
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("ddt");


    // Constructors

        //- Construct from Time and dictionary
        ddt
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ddt() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ddtTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
