/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::functionObjects::reference

Group
    grpFieldFunctionObjects

Description
    Calculates and outputs a field whose values are offset to a reference
    value obtained by sampling the field at a user-specified location.

    The field values are calculated using:

    \f[
        ref_c = s(f_{c}(t) - f_p + f_{off})
    \f]

    where
    \vartable
        ref_c    | field values at cell
        s        | optional scale factor (default = 1)
        f_{c}(t) | current field values at cell at this time
        f_p      | field value at position
        f_{off}  | offset field value (default = 0)
    \endvartable

Usage
    Example of function object specification to calculate the reference field:
    \verbatim
    pRef
    {
        type        reference;
        libs        ("libfieldFunctionObjects.so");
        ...
        field       p;
        result      pRef;
        position    (0 0 0);
        scale       1.2;
        offset      100000;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required | Default value
        type         | Type name: reference    | yes      |
        field        | Name of field           | yes      |
        result       | Name of result field    | no       | reference(\<field\>)
        position     | Position to sample      | no       | \<not used\>
        scale        | Scale value             | no       | 1
        offset       | Offset value            | no       | zero
        log          | Log to standard output  | no       | yes
    \endtable


See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    reference.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_reference_H
#define functionObjects_reference_H

#include "fieldExpression.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class reference Declaration
\*---------------------------------------------------------------------------*/

class reference
:
    public fieldExpression
{
    // Private data

        //- Local copy of dictionary used for construction
        dictionary localDict_;

        //- Sample location
        point position_;

        //- Flag to indicate that the position is set
        bool positionIsSet_;

        //- Sample cell
        label celli_;

        //- Interpolation scheme
        word interpolationScheme_;

        //- Scale factor
        scalar scale_;


    // Private Member Functions

        //- Calculate the reference field and return true if successful
        virtual bool calc();

        //- Utility function to calc a given type of field
        template<class Type>
        bool calcType();


public:

    //- Runtime type information
    TypeName("reference");


    // Constructors

        //- Construct from Time and dictionary
        reference
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reference();


    // Public Member Functions

        //- Read the fieldExpression data
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "referenceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
