/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::limitFields

Group
    grpFieldFunctionObjects

Description
    Limits fields to user-specified min and max bounds

Note
    For non-scalar field types, the user limit is used to create a
    scaling factor using the field magnitudes.

Usage
    Example of function object specification:

    \verbatim
    limitFields1
    {
        type        limitFields;
        libs        ("libfieldFunctionObjects.so");
        ...
        fields      (U);
        limit       max;
        max         100;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description                          | Required | Default
        type         | type name: limitFields               | yes |
        fields       | list of fields to process            | yes |
        limit        | bound to limit - see below           | yes |
        min          | min limit value                      | partly |
        max          | max limit value                      | partly |
    \endtable

    The \c limit entry can take the value:
    - \c min : specify a minimum value
    - \c max : specify a maximum value
    - \c both : specify a minimum value and a maximum value

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    limitFields.C
    limitFieldsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_limitFields_H
#define functionObjects_limitFields_H

#include "Enum.H"
#include "fvMeshFunctionObject.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class limitFields Declaration
\*---------------------------------------------------------------------------*/

class limitFields
:
    public fvMeshFunctionObject
{
public:

    // Public Enumerations

        enum limitType : unsigned
        {
            MIN = 0x1,          //!< limit by minimum value
            MAX = 0x2,          //!< limit by maximum value
            BOTH = (MIN | MAX)  //!< limit by both minimum and maximum values
        };


protected:

    // Protected Data

        //- Limit type names
        static const Enum<limitType> limitTypeNames_;

        //- Limiting type
        limitType limit_;

        //- Fields to limit
        volFieldSelection fieldSet_;

        //- Minimum limit
        scalar min_;

        //- Maximum limit
        scalar max_;


    // Protected Member Functions

        //- Limit a scalar field.
        //  \return true if field of this type was found.
        bool limitScalarField(const word& fieldName);

        //- Limit a field.
        //  \return true if field of this type was found.
        template<class Type>
        bool limitField(const word& fieldName);


        //- No copy construct
        limitFields(const limitFields&) = delete;

        //- No copy assignment
        void operator=(const limitFields&) = delete;


public:

    //- Runtime type information
    TypeName("limitFields");


    // Constructors

        //- Construct from Time and dictionary
        limitFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~limitFields() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the limitFields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "limitFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
