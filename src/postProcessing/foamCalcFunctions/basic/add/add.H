/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::calcTypes::add

Description
    Adds and field or value to base field.

    New field name specified by -resultName option, or automatically as:
        <baseFieldName>_plus_<addFieldName>
        <baseFieldName>_plus_value

    Example usage:
        add p -value 100000 -resultName pAbs
        add U -field U0

SourceFiles
    add.C
    writeAddField.C
    writeAddValue.C

\*---------------------------------------------------------------------------*/

#ifndef add_H
#define add_H

#include "calcType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace calcTypes
{

/*---------------------------------------------------------------------------*\
                          Class add Declaration
\*---------------------------------------------------------------------------*/

class add
:
    public calcType
{
public:

    enum calcTypes
    {
        FIELD,
        VALUE
    };


private:

    // Private data

        //- Name of base field (to add to)
        word baseFieldName_;

        //- Calc type as given by enumerations above
        calcTypes calcType_;

        //- Name of field to add
        word addFieldName_;

        //- String representation of value to add
        string addValueStr_;

        //- Name of result field
        word resultName_;


    // Private Member Functions

        // Output

            //- Calc and output field additions
            void writeAddFields
            (
                const Time& runTime,
                const fvMesh& mesh,
                const IOobject& baseFieldHeader
            );

            //- Calc and output field and value additions
            void writeAddValues
            (
                const Time& runTime,
                const fvMesh& mesh,
                const IOobject& baseFieldHeader
            );


        //- Disallow default bitwise copy construct
        add(const add&);

        //- Disallow default bitwise assignment
        void operator=(const add&);


protected:

    // Member Functions

        // Calculation routines

            //- Initialise - typically setting static variables,
            //  e.g. command line arguments
            virtual void init();

            //- Pre-time loop calculations
            virtual void preCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Time loop calculations
            virtual void calc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );


        // I-O

            //- Write add field
            template<class Type>
            void writeAddField
            (
                const IOobject& baseHeader,
                const IOobject& addHeader,
                const fvMesh& mesh,
                bool& processed
            );

            //- Write add value
            template<class Type>
            void writeAddValue
            (
                const IOobject& baseHeader,
                const string& valueStr,
                const fvMesh& mesh,
                bool& processed
            );


public:

    //- Runtime type information
    TypeName("add");


    // Constructors

        //- Construct null
        add();


    // Destructor

        virtual ~add();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "writeAddField.C"
#   include "writeAddValue.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
