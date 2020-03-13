/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::fieldValues::fieldValueDelta

Group
    grpFieldFunctionObjects

Description
    Computes a selected operation between two \c fieldValue function objects.

    The operation is applied to all results of each \c fieldValue object.
    Accordingly, each object must generate the same number and type of results.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    fieldValueDelta1
    {
        // Mandatory entries (unmodifiable)
        type            fieldValueDelta;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        operation       subtract;
        region1
        {
            ...
        }
        region2
        {
            ...
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: fieldValueDelta         | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      operation    | Operation type to apply to values  | word |  yes  | -
      region1      | Region1 properties                 | dict |  yes  | -
      region2      | Region2 properties                 | dict |  yes  | -
    \endtable

    Options for the \c operation entry:
    \plaintable
       add           | add
       subtract      | subtract
       min           | minimum
       max           | maximum
       average       | average
    \endplaintable

    The inherited entries are elaborated in:
     - \link fieldValue.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldValue
    - ExtendedCodeGuide::functionObjects::field::fieldValueDelta

SourceFiles
    fieldValueDelta.C
    fieldValueDeltaTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldValueDelta_H
#define functionObjects_fieldValueDelta_H

#include "stateFunctionObject.H"
#include "writeFile.H"
#include "fieldValue.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                       Class fieldValueDelta Declaration
\*---------------------------------------------------------------------------*/

class fieldValueDelta
:
    public fvMeshFunctionObject,
    public writeFile
{
public:
        //- Operation type enumeration
        enum operationType
        {
            opAdd,              //!< Add
            opSubtract,         //!< Subtract
            opMin,              //!< Minimum
            opMax,              //!< Maximum
            opAverage           //!< Average
        };

        //- Operation type names
        static const Enum<operationType> operationTypeNames_;


private:

    // Private Data

        //- Operation to apply to values
        operationType operation_;

        //- Field value region object 1
        autoPtr<fieldValue> region1Ptr_;

        //- Field value region object 2
        autoPtr<fieldValue> region2Ptr_;


    // Private Member Functions

        //- Templated function to apply the operation
        template<class Type>
        void applyOperation
        (
            const word& resultType,
            const word& name1,
            const word& name2,
            const word& entryName1,
            const word& entryName2,
            bool& found
        );


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Run-time type information
    TypeName("fieldValueDelta");


    // Constructors

        //- Construct from Time and dictionary
        fieldValueDelta
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        fieldValueDelta(const fieldValueDelta&) = delete;

        //- No copy assignment
        void operator=(const fieldValueDelta&) = delete;


    //- Destructor
    virtual ~fieldValueDelta() = default;


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldValueDeltaTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
