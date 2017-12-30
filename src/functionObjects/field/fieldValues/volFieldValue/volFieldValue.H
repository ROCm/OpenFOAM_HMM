/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::functionObjects::fieldValues::volFieldValue

Group
    grpFieldFunctionObjects

Description
    Provides a 'volRegion' specialization of the fieldValue function object.

    Given a list of user-specified fields and a 'volRegion', a number of
    operations can be performed, such as sums, averages and integrations.

Usage
    Example of function object specification:
    \verbatim
    volFieldValue1
    {
        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");

        log             true;
        writeControl    writeTime;
        writeFields     true;

        regionType      cellZone;
        name            c0;
        operation       volAverage;

        weightField     alpha1;

        fields
        (
            p
            U
        );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description              | Required    | Default value
        type         | Type name: volFieldValue            | yes   |
        log          | Write data to standard output       | no    | no
        writeFields  | Write the region field values       | yes   |
        regionType   | volRegion type: see below           | yes   |
        name         | Name of volRegion if required       | no    |
        operation    | Operation to perform                | yes   |
        weightField  | Name of field to apply weighting    | no    |
        fields       | List of fields to operate on        | yes   |
    \endtable

    Where \c regionType is defined by
    \plaintable
        cellZone     | requires a 'name' entry to specify the cellZone
        all          | all cells
    \endplaintable

    The \c operation is one of:
    \plaintable
       none                 | No operation
       min                  | Minimum
       max                  | Maximum
       sum                  | Sum
       sumMag               | Sum of component magnitudes
       average              | Ensemble average
       volAverage           | Volume weighted average
       volIntegrate         | Volume integral
       CoV                  | Coefficient of variation: standard deviation/mean
       weightedSum          | Weighted sum
       weightedAverage      | Weighted average
       weightedVolAverage   | Weighted volume average
       weightedVolIntegrate | Weighted volume integral
    \endplaintable

See also
    Foam::functionObjects::fieldValues::fieldValue
    Foam::functionObjects::volRegion
    Foam::functionObject

SourceFiles
    volFieldValue.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_volFieldValue_H
#define functionObjects_volFieldValue_H

#include "fieldValue.H"
#include "volRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                        Class volFieldValue Declaration
\*---------------------------------------------------------------------------*/

class volFieldValue
:
    public fieldValue,
    public volRegion
{

public:

    // Public data types

        //- Bitmask values for operation variants
        enum operationVariant
        {
            typeBase = 0,         //!< Base operation
            typeWeighted = 0x100, //!< Operation using weighting
        };

        //- Operation type enumeration
        enum operationType
        {
            // Normal operations

            opNone = 0,             //!< No operation
            opMin,                  //!< Minimum
            opMax,                  //!< Maximum
            opSum,                  //!< Sum
            opSumMag,               //!< Magnitude of sum
            opAverage,              //!< Average
            opVolAverage,           //!< Volume average
            opVolIntegrate,         //!< Volume integral
            opCoV,                  //!< Coefficient of variation

            // Weighted variants

            //! Weighted sum
            opWeightedSum = (opSum | typeWeighted),

            //! Weighted average
            opWeightedAverage = (opAverage | typeWeighted),

            //! Weighted volume average
            opWeightedVolAverage = (opVolAverage | typeWeighted),

            //! Weighted volume integral
            opWeightedVolIntegrate = (opVolIntegrate | typeWeighted),
        };

        //- Operation type names
        static const Enum<operationType> operationTypeNames_;


protected:

    // Protected data

        //- Operation to apply to values
        operationType operation_;

        //- Weight field name - only used for weighted modes
        word weightFieldName_;


    // Protected Member Functions

        //- True if the operation needs the cell-volume
        bool usesVol() const;

        //- True if the operation variant uses a weight-field
        bool usesWeight() const;

        //- True if operation variant uses a weight-field that is available.
        //  Checks for availability on any processor.
        inline bool canWeight(const scalarField& weightField) const;

        //- Initialise, e.g. cell addressing
        void initialise(const dictionary& dict);

        //- Return true if the field name is valid
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Insert field values into values list
        template<class Type>
        tmp<Field<Type>> getFieldValues
        (
            const word& fieldName,
            const bool mustGet = false
        ) const;

        //- Apply the 'operation' to the values
        template<class Type>
        Type processValues
        (
            const Field<Type>& values,
            const scalarField& V,
            const scalarField& weightField
        ) const;

        //- Helper function to output field values
        label writeAll
        (
            const scalarField& V,
            const scalarField& weightField
        );

        //- Templated helper function to output field values
        template<class Type>
        bool writeValues
        (
            const word& fieldName,
            const scalarField& V,
            const scalarField& weightField
        );

        //- Filter a field according to cellIds
        template<class Type>
        tmp<Field<Type>> filterField(const Field<Type>& field) const;


        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

public:

    //- Run-time type information
    TypeName("volFieldValue");


    // Constructors

        //- Construct from name, Time and dictionary
        volFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        volFieldValue
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~volFieldValue();


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volFieldValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
