/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::explicitSetValue

Description
    Explicit set values on fields.

    Sources described by:

    explicitSetValueCoeffs
    {
        fieldData            // field data - usage for multiple fields
        {
            k   30.7;
            epsilon  1.5;
        }
    }

SourceFiles
    explicitSetValue.C

\*---------------------------------------------------------------------------*/

#ifndef explicitSetValue_H
#define explicitSetValue_H

#include "cellSet.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "basicSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class explicitSetValue Declaration
\*---------------------------------------------------------------------------*/

class explicitSetValue
:
    public basicSource
{
    // Private data

        //- List of field types
        HashTable<scalar> scalarFields_;
        HashTable<vector> vectorFields_;

        //- Set value to field
        template<class Type>
        void setFieldValue(fvMatrix<Type>&, const Type&) const;

        //- Add field names and values to field table for types.
        template<class Type>
        void addField
        (
            HashTable<Type>& fields,
            const wordList& fieldTypes,
            const wordList& fieldNames,
            const dictionary& dict
        );


protected:

    // Protected functions

        //- Set the local field data
        void setFieldData(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("explicitSetValue");


    // Constructors

        //- Construct from components
        explicitSetValue
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<explicitSetValue> clone() const
        {
            notImplemented("autoPtr<explicitSetValue> clone() const");
            return autoPtr<explicitSetValue>(NULL);
        }


    // Member Functions

        // Check

            //- Return index of field name if found in fieldNames list
            virtual label applyToField(const word& fieldName) const;


        // Evaluation

            //- Set value on vector field
            virtual void setValue(fvMatrix<vector>& eqn, const label fieldI);

            //- Set value on scalar field
            virtual void setValue(fvMatrix<scalar>& eqn, const label fieldI);


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "explicitSetValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
