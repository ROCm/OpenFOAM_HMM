/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::DimensionedField

Description
    Field with dimensions and associated with geometry type GeoMesh which is
    used to size the field and a reference to it is maintained.

SourceFiles
    DimensionedFieldI.H
    DimensionedField.C
    DimensionedFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedField_H
#define DimensionedField_H

#include "regIOobject.H"
#include "Field.H"
#include "dimensionedType.H"
#include "orientedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class Type, class GeoMesh> class DimensionedField;

template<class Type, class GeoMesh>
Ostream& operator<<
(
    Ostream& os,
    const DimensionedField<Type, GeoMesh>& df
);

template<class Type, class GeoMesh>
Ostream& operator<<
(
    Ostream& os,
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
);


/*---------------------------------------------------------------------------*\
                      Class DimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class DimensionedField
:
    public regIOobject,
    public Field<Type>
{
public:

    // Public Typedefs

        //- Type of mesh on which this DimensionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of the field from which this DimensionedField is derived
        typedef Field<Type> FieldType;

        //- Component type of the elements of the field
        typedef typename Field<Type>::cmptType cmptType;


private:

    // Private data

        //- Reference to mesh
        const Mesh& mesh_;

        //- Dimension set for this field
        dimensionSet dimensions_;

        //- Oriented flag
        orientedType oriented_;


    // Private Member Functions

        //- Assert that non-zero field size == mesh size
        void checkFieldSize() const;

        void readIfPresent(const word& fieldDictEntry = "value");


public:

    //- Runtime type information
    TypeName("DimensionedField");


    // Static Member Functions

        //- Return a NullObjectRef DimensionedField
        inline static const DimensionedField<Type, GeoMesh>& null();


    // Constructors

        //- Construct from components
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const Field<Type>& field
        );

        //- Construct from components, transferring initial field content
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            Field<Type>&& field
        );

        //- Construct from components, transferring initial field content
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            List<Type>&& field
        );

        //- Construct from components, copy or transfer tmp content
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const tmp<Field<Type>>& tfield
        );

        //- Construct from components, setting the initial size and assigning
        //- the dimensions, but not initialising any field values.
        //  Used for temporary fields which are initialised after construction
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const bool checkIOFlags = true
        );

        //- Construct from components, setting the initial size and assigning
        //- both dimensions and values.
        //  The internal name for the dimensioned\<Type\> has no influence.
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const bool checkIOFlags = true
        );

        //- Construct from Istream.
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const word& fieldDictEntry = "value"
        );

        //- Construct from dictionary
        DimensionedField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dictionary& fieldDict,
            const word& fieldDictEntry = "value"
        );

        //- Copy construct
        DimensionedField(const DimensionedField<Type, GeoMesh>& df);

        //- Move construct
        DimensionedField(DimensionedField<Type, GeoMesh>&& df);

        //- Copy construct or reuse (move) as specified.
        DimensionedField(DimensionedField<Type, GeoMesh>& df, bool reuse);

        //- Construct from tmp\<DimensionedField\> deleting argument
        DimensionedField
        (
            const tmp<DimensionedField<Type, GeoMesh>>& tdf
        );

        //- Copy construct, resetting IO parameters
        DimensionedField
        (
            const IOobject& io,
            const DimensionedField<Type, GeoMesh>& df
        );

        //- Move construct, resetting IO parameters
        DimensionedField
        (
            const IOobject& io,
            DimensionedField<Type, GeoMesh>&& df
        );

        //- Copy or move construct, resetting IO parameters.
        DimensionedField
        (
            const IOobject& io,
            DimensionedField<Type, GeoMesh>& df,
            bool reuse
        );

        //- Construct from tmp\<DimensionedField\> deleting argument,
        //- resetting IO parameters.
        DimensionedField
        (
            const IOobject& io,
            const tmp<DimensionedField<Type, GeoMesh>>& tdf
        );

        //- Copy construct with a new name
        DimensionedField
        (
            const word& newName,
            const DimensionedField<Type, GeoMesh>& df
        );

        //- Move construct with a new name
        DimensionedField
        (
            const word& newName,
            DimensionedField<Type, GeoMesh>&& df
        );

        //- Copy or move construct, resetting name.
        DimensionedField
        (
            const word& newName,
            DimensionedField<Type, GeoMesh>& df,
            bool reuse
        );

        //- Construct with a new name from tmp\<DimensionedField\>
        DimensionedField
        (
            const word& newName,
            const tmp<DimensionedField<Type, GeoMesh>>& tdf
        );

        //- Clone
        tmp<DimensionedField<Type, GeoMesh>> clone() const;


    //- Destructor
    virtual ~DimensionedField() = default;


    // Member Functions

        void readField
        (
            const dictionary& fieldDict,
            const word& fieldDictEntry = "value"
        );

        //- Return mesh
        inline const Mesh& mesh() const;

        //- Return dimensions
        inline const dimensionSet& dimensions() const;

        //- Return non-const access to dimensions
        inline dimensionSet& dimensions();

        //- Return oriented type
        inline const orientedType& oriented() const noexcept;

        //- Return non-const access to the oriented type
        inline orientedType& oriented() noexcept;

        //- Set the oriented flag
        inline void setOriented(const bool oriented = true) noexcept;

        //- Return field
        inline const Field<Type>& field() const;

        //- Return field
        inline Field<Type>& field();

        //- Return a component field of the field
        tmp<DimensionedField<cmptType, GeoMesh>> component
        (
            const direction d
        ) const;

        //- Replace a component field of the field
        void replace
        (
            const direction d,
            const DimensionedField<cmptType, GeoMesh>& df
        );

        //- Replace a component field of the field
        void replace
        (
            const direction d,
            const tmp<DimensionedField<cmptType, GeoMesh>>& tdf
        );

        //- Return the field transpose (only defined for second rank tensors)
        tmp<DimensionedField<Type, GeoMesh>> T() const;

        //- Calculate and return arithmetic average
        dimensioned<Type> average() const;

        //- Calculate and return weighted average
        dimensioned<Type> weightedAverage
        (
            const DimensionedField<scalar, GeoMesh>& weightField
        ) const;

        //- Calculate and return weighted average
        dimensioned<Type> weightedAverage
        (
            const tmp<DimensionedField<scalar, GeoMesh>>& tweightField
        ) const;


        // Write

            bool writeData(Ostream& os, const word& fieldDictEntry) const;

            bool writeData(Ostream& os) const;


    // Member Operators

        void operator=(const DimensionedField<Type, GeoMesh>& df);
        void operator=(const tmp<DimensionedField<Type, GeoMesh>>& tdf);

        //- Assign dimensions and value.
        void operator=(const dimensioned<Type>& dt);

        void operator+=(const DimensionedField<Type, GeoMesh>& df);
        void operator+=(const tmp<DimensionedField<Type, GeoMesh>>& tdf);

        void operator-=(const DimensionedField<Type, GeoMesh>& df);
        void operator-=(const tmp<DimensionedField<Type, GeoMesh>>& tdf);

        void operator*=(const DimensionedField<scalar, GeoMesh>& df);
        void operator*=(const tmp<DimensionedField<scalar, GeoMesh>>& tdf);

        void operator/=(const DimensionedField<scalar, GeoMesh>& df);
        void operator/=(const tmp<DimensionedField<scalar, GeoMesh>>& tdf);

        void operator+=(const dimensioned<Type>& dt);
        void operator-=(const dimensioned<Type>& dt);

        void operator*=(const dimensioned<scalar>& dt);
        void operator/=(const dimensioned<scalar>& dt);


    // Ostream Operators

        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream& os,
            const DimensionedField<Type, GeoMesh>& df
        );

        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream& os,
            const tmp<DimensionedField<Type, GeoMesh>>& tdf
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldI.H"
#include "DimensionedFieldFunctions.H"

#ifdef NoRepository
    #include "DimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
