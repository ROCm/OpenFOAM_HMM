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
    Foam::Field

Description
    Generic templated field type.

SourceFiles
    FieldFunctions.H
    FieldFunctionsM.H
    FieldMapper.H
    FieldI.H
    FieldM.H
    Field.C
    FieldFunctions.C
    FieldFunctionsM.C

\*---------------------------------------------------------------------------*/

#ifndef Field_H
#define Field_H

#include "tmp.H"
#include "direction.H"
#include "labelList.H"
#include "scalarList.H"
#include "FieldBase.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class FieldMapper;
class dictionary;

template<class Type> class Field;
template<class Type> class SubField;

template<class Type>
Ostream& operator<<(Ostream&, const Field<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const tmp<Field<Type>>&);


/*---------------------------------------------------------------------------*\
                           Class Field Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Field
:
    public FieldBase,
    public List<Type>
{
public:

    //- Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    //- Declare type of subField
    typedef SubField<Type> subField;


    // Static Member Functions

        //- Return nullObject reference Field
        inline static const Field<Type>& null();


    // Constructors

        //- Default construct
        //  For temporary fields that are initialised after construction
        inline constexpr Field() noexcept;

        //- Construct given size
        //  For temporary fields that are initialised after construction
        inline explicit Field(const label len);

        //- Construct given size and initial value
        inline Field(const label len, const Type& val);

        //- Construct given size and initial values of zero
        inline Field(const label len, const Foam::zero);

        //- Construct with length=1, copying the value as the only content
        inline Field(const Foam::one, const Type& val);

        //- Construct with length=1, moving the value as the only content
        inline Field(const Foam::one, Type&& val);

        //- Construct with length=1, initializing content to zero
        inline Field(const Foam::one, const Foam::zero);

        //- Copy construct
        inline Field(const Field<Type>& fld);

        //- Copy construct from UList\<Type\>
        inline explicit Field(const UList<Type>& list);

        //- Copy construct from IndirectList
        template<class Addr>
        inline explicit Field(const IndirectListBase<Type, Addr>& list);

        //- Move construct from Field
        inline Field(Field<Type>&& fld);

        //- Move construct from List
        inline Field(List<Type>&& list);

        //- Move construct from DynamicList
        template<int SizeMin>
        inline Field(DynamicList<Type, SizeMin>&& list);

        //- Construct by 1 to 1 mapping from the given field
        Field
        (
            const UList<Type>& mapF,
            const labelUList& mapAddressing
        );

        //- Construct by 1 to 1 mapping from the given tmp field
        Field
        (
            const tmp<Field<Type>>& tmapF,
            const labelUList& mapAddressing
        );

        //- Construct by interpolative mapping from the given field
        Field
        (
            const UList<Type>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Construct by interpolative mapping from the given tmp field
        Field
        (
            const tmp<Field<Type>>& tmapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Construct by mapping from the given field
        Field
        (
            const UList<Type>& mapF,
            const FieldMapper& map,
            const bool applyFlip = true
        );

        //- Construct by mapping from the given field
        Field
        (
            const UList<Type>& mapF,
            const FieldMapper& map,
            const Type& defaultValue,
            const bool applyFlip = true
        );

        //- Construct by mapping from the given field
        Field
        (
            const UList<Type>& mapF,
            const FieldMapper& map,
            const UList<Type>& defaultValues,
            const bool applyFlip = true
        );

        //- Construct by mapping from the given tmp field
        Field
        (
            const tmp<Field<Type>>& tmapF,
            const FieldMapper& map,
            const bool applyFlip = true
        );

        //- Construct by mapping from the given tmp field.
        //- Uses supplied uniform value for unmapped items
        Field
        (
            const tmp<Field<Type>>& tmapF,
            const FieldMapper& map,
            const Type& defaultValue,
            const bool applyFlip = true
        );

        //- Construct by mapping from the given tmp field.
        //- Uses supplied values for unmapped items
        Field
        (
            const tmp<Field<Type>>& tmapF,
            const FieldMapper& map,
            const UList<Type>& defaultValues,
            const bool applyFlip = true
        );

        //- Copy construct or re-use as specified.
        inline Field(Field<Type>& fld, bool reuse);

        //- Copy or move construct from tmp
        inline Field(const tmp<Field<Type>>& tfld);

        //- Construct from Istream
        inline Field(Istream& is);

        //- Construct from a dictionary entry
        Field(const word& keyword, const dictionary& dict, const label len);

        //- Clone
        inline tmp<Field<Type>> clone() const;

        //- Return a pointer to a new Field created on freestore
        static autoPtr<Field<Type>> New(Istream& is)
        {
            return autoPtr<Field<Type>>(new Field<Type>(is));
        }

        //- Return a pointer to a new calculatedFvPatchFieldField created on
        //- freestore without setting patchField values
        template<class Type2>
        static tmp<Field<Type>> NewCalculatedType(const Field<Type2>& f)
        {
            return tmp<Field<Type>>::New(f.size());
        }


    // Member Functions

        //- 1 to 1 map from the given field
        void map
        (
            const UList<Type>& mapF,
            const labelUList& mapAddressing
        );

        //- 1 to 1 map from the given tmp field
        void map
        (
            const tmp<Field<Type>>& tmapF,
            const labelUList& mapAddressing
        );

        //- Interpolative map from the given field
        void map
        (
            const UList<Type>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Interpolative map from the given tmp field
        void map
        (
            const tmp<Field<Type>>& tmapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Map from the given field
        void map
        (
            const UList<Type>& mapF,
            const FieldMapper& map,
            const bool applyFlip = true
        );

        //- Map from the given tmp field
        void map
        (
            const tmp<Field<Type>>& tmapF,
            const FieldMapper& map,
            const bool applyFlip = true
        );

        //- Map from self
        void autoMap
        (
            const FieldMapper& map,
            const bool applyFlip = true
        );

        //- 1 to 1 reverse-map from the given field
        void rmap
        (
            const UList<Type>& mapF,
            const labelUList& mapAddressing
        );

        //- 1 to 1 reverse-map from the given tmp field
        void rmap
        (
            const tmp<Field<Type>>& tmapF,
            const labelUList& mapAddressing
        );

        //- Interpolative reverse map from the given field
        void rmap
        (
            const UList<Type>& mapF,
            const labelUList& mapAddressing,
            const UList<scalar>& weights
        );

        //- Interpolative reverse map from the given tmp field
        void rmap
        (
            const tmp<Field<Type>>& tmapF,
            const labelUList& mapAddressing,
            const UList<scalar>& weights
        );

        //- Negate this field (negative).
        //  Inverts the state for a bool field.
        void negate();

        //- Return a component field of the field
        tmp<Field<cmptType>> component(const direction) const;

        //- Replace a component field of the field
        void replace(const direction, const UList<cmptType>&);

        //- Replace a component field of the field
        void replace(const direction, const tmp<Field<cmptType>>&);

        //- Replace a component field of the field
        void replace(const direction, const cmptType&);

        template<class VSForm>
        VSForm block(const label start) const;

        //- Return the field transpose (only defined for second rank tensors)
        tmp<Field<Type>> T() const;

        //- Write the field as a dictionary entry
        void writeEntry(const word& keyword, Ostream& os) const;


    // Other Access

        //- Return SubField slice (non-const access) - no range checking
        SubField<Type> slice(const label pos, label len = -1);

        //- Return SubField slice (const access) - no range checking
        const SubField<Type> slice(const label pos, label len = -1) const;

        //- Return SubField slice (non-const access) - with range checking
        SubField<Type> slice(const labelRange& range);

        //- Return SubField slice (const access) - with range checking
        const SubField<Type> slice(const labelRange& range) const;


    // Member Operators

        //- Copy assignment
        void operator=(const Field<Type>&);
        void operator=(const tmp<Field<Type>>&);

        inline void operator=(const UList<Type>& rhs);
        inline void operator=(const SubField<Type>& rhs);

        //- Copy assign from IndirectList
        template<class Addr>
        inline void operator=(const IndirectListBase<Type, Addr>& rhs);

        //- Move assignment
        inline void operator=(Field<Type>&& rhs);
        inline void operator=(List<Type>&& rhs);

        template<int SizeMin>
        inline void operator=(DynamicList<Type, SizeMin>&& rhs);


        //- Value assignment
        inline void operator=(const Type& val);
        inline void operator=(const Foam::zero);

        template<class Form, class Cmpt, direction nCmpt>
        void operator=(const VectorSpace<Form,Cmpt,nCmpt>&);

        void operator+=(const UList<Type>&);
        void operator+=(const tmp<Field<Type>>&);

        void operator-=(const UList<Type>&);
        void operator-=(const tmp<Field<Type>>&);

        void operator*=(const UList<scalar>&);
        void operator*=(const tmp<Field<scalar>>&);

        void operator/=(const UList<scalar>&);
        void operator/=(const tmp<Field<scalar>>&);

        void operator+=(const Type&);
        void operator-=(const Type&);

        void operator*=(const scalar&);
        void operator/=(const scalar&);


    // IOstream operators

        friend Ostream& operator<< <Type>
        (Ostream&, const Field<Type>&);

        friend Ostream& operator<< <Type>
        (Ostream&, const tmp<Field<Type>>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FieldI.H"
#include "FieldFunctions.H"

#ifdef NoRepository
    #include "Field.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
