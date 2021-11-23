/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::exprTypeTraits

Description
    Simple type identifiers for polymorphic expression values.

    The definitions are similar to std::integral_constant in that they
    provide value, value_type (and name).

\*---------------------------------------------------------------------------*/

#ifndef expressions_exprTraits_H
#define expressions_exprTraits_H

// Regular field types
#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "tensor.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

//- An enumeration of known and expected expression value types.
//  Do not rely on the enumeration values for any direct coding.
//
//  \note NONE use used when initializing types, whereas INVALID is used
//  for unsupported types (never as a stored type).
//  This avoids false positives when testing.
//
//  Except NONE and INVALID, the enumerations will mostly not be used
//  directly, but through exprTypeTraits :: value

enum class valueTypeCode : unsigned char
{
    NONE = 0,               //!< No type, or default initialized type
    INVALID,                //!< Invalid/unknown/error type

    // Rank 0 types
    type_bool,              //!< Type is 'bool'
    type_label,             //!< Type is 'label'
    type_scalar,            //!< Type is 'scalar'

    // Rank 1 types
    type_vector,            //!< Type is 'vector'

    // Rank 2 types
    type_sphericalTensor,   //!< Type is 'sphericalTensor'
    type_symmTensor,        //!< Type is 'symmTensor'
    type_tensor             //!< Type is 'tensor'
};


// Global Functions

//- From string to valueTypeCode (if any)
valueTypeCode valueTypeCodeOf(const word& dataTypeName);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions


/*---------------------------------------------------------------------------*\
                       Class exprTypeTraits Declaration
\*---------------------------------------------------------------------------*/

// Generic enumerated traits is INVALID (unknown)
template<class Type>
struct exprTypeTraits
{
    typedef Type value_type;
    static constexpr const char* const name = "";
    static constexpr
        ::Foam::expressions::valueTypeCode value =
        ::Foam::expressions::valueTypeCode::INVALID;
};


#undef  defineExprTypeTraits
#define defineExprTypeTraits(Type, Name)                                      \
    template<>                                                                \
    struct exprTypeTraits<Type>                                               \
    {                                                                         \
        typedef Type value_type;                                              \
        static constexpr const char* const name = #Name;                      \
        static constexpr                                                      \
            ::Foam::expressions::valueTypeCode value =                        \
            ::Foam::expressions::valueTypeCode::type_##Name;                  \
    };


// Define with "name" to match regular pTraits typeName
defineExprTypeTraits(bool, bool);
defineExprTypeTraits(::Foam::label, label);
defineExprTypeTraits(::Foam::scalar, scalar);
defineExprTypeTraits(::Foam::vector, vector);
defineExprTypeTraits(::Foam::tensor, tensor);
defineExprTypeTraits(::Foam::sphericalTensor, sphericalTensor);
defineExprTypeTraits(::Foam::symmTensor, symmTensor);

#undef defineExprTypeTraits

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A word representation of a valueTypeCode. Empty for INVALID
word name(const expressions::valueTypeCode typeCode);


//- A word representation of a valueTypeCode. Empty for INVALID
template<>
struct nameOp<expressions::valueTypeCode>
{
    word operator()(const expressions::valueTypeCode typeCode) const
    {
        return Foam::name(typeCode);
    }
};

// No IOstream Operators for valueTypeCode at the moment (Nov 2021)

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
