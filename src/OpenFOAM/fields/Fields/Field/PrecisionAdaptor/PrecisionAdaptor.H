/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::PrecisionAdaptor

Description
    Conversion adaptor for Field that either wraps as a tmp reference
    or creates the necessary tmp and copies the values on construction
    and destruction.
    This provides automatic conversion between (scalar) types for use
    with linear solvers able to run mixed precision.

\*---------------------------------------------------------------------------*/

#ifndef PrecisionAdaptor_H
#define PrecisionAdaptor_H

#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- A const Field wrapper with possible data conversion
template<class Type, class InputType>
class ConstPrecisionAdaptor
:
    public tmp<Field<Type>>
{
public:

    // Constructors

        //- Construct from InputType
        ConstPrecisionAdaptor(const Field<InputType>& input)
        :
            tmp<Field<Type>>()
        {
            if (std::is_same<Type, InputType>::value)
            {
                // Set reference - cast for compiler to handle different types
                this->cref(reinterpret_cast<const Field<Type>&>(input));
            }
            else
            {
                this->reset(new Field<Type>(input.size()));
                std::copy(input.cbegin(), input.cend(), this->ref().begin());
            }
        }


    // Member Functions

        //- Return the field
        static const Field<Type>& get
        (
            const Field<InputType>& input,
            Field<Type>& dst
        )
        {
            if (std::is_same<Type, InputType>::value)
            {
                return reinterpret_cast<const Field<Type>&>(input);
            }
            else
            {
                dst.resize(input.size());
                std::copy(input.cbegin(), input.cend(), dst.begin());
                return dst;
            }
        }
};


//- A Field wrapper with possible data conversion
template<class Type, class InputType>
class PrecisionAdaptor
:
    public tmp<Field<Type>>
{
    // Private Data

        //- Reference to underlying field
        Field<InputType>& ref_;


public:

    // Constructors

        //- Construct from Field<InputType>, copying on input as required
        PrecisionAdaptor(Field<InputType>& input)
        :
            tmp<Field<Type>>(),
            ref_(input)
        {
            if (std::is_same<Type, InputType>::value)
            {
                // Set reference - cast for compiler to handle different types
                this->cref(reinterpret_cast<const Field<Type>&>(input));
            }
            else
            {
                this->reset(new Field<Type>(input.size()));
                std::copy(input.cbegin(), input.cend(), this->ref().begin());
            }
        }


    //- Destructor, copying on destroy
    ~PrecisionAdaptor()
    {
        if (this->isTmp())
        {
            const Field<Type>& store = this->cref();
            ref_.resize(store.size());
            std::copy(store.cbegin(), store.cend(), ref_.begin());
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
