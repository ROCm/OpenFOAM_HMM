/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Conversion adaptor for Field/List that either wrap the input as a
    reference, or creates a temporary pointer and copies the values
    on construction/destruction.

    This provides, for example, automatic conversion between types
    for linear solvers able to run mixed precision.

\*---------------------------------------------------------------------------*/

#ifndef PrecisionAdaptor_H
#define PrecisionAdaptor_H

#include "refPtr.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- A const Field/List wrapper with possible data conversion
template<class Type, class InputType, template<class> class Container = Field>
class ConstPrecisionAdaptor
:
    public refPtr<Container<Type>>
{
    // Private Member Functions

        //- Copy in field
        void copyInput(const Container<InputType>& input)
        {
            this->reset(new Container<Type>(input.size()));
            std::copy(input.cbegin(), input.cend(), this->ref().begin());
        }

        //- Construct from tmp Field, copy/move as required
        void moveInput(tmp<Container<InputType>>& input)
        {
            if (std::is_same<Type, InputType>::value)
            {
                auto& tinput = reinterpret_cast<tmp<Container<Type>>&>(input);

                if (tinput.is_pointer())
                {
                    // Acquire control of the managed pointer
                    this->reset(tinput.ptr());
                }
                else
                {
                    // Use const reference
                    this->cref(tinput.cref());
                }
            }
            else
            {
                this->copyInput(input.cref());
            }
            input.clear();
        }


public:

    //- The adapted field type
    typedef Container<Type> FieldType;


    // Constructors

        //- Construct from Container<InputType>, copying on input as required
        ConstPrecisionAdaptor(const Container<InputType>& input)
        :
            refPtr<Container<Type>>()
        {
            if (std::is_same<Type, InputType>::value)
            {
                // Use const reference directly
                this->cref(reinterpret_cast<const FieldType&>(input));
            }
            else
            {
                this->copyInput(input);
            }
        }


        //- Construct from tmp Container, copy/move as required
        ConstPrecisionAdaptor(tmp<Container<InputType>>&& input)
        :
            refPtr<Container<Type>>()
        {
            this->moveInput(input);
        }


        //- Construct from tmp Container, copy/move as required
        ConstPrecisionAdaptor(const tmp<Container<InputType>>& input)
        :
            refPtr<Container<Type>>()
        {
            this->moveInput(const_cast<tmp<Container<InputType>>&>(input));
        }


    // Member Functions

        // May need in the future: using refPtr<Container<Type>>::get;

        //- Return the field
        static const Container<Type>& get
        (
            const Container<InputType>& input,
            Container<Type>& dst
        )
        {
            if (std::is_same<Type, InputType>::value)
            {
                return reinterpret_cast<const FieldType&>(input);
            }
            else
            {
                dst.resize(input.size());
                std::copy(input.cbegin(), input.cend(), dst.begin());
                return dst;
            }
        }
};


//- A non-const Field/List wrapper with possible data conversion
template<class Type, class InputType, template<class> class Container = Field>
class PrecisionAdaptor
:
    public refPtr<Container<Type>>
{
    // Private Data

        //- Reference to underlying (input) data
        Container<InputType>& ref_;


    // Private Member Functions

        //- Copy in field
        void copyInput(const Container<InputType>& input, const bool copy)
        {
            this->reset(new Container<Type>(input.size()));
            if (copy)
            {
                std::copy(input.cbegin(), input.cend(), this->ref().begin());
            }
        }


public:

    //- The adapted field type
    typedef Container<Type> FieldType;


    // Constructors

        //- Construct from Container<InputType>, copying on input if required
        PrecisionAdaptor(Container<InputType>& input, const bool copy = true)
        :
            refPtr<Container<Type>>(),
            ref_(input)
        {
            if (std::is_same<Type, InputType>::value)
            {
                // Use non-const reference directly
                this->ref(reinterpret_cast<FieldType&>(ref_));
            }
            else
            {
                this->copyInput(input, copy);
            }
        }


    //- Destructor, copy back on destroy
    ~PrecisionAdaptor()
    {
        if (this->is_pointer())
        {
            const FieldType& store = this->cref();
            ref_.resize(store.size());  // extra safety
            std::copy(store.cbegin(), store.cend(), ref_.begin());
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
