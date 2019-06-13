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

#include "tmpNrc.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- A const Field wrapper with possible data conversion
template<class Type, class InputType, template<class> class Container = Field>
class ConstPrecisionAdaptor
:
    public tmpNrc<Container<Type>>
{
    // Private Member Functions

        //- Copy in field
        void copyInput(const Container<InputType>& input)
        {
            this->clear();

            Container<Type>* p = new Container<Type>(input.size());
            this->reset(p);
            std::copy(input.cbegin(), input.cend(), p->begin());
        }

        //- Construct from tmp Field, copy/move as required
        void moveInput(tmp<Container<InputType>>& input)
        {
            if (std::is_same<Type, InputType>::value)
            {
                auto& tinput = reinterpret_cast<tmp<Container<Type>>&>(input);

                if (tinput.isTmp())
                {
                    // Reset to tmp
                    this->reset(tinput.ptr());
                }
                else
                {
                    this->cref(tinput.cref());
                }
            }
            else
            {
                this->copyInput(input());
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
            tmpNrc<Container<Type>>()
        {
            if (std::is_same<Type, InputType>::value)
            {
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
            tmpNrc<Container<Type>>()
        {
            this->moveInput(input);
        }


        //- Construct from tmp Container, copy/move as required
        ConstPrecisionAdaptor(const tmp<Container<InputType>>& input)
        :
            tmpNrc<Container<Type>>()
        {
            this->moveInput(const_cast<tmp<Container<InputType>>&>(input));
        }


    // Member Functions

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


//- A Field wrapper with possible data conversion
template<class Type, class InputType, template<class> class Container = Field>
class PrecisionAdaptor
:
    public tmpNrc<Container<Type>>
{
    // Private Data

        //- Reference to underlying field
        Container<InputType>& ref_;


    // Private Member Functions

        //- Copy in field
        void copyInput(const Container<InputType>& input)
        {
            Container<Type>* p = new Container<Type>(input.size());
            this->reset(p);
            std::copy(input.cbegin(), input.cend(), p->begin());
        }


public:

    //- The adapted field type
    typedef Container<Type> FieldType;


    // Constructors

        //- Construct from Container<InputType>, copying on input as required
        PrecisionAdaptor(Container<InputType>& input)
        :
            tmpNrc<Container<Type>>(),
            ref_(input)
        {
            if (std::is_same<Type, InputType>::value)
            {
                this->cref(reinterpret_cast<const FieldType&>(input));
            }
            else
            {
                this->copyInput(input);
            }
        }


    //- Destructor, copying on destroy
    ~PrecisionAdaptor()
    {
        if (this->isTmp())
        {
            const FieldType& store = this->cref();
            ref_.resize(store.size());
            std::copy(store.cbegin(), store.cend(), ref_.begin());
        }
    }


    // Member Functions

        //- Allow modification without const-ref check
        FieldType& ref()
        {
            return this->constCast();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
