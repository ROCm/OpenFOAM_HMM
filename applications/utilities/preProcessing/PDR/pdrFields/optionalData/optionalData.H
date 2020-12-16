/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::optionalData

Description
    A simplified version of std::optional (c++17), with much simpler
    construction semantics.

\*---------------------------------------------------------------------------*/

#ifndef optionalData_H
#define optionalData_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class optionalData Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class optionalData
{
    // Private Data

        //- Validity of the value
        bool good_;

        //- The value
        T value_;


public:

    // Generated Methods

        //- Copy construct
        optionalData(const optionalData<T>&) = default;

        //- Move construct
        optionalData(optionalData<T>&&) = default;

        //- Copy assignment
        optionalData<T>& operator=(const optionalData<T>&) = default;

        //- Move assignment
        optionalData<T>& operator=(optionalData<T>&&) = default;


    // Constructors

        //- Default construct
        optionalData()
        :
            good_(false),
            value_()
        {}

        //- Copy construct from value
        optionalData(const T& val)
        :
            good_(true),
            value_(val)
        {}

        //- Move construct from value
        optionalData(T&& val)
        :
            good_(true),
            value_(std::move(val))
        {}


    // Member Functions

        //- True if it has a value
        bool has_value() const noexcept
        {
            return good_;
        }

        //- Access to the value
        T& value() noexcept
        {
            return value_;
        }

        //- Access to the value
        const T& value() const noexcept
        {
            return value_;
        }

        //- Return value or default
        const T& value_or(const T& deflt) const
        {
            return good_ ? value_ : deflt;
        }


    // Member Operators

        //- True if it has a value
        explicit operator bool() const noexcept
        {
            return good_;
        }

        //- Access the value
        const T& operator*() const noexcept
        {
            return value_;
        }

        //- Access the value
        T& operator*() noexcept
        {
            return value_;
        }

        //- Copy assignment from value
        void operator=(const T& val)
        {
            good_ = true;
            value_ = val;
        }

        //- Move assignment from value
        void operator=(T&& val)
        {
            good_ = true;
            value_ = std::move(val);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
