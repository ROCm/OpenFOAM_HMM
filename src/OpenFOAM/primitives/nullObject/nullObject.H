/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::NullObject

Description
    Singleton null-object class and instance.

    Its contents occupy enough space to also be reinterpreted
    as another class with a null pointer or zero long for its first
    member, with additional zero parameters for safe casting to List etc.

SourceFiles
    nullObject.C

\*---------------------------------------------------------------------------*/

#ifndef nullObject_H
#define nullObject_H

#include "labelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
class NullObject;

/*---------------------------------------------------------------------------*\
                         Class NullObject Declaration
\*---------------------------------------------------------------------------*/

class NullObject
{
    //- A %union of zero data types
    union zeros
    {
        void* ptr;
        unsigned long val;
    };


    // Private Data

        //- The zero data content
        zeros data_[4];


    // Constructors

        //- Private constructor for singleton only
        //  Could also rely on bit-wise zero initialization for union content
        NullObject()
        :
           data_{{nullptr}, {nullptr}, {nullptr}, {nullptr}}
        {}

        //- No copy construct
        NullObject(const NullObject&) = delete;

        //- No copy assignment
        void operator=(const NullObject&) = delete;


public:

    // Static Data

        //- A unique null object
        static const NullObject nullObject;


    // Member Functions

        //- A nullptr pointer content
        inline const void* pointer() const
        {
            return data_[0].ptr;
        }

        //- Zero valued integer content
        inline unsigned long value() const
        {
            return 0;
        }

        //- No elements
        inline bool empty() const
        {
            return true;
        }

        //- Zero elements
        inline label size() const
        {
            return 0;
        }

        //- No-op method (for HashTable replacement)
        inline const NullObject& toc() const
        {
            return *this;
        }

        //- No-op method (for HashTable replacement)
        inline const NullObject& sortedToc() const
        {
            return *this;
        }


    // Member Operators

        //- Swallow assignment (cf, std::ignore)
        template<class T>
        inline const NullObject& operator=(const T&) const
        {
            return *this;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Globals

//- Pointer to the unique nullObject
extern const NullObject* nullObjectPtr;


// IOstream Operators

//- Read from Istream consumes no content
inline Istream& operator>>(Istream& is, const NullObject&) noexcept
{
    return is;
}

//- Write to Ostream emits no content
inline Ostream& operator<<(Ostream& os, const NullObject&) noexcept
{
    return os;
}


// Global Functions

//- Pointer (of type T) to the nullObject
template<class T>
inline const T* NullObjectPtr()
{
    return reinterpret_cast<const T*>(nullObjectPtr);
}

//- Reference (of type T) to the nullObject
template<class T>
inline const T& NullObjectRef()
{
    return *reinterpret_cast<const T*>(nullObjectPtr);
}


//- True if ptr is a pointer (of type T) to the nullObject
template<class T>
inline bool isNull(const T* ptr)
{
    return ptr == NullObjectPtr<T>();
}

//- True if obj is a reference (of type T) to the nullObject
template<class T>
inline bool isNull(const T& obj)
{
    return &obj == NullObjectPtr<T>();
}


//- True if ptr is not a pointer (of type T) to the nullObject
template<class T>
inline bool notNull(const T* ptr)
{
    return ptr != NullObjectPtr<T>();
}

//- True if obj is not a reference (of type T) to the nullObject
template<class T>
inline bool notNull(const T& obj)
{
    return &obj != NullObjectPtr<T>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
