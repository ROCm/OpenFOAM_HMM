/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014 OpenFOAM Foundation
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
    member. There is an additional zero second parameter for safe
    casting to List etc.

SourceFiles
    nullObjectI.H
    nullObject.C

\*---------------------------------------------------------------------------*/

#ifndef nullObject_H
#define nullObject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
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
            return data_[0].val;
        }
};


// IOstream Operators

//- Read from Istream consumes no content.
inline Istream& operator>>(Istream& is, NullObject&)
{
    return is;
}

//- Write to Ostream emits no content.
inline Ostream& operator<<(Ostream& os, const NullObject&)
{
    return os;
}


//- Pointer to the unique nullObject
extern const NullObject* nullObjectPtr;


//- Reference to the nullObject of type T
template<class T>
inline const T& NullObjectRef();

//- Pointer to the nullObject of type T
template<class T>
inline const T* NullObjectPtr();


//- True if t is a reference to the nullObject of type T
template<class T>
inline bool isNull(const T& t);

//- True if t is not a reference to the nullObject of type T
template<class T>
inline bool notNull(const T& t);


//- True if t is a pointer to the nullObject of type T
template<class T>
inline bool isNull(const T* t);

//- True if t is not a pointer to the nullObject of type T
template<class T>
inline bool notNull(const T* t);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nullObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
