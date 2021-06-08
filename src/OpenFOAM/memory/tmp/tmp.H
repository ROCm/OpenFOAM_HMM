/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::tmp

Description
    A class for managing temporary objects.

    This is a combination of std::shared_ptr (with intrusive ref-counting)
    and a shared_ptr without ref-counting and null deleter.
    This allows the tmp to double as a pointer management and an indirect
    pointer to externally allocated objects.

SourceFiles
    tmpI.H

See also
    Foam::autoPtr
    Foam::refPtr
    Foam::refCount

\*---------------------------------------------------------------------------*/

#ifndef tmp_H
#define tmp_H

#include "refCount.H"
#include "word.H"
#include "stdFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class tmp Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class tmp
{
    // Private Data

        //- Object types
        enum refType
        {
            PTR,    //!< Managing a pointer (ref-counted)
            CREF,   //!< Using (const) reference to an object
            REF     //!< Using (non-const) reference to an object
        };

        //- The managed pointer or address of the object (reference)
        mutable T* ptr_;

        //- The type (managed pointer | object reference)
        mutable refType type_;


    // Private Member Operators

        //- Increment the ref-count for a managed pointer
        //- and check that it is not oversubscribed
        inline void incrCount();


public:

    // STL type definitions

        //- Type of object being managed or referenced
        typedef T element_type;

        //- Pointer to type of object being managed or referenced
        typedef T* pointer;


    //- Reference counter class
    typedef Foam::refCount refCount;


    // Factory Methods

        //- Construct tmp of T with forwarding arguments
        //  \param args list of arguments with which an instance of T
        //      will be constructed.
        //
        //  \note Similar to std::make_shared, but the overload for
        //      array types is not disabled.
        template<class... Args>
        inline static tmp<T> New(Args&&... args);

        //- Construct tmp from derived type with forwarding arguments
        //  \param args list of arguments with which an instance of U
        //      will be constructed.
        //
        //  \note Similar to New but for derived types
        template<class U, class... Args>
        inline static tmp<T> NewFrom(Args&&... args);


    // Constructors

        //- Default construct, no managed pointer.
        inline constexpr tmp() noexcept;

        //- Construct with no managed pointer.
        inline constexpr tmp(std::nullptr_t) noexcept;

        //- Construct, taking ownership of the pointer.
        inline explicit tmp(T* p);

        //- Construct for a const reference to an object.
        inline constexpr tmp(const T& obj) noexcept;

        //- Move construct, transferring ownership.
        //  Does not affect ref-count
        inline tmp(tmp<T>&& rhs) noexcept;

        //- Move construct, transferring ownership.
        //  Does not affect ref-count
        //  \note Non-standard definition - should be non-const
        inline tmp(const tmp<T>&& rhs) noexcept;

        //- Copy construct, incrementing ref-count of managed pointer.
        //  \note Non-standard definition - should be non-const
        inline tmp(const tmp<T>& rhs);

        //- Copy/move construct. Optionally reusing ref-counted pointer.
        inline tmp(const tmp<T>& rhs, bool reuse);


    //- Destructor: deletes managed pointer when the ref-count is 0
    inline ~tmp();


    // Member Functions

        //- The type-name, constructed from type-name of T
        inline static word typeName();


    // Query

        //- Deprecated(2020-07) True if a null managed pointer
        //
        //  \deprecated(2020-07) - use bool operator
        FOAM_DEPRECATED_FOR(2020-07, "bool operator")
        bool empty() const noexcept { return !ptr_; }

        //- True for non-null pointer/reference
        bool valid() const noexcept { return ptr_; }

        //- True if this is a managed pointer (not a reference)
        bool is_pointer() const noexcept { return type_ == PTR; }

        //- Identical to is_pointer()
        bool isTmp() const noexcept { return type_ == PTR; }

        //- True if this is a non-null managed pointer with a unique ref-count
        inline bool movable() const noexcept;


    // Access

        //- Return pointer without nullptr checking.
        T* get() noexcept { return ptr_; }

        //- Return const pointer without nullptr checking.
        const T* get() const noexcept { return ptr_; }

        //- Return const reference to the object or to the contents
        //- of a (non-null) managed pointer.
        //  Fatal for a null managed pointer
        inline const T& cref() const;

        //- Return non-const reference to the contents of a non-null
        //- managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline T& ref() const;

        //- Return non-const reference to the object or to the contents
        //- of a (non-null) managed pointer, with an additional const_cast.
        //  Fatal for a null pointer.
        inline T& constCast() const;


    // Edit

        //- Return managed pointer for reuse, or clone() the object reference.
        inline T* ptr() const;

        //- If object pointer points to valid object:
        //- delete object and set pointer to nullptr
        inline void clear() const noexcept;

        //- Delete managed temporary object and set to new given pointer
        inline void reset(T* p = nullptr) noexcept;

        //- Clear existing and transfer ownership.
        inline void reset(tmp<T>&& other) noexcept;

        //- Clear existing and set (const) reference
        inline void cref(const T& obj) noexcept;

        //- Clear existing and set (const) reference to pointer content.
        //  A null pointer is permitted (treated as a managed pointer).
        inline void cref(const T* p) noexcept;

        //- Clear existing and set to (non-const) reference
        inline void ref(T& obj) noexcept;

        //- Clear existing and set (non-const) reference to pointer content.
        //  A null pointer is permitted (treated as a managed pointer).
        inline void ref(T* p) noexcept;

        //- Swaps the managed object with other.
        inline void swap(tmp<T>& other) noexcept;


    // Member Operators

        //- Dereferences (const) pointer to the managed object.
        //  Fatal for a null managed pointer.
        inline const T* operator->() const;

        //- Dereferences (non-const) pointer to the managed object.
        //  Fatal for a null managed pointer or if the object is const.
        inline T* operator->();

        //- Return const reference to the object - same as cref() method.
        const T& operator()() const { return cref(); }

        //- Cast to underlying data type, using the cref() method.
        operator const T&() const { return cref(); }

        //- Non-null pointer/reference : valid()
        explicit operator bool() const noexcept { return ptr_; }

        //- Transfer ownership of the managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline void operator=(const tmp<T>& other);

        //- Clear existing and transfer ownership.
        inline void operator=(tmp<T>&& other) noexcept;

        //- Take ownership of the pointer.
        //  Fatal for a null pointer, or when the pointer is non-unique.
        inline void operator=(T* p);

        //- Reset via assignment from literal nullptr
        inline void operator=(std::nullptr_t) noexcept;
};


// Global Functions

//- Specialized Swap algorithm for tmp.
//  Swaps the pointers and types of lhs and rhs. Calls \c lhs.swap(rhs)
template<class T>
void Swap(tmp<T>& lhs, tmp<T>& rhs)
{
    lhs.swap(rhs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tmpI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
