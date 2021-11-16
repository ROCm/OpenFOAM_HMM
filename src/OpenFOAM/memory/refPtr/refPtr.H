/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::refPtr

Description
    A class for managing references or pointers (no reference counting)

SourceFiles
    refPtrI.H

See also
    Foam::autoPtr
    Foam::tmp

\*---------------------------------------------------------------------------*/

#ifndef refPtr_H
#define refPtr_H

#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class refPtr Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class refPtr
{
    // Private Data

        //- The stored reference type
        enum refType
        {
            PTR,        //!< A managed pointer (not ref-counted)
            CREF,       //!< A const reference to an object
            REF         //!< A non-const reference to an object
        };

        //- The managed pointer or address of the object (reference)
        mutable T* ptr_;

        //- The type (managed pointer | object reference)
        mutable refType type_;


public:

    // STL type definitions

        //- Type of object being managed or referenced
        typedef T element_type;

        //- Pointer to type of object being managed or referenced
        typedef T* pointer;


    //- Null reference counter class
    typedef Foam::refCount::zero refCount;


    // Factory Methods

        //- Construct refPtr of T with forwarding arguments
        //  \param args list of arguments with which an instance of T
        //      will be constructed.
        //
        //  \note Similar to std::make_shared, but the overload for
        //      array types is not disabled.
        template<class... Args>
        inline static refPtr<T> New(Args&&... args);

        //- Construct refPtr from derived type with forwarding arguments
        //  \param args list of arguments with which an instance of U
        //      will be constructed.
        //
        //  \note Similar to New but for derived types
        template<class U, class... Args>
        inline static refPtr<T> NewFrom(Args&&... args);


    // Constructors

        //- Construct with no managed pointer.
        inline constexpr refPtr() noexcept;

        //- Construct with no managed pointer (literal nullptr).
        inline constexpr refPtr(std::nullptr_t) noexcept;

        //- Construct, taking ownership of the pointer.
        inline explicit constexpr refPtr(T* p) noexcept;

        //- Move construct from autoPtr, transferring ownership.
        inline explicit refPtr(autoPtr<T>&& ptr) noexcept;

        //- Move construct from unique_ptr, transferring ownership.
        inline explicit refPtr(std::unique_ptr<T>&& ptr) noexcept;

        //- Construct for a const reference to an object.
        inline constexpr refPtr(const T& obj) noexcept;

        //- Move construct, transferring ownership.
        inline refPtr(refPtr<T>&& rhs) noexcept;

        //- Copy construct (shallow copy)
        inline refPtr(const refPtr<T>& rhs);

        //- Copy/move construct. Optionally reusing pointer.
        inline refPtr(const refPtr<T>& rhs, bool reuse);


    //- Destructor: deletes managed pointer
    inline ~refPtr();


    // Member Functions

        //- The type-name, constructed from type-name of T
        inline static word typeName();


    // Query

        //- True if pointer/reference is non-null.
        bool good() const noexcept { return bool(ptr_); }

        //- If the stored/referenced content is const
        bool is_const() const noexcept { return type_ == CREF; }

        //- True if this is a managed pointer (not a reference)
        bool is_pointer() const noexcept { return type_ == PTR; }

        //- True if this is a non-null managed pointer
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
        //  Fatal for a null managed pointer.
        inline T& constCast() const;

        //- Return a shallow copy as a wrapped reference, preserving the
        //- const/non-const status.
        inline refPtr<T> shallowClone() const noexcept;


    // Edit

        //- Release ownership and return the pointer.
        //- A no-op for reference objects (returns nullptr).
        //  \remark Method naming consistent with std::unique_ptr
        inline T* release() noexcept;

        //- Return managed pointer for reuse, or clone() the object reference.
        //  \remark Method naming consistent with Foam::tmp
        inline T* ptr() const;

        //- If object pointer points to valid object:
        //- delete object and set pointer to nullptr
        inline void clear() const noexcept;


        //- Clear existing and transfer ownership.
        inline void reset(refPtr<T>&& other) noexcept;

        //- Delete managed pointer and set to new given pointer
        inline void reset(T* p = nullptr) noexcept;

        //- Clear existing and transfer ownership from autoPtr.
        void reset(autoPtr<T>&& other) noexcept { reset(other.release()); }

        //- Clear existing and transfer ownership from unique_ptr
        void reset(std::unique_ptr<T>&& other) { reset(other.release()); }


        //- Clear existing and set (const) reference from other
        inline void cref(const refPtr<T>& other) noexcept;

        //- Clear existing and set (const) reference
        inline void cref(const T& obj) noexcept;

        //- Clear existing and set (const) reference to pointer content.
        //  The pointer can be null, which is handled like a clear().
        inline void cref(const T* p) noexcept;


        //- Clear existing and set (non-const) reference
        inline void ref(T& obj) noexcept;

        //- Clear existing and set (non-const) reference to pointer content.
        //  The pointer can be null, which is handled like a clear().
        inline void ref(T* p) noexcept;


        //- Swaps the managed object with other.
        inline void swap(refPtr<T>& other) noexcept;


    // Member Operators

        //- Return const reference to the object.
        //  Fatal if nothing is managed
        inline const T& operator*() const;

        //- Return reference to the managed object.
        //  Fatal if nothing is managed or if the object is const.
        inline T& operator*();

        //- Dereferences (const) pointer to the managed object.
        //  Fatal if nothing is managed.
        inline const T* operator->() const;

        //- Dereferences (non-const) pointer to the managed object.
        //  Fatal if nothing is managed or if the object is const.
        inline T* operator->();

        //- Return const reference to the object - same as cref() method.
        const T& operator()() const { return cref(); }


    // Casting

        //- True if pointer/reference is non-null. Same as good()
        explicit operator bool() const noexcept { return bool(ptr_); }

        //- Cast to underlying data type, using the cref() method.
        operator const T&() const { return cref(); }


    // Assignment

        //- Transfer ownership of the managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline void operator=(const refPtr<T>& other);

        //- Clear existing and transfer ownership.
        inline void operator=(refPtr<T>&& other) noexcept;

        //- Take ownership of the pointer.
        //  Fatal for a null pointer
        inline void operator=(T* p);

        //- Reset via assignment from literal nullptr
        inline void operator=(std::nullptr_t) noexcept;

        //- Conversion to tmp, releases pointer or shallow-copies reference
        inline operator tmp<T>();


    // Housekeeping

        //- Identical to good(), or bool operator
        bool valid() const noexcept { return bool(ptr_); }

        //- Deprecated(2020-07) True if a null managed pointer
        //
        //  \deprecated(2020-07) - use bool operator
        FOAM_DEPRECATED_FOR(2020-07, "bool operator")
        bool empty() const noexcept { return !ptr_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Functions

//- Specialized Swap algorithm for refPtr.
//  Swaps the pointers and types of lhs and rhs. Calls \c lhs.swap(rhs)
template<class T>
void Swap(refPtr<T>& lhs, refPtr<T>& rhs)
{
    lhs.swap(rhs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "refPtrI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
