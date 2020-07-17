/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::tmpNrc

Description
    A class for managing temporary objects without reference counting.

SourceFiles
    tmpNrcI.H

See also
    Foam::autoPtr
    Foam::tmp

\*---------------------------------------------------------------------------*/

#ifndef tmpNrc_H
#define tmpNrc_H

#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class tmpNrc Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class tmpNrc
{
    // Private Data

        //- Object types
        enum refType
        {
            PTR,    //!< Managing a pointer (not ref-counted)
            CREF    //!< Using (const) reference to an object
        };

        //- The managed pointer or address of the object (reference)
        mutable T* ptr_;

        //- The type (managed pointer | object reference)
        mutable refType type_;


public:

    // STL type definitions

        //- Type of object being managed
        typedef T element_type;

        //- Pointer to type of object being managed
        typedef T* pointer;


    //- Null reference counter class
    typedef Foam::refCount::zero refCount;


    // Factory Methods

        //- Construct tmpNrc of T with forwarding arguments
        //  \param args list of arguments with which an instance of T
        //      will be constructed.
        //
        //  \note Similar to std::make_shared, but the overload for
        //      array types is not disabled.
        template<class... Args>
        inline static tmpNrc<T> New(Args&&... args);

        //- Construct tmpNrc from derived type with forwarding arguments
        //  \param args list of arguments with which an instance of U
        //      will be constructed.
        //
        //  \note Similar to New but for derived types
        template<class U, class... Args>
        inline static tmpNrc<T> NewFrom(Args&&... args);


    // Constructors

        //- Default construct, no managed pointer.
        inline constexpr tmpNrc() noexcept;

        //- Construct with no managed pointer.
        inline constexpr tmpNrc(std::nullptr_t) noexcept;

        //- Construct, taking ownership of the pointer.
        inline explicit tmpNrc(T* p) noexcept;

        //- Construct for a const reference to an object.
        inline tmpNrc(const T& obj) noexcept;

        //- Move construct, transferring ownership.
        inline tmpNrc(tmpNrc<T>&& t) noexcept;

        //- Copy construct
        inline tmpNrc(const tmpNrc<T>& t);

        //- Copy construct. Optionally reusing pointer.
        inline tmpNrc(const tmpNrc<T>& t, bool reuse);


    //- Destructor: deletes managed pointer
    inline ~tmpNrc();


    // Member Functions

    // Query

        //- Deprecated(2020-07) True if a null managed pointer
        //
        //  \deprecated(2020-07) - use bool operator
        FOAM_DEPRECATED_FOR(2020-07, "bool operator")
        bool empty() const noexcept { return !ptr_; }

        //- True for non-null pointer/reference
        bool valid() const noexcept { return ptr_; }

        //- True if this is a managed pointer (not a reference)
        bool isTmp() const noexcept { return type_ == PTR; }

        //- True if this is a non-null managed pointer
        inline bool movable() const noexcept;

        //- Return type-name of the tmp, constructed from type-name of T
        inline word typeName() const;


    // Access

        //- Return pointer without nullptr checking.
        T* get() noexcept { return ptr_; }

        //- Return const pointer without nullptr checking.
        const T* get() const noexcept { return ptr_; }

        //- Return the const object reference or a const reference to the
        //- contents of a non-null managed pointer.
        //  Fatal for a null managed pointer
        inline const T& cref() const;

        //- Return non-const reference to the contents of a non-null
        //- managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline T& ref() const;

        //- Non-const dereference, even if the object is const.
        //  This is similar to ref(), but applies a const_cast to access
        //  const objects.
        //  Fatal for a null managed pointer.
        inline T& constCast() const;


    // Edit

        //- Return managed pointer for reuse, or clone() the const object
        //- reference.
        inline T* ptr() const;

        //- If object pointer points to valid object:
        //- delete object and set pointer to nullptr
        inline void clear() const noexcept;

        //- Delete managed temporary object and set to new given pointer
        inline void reset(T* p = nullptr) noexcept;

        //- Clear existing and transfer ownership.
        inline void reset(tmpNrc<T>&& other) noexcept;

        //- Delete managed temporary object and set to const reference
        inline void cref(const T& obj) noexcept;

        //- Swaps the managed object with other.
        inline void swap(tmpNrc<T>& other) noexcept;


    // Member Operators

        //- Return const reference to the object.
        //  Identical to cref() method.
        inline const T& operator()() const;

        //- Cast to underlying data type, using the cref() method.
        inline operator const T&() const;

        //- Dereferences (const) pointer to the managed object.
        //  Fatal for a null managed pointer.
        inline const T* operator->() const;

        //- Dereferences (non-const) pointer to the managed object.
        //  Fatal for a null managed pointer or if the object is const.
        inline T* operator->();

        //- Non-null pointer/reference : valid()
        explicit operator bool() const noexcept { return ptr_; }

        //- Transfer ownership of the managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline void operator=(const tmpNrc<T>& t);

        //- Clear existing and transfer ownership.
        inline void operator=(tmpNrc<T>&& other) noexcept;

        //- Take ownership of the pointer.
        //  Fatal for a null pointer
        inline void operator=(T* p);

        //- Reset via assignment from literal nullptr
        inline void operator=(std::nullptr_t) noexcept;

        //- Conversion to tmp - releases pointer or copies reference
        inline operator tmp<T>();
};


// Global Functions

//- Specializes the Swap algorithm for tmpNrc.
//  Swaps the pointers and types of lhs and rhs. Calls \c lhs.swap(rhs)
template<class T>
void Swap(tmpNrc<T>& lhs, tmpNrc<T>& rhs)
{
    lhs.swap(rhs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tmpNrcI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
