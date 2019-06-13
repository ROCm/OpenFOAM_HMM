/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef tmpNrc_H
#define tmpNrc_H

#include "refCount.H"
#include "word.H"
#include "tmp.H"
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class tmpNrc Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class tmpNrc
{
    // Private data

        //- Object types
        enum refType
        {
            PTR,    //!< Managing a pointer (not ref-counted)
            CREF    //!< Using a const-reference to an object
        };

        //- The managed pointer or the address of const-reference object
        mutable T* ptr_;

        //- The type (managed pointer | const-reference object)
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

        //- Construct with no managed pointer.
        inline constexpr tmpNrc() noexcept;

        //- Construct with no managed pointer.
        inline constexpr tmpNrc(std::nullptr_t) noexcept;

        //- Construct, taking ownership of the pointer.
        inline explicit tmpNrc(T* p) noexcept;

        //- Construct for a const reference to an object.
        inline tmpNrc(const T& obj) noexcept;

        //- Move construct
        inline tmpNrc(tmpNrc<T>&& t);

        //- Copy construct
        inline tmpNrc(const tmpNrc<T>& t);

        //- Construct copy, transferring content of temporary if required
        inline tmpNrc(const tmpNrc<T>& t, bool reuse);


    //- Destructor: deletes managed pointer
    inline ~tmpNrc();


    // Member Functions

    // Query

        //- True if this is a managed pointer (not a const reference)
        inline bool isTmp() const noexcept;

        //- True if this is a non-null managed pointer
        inline bool empty() const noexcept;

        //- True if this is a non-null managed pointer,
        //- or is a const object reference
        inline bool valid() const noexcept;

        //- True if this is a non-null managed pointer with a unique ref-count
        inline bool movable() const noexcept;

        //- Return type-name of the tmp, constructed from type-name of T
        inline word typeName() const;


    // Access

        //- Return pointer without nullptr checking.
        inline T* get() noexcept;

        //- Return const pointer without nullptr checking.
        inline const T* get() const noexcept;

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

        //- Release ownership of managed temporary object.
        //  After this call no object is managed.
        inline void reset() noexcept;

        //- Delete managed temporary object and set to new given pointer
        inline void reset(T* p) noexcept;

        //- Delete managed temporary object and set to const reference
        inline void cref(const T& obj) noexcept;

        //- Swaps the managed object with other tmpNrc.
        inline void swap(tmpNrc<T>& other) noexcept;


    // Member operators

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

        //- Take ownership of the pointer.
        //  Fatal for a null pointer, or when the pointer is non-unique.
        inline void operator=(T* p);

        //- Transfer ownership of the managed pointer.
        //  Fatal for a null managed pointer or if the object is const.
        inline void operator=(const tmpNrc<T>& t);

        //- Conversion to tmp
        inline operator tmp<T>();
};


// Global Functions

//- Specializes the Swap algorithm for tmp.
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
