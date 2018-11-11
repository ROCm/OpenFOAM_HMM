/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::autoPtr

Description
    Pointer management similar to std::unique_ptr, with some additional
    methods and type checking.

Note
    Parts of the interface now mirror std::unique_ptr, but since it pre-dates
    both C++11 and std::unique_ptr, it has some additional idiosyncrasies.
    The const-reference constructors and assignment operators
    actually use move semantics to allow their participation in
    default constructible, default assignable classes.

SourceFiles
    autoPtrI.H

\*---------------------------------------------------------------------------*/

#ifndef autoPtr_H
#define autoPtr_H

// Transitional features/misfeatures:
#define Foam_autoPtr_copyConstruct
#define Foam_autoPtr_copyAssign
#define Foam_autoPtr_castOperator

#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autoPtr Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class autoPtr
{
    //- Pointer to managed object
    T* ptr_;


public:

    // STL type definitions

        //- Type of object being managed
        typedef T element_type;

        //- Pointer to type of object being managed
        typedef T* pointer;


    // Factory Methods

        //- Construct autoPtr of T with forwarding arguments
        //  \param args list of arguments with which an instance of T
        //      will be constructed.
        //
        //  \note Similar to std::make_unique, but the overload for
        //      array types is not disabled.
        template<class... Args>
        inline static autoPtr<T> New(Args&&... args);

        //- Construct autoPtr from derived type with forwarding arguments
        //  \param args list of arguments with which an instance of U
        //      will be constructed.
        //
        //  \note Similar to New but for derived types.
        //  In the future check for is_convertible on the pointer types
        template<class U, class... Args>
        inline static autoPtr<T> NewFrom(Args&&... args);


    // Constructors

        //- Construct with no managed object.
        inline constexpr autoPtr() noexcept;

        //- Construct with no managed object (literal nullptr).
        inline constexpr autoPtr(std::nullptr_t) noexcept;

        //- Construct, taking ownership of the pointer.
        inline explicit autoPtr(T* p) noexcept;

        //- Move construct, transferring ownership.
        inline autoPtr(autoPtr<T>&& ap) noexcept;

        //- Move construct, transferring ownership from derived type.
        //  U must be derivable from T
        //  \note In the future check for is_convertible on the pointer types
        template<class U>
        inline explicit autoPtr(autoPtr<U>&& ap);

        //- A move construct disguised as a copy construct (transfers ownership)
        //  \remark This is a non-standard definition, and should ideally be
        //      marked as deleted - pending cleanup of code currently relying
        //      on this behaviour.
        #ifdef Foam_autoPtr_copyConstruct
        inline autoPtr(const autoPtr<T>& ap) noexcept;
        #else
        autoPtr(const autoPtr<T>& ap) = delete;
        #endif


    //- Destructs the managed object if such is present
    inline ~autoPtr() noexcept;


    // Member Functions

    // Check

        //- True if the managed pointer is null
        inline bool empty() const noexcept;

        //- True if the managed pointer is non-null
        inline bool valid() const noexcept;


    // Access

        //- Return pointer to managed object without nullptr checking.
        //  Pointer remains under autoPtr management.
        inline T* get() noexcept;

        //- Return const pointer to managed object without nullptr checking.
        //  Pointer remains under autoPtr management.
        inline const T* get() const noexcept;

        //- Return reference to the managed object without nullptr checking.
        //  When get() == nullptr, additional guards may be required to avoid
        //  inadvertent access to a nullptr.
        inline T& ref();


    // Edit

        //- Return pointer to the managed object and release ownership.
        inline T* release() noexcept;

        //- Return pointer to the managed object and release ownership.
        //- Identical behaviour to release().
        //  \note Provided for method naming consistent with Foam::tmp
        inline T* ptr() noexcept;

        //- Delete managed object and set pointer to nullptr
        inline void clear() noexcept;

        //- Delete managed object and set to new given pointer
        inline void reset(T* p = nullptr) noexcept;

        //- Delete managed object and set to new given pointer
        //  \remark This is a non-standard definition, but may provide better
        //      code documentation than a simple move assign would.
        inline void reset(autoPtr<T>&& other) noexcept;

        //- Delete managed object and set to new given pointer
        //- Identical behaviour to reset().
        //  \note Provided for backward compatibility - the older version
        //      enforced a run-time check (Fatal if pointer was already set)
        //      but this was rarely used.
        inline void set(T* p) noexcept;

        //- Swaps the managed object with other autoPtr.
        inline void swap(autoPtr<T>& other) noexcept;


    // Other

        //- Construct copy by invoking clone on underlying managed object
        //  \param args list of arguments for clone
        template<class... Args>
        inline autoPtr<T> clone(Args&&... args) const;


    // Member Operators

        //- Return reference to the managed object.
        //  Fatal error if no pointer is managed
        inline T& operator*();

        //- Return const reference to the object.
        //  Fatal error if no pointer is managed
        inline const T& operator*() const;

        //- Dereferences (non-const) pointer to the managed object
        //  Fatal error if no pointer is managed
        inline T* operator->();

        //- Dereferences (const) pointer to the managed object
        //  Fatal error if no pointer is managed
        inline const T* operator->() const;

        //- Return reference to the object data.
        //  Fatal error if no pointer is managed
        inline T& operator()();

        //- Return const reference to the object data
        //  Fatal error if no pointer is managed
        inline const T& operator()() const;

        //- Automatic cast conversion to the underlying type reference
        //  Fatal error if no pointer is managed
        #ifdef Foam_autoPtr_castOperator
        inline operator const T&() const { return operator*(); }
        #else
        operator const T&() const = delete;
        #endif

        //- True if the managed pointer is non-null
        explicit inline operator bool() const noexcept;

        //- Transfer object ownership from parameter
        inline void operator=(autoPtr<T>&& ap) noexcept;

        //- Transfer object ownership from parameter
        template<class U>
        inline void operator=(autoPtr<U>&& ap) noexcept;

        #ifdef Foam_autoPtr_copyAssign
        //- A move assignment disguised as a copy assignment
        //  \remark Non-standard definition - should just be movable
        inline void operator=(const autoPtr<T>& ap) noexcept;
        #else
        void operator=(const autoPtr<T>& ap) = delete;
        #endif

        //- Allow reset via assignment from literal nullptr
        inline void operator=(std::nullptr_t) noexcept;


    // Housekeeping

        //- No copy assignment from plain pointer
        //  \deprecated(2018-02) Convenient, but uncontrolled access
        void operator=(T* p) = delete;
};


// Global Functions

//- Specializes the Swap algorithm for autoPtr.
//  Swaps the pointers of lhs and rhs. Calls \c lhs.swap(rhs)
template<class T>
void Swap(autoPtr<T>& lhs, autoPtr<T>& rhs)
{
    lhs.swap(rhs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtrI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
