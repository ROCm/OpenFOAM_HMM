/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::Xfer

Description
    A simple container for copying or transferring objects of type \<T\>.

    The wrapped object of type \<T\> must implement a transfer() method and
    an operator=() copy method.

\deprecated(2018-03) This class is an artifact from pre-C++11 code.

Note
    This class is an artifact from pre-C++11 code, where it was used
    as a workaround for the lack of movable references (rvalue).
    The interfaces previously using Xfer to reclaim memory now use
    movable references directly or, in the rare case, an autoPtr.

See also
    xferCopy, xferCopyTo, xferMove, xferMoveTo

\*---------------------------------------------------------------------------*/

#ifndef Xfer_H
#define Xfer_H

#include "nullObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Xfer Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Xfer
{
    // Private data

        //- Pointer to underlying datatype
        mutable T* ptr_;


public:

    // STL type definition similar to std::shared_ptr

        //- Type of object being managed
        typedef T element_type;


    // Constructors

        //- Construct with no managed object.
        inline constexpr Xfer() noexcept : ptr_(nullptr)
        {}

        //- Store object pointer and manage its deletion
        //  Can also be used later to transfer by assignment
        inline explicit Xfer(T* p) noexcept : ptr_(p)
        {}

        //- Move construct, transferring ownership
        inline Xfer(const Xfer<T>& xf) noexcept : ptr_(xf.ptr_)
        {
            xf.ptr_ = nullptr;
        }


    //- Destructor
    inline ~Xfer() noexcept
    {
        delete ptr_;
        ptr_ = nullptr;
    }


    // Static Member Functions

        //- Return a null object reference
        inline static const Xfer<T>& null()
        {
            return NullObjectRef<Xfer<T>>();
        }


    // Member Functions

        //- Test for valid pointer
        inline bool valid() const noexcept { return ptr_; }

        //- Pointer to the underlying object.
        inline T* get() const noexcept { return ptr_; }

        //- Rvalue reference to the underlying object.
        inline T&& rvalue() const { return std::move(*ptr_); }

        //- Swaps the managed objects
        inline void swap(Xfer<T>& other) noexcept
        {
            T* p = ptr_;
            ptr_ = other.ptr_;
            other.ptr_ = p;
        }


    // Member Operators

        //- Reference to the underlying object
        inline T& operator*() const { return *ptr_; }

        //- Pointer to the underlying object
        inline T* operator->() const { return ptr_; }

        //- Reference to the underlying object
        inline T& operator()() const { return *ptr_; }

        //- Move assignment, transferring ownership
        inline void operator=(const Xfer<T>& rhs) noexcept
        {
            if (this != &rhs)
            {
                // clear and swap
                delete ptr_;
                ptr_ = rhs.ptr_;
                rhs.ptr_ = nullptr;
            }
        }
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //


//- Copy construct contents of the \a obj
//  \deprecated - Xfer should not be used in new code
//  \sa xferCopyTo, xferMove, xferMoveTo and Foam::Xfer
template<class T>
inline Xfer<T> xferCopy(const T& obj)
{
    return Foam::Xfer<T>(new T(obj));
}


//- Transfer construct contents of the \a obj
//  \deprecated - Xfer should not be used in new code
//
//  \sa xferCopy, xferCopyTo, xferMoveTo and Foam::Xfer
template<class T>
inline Xfer<T> xferMove(T& obj)
{
    T* ptr = new T;
    ptr->transfer(obj);
    return Foam::Xfer<T>(ptr);
}


//- Copy construct contents of the \a obj from dissimilar type
//  \deprecated - Xfer should not be used in new code
//
//  \sa xferCopy, xferMove, xferMoveTo and Foam::Xfer
template<class T, class From>
inline Xfer<T> xferCopyTo(const From& obj)
{
    return Foam::Xfer<T>(new T(obj));
}


//- Transfer construct contents of the \a obj from dissimilar type
//  \deprecated - Xfer should not be used in new code
//
//  \par Example Use
//  \code
//      DynamicList<label> dynLst;
//      ...
//      labelList plainLst( xferMoveTo<labelList>(dynLst) );
//  \endcode
//
//  \sa xferCopy, xferCopyTo, xferMove and Foam::Xfer
template<class T, class From>
inline Xfer<T> xferMoveTo(From& obj)
{
    T* ptr = new T;
    ptr->transfer(obj);
    return Foam::Xfer<T>(ptr);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
