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

    Since it is decided upon construction of the Xfer object whether the
    parameter is to be copied or transferred, the contents of the resulting
    Xfer object can be transferred unconditionally. This greatly simplifies
    defining constructors or methods in other classes with mixed
    transfer/copy semantics without requiring 2^N different versions.

    When transferring between dissimilar types, the xferCopyTo() and
    xferMoveTo() functions can prove useful. An example is transferring
    from a DynamicList to a List.

    \code
        List<label> list1;
        DynamicList<label> list2;
        ...

        SomeClass obj(xferCopy(list1), xferMoveTo<labelList>(list1));
    \endcode

See also
    xferCopy, xferCopyTo, xferMove, xferMoveTo

SourceFiles
    XferI.H

\*---------------------------------------------------------------------------*/

#ifndef Xfer_H
#define Xfer_H

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

        //- Store object pointer and manage its deletion
        //  Can also be used later to transfer by assignment
        inline explicit Xfer(T* p = nullptr);

        //- Construct, transferring its contents
        inline Xfer(const Xfer<T>& xf);


    //- Destructor
    inline ~Xfer();


    // Static Member Functions

        //- Return a null object reference
        inline static const Xfer<T>& null();

        //- Construct new Xfer container with forwarding arguments to T
        template<class... Args>
        inline static Xfer<T> New(Args&&... args);


    // Member Functions

        //- Test for valid pointer and not a nullObject reference.
        inline bool valid() const;

        //- Pointer to the underlying object.
        inline T* get() const;

        //- Rvalue reference to the underlying object.
        inline T&& rvalue() const;

        //- Swaps the managed objects
        inline void swap(Xfer<T>& other) noexcept;


    // Member Operators

        //- Transfer contents into the object
        inline void operator=(const Xfer<T>& xf);

        //- Reference to the underlying object
        inline T& operator*() const;

        //- Pointer to the underlying object
        inline T* operator->() const;

        //- Reference to the underlying object
        inline T& operator()() const;

};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //


//- Copy construct contents of the \a obj
//
//  \sa xferCopyTo, xferMove, xferMoveTo and Foam::Xfer
template<class T>
inline Xfer<T> xferCopy(const T& obj);

//- Transfer construct contents of the \a obj
//
//  \sa xferCopy, xferCopyTo, xferMoveTo and Foam::Xfer
template<class T>
inline Xfer<T> xferMove(T& obj);


//- Copy construct contents of the \a obj from dissimilar type
//
//  \sa xferCopy, xferMove, xferMoveTo and Foam::Xfer
template<class T, class From>
inline Xfer<T> xferCopyTo(const From& obj);


//- Transfer construct contents of the \a obj from dissimilar type
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
inline Xfer<T> xferMoveTo(From& obj);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "XferI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
