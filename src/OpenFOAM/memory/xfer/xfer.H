/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::xfer

Description
    A simple container for copying or transferring objects of type \<T\>.

    The wrapped object of type \<T\> must implement a transfer() method and
    an operator=() copy method.

    Since it is decided upon construction of the xfer object whether the
    parameter is to be copied or transferred, the contents of the resulting
    xfer object can be transferred unconditionally. This greatly simplifies
    defining constructors or methods in other classes with mixed
    transfer/copy semantics without requiring 2^N different versions.

    When transferring between dissimilar types, the xferCopyTo() and
    xferMoveTo() functions can prove useful. An example is transferring
    from a DynamicList to a List.  Since the
    List\<T\>::transfer(List\<T\>&) method could result in some allocated
    memory becoming inaccessible, the xferMoveTo() function can be used to
    invoke the correct List\<T\>::transfer(DynamicList\<T\>&) method.

    @code
        DynamicList<label> dynLst;
        ...
        labelList plainLst( xferMoveTo<labelList>(dynLst) );
    @endcode

SeeAlso
    xferCopy, xferCopyTo, xferMove, xferMoveTo, xferTmp

SourceFiles
    xferI.H

\*---------------------------------------------------------------------------*/

#ifndef xfer_H
#define xfer_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class T> class tmp;

/*---------------------------------------------------------------------------*\
                           Class xfer Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class xfer
{
    // Private data

        //- Pointer to underlying datatype
        mutable T* ptr_;

public:

    // Constructors

        //- Store object pointer and manage its deletion
        //  Can also be used later to transfer by assignment
        inline explicit xfer(T* = 0);

        //- Construct by copying or by transferring the parameter contents
        inline explicit xfer(T&, bool allowTransfer=false);

        //- Construct by copying the parameter contents
        inline explicit xfer(const T&);

        //- Construct by transferring the contents
        inline xfer(const xfer<T>&);

    // Destructor

        inline ~xfer();

    // Member Functions

        //- Return a null object reference
        inline static const xfer<T>& null();

    // Member Operators

        //- Transfer the contents into the object
        inline void operator=(T&);

        //- Transfer the contents into the object
        inline void operator=(const xfer<T>&);

        //- Reference to the underlying datatype
        inline T& operator()() const;

        //- Pointer to the underlying datatype
        inline T* operator->() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/**
 * Construct by copying the contents of the @a arg
 *
 * @sa xferCopyTo, xferMove, xferMoveTo, xferTmp and Foam::xfer
*/
template<class T>
inline xfer<T> xferCopy(const T&);

/**
 * Construct by transferring the contents of the @a arg
 *
 * @sa xferCopy, xferCopyTo, xferMoveTo, xferTmp and Foam::xfer
*/
template<class T>
inline xfer<T> xferMove(T&);


/**
 * Construct by transferring the contents of the @a arg
 *
 * @sa xferCopy, xferCopyTo, xferMove, xferMoveTo and Foam::xfer
*/
template<class T>
inline xfer<T> xferTmp(Foam::tmp<T>&);


/**
 * Construct by copying the contents of the @a arg
 * between dissimilar types
 *
 * @sa xferCopy, xferMove, xferMoveTo, xferTmp and Foam::xfer
*/
template<class To, class From>
inline xfer<To> xferCopyTo(const From&);


/**
 * Construct by transferring the contents of the @a arg
 * between dissimilar types
 *
 * @par Example Use
 * @code
 *     DynamicList<label> dynLst;
 *     ...
 *     labelList plainLst( xferMoveTo<labelList>(dynLst) );
 * @endcode
 *
 * @sa xferCopy, xferCopyTo, xferMove, xferTmp and Foam::xfer
*/
template<class To, class From>
inline xfer<To> xferMoveTo(From&);


} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "xferI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
