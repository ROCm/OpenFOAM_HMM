/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    A simple container that can be used to copy or transfer the contents
    of objects of type \<T\>.

    Since it is decided upon construction of the xfer object whether the
    parameter is to be copied or transferred, the contents of the resulting
    object can be transferred unconditionally.

    This greatly simplifies defining the constructors for other classes
    with mixed transfer/copy semantics.

    The wrapped object of type \<T\> must implement a transfer() method and
    an operator=() copy method.

Note
    The macros xferCopy(T,arg) and xferMove(T,arg) can be used as
    workarounds for passing temporaries to copy-constructors.

SourceFiles
    xferI.H

\*---------------------------------------------------------------------------*/

#ifndef xfer_H
#define xfer_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class xfer Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class xfer
{
    // Private data

        //- Pointer to temporary object
        mutable T* ptr_;

public:

    // Constructors

        //- Store object pointer and manage its deletion
        //  Can also be used later to transfer by assignment
        inline explicit xfer(T* = 0);

        //- Construct by copying or by transferring the parameter contents
        inline xfer(T&, bool allowTransfer=false);

        //- Construct by copying the parameter contents
        inline xfer(const T&);

        //- Construct by transferring the contents
        inline xfer(const xfer<T>&);

    // Destructor

        inline ~xfer();

    // Member Operators

        //- Transfer the contents into the object
        inline void operator=(T&);

        //- Transfer the contents into the object
        inline void operator=(const xfer<T>&);

        //- Return a non-const reference to const object
        //  @sa xferCopy, xferMove macros alternatives for copy-constructors
        inline const xfer<T>& operator()() const;

        //- Reference to the underlying datatype
        inline T& operator*() const;

        //- Pointer to the underlying datatype
        inline T* operator->() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

/**
 * @def xferCopy(T,arg)
 * Construct by copying the contents of the @a arg
 * and return a const reference to an xfer of type \<T\>
 *
 * Useful for constructors where the argument is temporary.
 * This is a workaround for a template resolution issue.
 *
 * @sa xferMove and Foam::xfer
*/
#define xferCopy(T,arg) \
    (static_cast<const Foam::xfer< T >&>(Foam::xfer< T >(arg)()))

/**
 * @def xferMove(T,arg)
 * Construct by transferring the contents of the @a arg
 * and return a const reference to an xfer of type \<T\>
 *
 * Useful for constructors where the argument is temporary.
 * This is a workaround for a template resolution issue.
 *
 * @par Example Use
 * @code
 *     List<label> a;
 *     ...
 *     List<label> b(xferMove(List<label>, a));
 * @endcode
 *
 * @sa xferCopy and Foam::xfer
*/
#define xferMove(T,arg) \
    (static_cast<const Foam::xfer< T >&>(Foam::xfer< T >(arg, true)()))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "xferI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
