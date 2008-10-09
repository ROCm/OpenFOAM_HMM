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
    A simple container that can be used to transfer the contents of objects
    of type \<T\> of rather than copying them.

    The wrapped object of type \<T\> must implement a transfer() method.

Note
    The macro xferTmp(T,arg) can be used as a workaround for passing
    temporaries to copy-constructors.

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

        //- Null constructor, transfer later by assignment
        inline xfer();

        //- Construct by transferring the parameter contents into the object
        inline xfer(T&);

        //- Construct by transferring the parameter contents into the object
        inline xfer(const xfer<T>&);

    // Destructor

        inline ~xfer();

    // Member Operators

        //- Transfer the contents into the object
        inline void operator=(T&);

        //- Transfer the contents into the object
        inline void operator=(const xfer<T>&);

        //- Return a non-const reference to const object
        //  @sa xferTmp macro for an alternative for copy-constructors
        inline const xfer<T>& operator()() const;

        //- Reference to the underlying datatype
        inline T& operator*() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

/**
 * @def xferTmp(T,arg)
 * Construct a temporary and return a const reference to an xfer of
 * type \<T\> and return a const reference.
 *
 * Useful for copy-constructors where the argument is temporary.
 * This is a workaround for a template resolution problem.
 *
 * @par Example Use
 * @code
 *     List<label> a;
 *     ...
 *     List<label> b(xferTmp(List<label>, a));
 *
 * @endcode
 * @sa xfer class
*/

#define xferTmp(T,arg) \
    (static_cast<const Foam::xfer< T >&>(Foam::xfer< T >(arg)()))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "xferI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
