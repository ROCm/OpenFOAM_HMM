/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    An auto-pointer similar to the STL auto_ptr but with automatic casting
    to a reference to the type and with pointer allocation checking on access.

SourceFiles
    autoPtrI.H

\*---------------------------------------------------------------------------*/

#ifndef autoPtr_H
#define autoPtr_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autoPtr Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class autoPtr
{
    // Public data

        //- Pointer to object
        mutable T* ptr_;


public:

    typedef T Type;


    // Constructors

        //- Construct null with nullptr
        inline explicit autoPtr();

        //- Store object pointer
        inline explicit autoPtr(T* p);

        //- Construct as copy by transferring pointer to this autoPtr and
        //- setting the arguments pointer to nullptr
        inline autoPtr(const autoPtr<T>& ap);

        //- Construct either by transferring pointer or cloning.
        //  Should only be called with type that supports cloning
        inline autoPtr(const autoPtr<T>& ap, const bool reuse);


    //- Destructor, delete object if pointer is not nullptr
    inline ~autoPtr();


    // Member Functions

        // Check

            //- Return true if the autoPtr is empty (ie, no pointer set)
            inline bool empty() const;

            //- Return true if the autoPtr valid (ie, the pointer is set)
            inline bool valid() const;


        // Edit

            //- Release ownership of the object pointer and return the pointer
            inline T* ptr();

            //- Set pointer to that given.
            //  If object pointer already set issue a FatalError
            inline void set(T* p);

            //- If object pointer already set, delete object and set to given
            //- pointer
            inline void reset(T* p = nullptr);

            //- Delete object (if the pointer is valid)
            //- and set pointer to nullptr
            inline void clear();


        // Access

            //- Return the pointer, without nullptr checking.
            //  Pointer remains under autoPtr management.
            inline T* rawPtr();

            //- Const access to the pointer, without nullptr checking.
            //  Pointer remains under autoPtr management.
            inline const T* rawPtr() const;

            //- Return the reference, without nullptr checking.
            inline T& rawRef();

            //- Return the const reference, without nullptr checking.
            inline const T& rawRef() const;


        // Member operators

            //- Return reference to the object data.
            //  Fatal if no pointer is allocated.
            inline T& operator()();

            //- Return const reference to the object data
            //  Fatal if no pointer is allocated.
            inline const T& operator()() const;

            //- Const cast to the underlying type reference
            //  Fatal if no pointer is allocated.
            inline operator const T&() const;

            //- Dereferences (non-const) pointer to the managed object
            //  Fatal if no pointer is allocated.
            inline T* operator->();

            //- Dereferences (const) pointer to the managed object
            //  Fatal if no pointer is allocated.
            inline const T* operator->() const;

            //- Take over the object pointer from parameter
            inline void operator=(T* p);

            //- Take over the object pointer from parameter
            inline void operator=(const autoPtr<T>& ap);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtrI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
