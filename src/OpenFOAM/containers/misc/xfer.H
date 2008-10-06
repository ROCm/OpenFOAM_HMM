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
    A simple container for objects of type \<T\> that can be used for
    transferring the contents rather than copying them.

    The wrapped object of type \<T\> must implement a transfer() method.

SourceFiles
    xferI.H

\*---------------------------------------------------------------------------*/

#ifndef xfer_H
#define xfer_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
template<class T> class xfer;

template<class T> Istream& operator>>(Istream&, xfer<T>&);
template<class T> Ostream& operator<<(Ostream&, const xfer<T>&);


/*---------------------------------------------------------------------------*\
                         Class xfer Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class xfer
{
    // Private data

        //- Reference to the underlying datatype
        T& ref_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        xfer(const xfer<T>&);

        //- Disallow default bitwise assignment
        void operator=(const xfer<T>&);

public:

    // Constructors

        //- Construct null, delaying the transfer to later
        inline xfer();

        //- Construct by transferring the parameter into this object
        inline xfer(T&);

        //- Construct by transferring the parameter into this object
        inline xfer(xfer<T>&);

    // Destructor

        inline ~xfer();

    // Member Functions

        //- Access to the underlying datatype
        inline T& ref();

        //- Transfer the contents of the parameter into this object
        inline void transfer(T&);

        //- Transfer the contents of the parameter into this object
        inline void transfer(xfer<T>&);

        //- Transfer the contents from this object to the parameter
        inline void yield(T&);

        //- Transfer the contents from this object to the parameter
        inline void yield(xfer<T>&);

    // Member Operators

        //- Transfer the contents into this object
        inline void operator=(T&);

        //- Transfer the contents into this object
        inline void operator=(xfer<T>&);

        //- Return a non-const reference to const object
        //  Useful for read-constructors where the argument is temporary
        inline xfer<T>& operator()() const;


    // Istream operator

        //- Read from Istream, discarding current contents
        friend Istream& operator>>
        #ifndef __CINT__
        <T>
        #endif
        (Istream&, xfer<T>&);

    // Ostream operator

        //- Write to Ostream
        friend Ostream& operator<<
        #ifndef __CINT__
        <T>
        #endif
        (Ostream&, const xfer<T>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "xferI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
