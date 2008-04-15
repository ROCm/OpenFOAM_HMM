/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    FoamXString

Description

SourceFiles
    FoamXString.C

\*---------------------------------------------------------------------------*/

#ifndef FoamXString_H
#define FoamXString_H

// Foam header files.
#include "Istream.H"
#include "Ostream.H"
#include "token.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                           Class FoamXString Declaration
\*---------------------------------------------------------------------------*/

class FoamXString
:
    public CORBA::String_var
{
    // Private data

        // Additional information to allow mapping between CORBA and Foam types.

            //- Foam type  - WORD or STRING.
            Foam::token::tokenType tokenType_;

            //- FoamX type - Type_Char, Type_Word or Type_String.
            FoamXServer::FoamXType foamXType_;
  

public:

    // Constructors

        FoamXString();

        FoamXString(const FoamXString&);

        FoamXString(const Foam::string&);

        FoamXString(const Foam::word&);

        FoamXString(const char* str);

        FoamXString(Foam::Istream& is);


    // Destructor

        ~FoamXString();


    // Member Functions

        // I/O methods.

            void write(Foam::Ostream& os) const;

            void read(Foam::Istream& is);


    // Member Operators

        // Assignment operators.

            FoamXString& operator=(const FoamXString&);

            FoamXString& operator=(const Foam::string&);

            FoamXString& operator=(const Foam::word&);

            FoamXString& operator=(const char*);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
