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
    FoamXTypes

Description

\*---------------------------------------------------------------------------*/

#ifndef FoamXTypes_H
#define FoamXTypes_H

// Foam header files
#include "HashTable.H"
#include "wordList.H"

// FoamX header files
#include "FoamX.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                       Class FoamXTypes Declaration
\*---------------------------------------------------------------------------*/

class FoamXTypes
{
    // Private data

        //- The names of the types
        static Foam::wordList typeNames_;

        //- Lookup table of types
        static Foam::HashTable<FoamXServer::FoamXType> types_;

        //- Return this name when the type is invalid
        static Foam::word invalidTypeName_;


public:

    // Constructors

        //- Construct null
        FoamXTypes();


    // Member Functions

        static bool found(const Foam::word& typeName);
        static const Foam::word& typeName(FoamXServer::FoamXType fxType);
        static FoamXServer::FoamXType lookupType(const Foam::word& typeName);
        static bool isNumber(FoamXServer::FoamXType type);
        static bool isPrimitive(FoamXServer::FoamXType type);
        static bool isPrimitive(const Foam::word& typeName);
        static bool isCompound(FoamXServer::FoamXType type);
        static bool isCompound(const Foam::word& typeName);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
