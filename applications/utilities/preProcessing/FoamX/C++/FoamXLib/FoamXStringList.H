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
    FoamXStringList.C

\*---------------------------------------------------------------------------*/

#ifndef FoamXStringList_H
#define FoamXStringList_H

// Foam header files.
#include "stringList.H"

// FoamX header files.
#include "FoamX.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                           Class FoamXStringList Declaration
\*---------------------------------------------------------------------------*/

class FoamXStringList
:
    public FoamXServer::StringList
{

public:

    // Constructors

        FoamXStringList();

        FoamXStringList(const FoamXStringList&);

        FoamXStringList(const Foam::stringList&);

        FoamXStringList(Foam::Istream& is);


    // Member Functions

        // Utility functions.

            Foam::label size() const
            {
                return length();
            }

            void append(const char*);

            //- Returns -1 if string not found on list, or the index otherwise.
            int find(const char*);

            bool found(const char* name)
            {
                return (find(name) != -1);
            }

            bool remove(const char*);


        // Input/Output.

            void write(Foam::Ostream& os) const;

            void read(Foam::Istream& is);


    // Member Operators

        FoamXStringList& operator=(const FoamXServer::StringList&);
        FoamXStringList& operator=(const Foam::stringList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
