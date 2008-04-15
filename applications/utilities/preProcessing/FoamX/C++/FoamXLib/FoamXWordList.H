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
    FoamXWordList.C

\*---------------------------------------------------------------------------*/

#ifndef FoamXWordList_H
#define FoamXWordList_H

// Foam header files.
#include "wordList.H"

// FoamX header files.
#include "FoamXStringList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                           Class FoamXWordList Declaration
\*---------------------------------------------------------------------------*/

class FoamXWordList
:
    public FoamXStringList
{

public:

    // Constructors

        FoamXWordList();

        FoamXWordList(const Foam::wordList&);

        FoamXWordList(Foam::Istream& is);


    // Member Functions

        // Input/Output.

            void write(Foam::Ostream& os) const;

            void read(Foam::Istream& is);


    // Member Operators

        FoamXWordList& operator=(const FoamXServer::StringList&);
        FoamXWordList& operator=(const Foam::wordList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
