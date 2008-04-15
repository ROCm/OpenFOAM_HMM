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
    FoamXFileNameList.C

\*---------------------------------------------------------------------------*/

#ifndef FoamXFileNameList_H
#define FoamXFileNameList_H

// Foam header files.
#include "fileNameList.H"

// FoamX header files.
#include "FoamXStringList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                      Class FoamXFileNameList Declaration
\*---------------------------------------------------------------------------*/

class FoamXFileNameList
:
    public FoamXStringList
{

public:

    // Constructors

        FoamXFileNameList();

        FoamXFileNameList(const Foam::fileNameList&);

        FoamXFileNameList(Foam::Istream& is);


    // Member Operators

        FoamXFileNameList& operator=(const FoamXServer::StringList&);
        FoamXFileNameList& operator=(const Foam::fileNameList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
