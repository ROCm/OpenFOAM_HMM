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

\*---------------------------------------------------------------------------*/

// FoamX header files.
#include "FoamX.H"
#include "FoamXFileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXFileNameList::FoamXFileNameList()
{}

FoamX::FoamXFileNameList::FoamXFileNameList(const Foam::fileNameList& list)
{
    // Reset list length.
    length(list.size());

    // Copy list elements.
    forAll(list, i)
    {
        (*this)[i] = (const char*)list[i].c_str();
    }
}

FoamX::FoamXFileNameList::FoamXFileNameList(Foam::Istream& is)
:
    FoamXStringList(is)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXFileNameList& FoamX::FoamXFileNameList::operator=
(
    const FoamXServer::StringList& list
)
{
    // Copy string values. Need to call base classes assignment operator.
    FoamXServer::StringList::operator=(list);

    return *this;
}

FoamX::FoamXFileNameList& FoamX::FoamXFileNameList::operator=
(
    const Foam::fileNameList& list
)
{
    length(list.size());

    forAll(list, i)
    {
        operator[](i) = list[i].c_str();
    }

    return *this;
}


// ************************************************************************* //
