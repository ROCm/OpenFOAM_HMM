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
#include "FoamXWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXWordList::FoamXWordList()
{}

FoamX::FoamXWordList::FoamXWordList(const Foam::wordList& list)
{
    // Reset list length.
    length(list.size());

    // Copy list elements.
    forAll(list, i)
    {
        (*this)[i] = list[i].c_str();
    }
}

FoamX::FoamXWordList::FoamXWordList(Foam::Istream& is)
{
    read(is);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXWordList& FoamX::FoamXWordList::operator=
(
    const FoamXServer::StringList& list
)
{
    // Copy string values. Need to call base classes assignment operator.
    FoamXServer::StringList::operator=(list);

    return *this;
}

FoamX::FoamXWordList& FoamX::FoamXWordList::operator=
(
    const Foam::wordList& list
)
{
    length(list.size());

    forAll(list, i)
    {
        operator[](i) = list[i].c_str();
    }

    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::FoamXWordList::read(Foam::Istream& is)
{
    operator=(static_cast<const Foam::wordList&>(Foam::wordList(is)));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::FoamXWordList::write(Foam::Ostream& os) const
{
    //os << (int)length();

    os << Foam::token::BEGIN_LIST;

    for (unsigned int i = 0; i < length(); i++)
    {
        os << Foam::word(operator[](i));
    }

    os << Foam::token::END_LIST;
}


// ************************************************************************* //
