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
#include "FoamXStringList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXStringList::FoamXStringList()
{}

FoamX::FoamXStringList::FoamXStringList(const FoamX::FoamXStringList& list)
:
    StringSequenceTmpl<CORBA::String_var>()
{
    // Copy string values. Need to call base classes assignment operator.
    FoamXServer::StringList::operator=(list);
}

FoamX::FoamXStringList::FoamXStringList(const Foam::stringList& list)
{
    // Reset list length.
    length(list.size());

    // Copy list elements.
    forAll(list, i)
    {
        (*this)[i] = list[i].c_str();
    }
}

FoamX::FoamXStringList::FoamXStringList(Foam::Istream& is)
{
    read(is);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXStringList& FoamX::FoamXStringList::operator=
(
    const FoamXServer::StringList& list
)
{
    // Copy string values. Need to call base classes assignment operator.
    FoamXServer::StringList::operator=(list);

    return *this;
}

FoamX::FoamXStringList& FoamX::FoamXStringList::operator=
(
    const Foam::stringList& list
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

void FoamX::FoamXStringList::append(const char* str)
{
    // Append string to list.
    length(length() + 1);
    (*this)[length() - 1] = str;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int FoamX::FoamXStringList::find(const char* str)
{
    int ret = -1;

    for (unsigned int i = 0; i <length(); i++)
    {
        if (strcmp((*this)[i], str) == 0)
        {
            ret = i;
            break;
        }
    }

    return ret;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool FoamX::FoamXStringList::remove(const char* str)
{
    bool bRet = false;

    // See if this list contains the specified entry.
    int index = find(str);

    // Remove name from list.
    if (index != -1)
    {
        FoamXServer::StringList* newList = new FoamXServer::StringList();
        newList->length(length() - 1);
        for (unsigned int i = 0, j = 0; i <length(); i++)
        {
            if (i != (unsigned int)(index))
            {
                (*newList)[j++] = (*this)[i];
            }
        }
        FoamXServer::StringList::operator=(*newList);
        delete newList;
        bRet = true;
    }

    return bRet;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::FoamXStringList::read(Foam::Istream& is)
{
    operator=(static_cast<const Foam::stringList&>(Foam::stringList(is)));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::FoamXStringList::write(Foam::Ostream& os) const
{
    //os << Foam::label(length());

    os << Foam::token::BEGIN_LIST;

    for (unsigned int i = 0; i < length(); i++)
    {
        os << Foam::string(operator[](i));
    }

    os << Foam::token::END_LIST;
}


// ************************************************************************* //
