/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Application

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "IOstreams.H"
#include "IStringStream.H"
#include "scalar.H"
#include "vector.H"
#include "ListOps.H"
#include "List.H"
#include "PackedBoolList.H"
#include <bitset>

using namespace Foam;

template <int nBits>
void printPackedList(const PackedList<nBits>& L)
{
    const List<unsigned int>& stor = L.storage();

    cout<< "PackedList<" << nBits << ">"
        << " max_bits:"  << L.max_bits()
        << " max_value:" << L.max_value()
        << " packing:"   << L.packing() << nl;

    cout<< "values:  " << L.size() << "/" << L.capacity() << " ( ";
    forAll(L, i)
    {
        cout<< L[i] << ' ';
    }
    cout<< ")\n\n";

    cout<< "storage: " << stor.size() << "( ";
    forAll(stor, i)
    {
        cout<< std::bitset<32>(stor[i]) << ' ';
    }
    cout<< ")\n" << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{    
    cout<< "PackedList::max_bits() = " << PackedList<0>::max_bits() << nl;

    PackedList<3> list1(5,1);

    printPackedList(list1);

    list1 = 2;
    printPackedList(list1);

    list1.resize(6, 3);
    printPackedList(list1);

    list1 = false;
    printPackedList(list1);

    list1 = true;
    printPackedList(list1);

    list1.resize(12);
    printPackedList(list1);

    list1.resize(25, list1.max_value());
    printPackedList(list1);

    list1.resize(8);
    printPackedList(list1);

    return 0;
}


// ************************************************************************* //
