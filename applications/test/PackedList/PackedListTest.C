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

#include "IOstreams.H"
#include "PackedBoolList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    bool changed;
    Info<< "PackedList max_bits() = " << PackedList<0>::max_bits() << nl;

    Info<< "\ntest allocation with value\n";
    PackedList<3> list1(5,1);
    list1.print(Info);

    Info<< "\ntest assign uniform value\n";
    list1 = 4;
    list1.print(Info);

    Info<< "\ntest resize with value (without reallocation)\n";
    list1.resize(8, list1.max_value());
    list1.print(Info);

    Info<< "\ntest set() function\n";
    list1.set(1, 5);
    list1.print(Info);

    Info<< "\ntest assign bool\n";
    list1 = false;
    list1.print(Info);

    Info<< "\ntest assign bool\n";
    list1 = true;
    list1.print(Info);

    Info<< "\ntest resize without value (with reallocation)\n";
    list1.resize(12);
    list1.print(Info);

    Info<< "\ntest resize with value (with reallocation)\n";
    list1.resize(25, list1.max_value());
    list1.print(Info);

    Info<< "\ntest resize smaller (should not touch allocation)\n";
    list1.resize(8);
    list1.print(Info);

    Info<< "\ntest append() operation\n";
    list1.append(2);
    list1.append(3);
    list1.append(4);
    list1.print(Info);

    Info<< "\ntest reserve() operation\n";
    list1.reserve(32);
    list1.print(Info);

    Info<< "\ntest shrink() operation\n";
    list1.shrink();
    list1.print(Info);

    Info<< "\ntest setCapacity() operation\n";
    list1.setCapacity(15);
    list1.print(Info);

    Info<< "\ntest setCapacity() operation\n";
    list1.setCapacity(100);
    list1.print(Info);

    Info<< "\ntest operator[] assignment\n";
    list1[16] = 5;
    list1.print(Info);

    Info<< "\ntest operator[] assignment with auto-vivify\n";
    list1[36] = list1.max_value();
    list1.print(Info);

    Info<< "\ntest setCapacity smaller\n";
    list1.setCapacity(24);
    list1.print(Info);

    // add in some misc values
    list1[31] = 1;
    list1[32] = 2;
    list1[33] = 3;

    Info<< "\ntest iterator\n";
    PackedList<3>::iterator iter = list1.begin();
    Info<< "begin():";
    iter.print(Info) << "\n";

    Info<< "\ntest iterator operator=\n";
    changed = (iter = 5);

    Info<< "iterator:" << iter() << "\n";
    Info<< "changed:" << changed << "\n";
    changed = (iter = 5);
    Info<< "changed:" << changed << "\n";
    list1.print(Info);

    Info<< "\ntest get() method\n";
    Info<< "get(10):" << list1.get(10) << " and list[10]:" << list1[10] << "\n";
    list1.print(Info);

    Info<< "\ntest iterator indexing\n";
    Info<< "cend() ";
    list1.cend().print(Info) << "\n";

    for 
    (
        PackedList<3>::const_iterator cit = list1[30];
        cit != list1.cend();
        ++cit)
    {
        cit.print(Info);
    }

    Info<< "\ntest operator[] auto-vivify\n";
    const unsigned int val = list1[45];

    Info<< "list[45]:" << val << "\n";
    list1[45] = list1.max_value();
    Info<< "list[45]:" << list1[45] << "\n";
    list1[49] = list1.max_value();
    list1.print(Info);


    Info<< "\ntest copy constructor + append\n";
    PackedList<3> list2(list1);
    list2.append(4);
    Info<< "source list:\n";
    list1.print(Info);
    Info<< "destination list:\n";
    list2.print(Info);

    Info<< "\ntest pattern that fills all bits\n";
    PackedList<4> list3(8, 8);
    
    label pos = list3.size() - 1;

    list3[pos--] = list3.max_value();
    list3[pos--] = 0;
    list3[pos--] = list3.max_value();
    list3.print(Info);

    Info<< "removed final value: " << list3.remove() << endl;
    list3.print(Info);

    Info<< "\n\nDone.\n";

    return 0;
}


// ************************************************************************* //
