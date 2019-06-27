/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include <memory>
#include <iostream>
#include "autoPtr.H"
#include "HashPtrTable.H"

using namespace Foam;

template<class T>
void printTable(const HashPtrTable<T>& table)
{
    Info<< table.size() << nl << "(" << nl;

    forAllConstIters(table, iter)
    {
        const T* ptr = iter.val();
        Info<< iter.key() << " = ";
        if (ptr)
        {
            Info<< *ptr << " (" << uintptr_t(ptr) << ")";
        }
        else
        {
            Info<< "nullptr";
        }
        Info<< nl;
    }

    Info<< ")" << endl;

    // Iterate across values, with for-range
    Info<< "values (";
    for (const auto& ptr : table)
    {
        Info<< ' ';
        if (ptr)
        {
            Info<< *ptr;
        }
        else
        {
            Info<< "nullptr";
        }
    }
    Info<< " )" << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    HashPtrTable<double> myTable;
    myTable.set("abc", new double(42.1));
    myTable.set("def", nullptr);
    myTable.set("pi", new double(3.14159));
    myTable.set("natlog", new double(2.718282));
    myTable.insert("sqrt2", autoPtr<double>::New(1.414214));

    HashTable<std::unique_ptr<double>, word, string::hash> myTable1;

    myTable1.set("abc", std::unique_ptr<double>(new double(42.1)));
    myTable1.set("pi", std::unique_ptr<double>(new double(3.14159)));
    myTable1.set("natlog", std::unique_ptr<double>(new double(2.718282)));

    HashTable<autoPtr<double>, word, string::hash> myTable2;

    myTable2.set("abc", autoPtr<double>(new double(42.1)));
    myTable2.set("pi", autoPtr<double>(new double(3.14159)));
    myTable2.set("natlog", autoPtr<double>(new double(2.718282)));
    myTable2.insert("sqrt2", autoPtr<double>::New(1.414214));

    // Info<< myTable << endl;
    printTable(myTable);
    Info<< myTable2 << nl;

    auto iter2 = myTable2.find("pi");

    Info<<"PI: " << **iter2 << nl;


    HashPtrTable<double> copy(myTable);

    // Info<< copy << endl;
    printTable(copy);
    Info<< copy << endl;

    Info<<"\nerase some existing and non-existing entries" << nl;

    auto iter = myTable.find("pi");
    myTable.erase(iter);

    iter = myTable.find("unknownKey");
    myTable.erase(iter);

    myTable.erase("abc");
    myTable.erase("unknownKey");

    printTable(myTable);

    HashPtrTable<double> moved(std::move(copy));

    Info<< nl << "test movable" << nl;
    Info<<"input:" << nl;
    printTable(copy);

    Info<<"output:" << nl;
    printTable(moved);

    HashPtrTable<double> other;

    Info<<"move assign" << nl;

    other = std::move(moved);
    printTable(other);

    Info<<"old" << nl;
    printTable(moved);

    return 0;
}

// ************************************************************************* //
