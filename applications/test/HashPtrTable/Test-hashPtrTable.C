/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include <iostream>
#include "HashPtrTable.H"

using namespace Foam;

template<class T>
void printTable(const HashPtrTable<T>& table)
{
    Info<< table.size() << nl << "(" << nl;

    for
    (
        typename HashPtrTable<T>::const_iterator iter = table.cbegin();
        iter != table.cend();
        ++iter
    )
    {
        const T* ptr = *iter;
        Info<< iter.key() << " = ";
        if (ptr)
        {
            Info<< *ptr;
        }
        else
        {
            Info<< "nullptr";
        }
        Info<< nl;
    }

    Info<< ")" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    HashPtrTable<double> myTable;
    myTable.insert("abc", new double(42.1));
    myTable.insert("def", nullptr);
    myTable.insert("ghi", new double(3.14159));

    // Info<< myTable << endl;
    printTable(myTable);

    HashPtrTable<double> copy(myTable);

    // Info<< copy << endl;
    printTable(copy);
    Info<< copy << endl;

    return 0;
}


// ************************************************************************* //
