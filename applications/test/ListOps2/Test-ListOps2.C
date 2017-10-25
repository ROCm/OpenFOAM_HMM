/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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

Application
    Test-ListOps2

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "List.H"
#include "FixedList.H"
#include "DynamicList.H"
#include "SubList.H"
#include "ListOps.H"
#include "FlatOutput.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ListType>
void testFind(const ListType& list)
{
    Info<< nl << "list: " << flatOutput(list) << nl << endl;

    for (auto val : { 20, 35, 6, 13, 12, -12})
    {
        Info<< "lookup " << val
            << " found " << list.found(val)
            << " index " << list.find(val) << nl;
    }
}


template<class ListType>
void testMoving(ListType& list)
{
    Info<< nl << "list: " << flatOutput(list) << nl << endl;

    {
        const label i = 3;
        list.swapFirst(i);
        Info<<"swapFirst: " << i << " = " << flatOutput(list) << nl;
    }

    {
        const label i = 6;
        list.moveFirst(i);
        Info<<"moveFirst: " << i << " = " << flatOutput(list) << nl;
    }

    {
        const label i = 6;
        list.moveLast(i);
        Info<<"moveLast: " << i  << " = " << flatOutput(list) << nl;
    }

    {
        const label i = 8;
        list.swapLast(i);
        Info<<"swapLast: " << i << " = " << flatOutput(list) << nl;
    }
}


// Main program:

int main(int argc, char *argv[])
{
    List<label> list1(identity(15));
    shuffle(list1);

    FixedList<label, 15> list2(list1);
    inplaceReverseList(list2);

    DynamicList<label> list3(list1);
    inplaceReverseList(list3);

    testFind(list1);
    testFind(list2);
    testFind(list3);

    testMoving(list1);
    testMoving(list2);
    testMoving(list3);

    // Test remove
    {
        auto& list = list3;

        Info<< nl << "list: " << flatOutput(list) << nl << endl;

        list.remove();
        Info<<"remove = " << flatOutput(list) << nl;

        {
            const label i = 6;
            list.remove(i);
            Info<<"rmSwap: " << i << " = " << flatOutput(list) << nl;
        }

        {
            const label i = 3;
            list.remove(i);
            Info<<"rmSwap: " << i << " = " << flatOutput(list) << nl;
        }

        {
            const label i = 8;
            list.remove(i);
            Info<<"rmMove: " << i << " = " << flatOutput(list) << nl;
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
