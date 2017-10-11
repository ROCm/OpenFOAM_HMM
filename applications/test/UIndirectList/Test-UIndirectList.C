/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "UIndirectList.H"
#include "DynamicList.H"
#include "IOstreams.H"
#include "ListOps.H"
#include "labelIndList.H"

using namespace Foam;

template<class ListType>
void printInfo(const ListType& lst)
{
    Info<< "addr: " << flatOutput(lst.addressing()) << nl
        << "list: " << flatOutput(lst) << nl
        << endl;
}

template<class T, class ListType>
void testFind(const T& val, const ListType& lst)
{
    Info<< nl
        << "Search for "<< val << " in " << flatOutput(lst) << nl
        <<" found() = " << lst.found(val)
        <<" find() = " << lst.find(val)
        <<" rfind() = " << lst.rfind(val)
        <<" find(2) = " << lst.find(val, 2)
        <<" rfind(2) = " << lst.rfind(val, 2)
        <<" findIndex = " << findIndex(lst, val) << nl
        << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<label> completeList(20);

    forAll(completeList, i)
    {
        completeList[i] = 10*i;
    }

    Info<< "raw : " << flatOutput(completeList) << nl << endl;

    List<label> addresses{1, 0, 3, 7, 4, 8, 5, 1, 0, 3, 7, 4, 8, 5, };

    labelUIndList idl1(completeList, addresses);

    printInfo(idl1);

    for (const label val : { 10, 30, 40, 50, 90, 80, 120 } )
    {
        testFind(val, idl1);
    }

    Info<< flatOutput(idl1) << nl;

    idl1[1] = -666;

    Info<< "idl1[1] changed: " << flatOutput(idl1) << endl;

    idl1 = -999;

    Info<< "idl1 changed: " << flatOutput(idl1) << endl;

    labelUIndList idl2(idl1);

    Info<< "idl2: " << flatOutput(idl2) << endl;

    {
        List<label> ident(idl1.size());

        forAll(ident, i)
        {
            ident[i] = ident.size() - i;
        }
        idl1 = ident;
    }

    Info<< "idl1 assigned from UList: " << flatOutput(idl1) << endl;

    // test List operations

    List<label> flatList(labelUIndList(completeList, addresses));
    Info<< "List construct from UIndirectList: " << flatOutput(flatList) << nl;

    flatList = labelUIndList(completeList, addresses);
    Info<< "List assign from UIndirectList: " << flatOutput(flatList) << nl;

    flatList.append(labelUIndList(completeList, addresses));
    Info<< "List::append(UIndirectList): " << flatOutput(flatList) << nl;


    DynamicList<label> dynList(labelUIndList(completeList, addresses));
    Info<< "DynamicList construct from UIndirectList: " << flatOutput(dynList)
        << nl;

    dynList.append(labelUIndList(completeList, addresses));
    Info<< "DynamicList::append(UIndirectList): " << flatOutput(dynList) << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
