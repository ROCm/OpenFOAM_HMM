/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "DynamicList.H"
#include "IOstreams.H"
#include "FlatOutput.H"
#include "ListOps.H"
#include "labelRange.H"
#include "labelIndList.H"

using namespace Foam;

template<class T>
void printInfo
(
    const word& tag,
    const UList<T>& lst,
    const bool showSize = false
)
{
    Info<< "<" << tag;
    if (showSize)
    {
        Info<< " size=\"" << lst.size() << "\"";
        if (lst.cdata())
        {
            Info<< " ptr=\"" << long(lst.cdata()) << "\"";
        }
        else
        {
            Info<< " ptr=\"nullptr\"";
        }
    }
    Info<< ">" << nl << flatOutput(lst) << nl << "</" << tag << ">" << endl;
}


template<class T, int SizeMin>
void printInfo
(
    const word& tag,
    const DynamicList<T, SizeMin>& lst,
    const bool showSize = false
)
{
    Info<< "<" << tag;
    if (showSize)
    {
        Info<< " size=\"" << lst.size()
            << "\" capacity=\"" << lst.capacity() << "\"";
        if (lst.cdata())
        {
            Info<< " ptr=\"" << long(lst.cdata()) << "\"";
        }
        else
        {
            Info<< " ptr=\"nullptr\"";
        }
    }
    Info<< ">" << nl << flatOutput(lst) << nl << "</" << tag << ">" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<DynamicList<label>> ldl(2);

    ldl[0](0) = 0;
    ldl[0](2) = 2;
    ldl[0](3) = 3;
    ldl[0](1) = 1;

    ldl[0].setCapacity(5);    // increase allocated size
    ldl[1].setCapacity(10);   // increase allocated size
    ldl[0].reserve(15);       // should increase allocated size
    ldl[1].reserve(5);        // should not decrease allocated size
    ldl[1](3) = 2;            // allocates space and sets value

    // this works without a segfault, but doesn't change the list size
    ldl[0][4] = 4;

    ldl[1] = 3;

    Info<< "<ldl>" << flatOutput(ldl) << "</ldl>" << nl << "sizes: ";
    forAll(ldl, i)
    {
        Info<< " " << ldl[i].size() << "/" << ldl[i].capacity();
    }
    Info<< endl;

    List<List<label>> ll(2);
    ll[0].transfer(ldl[0]);
    ll[1].transfer(ldl[1].shrink());

    Info<< "<ldl>" << flatOutput(ldl) << "</ldl>" << nl << "sizes: ";
    forAll(ldl, i)
    {
        Info<< " " << ldl[i].size() << "/" << ldl[i].capacity();
    }
    Info<< endl;

    Info<< "<ll>" << ll << "</ll>" << nl << endl;


    // test the transfer between DynamicLists
    DynamicList<label> dlA
    {
        0, 1, 2, 3, 4
    };
    dlA.append({ 5, 6 });
    dlA = { 1, 2, 4 };

    DynamicList<label> dlB;

    dlA.setCapacity(10);

    Info<< "<dlA>" << flatOutput(dlA) << "</dlA>" << nl << "sizes: "
        << " " << dlA.size() << "/" << dlA.capacity() << endl;

    dlB.transfer(dlA);

    // provokes memory error if previous transfer did not maintain
    // the correct allocated space
    dlB[6] = 6;

    Info<< "Transferred to dlB" << endl;
    Info<< "<dlA>" << flatOutput(dlA) << "</dlA>" << nl << "sizes: "
        << " " << dlA.size() << "/" << dlA.capacity() << endl;
    Info<< "<dlB>" << flatOutput(dlB) << "</dlB>" << nl << "sizes: "
        << " " << dlB.size() << "/" << dlB.capacity() << endl;

    // try with a normal list:
    List<label> lstA;
    lstA.transfer(dlB);
    Info<< "Transferred to normal list" << endl;
    printInfo("lstA", lstA, true);
    printInfo("dlB", dlB, true);

    // Copy back and append a few time
    for (label i=0; i < 3; i++)
    {
        dlB.append(lstA);
    }

    Info<< "appended list a few times" << endl;
    printInfo("dlB", dlB, true);

    // assign the list (should maintain allocated space)
    dlB = lstA;
    Info<< "assigned list" << endl;
    printInfo("dlB", dlB, true);

    // Copy back and append a few time
    for (label i=0; i < 3; i++)
    {
        dlB.append(lstA);
    }


    // check allocation granularity
    DynamicList<label> dlC;

    printInfo("dlC", dlC, true);

    dlC.reserve(dlB.size());
    dlC = dlB;

    printInfo("dlC", dlC, true);

    List<label> lstB(dlC.xfer());

    Info<< "Transferred to normal list via the xfer() method" << endl;
    printInfo("lstB", lstB, true);
    printInfo("dlC", dlC, true);

    DynamicList<label> dlD(lstB.xfer());

    Info<< "Transfer construct from normal list" << endl;
    printInfo("lstB", lstB, true);
    printInfo("dlD", dlD, true);

    DynamicList<label,10> dlE1(10);
    DynamicList<label> dlE2(dlE1);   // construct dissimilar

    printInfo("dlE1", dlE1, true);
    printInfo("dlE2", dlE2, true);

    for (label elemI=0; elemI < 5; ++elemI)
    {
        dlE1.append(4 - elemI);
        dlE2.append(elemI);
    }

    printInfo("dlE2", dlE2, true);

    DynamicList<label> dlE3(dlE2);   // construct identical
    printInfo("dlE3", dlE3, true);

    dlE3 = dlE1;   // assign dissimilar
    printInfo("dlE3", dlE3, true);

    dlE3 = dlE2;   // assign identical
    printInfo("dlE3", dlE3, true);

    DynamicList<label> dlE4(reorder(identity(dlE3.size()), dlE3));
    printInfo("dlE4", dlE4, true);

    printInfo("dlE3", dlE3, true);


    {
        DynamicList<label> addr(10);
        addr.append(3);
        addr.append(1);
        addr.append(2);

        forAll(dlE2, i)
        {
            dlE2[i] *= 10;
        }

        labelUIndList uil
        (
            dlE2, addr
        );
        Info<< "use UIndirectList " << uil << " remapped from " << dlE2 << endl;
        dlE4 = uil;
        printInfo("dlE4", dlE4, true);
    }

    {
        Info<< nl << "Test moving:" << nl;

        labelList input1 = identity(15);
        labelList input2 = identity(15);
        inplaceReverseList(input2);

        DynamicList<label> list1(std::move(input1));
        DynamicList<label> list2;

        Info<< "move construct:" << nl
            << "input: " << flatOutput(input1) << nl
            << "list:  " << flatOutput(list1) << endl;

        list1 = std::move(input2);

        Info<< "move assignment:" << nl
            << "input: " << flatOutput(input2) << nl
            << "list:  " << flatOutput(list1) << endl;

        list2 = std::move(list1);
        Info<< "list in:  " << flatOutput(list1) << nl
            << "list out: " << flatOutput(list2) << endl;

        input2 = std::move(identity(15));
        list2 = std::move(input2);
        Info<< "list in:  " << flatOutput(input2) << nl
            << "list out: " << flatOutput(list2) << endl;

        input1 = identity(15);
        input2 = identity(15);
        inplaceReverseList(input2);

        Info<< "test move-append with "
            << flatOutput(input1) << " and " << flatOutput(input2) << endl;

        list2.append(std::move(list1));
        list2.append(std::move(input1));
        list2.append(std::move(input2));

        Info<< "result: " << flatOutput(list2) << nl
            << "inputs: " << flatOutput(list1) << " / "
            << flatOutput(input1) << " / "
            << flatOutput(input2) << nl;

        Info<< "test move dissimilar sizing:" << nl;
        list1 = list2;
        list1.reserve(100);

        // DynamicList<label,1000> list3; // (std::move(list1));
        DynamicList<label,1000> list3(std::move(list1));
        Info<< "orig: " << flatOutput(list1) << nl;

        // list3.swap(list1);
        // list3 = std::move(list1);

        printInfo("input",  list1, true);
        printInfo("output", list3, true);

        input1 = list2;

        Info<< nl << "test subset/remove with "
            << flatOutput(input1) << endl;

        for
        (
            const labelRange range :
            {
                labelRange(-10, 8),  // invalid range
                labelRange(40, 18),  // trailing portion
                labelRange(-5, 10),  // leading portion
                labelRange(10, 8),   // mid-portion
                labelRange(0, input1.size()), // everything
            }
        )
        {
            list1 = input1;
            list2 = input1;

            list1.remove(range);
            list2.subset(range);

            Info<< "input = " << flatOutput(input1) << nl
                << "remove " << range << " = " << flatOutput(list1) << nl
                << "subset " << range << " = " << flatOutput(list2) << nl;
        }
    }

    Info<< "\nEnd\n";

    return 0;
}


// ************************************************************************* //
