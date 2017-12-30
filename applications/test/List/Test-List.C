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

Application
    Test-List

Description
    Simple tests and examples of use of List

See also
    Foam::List

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "wordReList.H"

#include "IOstreams.H"
#include "StringStream.H"
#include "scalar.H"
#include "vector.H"

#include "labelRange.H"
#include "scalarList.H"
#include "ListOps.H"
#include "SubList.H"

#include <list>
#include <numeric>

using namespace Foam;

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
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("reList", "reList");
    argList::addOption("wordList", "wordList");
    argList::addOption("stringList", "stringList");
    argList::addOption("float", "xx");
    argList::addBoolOption("flag");

    #include "setRootCase.H"

    {
        List<label> ident(15);
        std::iota(ident.begin(), ident.end(), 0);

        Info<<"Ident:";
        forAllConstIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<<"reverse:";
        forAllReverseIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<<"const reverse:";
        forAllConstReverseIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;
    }


    if (false)
    {
        labelList intlist(IStringStream("(0 1 2)")());
        Info<<"construct from Istream: " << intlist << endl;

        IStringStream("(3 4 5)")() >> static_cast<labelUList&>(intlist);
        Info<<"is >>: " << intlist << endl;

        IStringStream("(6 7 8)")() >> intlist;
        Info<<"is >>: " << intlist << endl;
    }

    List<vector> list1(IStringStream("1 ((0 1 2))")());
    Info<< "list1: " << list1 << endl;

    List<vector> list2
    {
        vector(0, 1, 2),
        vector(3, 4, 5),
        vector(6, 7, 8),
        vector(0, 1, 2),
        vector(3, 4, 5),
        vector(6, 7, 8),
    };
    Info<< "list2: " << list2 << endl;

    Info<< "forAllConstIters(list2): ";
    forAllConstIters(list2, iter) { Info<< " " << *iter; }
    Info<< endl;

    Info<< "forAllConstReverseIters(list2): ";
    forAllConstReverseIters(list2, iter) { Info<< " " << *iter; }
    Info<< endl;

    Info<< "forAllConstIters(list2): ";
    forAllIters(list2, iter) { *iter *= 2; Info<< " " << *iter; }
    Info<< endl;

    Info<< "forAllReverseIters(list2): ";
    forAllReverseIters(list2, iter) { *iter *= 0.5; Info<< " " << *iter; }
    Info<< endl;

    list1.append(list2);
    Info<< "list1.append(list2): " << list1 << endl;

    for (const vector& val : { vector(3, 4, 5), vector(10,11, 12)} )
    {
        testFind(val, list2);
    }

    list2.setSize(10, vector(1, 2, 3));
    Info<< "list2: " << list2 << endl;

    List<vector> list3(list2.xfer());
    Info<< "Transferred via the xfer() method" << endl;
    Info<< "list2: " << list2 << nl
        << "list3: " << list3 << endl;

    List<vector> list4
    {
        vector(0, 1, 2),
        vector(3, 4, 5),
        vector(6, 7, 8)
    };
    Info<< "list4: " << list4 << endl;

    List<vector> list5
    {
        {5, 3, 1},
        {10, 2, 2},
        {8, 1, 0}
    };
    Info<< "list5: " << list5 << endl;
    list5 =
    {
        {8, 1, 0},
        {5, 3, 1},
        {10, 2, 2}

    };
    Info<< "list5: " << list5 << endl;

    list4.swap(list5);
    Info<< "Swapped via the swap() method" << endl;
    Info<< "list4: " << list4 << nl
        << "list5: " << list5 << endl;

    List<vector> list6(list4.begin(), list4.end());
    Info<< "list6: " << list6 << endl;

    // Subset
    const labelList map{0, 2};
    List<vector> subList3(list3, map);
    Info<< "Elements " << map << " out of " << list3
        << " => " << subList3 << endl;

    // test flattened output
    {
        Info<< nl;

        labelList longLabelList = identity(15);

        // This does not work:
        // scalarList slist = identity(15);
        //
        // More writing, but does work:
        scalarList slist
        (
            labelRange::null.begin(),
            labelRange::identity(15).end()
        );

        Info<<"scalar identity:" << flatOutput(slist) << endl;

        Info<< "labels (contiguous=" << contiguous<label>() << ")" << nl;

        Info<< "normal: " << longLabelList << nl;
        Info<< "flatOutput: " << flatOutput(longLabelList) << nl;
        // Info<< "flatOutput(14): " << flatOutput(longLabelList, 14) << nl;

        stringList longStringList(12);
        forAll(longStringList, i)
        {
            longStringList[i].resize(3, 'a' + i);
        }

        Info<< "string (contiguous=" << contiguous<string>() << ")" << nl;

        Info<< "normal: " << longStringList << nl;
        Info<< "flatOutput: " << flatOutput(longStringList) << nl;
        // contiguous longStringList[i].resize(3, 'a' + i);
    }

    // test SubList and labelRange
    {
        Info<< nl;
        labelList longLabelList = identity(25);
        reverse(longLabelList);

        FixedList<label, 6> fixedLabelList{0,1,2,3,4,5};
        const labelList constLabelList = identity(25);

        Info<< "full-list: " << flatOutput(longLabelList) << nl;

        labelRange range1(-15, 25);
        Info<<"sub range:" << range1 << "=";
        Info<< SubList<label>(longLabelList, range1) << nl;

        labelRange range2(7, 8);
        Info<<"sub range:" << range2 << "=";
        Info<< SubList<label>(longLabelList, range2) << nl;

        // labelRange range2(7, 8);
        Info<<"use range " << range2 << " to set value";
        SubList<label>(longLabelList, range2) = -15;
        Info<< "=> " << flatOutput(longLabelList) << nl;

        // This syntax looks even nicer:

        // GOOD: does not compile
        // > constLabelList[labelRange(23,5)] = 5;

        // Check correct overlaps
        longLabelList[labelRange(-10, 12)] = 200;
        longLabelList[{18,3}] = 100;
        longLabelList[{23,3}] = 400;
        // and complete misses
        longLabelList[{500,50}] = 100;

        // labelRange automatically suppresses -ve size -> nop
        longLabelList[{5,-5}] = 42;
        longLabelList[{21,100}] = 42;

        //Good: does not compile
        //> longLabelList[labelRange(20,50)] = constLabelList;

        //Good: does not compile
        // longLabelList[labelRange(20,50)] = fixedLabelList;

        Info<< "updated: " << constLabelList[labelRange(23,5)] << nl;
        Info<< "updated: " << flatOutput(longLabelList) << nl;

        //Nope: sort(longLabelList[labelRange(18,5)]);
        {
            // Instead
            UList<label> sub = longLabelList[labelRange(0, 8)];
            sort(sub);
        }
        Info<< "sub-sorted: " << flatOutput(longLabelList) << nl;

        // construct from a label-range
        labelRange range(25,15);

        labelList ident(range.begin(), range.end());
        Info<<"range-list (label)=" << ident << nl;

        List<scalar> sident(range.begin(), range.end());
        Info<<"range-list (scalar)=" << sident << nl;

        // Sub-ranges also work
        List<scalar> sident2(range(3), range(10));
        Info<<"range-list (scalar)=" << sident2 << nl;

        // VERY BAD IDEA: List<scalar> sident3(range(10), range(3));

        // This doesn't work, and don't know what it should do anyhow
        // List<vector> vident(range.begin(), range.end());
        // Info<<"range-list (vector)=" << vident << nl;

        // Even weird things like this
        List<scalar> sident4
        (
            labelRange().begin(),
            labelRange::identity(8).end()
        );
        Info<<"range-list (scalar)=" << sident4 << nl;
    }

    wordReList reLst;
    wordList wLst;
    stringList sLst;

    scalar xxx(-1);

    if (args.optionFound("flag"))
    {
        Info<<"-flag:" << args["flag"] << endl;
    }

    if (args.optionReadIfPresent<scalar>("float", xxx))
    {
        Info<<"read float " << xxx << endl;
    }

    if (args.optionFound("reList"))
    {
        reLst = args.optionReadList<wordRe>("reList");
    }

    if (args.optionFound("wordList"))
    {
        wLst = args.optionReadList<word>("wordList");
    }

    if (args.optionFound("stringList"))
    {
        sLst = args.optionReadList<string>("stringList");
    }

    Info<< nl
        << "-reList:     " << flatOutput(reLst) << nl
        << "-wordList:   " << flatOutput(wLst)  << nl
        << "-stringList: " << flatOutput(sLst)  << endl;

    return 0;
}

// ************************************************************************* //
