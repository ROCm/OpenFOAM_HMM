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
#include "IStringStream.H"
#include "scalar.H"
#include "vector.H"

#include "labelRange.H"
#include "ListOps.H"
#include "SubList.H"

#include <list>

using namespace Foam;


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
        vector(6, 7, 8)
    };
    Info<< "list2: " << list2 << endl;

    list1.append(list2);
    Info<< "list1.append(list2): " << list1 << endl;

    Info<< findIndex(list2, vector(3, 4, 5)) << endl;

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

        Info<< "labels (contiguous=" << contiguous<label>() << ")" << nl;

        Info<< "normal: " << longLabelList << nl;
        Info<< "flatOutput: " << flatOutput(longLabelList) << nl;
        // Info<< "flatOutput(14): " << flatOutput(longLabelList, 14) << nl;
        // Info<< "flatOutput(15): " << flatOutput(longLabelList, 15) << nl;

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
