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

#include "hashedWordList.H"
#include "nil.H"
#include "HashSet.H"
#include "Map.H"
#include "labelPairHashes.H"
#include "FlatOutput.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    hashedWordList words
    {
        "abc",
        "def",
        "ghi"
    };
    words = { "def", "ghi", "xy", "all", "end", "all" };

    wordHashSet setA
    {
        "xx",
        "yy",
        "zz"
    };

    setA = { "kjhk", "kjhk2", "abced" };

    HashTable<label, word> tableA
    {
        { "value1", 1 },
        { "value2", 2 },
        { "value3", 3 }
    };

    HashTable<nil> tableB;
    tableB.insert("value4", nil());
    tableB.insert("value5", nil());
    tableB.insert("value6", nil());

    Info<< "tableA keys: "; tableA.writeKeys(Info) << endl;

    Info<< "tableB content: " << tableB << endl;

    auto keyIterPair = tableA.keys();
    for (const auto& i : keyIterPair)
    {
        Info<<" keys: " << i << endl;
    }

    Map<label> mapA
    {
        { 1, 1 },
        { 2, 2 },
        { 3, 3 }
    };
    mapA.set(4, 4);

    Info<< "hashedWordList: " << words << nl
        << "with lookup: "  << words.lookup() << endl;

    words.sort();
    Info<< "hashedWordList: " << words << nl
        << "with lookup: "  << words.lookup() << endl;

    words.uniq();
    Info<< "hashedWordList: " << words << nl
        << "with lookup: "  << words.lookup() << endl;

    {
        List<word> input = { "def", "ghi", "xy", "all", "end", "all", "def" };
        hashedWordList words1(input, true);

        Info<< "input word list: " << input << nl
            << "without dup: "  << words1 << endl;

        Info<< "from wordHashSet: " << hashedWordList(setA) << endl;
        Info<< "from HashTable: " << hashedWordList(tableA) << endl;
        Info<< "from HashTable: " << hashedWordList(tableB) << endl;

        // even this works
        Info<< "from hashSet: "
            << hashedWordList
               (
                   wordHashSet(setA)
                 | wordHashSet(tableA) | wordHashSet(tableB)
               ) << endl;
    }

    Info<< "wordHashSet: "    << setA << endl;
    Info<< "Table-HashSet: "  << tableA << endl;
    Info<< "Map<label>: "     << mapA << endl;

    Info<< "create from HashSet: ";
    Info<< wordHashSet(setA) << endl;
    Info<< "create from HashTable<T>: ";
    Info<< wordHashSet(tableA) << endl;
    Info<< "create from HashTable<nil>: ";
    Info<< wordHashSet(tableB) << endl;

    Info<< "create from Map<label>: ";
    Info<< labelHashSet(mapA) << endl;

    Info<<"combined toc: "
        << (wordHashSet(setA) | wordHashSet(tableA) | wordHashSet(tableB))
        << nl;

    labelHashSet setB
    {
        1, 11, 42
    };

    setB = FixedList<label, 4>({1, 2, 3, 4});
    setB = {1, 2, 4};
    setB = List<label>({1, 2, 4});
    Info<< "setB : " << setB << endl;

    labelPair pair(12, 15);
    setB.set(pair);

    Info<< "setB : " << setB << endl;
    setB.unset(pair);


    labelHashSet setC(1);
    setC.insert(2008);
    setC.insert(1984);

    Info<< "setC : " << setC << endl;

    labelHashSet setD(1);
    setD.insert({11, 100, 49, 36, 2008});

    Info<< "setD : " << flatOutput(setD) << endl;

    Info<< "setB == setC: " << (setB == setC) << endl;
    Info<< "setC != setD: " << (setC != setD) << endl;

    // test operations
    setB += setC;
    Info<< "setB += setC : " << setB << endl;

    setB &= setD;
    Info<< "setB &= setD : " << setB << endl;

    Info<< "setB : " << setB << endl;
    Info<< "setC : " << setC << endl;
    Info<< "setD : " << setD << endl;
    Info<< "setB ^ setC ^ setD : " << (setB ^ setC ^ setD) << endl;

    // test operator[]

    Info<< "setD : " << setD << endl;
    if (setD[0])
    {
        Info<< "setD has 0" << endl;
    }
    else
    {
        Info<< "setD has no 0" << endl;
    }


    if (setD[11])
    {
        Info<< "setD has 11" << endl;
    }
    else
    {
        Info<< "setD has no 11" << endl;
    }

    Info<< "setB : " << flatOutput(setB) << endl;
    Info<< "setD : " << flatOutput(setD) << endl;

    setD -= setB;
    Info<< "setD -= setB : " << flatOutput(setD) << endl;

    // This should not work (yet?)
    // setD[12] = true;

    List<label> someLst(10);
    forAll(someLst, elemI)
    {
        someLst[elemI] = elemI*elemI;
    }

    label added = setD.set(someLst);
    Info<< "added " << added << " from " << someLst.size() << endl;
    Info<< "setD : " << flatOutput(setD) << endl;

    Info<< "setD for-range()" << nl;

    for (auto i : setD)
    {
        Info << i << endl;
    }

    Info<< nl << "Test swapping, moving etc." << nl;
    setA.insert({ "some", "more", "entries" });

    Info<< "input" << nl;
    Info<< "setA: " << setA << nl;

    wordHashSet setA1(std::move(setA));

    Info<< "move construct" << nl;
    Info<< "setA: " << setA << nl
        << "setA1: " << setA1 << nl;


    wordHashSet setB1;
    Info<< "move assign" << nl;
    setB1 = std::move(setA1);

    Info<< "setA1: " << setA1 << nl
        << "setB1: " << setB1 << nl;

    setB1.swap(setA1);

    Info<< "setA1: " << setA1 << nl
        << "setB1: " << setB1 << nl;

    return 0;
}


// ************************************************************************* //
