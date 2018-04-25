/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Test-bitSet2

Description
    Test bitSet functionality

\*---------------------------------------------------------------------------*/

#include "uLabel.H"
#include "boolList.H"
#include "DynamicList.H"
#include "IOstreams.H"
#include "ITstream.H"
#include "StringStream.H"
#include "bitSet.H"
#include "FlatOutput.H"

using namespace Foam;

inline Ostream& report
(
    const bitSet& bitset,
    bool showBits = false,
    bool debugOutput = false
)
{
    Info<< "size=" << bitset.size() << "/" << bitset.capacity()
        << " count=" << bitset.count()
        << " all:" << bitset.all()
        << " any:" << bitset.any()
        << " none:" << bitset.none() << nl;

    Info<< "values: " << flatOutput(bitset) << nl;
    if (showBits)
    {
        bitset.printBits(Info, debugOutput) << nl;
    }

    return Info;
}


template<class UIntType>
std::string toString(UIntType value, char off='.', char on='1')
{
    std::string str(std::numeric_limits<UIntType>::digits, off);

    unsigned n = 0;

    // Starting from most significant bit - makes for easy reading.
    for
    (
        unsigned test = (1u << (std::numeric_limits<UIntType>::digits-1));
        test;
        test >>= 1u
    )
    {
        str[n++] = ((value & test) ? on : off);
    }

    return str;
}


inline bool compare
(
    const bitSet& bitset,
    const std::string& expected
)
{
    const List<unsigned int>& store = bitset.storage();

    std::string has;

    for (label blocki=0; blocki < bitset.nBlocks(); ++blocki)
    {
        has += toString(store[blocki]);
    }

    if (has == expected)
    {
        Info<< "pass: " << has << nl;
        return true;
    }

    Info<< "fail:   " << has << nl;
    Info<< "expect: " << expected << nl;

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    bitSet list1(22);
    // Set every third one on
    forAll(list1, i)
    {
        list1[i] = !(i % 3);
    }
    Info<< "\nalternating bit pattern\n";
    compare(list1, "..........1..1..1..1..1..1..1..1");

    list1.unset(labelRange(13, 20));  // In range

    Info<< "\nafter clear [13,..]\n";
    compare(list1, "...................1..1..1..1..1");

    report(list1, true);

    list1.unset(labelRange(40, 20));  // out of range
    Info<< "\nafter clear [40,..]\n";
    compare(list1, "...................1..1..1..1..1");

    report(list1, true);


    Info<< "first: " << list1.find_first()
        << " last: " << list1.find_last() << endl;

    Info<< "iterate through:";
    for (const label idx : list1)
    {
        Info<<" " << idx;
    }
    Info<< nl;

    Info<< "\nalternating bit pattern\n";
    report(list1, true);

    bitSet list2 = ~list1;

    Info<< "\nflipped bit pattern\n";
    report(list2, true);

    // set every other on
    forAll(list2, i)
    {
        list2[i] = !(i % 2);
    }

    Info<< "\nstarting pattern\n";
    report(list2, true);

    list2.resize(28, false);
    list2.resize(34, true);
    list2.resize(40, false);
    for (label i=0; i < 4; ++i)
    {
        list2[i] = true;
    }

    Info<< "\nresized with false, [28,34) true + 6 false, bottom 4 bits true\n";

    compare
    (
        list2,
        "1111.......1.1.1.1.1.1.1.1.11111"
        "..............................11"
    );


    report(list2, true);

    labelList list2Labels = list2.toc();

    Info<< "\noperator|\n";

    list1.printBits(Info);
    list2.printBits(Info);
    Info<< "==\n";
    (list1 | list2).printBits(Info);

    Info<< "\noperator& : does trim\n";
    report((list1 & list2), true);

    Info<< "\noperator^\n";
    report((list1 ^ list2), true);

    Info<< "\noperator|=\n";
    {
        bitSet list3 = list1;
        report((list3 |= list2), true);
    }

    Info<< "\noperator&=\n";
    {
        bitSet list3 = list1;
        report((list3 &= list2), true);
    }

    Info<< "\noperator^=\n";
    {
        bitSet list3 = list1;
        report((list3 ^= list2), true);
    }

    Info<< "\noperator-=\n";
    {
        bitSet list3 = list1;
        report((list3 -= list2), true);
    }

    bitSet list4
    (
        ITstream
        (
            "input",
            "(1 n 1 n 1 n 1 1 off 0 0 f f 0 y yes y true y false on t)"
        )()
    );

    Info<< "\ntest Istream constructor\n";

    report(list4, true);

    Info<< "\nclear/assign from labelList\n";
    list4.clear();
    list4.setMany(labelList{0, 1, 2, 3, 12, 13, 14, 19, 20, 21});

    report(list4, true);

    // Not yet:
    // bitSet list5{0, 1, 2, 3, 12, 13, 14, 19, 20, 21};
    // list5.printInfo(Info, true);
    // Info<< list5 << " indices: " << list5.toc() << nl;

    Info<< "\nassign from indices\n";
    list4.read
    (
        IStringStream
        (
            "{0 1 2 3 12 13 14 19 20 21}"
        )()
    );

    report(list4, true);
    compare(list4, "..........111....111........1111");

    list4.set(labelRange(28, 6));  // extends size

    Info<<"extended\n";
    compare
    (
        list4,
        "1111......111....111........1111"
        "..............................11"
    );

    list4.set(labelRange(40, 6));  // extends size
    Info<<"extended\n";
    compare
    (
        list4,
        "1111......111....111........1111"
        "..................111111......11"
    );

    list4.unset(labelRange(14, 19));
    Info<<"cleared [14,33)\n";
    compare
    (
        list4,
        "..................11........1111"
        "..................111111......1."
    );


    // Construct from labelUList, labelUIndList
    {
        DynamicList<label> indices({10, 50, 300});

        Info<< "set: " << flatOutput(indices) << endl;

        bitSet bools1(indices);

        Info<< "used: " << bools1.size() << " "
            << flatOutput(bools1.toc()) << endl;
    }

    return 0;
}


// ************************************************************************* //
