/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "uLabel.H"
#include "IOstreams.H"
#include "PackedBoolList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    PackedBoolList list1(20);
    // set every other one on
    forAll(list1, i)
    {
        list1[i] = i % 2;
    }

    Info<< "\nalternative bit pattern\n";
    list1.print(Info, true);

    PackedBoolList list2 = ~list1;

    Info<< "\ncomplementary bit pattern\n";
    list2.print(Info, true);

    list2.resize(24, true);
    list2.resize(28, false);
    for (label i=0; i < 4; ++i)
    {
        list2[i] = true;
    }

    Info<< "\nresized with 4 true + 4 false, bottom 4 bits true\n";
    list2.print(Info, true);

    labelList list2Labels = list2.used();

    Info<< "\noperator|\n";
    (list1 | list2).print(Info, true);

    Info<< "\noperator& : does trim\n";
    (list1 & list2).print(Info, true);

    Info<< "\noperator^\n";
    (list1 ^ list2).print(Info, true);


    Info<< "\noperator|=\n";
    {
        PackedBoolList list3 = list1;
        (list3 |= list2).print(Info, true);
    }

    Info<< "\noperator|= with UList<label>\n";
    {
        PackedBoolList list3 = list1;
        (list3 |= list2Labels).print(Info, true);
    }

    Info<< "\noperator&=\n";
    {
        PackedBoolList list3 = list1;
        (list3 &= list2).print(Info, true);
    }

    Info<< "\noperator+=\n";
    {
        PackedBoolList list3 = list1;
        (list3 += list2).print(Info, true);
    }

    Info<< "\noperator+= with UList<label>\n";
    {
        PackedBoolList list3 = list1;
        (list3 += list2Labels).print(Info, true);
    }

    Info<< "\noperator-=\n";
    {
        PackedBoolList list3 = list1;
        (list3 -= list2).print(Info, true);
    }

    Info<< "\noperator-= with UList<label>\n";
    {
        PackedBoolList list3 = list1;
        (list3 -= list2Labels).print(Info, true);
    }

    return 0;
}


// ************************************************************************* //
