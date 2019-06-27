/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    Test slice range
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelList.H"
#include "FixedList.H"
#include "sliceRange.H"
#include "SliceList.H"
#include "IndirectList.H"
#include "Random.H"

using namespace Foam;


typedef FixedList<label, 3> sliceCoeffs;

void printInfo(const sliceCoeffs& coeffs)
{
    sliceRange range(coeffs);

    Info<< nl
        << "coeffs: " << coeffs << nl
        << "range: " << range << nl
        << "first: " << range.first() << nl
        << "*begin  " << *range.begin() << nl
        << "last:  " << range.last() << nl
        << "*end    " << *range.end() << nl
        << "values: " << flatOutput(range.labels()) << nl;

    Info<< "for  :";
    for (const label val : range)
    {
        Info<< ' ' << val;
    }
    Info<< nl;

    if (range.empty())
    {
        Info<< "empty"<< nl;
        return;
    }

    Info<< "mid-point [" << (range.size()/2) << "] = "
        << range[range.size()/2] << nl;

    const auto endIter = range.cend();

    for (const label i : {label(-1), (range.size()/2), range.size()})
    {
        const auto iter = range.at(i);

        Info<< "at(" << i << ") = " << *iter;

        if (iter == endIter)
        {
            Info<< " (out-of-range)";
        }

        Info<< nl;
    }

    // Copy construct
    sliceRange range2(range);

    // Copy assign
    range2 = range;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();

    argList args(argc, argv, false, true);


    for
    (
        sliceCoeffs coeffs :
        {
            sliceCoeffs{25, 8, 2},
            sliceCoeffs{15, 5, 3},
            sliceCoeffs{15, -5, 2},
        }
    )
    {
        printInfo(coeffs);
    }


    // Generator
    {
        sliceRange range(25, 8, 3);

        Info<< "Generator for " << range << nl;

        auto gen = range.generator();

        for (label i=0; i < 10; ++i)
        {
            Info<< "    " << gen() << nl;
        }
    }

    // Sliced lists
    {
        List<scalar> list1(100);

        Random rnd(1234);

        for (scalar& val : list1)
        {
            val = 100 * rnd.sample01<scalar>();
        }

        Info<< nl << "Random list: " << flatOutput(list1) << nl;

        SliceList<scalar> slice1(list1, sliceRange(0, 15, 3));

        Info<< nl << "slicing with: " << slice1.addressing() << nl;

        Info<< nl << "Sliced list: " << flatOutput(slice1) << nl;

        for (scalar& val : slice1)
        {
            val = -val;
        }

        // Changed list via slice:
        Info<< nl << "Changed via slice: " << flatOutput(list1) << nl;

        // Some indirect list

        IndirectList<scalar> indlist
        (
            list1,
            identity(slice1.size(), list1.size()-slice1.size())
        );

        Info<< nl << "Indirect slice: " << flatOutput(indlist) << nl;

        indlist = 1000;
        Info<< nl << "zeroed slice: " << flatOutput(indlist) << nl;

        slice1 = indlist;

        Info<< nl << "self-copy: " << flatOutput(list1) << nl;

        slice1 = 100000;

        Info<< nl << "set values: " << flatOutput(slice1) << nl
            << " = " << flatOutput(list1) << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
