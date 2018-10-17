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

Application
    Tuple2Test

Description

\*---------------------------------------------------------------------------*/

#include "Tuple2.H"
#include "label.H"
#include "scalar.H"
#include "List.H"
#include "ops.H"
#include <functional>

using namespace Foam;

// Test for special comparison operation using compareOp
// Normal sort on label, reverse sort on scalar
struct special1
{
    typedef Tuple2<label, scalar> type;

    bool operator()(const type& a, const type& b) const
    {
        int val = compareOp<label>()(a.first(), b.first());
        return (val == 0) ? (b.second() < a.second()) : (val < 0);
    }
};


// Test for special comparison operation using compareOp
// Normal sort on scalar, reverse sort on label
struct special2
{
    typedef Tuple2<label, scalar> type;

    bool operator()(const type& a, const type& b) const
    {
        scalar val = compareOp<scalar>()(a.second(), b.second());
        return (val == 0) ? (b.first() < a.first()) : (val < 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    typedef Tuple2<label, scalar> indexedScalar;

    indexedScalar t2(1, 3.2);

    Info<< "tuple: "
        << t2 << " "
        << t2.first() << " " << t2.second() << nl;

    // As list. Generated so that we have duplicate indices
    List<indexedScalar> list1(3*4);
    for (label i = 0; i < 4; ++i)
    {
        const label j = (i+1);
        const label idx = ((i % 2) ? -1 : 1) * (j);

        list1[i]   = indexedScalar(idx, (j*j));
        list1[i+4] = indexedScalar(idx, 2*j);    // duplicate index
        list1[i+8] = indexedScalar(idx+12, 2*j); // duplicate value
    }

    Info<< "Unsorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, std::less<indexedScalar>());

    Info<< "sorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, std::greater<indexedScalar>());

    Info<< "reverse sorted tuples:" << nl << list1 << nl;

    Foam::sort(list1, special1());

    Info<< "special sorted tuples - sort on index, reverse on value:"
        << nl << list1 << nl;

    Foam::sort(list1, special2());

    Info<< "special sorted tuples - sort on value, reverse on index:"
        << nl << list1 << nl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
