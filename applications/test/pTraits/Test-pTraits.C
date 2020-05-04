/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "IOstreams.H"
#include "pTraits.H"
#include "vector.H"
#include "tensor.H"
#include "uLabel.H"

#include <type_traits>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

template<class T>
void printTraits()
{
    Info<< pTraits<T>::typeName
        << ": zero=" << pTraits<T>::zero
        << " one=" << pTraits<T>::one
        << " integral=" << std::is_integral<T>::value
        << " floating=" << std::is_floating_point<T>::value
        << endl;
}


template<class T>
void printTraits(const pTraits<T>& p)
{
    Info<< p.typeName << " == " << p << endl;
}


#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#pragma GCC diagnostic warning "-Wuninitialized"

int main()
{
    printTraits<bool>();
    printTraits<label>();
    printTraits<scalar>();
    printTraits<vector>();
    printTraits<tensor>();

    {
        pTraits<bool> b(true);
        printTraits(b);
    }

    {
        pTraits<label> l(100);
        printTraits(l);
    }

    printTraits(pTraits<scalar>(3.14159));

    label abc;
    Info<< "uninitialized primitive:"<< abc << endl;

    label def = label();
    Info<< "initialized primitive:"<< def << endl;

    Info<< nl << "some interesting label limits:" << nl;
    std::cout<< "sizeof = " << sizeof(label) << nl;
    std::cout<< "min = " << pTraits<label>::min << nl;
    std::cout<< "max = " << pTraits<label>::max << nl;
    std::cout<< "umax = " << pTraits<uLabel>::max << nl;

    std::cout<< "max_2 = " << pTraits<label>::max/2 << " <=> "
        << (1L << (sizeof(label)*8-2)) << nl;

    std::cout<< "max_4 = " << pTraits<label>::max/4 << " <=> "
        << (1L << (sizeof(label)*8-3)) << nl;

    std::cout<< "max_8 = " << pTraits<label>::max/8 << " <=> "
        << (1L << (sizeof(label)*8-4)) << nl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
