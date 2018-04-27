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
    Test-contiguous

Description
    Simple test of contiguous data

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "wordRes.H"
#include "contiguous.H"

#include "IOstreams.H"
#include "scalar.H"
#include "vector.H"

#include "labelRange.H"
#include "scalarList.H"
#include "HashOps.H"
#include "FixedList.H"
#include "Pair.H"

namespace Foam
{

// Wrong, but interesting to test
// template<> struct contiguous<Pair<word>> : std::true_type {};

} // End namespace Foam

using namespace Foam;


template<class T>
void printContiguous()
{
    Info<<"contiguous " << typeid(T).name() << " () = "
        << contiguous<T>()
        // << " value = " << contiguous<T>::value
        << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();

    #include "setRootCase.H"

    printContiguous<label>();
    printContiguous<double>();
    printContiguous<FixedList<int, 2>>();
    printContiguous<FixedList<int, 3>>();
    printContiguous<Pair<long>>();

    printContiguous<FixedList<word, 2>>();
    printContiguous<Pair<word>>();

    printContiguous<FixedList<FixedList<int, 2>, 2>>();

    return 0;
}

// ************************************************************************* //
