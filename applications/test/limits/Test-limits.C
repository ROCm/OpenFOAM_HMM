/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Print max limits.

\*---------------------------------------------------------------------------*/

#include <limits>
#include "int.H"
#include "uint.H"
#include "string.H"
#include "scalar.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    //NONE Info<<"int16:" << pTraits<int16_t>::max << nl;
    Info<< "=max=" << nl;
    Info<< "uint8:" << pTraits<uint8_t>::max << nl;
    Info<< "int16:" << std::numeric_limits<int16_t>::max() << nl;
    Info<< "int32:" << pTraits<int32_t>::max << nl;
    Info<< "uint32:" << pTraits<uint32_t>::max << nl;
    Info<< "int64:" << pTraits<int64_t>::max << nl;
    Info<< "uint64:" << pTraits<uint64_t>::max << nl;

    Info<< nl;

    cout<< "int16:" << std::numeric_limits<int16_t>::max() << nl;
    cout<< "int32:" << pTraits<int32_t>::max << nl;
    cout<< "uint32:" << pTraits<uint32_t>::max << nl;
    cout<< "int64:" << pTraits<int64_t>::max << nl;
    cout<< "uint64:" << pTraits<uint64_t>::max << nl;

    Info<< nl << "=digits=" << nl;

    cout<< "int16:" << std::numeric_limits<int16_t>::digits << nl;
    cout<< "int32:" << std::numeric_limits<int32_t>::digits << nl;
    cout<< "uint32:" << std::numeric_limits<uint32_t>::digits << nl;
    cout<< "int64:" << std::numeric_limits<int64_t>::digits << nl;
    cout<< "uint64:" << std::numeric_limits<uint64_t>::digits << nl;
    cout<< "float:"  << std::numeric_limits<float>::digits << nl;
    cout<< "double:" << std::numeric_limits<double>::digits << nl;
    cout<< "long double:" << std::numeric_limits<long double>::digits << nl;

    Info<< nl << "=float=" << nl;

    cout<< "max:" << std::numeric_limits<float>::max()
        << " VGREAT:" << floatScalarVGREAT << nl;
    cout<< "min:" << std::numeric_limits<float>::min()
        << " VSMALL:" << floatScalarVSMALL << nl;
    cout<< "epsilon:" << std::numeric_limits<float>::epsilon()
        << " SMALL:" << floatScalarSMALL << nl;
    cout<< "1/epsilon:" << 1.0f/std::numeric_limits<float>::epsilon()
        << " GREAT:" << floatScalarGREAT << nl;

    Info<< nl << "=double=" << nl;

    cout<< "max:" << std::numeric_limits<double>::max()
        << " VGREAT:" << doubleScalarVGREAT << nl;
    cout<< "min:" << std::numeric_limits<double>::min()
        << " VSMALL:" << doubleScalarVSMALL << nl;
    cout<< "epsilon:" << std::numeric_limits<double>::epsilon()
        << " SMALL:" << doubleScalarSMALL << nl;
    cout<< "1/epsilon:" << 1.0f/std::numeric_limits<double>::epsilon()
        << " GREAT:" << doubleScalarGREAT << nl;

    Info << "---\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
