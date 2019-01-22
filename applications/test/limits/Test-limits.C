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

Description
    Print max limits.

\*---------------------------------------------------------------------------*/

#include <limits>
#include "int.H"
#include "uint.H"
#include "string.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    //NONE Info<<"int16:" << pTraits<int16_t>::max << nl;
    Info<<"=max=" << nl;
    Info<<"uint8:" << std::numeric_limits<uint8_t>::max() << nl;
    Info<<"int16:" << std::numeric_limits<int16_t>::max() << nl;
    Info<<"int32:" << pTraits<int32_t>::max << nl;
    Info<<"uint32:" << pTraits<uint32_t>::max << nl;
    Info<<"int64:" << pTraits<int64_t>::max << nl;
    Info<<"uint64:" << pTraits<uint64_t>::max << nl;

    Info<< nl;

    cout<<"int16:" << std::numeric_limits<int16_t>::max() << nl;
    cout<<"int32:" << pTraits<int32_t>::max << nl;
    cout<<"uint32:" << pTraits<uint32_t>::max << nl;
    cout<<"int64:" << pTraits<int64_t>::max << nl;
    cout<<"uint64:" << pTraits<uint64_t>::max << nl;

    Info<< nl << "=digits=" << nl;

    cout<<"int16:" << std::numeric_limits<int16_t>::digits << nl;
    cout<<"int32:" << std::numeric_limits<int32_t>::digits << nl;
    cout<<"uint32:" << std::numeric_limits<uint32_t>::digits << nl;
    cout<<"int64:" << std::numeric_limits<int64_t>::digits << nl;
    cout<<"uint64:" << std::numeric_limits<uint64_t>::digits << nl;
    cout<<"float:"  << std::numeric_limits<float>::digits << nl;
    cout<<"double:" << std::numeric_limits<double>::digits << nl;
    cout<<"long double:" << std::numeric_limits<long double>::digits << nl;

    Info << "---\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
