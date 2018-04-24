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
    Test the sizeof for basic types. Can be compiled and run without
    any OpenFOAM libraries.

\*---------------------------------------------------------------------------*/

#include <cstdint>
#include <climits>
#include <cstdlib>
#include <limits>
#include <iostream>

// Can also compile without OpenFOAM
#ifdef WM_LABEL_SIZE
    #include "IOstreams.H"
#endif

template<class T>
void print(const char* msg)
{
    std::cout<< msg << ' ' << sizeof(T) << '\n';
}


#ifdef WM_LABEL_SIZE
template<class T>
void printMax(const char* msg)
{
    std::cout<< msg << ' ' << sizeof(T)
        << " max "
        << std::numeric_limits<T>::max() << '\n';

    Foam::Info<< msg << ' ' << sizeof(T)
        << " max "
        << long(std::numeric_limits<T>::max()) << '\n';
}
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    std::cout<<"machine sizes\n---\n\n";

    print<mode_t>("mode_t");
    print<short>("short");
    print<int>("int");
    print<long>("long");
    print<long long>("long long");
    print<unsigned short>("ushort");
    print<unsigned int>("uint");
    print<unsigned long>("ulong");
    print<unsigned long long>("ulong-long");
    print<int16_t>("int16");
    print<int32_t>("int32");
    print<int64_t>("int64");
    print<uint16_t>("uint16");
    print<uint32_t>("uint32");
    print<uint64_t>("uint64");
    print<float>("float");
    print<double>("double");
    print<std::string>("std::string");
    print<std::string::size_type>("std::string::size_type");

    #ifdef WM_LABEL_SIZE
    std::cout<<"\nmax values\n---\n\n";

    printMax<mode_t>("mode_t");
    Foam::Info<< "mode_t 0777: " << mode_t(0777) << '\n';

    printMax<short>("short");
    printMax<int>("int");
    printMax<long>("long");
    printMax<unsigned short>("ushort");
    printMax<unsigned int>("uint");
    printMax<unsigned long>("ulong");
    printMax<float>("float");
    printMax<double>("double");
    #endif

    std::cout << "\n---\nEnd\n\n";

    return 0;
}


// ************************************************************************* //
