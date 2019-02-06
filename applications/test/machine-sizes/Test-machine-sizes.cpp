/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011 OpenFOAM Foundation
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
    Test the sizeof for basic types.
    Can be compiled and run without any OpenFOAM libraries.

        g++ -std=c++11 -oTest-machine-sizes Test-machine-sizes.cpp

\*---------------------------------------------------------------------------*/

#include <cstdint>
#include <climits>
#include <cstdlib>
#include <iostream>

template<class T>
void print(const char* msg)
{
    std::cout<< msg << ' ' << sizeof(T) << '\n';
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    std::cout<<"machine sizes\n---\n\n";

    print<short>("short");
    print<int>("int");
    print<long>("long");
    print<long long>("long long");
    print<float>("float");
    print<double>("double");
    print<long double>("long double");
    print<std::string>("std::string");
    print<std::string::size_type>("std::string::size_type");

    std::cout << "\n---\nEnd\n\n";

    return 0;
}


// ************************************************************************* //
