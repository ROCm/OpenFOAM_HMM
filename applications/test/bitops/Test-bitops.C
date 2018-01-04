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
    Test some bit-operations.

\*---------------------------------------------------------------------------*/

#include "bool.H"
#include "IOstreams.H"

#include <algorithm>
#include <type_traits>
#include <limits>

template<typename UIntType, UIntType v, unsigned int n>
struct set_lower_bits
:
    public std::integral_constant
    <
        UIntType,
        v | (v >> n) | set_lower_bits<UIntType, v | (v >> n), (n >> 1)>::value
    >
{};


template<typename UIntType, UIntType v>
struct set_lower_bits<UIntType, v, 1>
:
    public std::integral_constant<UIntType, v | (v >> 1)>
{};


template<typename UIntType, UIntType v>
struct pow2ceil
:
    public std::integral_constant
    <
        UIntType,
        set_lower_bits
        <
            UIntType,
            v - 1,
            (std::numeric_limits<UIntType>::digits >> 1)
        >::value + 1
    >
{};


// No quite
template<size_t N>
struct topbit : std::integral_constant<size_t, topbit<(N >> 1)>{} + 1> {};

template<> struct topbit<1> : std::integral_constant<size_t,1> {};
template<> struct topbit<0> : std::integral_constant<size_t,0> {};

using namespace Foam;


template<size_t N>
void printTopbit()
{
    std::cout
        <<"topbit<" << N << "> : " << topbit<N>::value << '\n';
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<<"pow2ceil<10>: "
        << pow2ceil<unsigned, 10>::value << nl;

    Info<<"pow2ceil<16>: "
        << pow2ceil<unsigned, 16>::value << nl;

    Info<<"pow2ceil<369>: "
        << pow2ceil<unsigned, 369>::value << nl;

    printTopbit<0>();
    printTopbit<1>();
    printTopbit<2>();
    printTopbit<3>();
    printTopbit<4>();
    printTopbit<5>();
    printTopbit<6>();
    printTopbit<7>();
    printTopbit<8>();
    printTopbit<9>();
    printTopbit<10>();
    printTopbit<11>();
    printTopbit<15>();
    printTopbit<16>();
    printTopbit<17>();
    printTopbit<18>();

    Info << "---\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
