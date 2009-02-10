/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    testSHA1

Description


\*---------------------------------------------------------------------------*/

#include "SHA1.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char * argv[])
{
    SHA1 sha;
    SHA1::Digest shaDig;

    std::string str("The quick brown fox jumps over the lazy dog");
    Info<< shaDig << nl;
    Info<< SHA1("The quick brown fox jumps over the lazy dog") << nl;

    sha.append("The quick brown fox jumps over the lazy dog");
    Info<< sha << nl;

    sha.clear();
    sha.append("The quick brown fox jumps over the lazy dog");
    shaDig = sha;

    sha.append("\n");
    Info<< sha << nl;
    Info<< shaDig << nl;

    if (sha == shaDig)
    {
        Info<<"SHA1 digests are identical\n";
    }
    else
    {
        Info<<"SHA1 digests are different\n";
    }
    Info<<"lhs:" << sha << " rhs:" << shaDig << endl;
    
    // start over:
    sha.clear();
    sha.append(str);

    SHA1::Digest shaDig_A = sha;

    SHA1 sha_A = sha;

    sha.append("\n");

    Info<< "digest1: " << sha_A << nl;
    Info<< "digest2: " << sha << nl;

    
    return 0;
}
