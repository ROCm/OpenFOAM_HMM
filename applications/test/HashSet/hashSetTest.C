/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "HashSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    HashSet<string> setA(0);

    setA.insert("kjhk");
    setA.insert("kjhk2");

    Info<< setA << endl;

    labelHashSet setB(1);
    setB.insert(11);
    setB.insert(42);

    Info<< "setB : " << setB << endl;

    labelHashSet setC(1);
    setC.insert(2008);
    setC.insert(1984);

    Info<< "setC : " << setC << endl;

    labelHashSet setD(1);
    setD.insert(11);
    setD.insert(100);
    setD.insert(2008);

    Info<< "setD : " << setD << endl;

    Info<< "setB == setC: " << (setB == setC) << endl;
    Info<< "setC != setD: " << (setC != setD) << endl;

    // test operations
    setB += setC;
    Info<< "setB += setC : " << setB << endl;

    setB &= setD;
    Info<< "setB &= setD : " << setB << endl;

    setB += setC;
    setB -= setD;
    Info<< "setB += setC -= setD : " << setB << endl;

    return 0;
}


// ************************************************************************* //
