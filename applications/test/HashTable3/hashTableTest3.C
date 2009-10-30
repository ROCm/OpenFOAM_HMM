/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

    Test speeds for some HashTable operations

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "HashTable.H"
#include "HashPtrTable.H"
#include "Map.H"
#include "StaticHashTable.H"
#include "HashTbl.H"
#include "cpuTime.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    const label nLoops = 30;
    const label nBase  = 100000;
    const label nSize  = nLoops * nBase;

    cpuTime timer;

    // ie, a
    // Map<label> map(2 * nSize);
    // HashTable<label, label, Hash<label> > map(2 * nSize);
    // StaticHashTable<label, label, Hash<label> > map(2 * nSize);
    HashTbl<label, label, Hash<label> > map(2 * nSize);

    Info<< "Constructed map of size: " << nSize
        << "  " << timer.cpuTimeIncrement() << " s\n\n";

    for (label i = 0; i < nSize; i++)
    {
        map.insert(i, i);
    }
    Info<< "Inserted " << nSize << " elements: "
        << timer.cpuTimeIncrement() << " s\n\n";


    label elemI = 0;
    for (label iLoop = 0; iLoop < nLoops; iLoop++)
    {
        for (label i = 0; i < nBase; i++)
        {
            map.erase(elemI++);
        }
        Info<< "loop " << iLoop << " - Erased " << nBase << " elements: "
            << "  " << timer.cpuTimeIncrement() << " s\n";
    }

    return 0;
}

// ************************************************************************* //
