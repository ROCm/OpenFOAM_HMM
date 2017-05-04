/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "StaticHashTable.H"
#include "uLabel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(StaticHashTableCore, 0);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::StaticHashTableCore::canonicalSize(const label size)
{
    if (size < 1)
    {
        return 0;
    }

    // Enforce power of two - makes for a vey fast modulus etc.
    // The value '8' is some arbitrary lower limit.
    // If the hash table is too small, there will be many table collisions!

    const uLabel unsigned_size = size;
    uLabel powerOfTwo = 8;

    if (size < powerOfTwo)
    {
        return powerOfTwo;
    }
    else if (unsigned_size & (unsigned_size-1))  // <- Modulus of i^2
    {
        // Determine power-of-two. Brute-force is fast enough.
        while (powerOfTwo < unsigned_size)
        {
            powerOfTwo <<= 1;
        }

        return powerOfTwo;
    }
    else
    {
        return unsigned_size;
    }
}


// ************************************************************************* //
