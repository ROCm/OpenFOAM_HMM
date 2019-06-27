/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "PDRblock.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::PDRblock::location::findCell(const scalar p) const
{
    if (scalarList::empty())
    {
        return -1;
    }
    else if (equal(p, first()))
    {
        return 0;
    }
    else if (equal(p, last()))
    {
        return nCells()-1;
    }
    else if (p < first() || p > last())
    {
        // The point is out-of-bounds
        return -1;
    }

    // Binary search, finds lower index and thus corresponds to the
    // cell in which the point is found
    return findLower(*this, p);
}


Foam::label Foam::PDRblock::location::findIndex
(
    const scalar p,
    const scalar tol
) const
{
    if (scalarList::empty())
    {
        return -1;
    }
    else if (Foam::mag(p - first()) <= tol)
    {
        return 0;
    }
    else if (Foam::mag(p - last()) <= tol)
    {
        return scalarList::size()-1;
    }
    else if (p < first() || p > last())
    {
        // The point is out-of-bounds
        return -1;
    }

    // Linear search
    label i = 0;
    scalar delta = GREAT;

    for (const scalar& val : *this)
    {
        const scalar diff = mag(p - val);

        if (diff <= tol)
        {
            return i;
        }
        else if (delta < diff)
        {
            // Moving further away
            break;
        }

        delta = diff;
        ++i;
    }

    // This point is within bounds, but not near a grid-point
    return -2;
}


// ************************************************************************* //
