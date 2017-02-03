/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "boundBox.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned Size>
Foam::boundBox::boundBox
(
    const UList<point>& points,
    const FixedList<label, Size>& indices,
    bool doReduce
)
:
    min_(invertedBox.min()),
    max_(invertedBox.max())
{
    add(points, indices);

    if (doReduce)
    {
        reduce();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<unsigned Size>
void Foam::boundBox::add
(
    const FixedList<point, Size>& points
)
{
    // a FixedList is never empty
    for (unsigned i=0; i < Size; ++i)
    {
        add(points[i]);
    }
}


template<unsigned Size>
void Foam::boundBox::add
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
)
{
    // points may be empty, but a FixedList is never empty
    if (!points.empty())
    {
        for (unsigned i=0; i < Size; ++i)
        {
            add(points[indices[i]]);
        }
    }
}


template<unsigned Size>
bool Foam::boundBox::contains
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
) const
{
    // points may be empty, but a FixedList is never empty
    if (points.empty())
    {
        return false;
    }

    forAll(indices, i)
    {
        if (!contains(points[indices[i]]))
        {
            return false;
        }
    }

    return true;
}


template<unsigned Size>
bool Foam::boundBox::containsAny
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
) const
{
    // points may be empty, but a FixedList is never empty
    if (points.empty())
    {
        return false;
    }

    forAll(indices, i)
    {
        if (contains(points[indices[i]]))
        {
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
