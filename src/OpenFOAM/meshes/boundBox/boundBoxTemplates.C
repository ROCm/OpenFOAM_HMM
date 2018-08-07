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
    boundBox()
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
    for (const point& p : points)
    {
        add(p);
    }
}


template<unsigned Size>
void Foam::boundBox::add
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
)
{
    const label len = points.size();

    // Skip if points is empty
    if (len)
    {
        for (const label pointi : indices)
        {
            if (pointi >= 0 && pointi < len)
            {
                add(points[pointi]);
            }
        }
    }
}


template<class IntContainer>
void Foam::boundBox::add
(
    const UList<point>& points,
    const IntContainer& indices
)
{
    const label len = points.size();

    // Skip if points is empty
    if (len)
    {
        for (const label pointi : indices)
        {
            if (pointi >= 0 && pointi < len)
            {
                add(points[pointi]);
            }
        }
    }
}


template<unsigned Size>
inline bool Foam::boundBox::contains
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
) const
{
    const label len = points.size();

    if (!len)
    {
        return true;
    }

    for (const label pointi : indices)
    {
        if (pointi >= 0 && pointi < len)
        {
            if (!contains(points[pointi]))
            {
                return false;
            }
        }
    }

    return true;
}


template<class IntContainer>
inline bool Foam::boundBox::contains
(
    const UList<point>& points,
    const IntContainer& indices
) const
{
    const label len = points.size();

    if (!len)
    {
        return true;
    }

    for (const label pointi : indices)
    {
        if (pointi >= 0 && pointi < len)
        {
            if (!contains(points[pointi]))
            {
                return false;
            }
        }
    }

    return true;
}


template<unsigned Size>
inline bool Foam::boundBox::containsAny
(
    const UList<point>& points,
    const FixedList<label, Size>& indices
) const
{
    const label len = points.size();

    if (!len)
    {
        return true;
    }

    label failed = 0;

    for (const label pointi : indices)
    {
        if (pointi >= 0 && pointi < len)
        {
            if (contains(points[pointi]))
            {
                return true;
            }

            ++failed;
        }
    }

    return !failed;
}


template<class IntContainer>
inline bool Foam::boundBox::containsAny
(
    const UList<point>& points,
    const IntContainer& indices
) const
{
    const label len = points.size();

    if (!len)
    {
        return true;
    }

    label failed = 0;

    for (const label pointi : indices)
    {
        if (pointi >= 0 && pointi < len)
        {
            if (contains(points[pointi]))
            {
                return true;
            }

            ++failed;
        }
    }

    return !failed;
}


// ************************************************************************* //
