/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::PackedBoolList::PackedBoolList(const label size)
:
    PackedList<1>(size)
{}


inline Foam::PackedBoolList::PackedBoolList
(
    const label size,
    const bool val
)
:
    PackedList<1>(size, (val ? 1u : 0u))
{}


inline Foam::PackedBoolList::PackedBoolList(const PackedBoolList& lst)
:
    PackedList<1>(lst)
{}


inline Foam::PackedBoolList::PackedBoolList(PackedBoolList&& lst)
:
    PackedList<1>()
{
    transfer(lst);
}


inline Foam::PackedBoolList::PackedBoolList(const UList<bool>& lst)
:
    PackedList<1>(lst.size())
{
    // Set according to indices that are true
    const label len = lst.size();

    for (label i = 0; i < len; ++i)
    {
        if (lst[i])
        {
            this->set(i, 1u);
        }
    }
}


inline Foam::PackedBoolList::PackedBoolList(const labelUList& indices)
:
    PackedBoolList((indices.size() ? indices.last() : 0), indices)
{}


inline Foam::PackedBoolList::PackedBoolList(const labelUIndList& indices)
:
    PackedBoolList((indices.size() ? indices.last() : 0), indices)
{}


inline Foam::PackedBoolList::PackedBoolList
(
    const label size,
    const labelUList& indices
)
:
    PackedList<1>(size)
{
    setMany(indices);
}


inline Foam::PackedBoolList::PackedBoolList
(
    const label size,
    const labelUIndList& indices
)
:
    PackedList<1>(size)
{
    setMany(indices);
}


inline Foam::autoPtr<Foam::PackedBoolList>
Foam::PackedBoolList::clone() const
{
    return autoPtr<PackedBoolList>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::PackedBoolList::test(const label pos) const
{
    return get(pos);
}


inline void Foam::PackedBoolList::transfer(PackedBoolList& lst)
{
    PackedList<1>::transfer(static_cast<PackedList<1>&>(lst));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::PackedBoolList::operator=(const bool val)
{
    PackedList<1>::operator=(val);
}


inline void Foam::PackedBoolList::operator=(const PackedBoolList& lst)
{
    PackedList<1>::operator=(lst);
}


inline void Foam::PackedBoolList::operator=(PackedBoolList&& lst)
{
    transfer(lst);
}


// ************************************************************************* //
