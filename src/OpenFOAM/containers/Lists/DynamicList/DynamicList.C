/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "DynamicList.H"
#include "labelRange.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::label Foam::DynamicList<T, SizeMin>::removeElements
(
    const labelRange& slice
)
{
    if (!slice.size())
    {
        // Noop
        return 0;
    }
    else if (slice.after() >= this->size())
    {
        // Remove tail
        this->resize(slice.first());
    }
    else
    {
        // Copy (swap) down
        label j = slice.first();
        const label len = this->size();

        for (label i = slice.after(); i < len; ++i, ++j)
        {
            Foam::Swap(this->operator[](i), this->operator[](j));
        }

        resize(this->size() - slice.size());
    }

    return slice.size();
}


template<class T, int SizeMin>
Foam::label Foam::DynamicList<T, SizeMin>::subsetElements
(
    const labelRange& slice
)
{
    if (slice.first() > 0)
    {
        // Copy (swap) down
        label j = slice.first();
        const label len = slice.size();

        for (label i = 0; i < len; ++i, ++j)
        {
            Foam::Swap(this->operator[](i), this->operator[](j));
        }
    }

    // Don't need min size, since slice size was already checked before
    resize(slice.size());
    return this->size();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, int SizeMin>
Foam::DynamicList<T, SizeMin>::DynamicList(Istream& is)
:
    List<T>(is),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DynamicList<T, SizeMin>& lst
)
{
    os << static_cast<const List<T>&>(lst);
    return os;
}


template<class T, int SizeMin>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    DynamicList<T, SizeMin>& lst
)
{
    is >> static_cast<List<T>&>(lst);
    lst.capacity_ = lst.List<T>::size();

    return is;
}


// ************************************************************************* //
