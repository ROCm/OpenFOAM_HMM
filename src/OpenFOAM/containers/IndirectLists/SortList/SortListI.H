/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::SortList<T>::SortList(const UList<T>& values)
:
    IndirectList<T>(values, identity(values.size()))
{
    sort();
}


template<class T>
template<class Compare>
inline Foam::SortList<T>::SortList(const UList<T>& values, const Compare& comp)
:
    IndirectList<T>(values, identity(values.size()))
{
    sort<Compare>(comp);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::labelUList& Foam::SortList<T>::indices() const
{
    return this->addressing();
}


template<class T>
inline Foam::labelList& Foam::SortList<T>::indices()
{
    return this->addressing();
}


template<class T>
inline void Foam::SortList<T>::reverse()
{
    Foam::reverse(this->indices());
}


template<class T>
inline void Foam::SortList<T>::reset()
{
    const label len = this->values().size();

    labelList& addr = this->indices();
    addr.resize(len);
    ListOps::identity(addr);
}


template<class T>
template<class Compare>
inline void Foam::SortList<T>::sort(const Compare& comp)
{
    UList<T>& vals = this->values();
    labelList& addr = this->indices();

    if (addr.size() != vals.size())
    {
        addr.resize(vals.size());
        ListOps::identity(addr);
    }

    std::stable_sort
    (
        addr.begin(),
        addr.end(),
        [&](label a, label b) -> bool { return comp(vals[a], vals[b]); }
    );
}


template<class T>
inline void Foam::SortList<T>::sort()
{
    UList<T>& vals = this->values();
    labelList& addr = this->indices();

    if (addr.size() != vals.size())
    {
        addr.resize(vals.size());
        ListOps::identity(addr);
    }

    // Forward sort of indices
    std::stable_sort
    (
        addr.begin(),
        addr.end(),
        [&](label a, label b) -> bool { return vals[a] < vals[b]; }
    );
}


template<class T>
inline void Foam::SortList<T>::reverseSort()
{
    UList<T>& vals = this->values();
    labelList& addr = this->indices();

    if (addr.size() != vals.size())
    {
        addr.resize(vals.size());
        ListOps::identity(addr);
    }

    // Reverse sort of indices
    std::stable_sort
    (
        addr.begin(),
        addr.end(),
        [&](label a, label b) -> bool { return vals[b] < vals[a]; }
    );
}


// ************************************************************************* //
