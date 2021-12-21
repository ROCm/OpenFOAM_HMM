/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "PtrListDetail.H"
#include <utility>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::label Foam::Detail::PtrListDetail<T>::count() const
{
    label ngood = 0;

    for (const T* ptr : *this)
    {
        if (ptr)
        {
            ++ngood;
        }
    }

    return ngood;
}


template<class T>
Foam::label Foam::Detail::PtrListDetail<T>::findNull() const
{
    label idx = 0;

    for (const T* ptr : *this)
    {
        if (!ptr)
        {
            return idx;
        }

        ++idx;
    }

    return -1;
}


template<class T>
void Foam::Detail::PtrListDetail<T>::setNull()
{
    List<T*>& ptrs = *this;
    const label len = ptrs.size();

    for (label i=0; i<len; ++i)
    {
        ptrs[i] = nullptr;
    }
}


template<class T>
void Foam::Detail::PtrListDetail<T>::free()
{
    List<T*>& ptrs = *this;
    const label len = ptrs.size();

    for (label i=0; i<len; ++i)
    {
        T* ptr = ptrs[i];

        if (ptr)
        {
            delete ptr;
        }

        ptrs[i] = nullptr;
    }
}


template<class T>
template<class... Args>
Foam::Detail::PtrListDetail<T>
Foam::Detail::PtrListDetail<T>::clone(Args&&... args) const
{
    const List<T*>& ptrs = *this;
    const label len = ptrs.size();

    PtrListDetail<T> cloned(len);

    for (label i=0; i<len; ++i)
    {
        const T* ptr = ptrs[i];

        if (ptr)
        {
            cloned[i] = ptr->clone(std::forward<Args>(args)...).ptr();
        }
    }

    return cloned;
}


// ************************************************************************* //
