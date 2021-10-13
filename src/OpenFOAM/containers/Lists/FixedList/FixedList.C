/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "FixedList.H"
#include "ListLoopM.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T, unsigned N>
std::streamsize Foam::FixedList<T, N>::byteSize()
{
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << abort(FatalError);
    }
    return FixedList<T, N>::size_bytes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned N>
Foam::label Foam::FixedList<T, N>::find(const T& val, label pos) const
{
    if (pos >= 0)
    {
        List_CONST_ACCESS(T, *this, list);

        while (pos < label(N))
        {
            if (list[pos] == val)
            {
                return pos;
            }

            ++pos;
        }
    }

    return -1;
}


template<class T, unsigned N>
Foam::label Foam::FixedList<T, N>::rfind(const T& val, label pos) const
{
    // pos == -1 has same meaning as std::string::npos - search from end
    if (pos < 0 || pos >= label(N))
    {
        pos = label(N)-1;
    }

    List_CONST_ACCESS(T, *this, list);

    while (pos >= 0)
    {
        if (list[pos] == val)
        {
            return pos;
        }

        --pos;
    }

    return -1;
}


template<class T, unsigned N>
void Foam::FixedList<T, N>::moveFirst(const label i)
{
    checkIndex(i);

    for (label lower = 0; lower < i; ++lower)
    {
        Foam::Swap(v_[lower], v_[i]);
    }
}


template<class T, unsigned N>
void Foam::FixedList<T, N>::moveLast(const label i)
{
    checkIndex(i);

    for (label upper = label(N-1); upper > i; --upper)
    {
        Foam::Swap(v_[i], v_[upper]);
    }
}


template<class T, unsigned N>
void Foam::FixedList<T, N>::swapFirst(const label i)
{
    checkIndex(i);

    if (i > 0)
    {
        Foam::Swap(v_[0], v_[i]);
    }
}


template<class T, unsigned N>
void Foam::FixedList<T, N>::swapLast(const label i)
{
    checkIndex(i);

    const label upper = label(N-1);

    if (i < upper)
    {
        Foam::Swap(v_[i], v_[upper]);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator==(const FixedList<T, N>& list) const
{
    List_CONST_ACCESS(T, *this, lhs);
    List_CONST_ACCESS(T, (list), rhs);

    // List sizes are identical by definition (template parameter)
    for (unsigned i = 0; i < N; ++i)
    {
        if (!(lhs[i] == rhs[i]))
        {
            return false;
        }
    }

    // Contents appear to be identical.
    return true;
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator<(const FixedList<T, N>& list) const
{
    List_CONST_ACCESS(T, *this, lhs);
    List_CONST_ACCESS(T, (list), rhs);

    // List sizes are identical by definition (template parameter)
    for (unsigned i=0; i<N; ++i)
    {
        if (lhs[i] < rhs[i])
        {
            return true;
        }
        else if (rhs[i] < lhs[i])
        {
            return false;
        }
    }

    // Contents appear to be identical.
    return false;
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator!=(const FixedList<T, N>& list) const
{
    return !operator==(list);
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator>(const FixedList<T, N>& list) const
{
    return list.operator<(*this);
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator<=(const FixedList<T, N>& list) const
{
    return !list.operator<(*this);
}


template<class T, unsigned N>
bool Foam::FixedList<T, N>::operator>=(const FixedList<T, N>& list) const
{
    return !operator<(list);
}


// ************************************************************************* //
