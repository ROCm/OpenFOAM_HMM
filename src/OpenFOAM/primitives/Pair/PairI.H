/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "Swap.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
inline int Foam::Pair<T>::compare(const Pair<T>& a, const Pair<T>& b)
{
    if (a.first() == b.first() && a.second() == b.second())
    {
        return 1;
    }
    if (a.first() == b.second() && a.second() == b.first())
    {
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::Pair<T>::Pair()
{}


template<class T>
inline Foam::Pair<T>::Pair(const T& f, const T& s)
{
    first()  = f;
    second() = s;
}


template<class T>
inline Foam::Pair<T>::Pair(const FixedList<T, 2>& lst)
:
    FixedList<T, 2>(lst)
{}


template<class T>
inline Foam::Pair<T>::Pair(const T& f, const T& s, const bool doSort)
{
    if (doSort && f < s)
    {
        first()  = s;
        second() = f;
    }
    else
    {
        first()  = f;
        second() = s;
    }
}


template<class T>
inline Foam::Pair<T>::Pair(const FixedList<T, 2>& lst, const bool doSort)
:
    Pair<T>(lst.first(), lst.last(), doSort)
{}



template<class T>
inline Foam::Pair<T>::Pair(Istream& is)
:
    FixedList<T, 2>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::Pair<T>::second() const
{
    return last();
}


template<class T>
inline T& Foam::Pair<T>::second()
{
    return last();
}


template<class T>
inline const T& Foam::Pair<T>::other(const T& a) const
{
    if (first() == second())
    {
        FatalErrorInFunction
            << "Call to other only valid for Pair with differing"
            << " elements:" << *this << abort(FatalError);
    }
    else if (a == first())
    {
        return second();
    }
    else
    {
        if (a != second())
        {
            FatalErrorInFunction
                << "Pair " << *this
                << " does not contain " << a << abort(FatalError);
        }
        return first();
    }
}


template<class T>
inline bool Foam::Pair<T>::sorted() const
{
    return (first() < second());
}


template<class T>
inline void Foam::Pair<T>::flip()
{
    Foam::Swap(first(), second());
}


template<class T>
inline void Foam::Pair<T>::sort()
{
    if (second() < first())
    {
        flip();
    }
}


// ************************************************************************* //
