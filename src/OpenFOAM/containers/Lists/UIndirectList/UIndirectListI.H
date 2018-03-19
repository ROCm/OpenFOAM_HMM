/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
inline bool Foam::UIndirectList<T>::uniform() const
{
    const label len = this->size();

    if (len > 1)
    {
        const T& val = (*this)[0];

        for (label i=1; i<len; ++i)
        {
            if (val != (*this)[i])
            {
                return false;
            }
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::UIndirectList<T>::UIndirectList
(
    const UList<T>& values,
    const labelUList& addr
)
:
    values_(const_cast<UList<T>&>(values)),
    addressing_(addr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::UIndirectList<T>::size() const
{
    return addressing_.size();
}


template<class T>
inline bool Foam::UIndirectList<T>::empty() const
{
    return addressing_.empty();
}


template<class T>
inline T& Foam::UIndirectList<T>::first()
{
    return values_[addressing_.first()];
}


template<class T>
inline const T& Foam::UIndirectList<T>::first() const
{
    return values_[addressing_.first()];
}


template<class T>
inline T& Foam::UIndirectList<T>::last()
{
    return values_[addressing_.last()];
}


template<class T>
inline const T& Foam::UIndirectList<T>::last() const
{
    return values_[addressing_.last()];
}


template<class T>
inline const Foam::UList<T>& Foam::UIndirectList<T>::completeList() const
{
    return values_;
}


template<class T>
inline const Foam::List<Foam::label>& Foam::UIndirectList<T>::addressing() const
{
    return addressing_;
}


template<class T>
inline bool Foam::UIndirectList<T>::found
(
    const T& val,
    const label pos
) const
{
    return this->find(val, pos) != -1;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::List<T> Foam::UIndirectList<T>::operator()() const
{
    List<T> result(size());

    forAll(*this, i)
    {
        result[i] = operator[](i);
    }

    return result;
}


template<class T>
inline T& Foam::UIndirectList<T>::operator[](const label i)
{
    return values_[addressing_[i]];
}


template<class T>
inline const T& Foam::UIndirectList<T>::operator[](const label i) const
{
    return values_[addressing_[i]];
}


template<class T>
inline void Foam::UIndirectList<T>::operator=(const UList<T>& rhs)
{
    const label len = addressing_.size();

    if (len != rhs.size())
    {
        FatalErrorInFunction
            << "Addressing and list of addressed elements "
               "have different sizes: " << len << " " << rhs.size()
            << abort(FatalError);
    }

    for (label i = 0; i < len; ++i)
    {
        values_[addressing_[i]] = rhs[i];
    }
}


template<class T>
inline void Foam::UIndirectList<T>::operator=(const UIndirectList<T>& rhs)
{
    const label len = addressing_.size();

    if (len != rhs.size())
    {
        FatalErrorInFunction
            << "Addressing and list of addressed elements "
               "have different sizes: " << len << " " << rhs.size()
            << abort(FatalError);
    }

    for (label i = 0; i < len; ++i)
    {
        values_[addressing_[i]] = rhs[i];
    }
}


template<class T>
inline void Foam::UIndirectList<T>::operator=(const T& val)
{
    for (const label idx : addressing_)
    {
        values_[idx] = val;
    }
}


// ************************************************************************* //
