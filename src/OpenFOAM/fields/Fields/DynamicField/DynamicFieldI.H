/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
template<class ListType>
inline void Foam::DynamicField<T, SizeMin>::assignDynList
(
    const ListType& list
)
{
    const label newLen = list.size();

    if (newLen <= capacity_)
    {
        // Can copy w/o reallocating - adjust addressable size accordingly.
        List<T>::setAddressableSize(newLen);
        List<T>::operator=(list);
    }
    else
    {
        // Ensure list size consistency prior to copying.
        List<T>::setAddressableSize(capacity_);

        List<T>::operator=(list);
        capacity_ = List<T>::size();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline constexpr Foam::DynamicField<T, SizeMin>::DynamicField() noexcept
:
    Field<T>(),
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField(const label len)
:
    Field<T>(),
    capacity_(0)
{
    reserve(len);
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const label len,
    const T& val
)
:
    Field<T>(len, val),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const label len,
    const Foam::zero
)
:
    Field<T>(len, Zero),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const DynamicField<T, SizeMin>& list
)
:
    Field<T>(list),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const DynamicField<T, AnySizeMin>& list
)
:
    Field<T>(list),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& list
)
:
    Field<T>(list),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
template<class Addr>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const IndirectListBase<T, Addr>& list
)
:
    Field<T>(list),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    List<T>&& content
)
:
    Field<T>(std::move(content)),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    DynamicField<T, SizeMin>&& content
)
:
    Field<T>(),
    capacity_(0)
{
    transfer(content);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    DynamicField<T, AnySizeMin>&& content
)
:
    Field<T>(),
    capacity_(0)
{
    transfer(content);
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& mapF,
    const labelUList& mapAddressing
)
:
    Field<T>(mapF, mapAddressing),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& mapF,
    const labelListList& mapAddressing,
    const scalarListList& weights
)
:
    Field<T>(mapF, mapAddressing, weights),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& mapF,
    const FieldMapper& map
)
:
    Field<T>(mapF, map),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField(Istream& is)
:
    Field<T>(is),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::tmp<Foam::DynamicField<T, SizeMin>>
Foam::DynamicField<T, SizeMin>::clone() const
{
    return tmp<DynamicField<T, SizeMin>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::label Foam::DynamicField<T, SizeMin>::capacity() const noexcept
{
    return capacity_;
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::setCapacity
(
    const label newCapacity
)
{
    label currLen = List<T>::size();
    capacity_ = newCapacity;

    if (currLen > capacity_)
    {
        // Truncate addressed sizes too
        currLen = capacity_;
    }

    List<T>::resize(capacity_);
    List<T>::setAddressableSize(currLen);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::reserve
(
    const label len
)
{
    if (capacity_ < len)
    {
        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(len, label(2 * capacity_)));

        // Adjust allocated size, leave addressed size untouched
        const label currLen = List<T>::size();
        List<T>::resize(capacity_);
        List<T>::setAddressableSize(currLen);
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::resize
(
    const label newLen
)
{
    if (capacity_ < newLen)
    {
        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(newLen, label(2 * capacity_)));

        List<T>::resize(capacity_);
    }

    // Adjust addressed size
    List<T>::setAddressableSize(newLen);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::resize
(
    const label newLen,
    const T& val
)
{
    label currLen = List<T>::size();
    resize(newLen);

    // Fill newly exposed with constant value
    while (currLen < newLen)
    {
        this->operator[](currLen++) = val;
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::clear() noexcept
{
    List<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::clearStorage()
{
    List<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicField<T, SizeMin>::expandStorage() noexcept
{
    const label currLen = List<T>::size();

    // Allow addressing into the entire list
    List<T>::setAddressableSize(capacity_);

    return currLen;
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>&
Foam::DynamicField<T, SizeMin>::shrink()
{
    const label currLen = Field<T>::size();
    if (currLen < capacity_)
    {
        // Use the full list when resizing
        List<T>::setAddressableSize(capacity_);

        // Capacity and size are identical
        capacity_ = currLen;
        List<T>::resize(currLen);
        // Redundant: List<T>::setAddressableSize(currLen);
    }
    return *this;
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::swap
(
    DynamicField<T, AnySizeMin>& other
)
{
    // Cannot compare 'this' for different types, so use cdata()
    if (this->cdata() == other.cdata())
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UList<T>::swap(other);

    // Swap capacity
    std::swap(this->capacity_, other.capacity_);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::transfer(List<T>& list)
{
    // Take over storage, clear addressing for list
    capacity_ = list.size();
    Field<T>::transfer(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::transfer
(
    DynamicList<T, AnySizeMin>& list
)
{
    // Cannot compare 'this' for different types, so use cdata()
    if (this->cdata() == list.cdata())
    {
        return;  // Self-assignment is a no-op
    }

    // Take over storage as-is (without shrink, without using SizeMin)
    // clear addressing and storage for old list.
    capacity_ = list.capacity();

    Field<T>::transfer(static_cast<Field<T>&>(list));
    list.clearStorage();  // Ensure capacity=0
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::transfer
(
    DynamicField<T, AnySizeMin>& list
)
{
    // Cannot compare 'this' for different types, so use cdata()
    if (this->cdata() == list.cdata())
    {
        return;  // Self-assignment is a no-op
    }

    // Take over storage as-is (without shrink, without using SizeMin)
    // clear addressing and storage for old list.
    capacity_ = list.capacity();

    Field<T>::transfer(static_cast<Field<T>&>(list));
    list.clearStorage();  // Ensure capacity=0
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>&
Foam::DynamicField<T, SizeMin>::append
(
    const T& val
)
{
    const label idx = List<T>::size();
    resize(idx + 1);

    this->operator[](idx) = val;  // copy element
    return *this;
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>&
Foam::DynamicField<T, SizeMin>::append
(
    T&& val
)
{
    const label idx = List<T>::size();
    resize(idx + 1);

    this->operator[](idx) = std::move(val);  // move assign element
    return *this;
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>&
Foam::DynamicField<T, SizeMin>::append
(
    const UList<T>& list
)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "Attempted appending to self"
            << abort(FatalError);
    }

    label idx = List<T>::size();
    resize(idx + list.size());

    for (const T& val : list)
    {
        this->operator[](idx++) = val;  // copy element
    }
    return *this;
}


template<class T, int SizeMin>
inline T Foam::DynamicField<T, SizeMin>::remove()
{
    // Location of last element and simultaneously the new size
    const label idx = List<T>::size() - 1;

    if (idx < 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    const T& val = List<T>::operator[](idx);

    List<T>::setAddressableSize(idx);

    return val;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline T& Foam::DynamicField<T, SizeMin>::operator()
(
    const label i
)
{
    if (i >= Field<T>::size())
    {
        resize(i + 1);
    }

    return this->operator[](i);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const T& val
)
{
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const Foam::zero
)
{
    UList<T>::operator=(Zero);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const UList<T>& list
)
{
    assignDynList(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const DynamicField<T, SizeMin>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    assignDynList(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    List<T>&& list
)
{
    transfer(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    DynamicField<T, SizeMin>&& list
)
{
    transfer(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    DynamicField<T, AnySizeMin>&& list
)
{
    transfer(list);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    DynamicField<T, SizeMin>& rhs
)
{
    // Use entire storage - ie, resize(capacity())
    (void) rhs.expandStorage();

    is >> static_cast<Field<T>&>(rhs);
    rhs.capacity_ = rhs.Field<T>::size();

    return is;
}


template<class T, int SizeMin>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DynamicField<T, SizeMin>& rhs
)
{
    os << static_cast<const Field<T>&>(rhs);
    return os;
}


// ************************************************************************* //
