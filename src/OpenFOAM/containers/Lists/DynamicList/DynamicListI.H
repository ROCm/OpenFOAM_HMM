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

#include "FixedList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
template<class ListType>
inline void Foam::DynamicList<T, SizeMin>::assignDynList
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
inline constexpr Foam::DynamicList<T, SizeMin>::DynamicList() noexcept
:
    List<T>(),
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList(const label len)
:
    List<T>(),
    capacity_(0)
{
    reserve(len);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const label len,
    const T& val
)
:
    List<T>(len, val),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const label len,
    const Foam::zero
)
:
    List<T>(len, Zero),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const DynamicList<T, SizeMin>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const DynamicList<T, AnySizeMin>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const UList<T>& list
)
:
    List<T>(list),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<unsigned N>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const FixedList<T, N>& list
)
:
    capacity_(0)
{
    this->operator=(list);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    std::initializer_list<T> lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<class Addr>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const IndirectListBase<T, Addr>& lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicList<T, SizeMin>&& lst
)
:
    capacity_(0)
{
    transfer(lst);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    DynamicList<T, AnySizeMin>&& lst
)
:
    capacity_(0)
{
    transfer(lst);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    List<T>&& lst
)
:
    capacity_(0)
{
    transfer(lst);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList(Istream& is)
:
    List<T>(),
    capacity_(0)
{
    this->readList(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::capacity() const noexcept
{
    return capacity_;
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setCapacity
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
inline void Foam::DynamicList<T, SizeMin>::reserve
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
inline void Foam::DynamicList<T, SizeMin>::resize
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
inline void Foam::DynamicList<T, SizeMin>::resize
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
inline void Foam::DynamicList<T, SizeMin>::clear() noexcept
{
    List<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::clearStorage()
{
    List<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::expandStorage() noexcept
{
    const label currLen = List<T>::size();

    // Allow addressing into the entire list
    List<T>::setAddressableSize(capacity_);

    return currLen;
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::shrink()
{
    const label currLen = List<T>::size();
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
inline void Foam::DynamicList<T, SizeMin>::swap
(
    DynamicList<T, AnySizeMin>& other
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
inline void
Foam::DynamicList<T, SizeMin>::transfer(List<T>& list)
{
    // Take over storage, clear addressing for list
    capacity_ = list.size();
    List<T>::transfer(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer
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
    // clear addressing and storage for old lst.
    capacity_ = list.capacity();

    List<T>::transfer(static_cast<List<T>&>(list));
    list.clearStorage();  // Ensure capacity=0
}


template<class T, int SizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer
(
    SortableList<T>& list
)
{
    list.shrink();  // Shrink away sort indices
    capacity_ = list.size(); // Capacity after transfer == list size
    List<T>::transfer(list);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
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
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
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
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    const UList<T>& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "Attempted appending to self"
            << abort(FatalError);
    }

    label idx = List<T>::size();
    resize(idx + lst.size());

    for (const T& val : lst)
    {
        this->operator[](idx++) = val;  // copy element
    }
    return *this;
}


template<class T, int SizeMin>
template<unsigned N>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    const FixedList<T, N>& lst
)
{
    label idx = List<T>::size();
    resize(idx + lst.size());

    for (const T& val : lst)
    {
        this->operator[](idx++) = val;  // copy element
    }
    return *this;
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    std::initializer_list<T> lst
)
{
    label idx = List<T>::size();
    resize(idx + lst.size());

    for (const T& val : lst)
    {
        this->operator[](idx++) = val;  // copy element
    }
    return *this;
}


template<class T, int SizeMin>
template<class Addr>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    const IndirectListBase<T, Addr>& lst
)
{
    label idx = List<T>::size();
    const label n = lst.size();

    resize(idx + n);

    for (label i=0; i<n; ++i)
    {
        this->operator[](idx++) = lst[i];  // copy element
    }
    return *this;
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    List<T>&& list
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

    for (T& val : list)
    {
        Foam::Swap(this->operator[](idx++), val);  // moved content
    }

    list.clear();
    return *this;
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    DynamicList<T, SizeMin>&& lst
)
{
    append(std::move(static_cast<List<T>&>(lst)));
    lst.clearStorage();  // Ensure capacity=0
    return *this;
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    DynamicList<T, AnySizeMin>&& lst
)
{
    append(std::move(static_cast<List<T>&>(lst)));
    lst.clearStorage();  // Ensure capacity=0
    return *this;
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    SortableList<T>&& lst
)
{
    lst.shrink();  // Shrink away sort indices
    append(std::move(static_cast<List<T>&>(lst)));
    return *this;
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::appendUniq(const T& val)
{
    if (this->found(val))
    {
        return 0;
    }
    else
    {
        this->append(val);
        return 1;  // Increased list length by one
    }
}


template<class T, int SizeMin>
inline T Foam::DynamicList<T, SizeMin>::remove()
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


template<class T, int SizeMin>
inline T Foam::DynamicList<T, SizeMin>::remove
(
    const label idx,
    const bool fast
)
{
    if (fast)
    {
        // Simply swap idx <=> last
        this->swapLast(idx);
    }
    else
    {
        // Move element to the end and move everything down
        this->moveLast(idx);
    }

    // Element to remove is now at the end
    return this->remove();
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::remove
(
    const labelRange& range
)
{
    return this->removeElements(this->validateRange(range));
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::remove
(
    std::initializer_list<label> start_size
)
{
    return this->removeElements(this->validateRange(start_size));
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::subset
(
    const labelRange& range
)
{
    return this->subsetElements(this->validateRange(range));
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::subset
(
    std::initializer_list<label> start_size
)
{
    return this->subsetElements(this->validateRange(start_size));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline T& Foam::DynamicList<T, SizeMin>::operator()
(
    const label i
)
{
    if (i >= List<T>::size())
    {
        resize(i + 1);
    }

    return this->operator[](i);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const T& val
)
{
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const Foam::zero
)
{
    UList<T>::operator=(Zero);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const UList<T>& lst
)
{
    assignDynList(lst);
}


template<class T, int SizeMin>
template<unsigned N>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const FixedList<T, N>& lst
)
{
    const label n = lst.size();

    resize(n);

    for (label i=0; i<n; ++i)
    {
        this->operator[](i) = lst[i];  // copy element
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const DynamicList<T, SizeMin>& lst
)
{
    if (this == &lst)
    {
        return;  // Self-assignment is a no-op
    }

    assignDynList(lst);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const DynamicList<T, AnySizeMin>& list
)
{
    // Cannot compare 'this' for different types, so use cdata()
    if (this->cdata() == list.cdata())
    {
        return;  // Self-assignment is a no-op
    }

    assignDynList(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    std::initializer_list<T> lst
)
{
    assignDynList(lst);
}


template<class T, int SizeMin>
template<class Addr>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const IndirectListBase<T, Addr>& lst
)
{
    assignDynList(lst);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    List<T>&& lst
)
{
    clear();
    transfer(lst);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    DynamicList<T, SizeMin>&& lst
)
{
    if (this == &lst)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    transfer(lst);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    DynamicList<T, AnySizeMin>&& list
)
{
    // Cannot compare 'this' for different types, so use cdata()
    if (this->cdata() == list.cdata())
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    transfer(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    SortableList<T>&& lst
)
{
    clear();
    transfer(lst);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::Istream& Foam::DynamicList<T, SizeMin>::readList
(
    Istream& is
)
{
    DynamicList<T, SizeMin>& list = *this;

    // Use entire storage - ie, resize(capacity())
    (void) list.expandStorage();

    static_cast<List<T>&>(list).readList(is);
    list.capacity_ = list.size();

    return is;
}


template<class T, int SizeMin>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    DynamicList<T, SizeMin>& list
)
{
    return list.readList(is);
}


template<class T, int SizeMin>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DynamicList<T, SizeMin>& list
)
{
    os << static_cast<const List<T>&>(list);
    return os;
}


// ************************************************************************* //
