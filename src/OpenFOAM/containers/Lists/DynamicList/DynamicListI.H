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

#include "FixedList.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
template<class ListType>
inline void Foam::DynamicList<T, SizeMin>::assignDynList
(
    const ListType& lst
)
{
    const label newSize = lst.size();

    if (capacity_ >= newSize)
    {
        // Can copy w/o reallocating - adjust addressable size accordingly.
        List<T>::size(newSize);
        List<T>::operator=(lst);
    }
    else
    {
        // Ensure list size consistency prior to copying.
        List<T>::size(capacity_);

        List<T>::operator=(lst);
        capacity_ = List<T>::size();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList()
:
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList(const label nElem)
:
    capacity_(0)
{
    reserve(nElem);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const label nElem,
    const T& val
)
:
    List<T>(nElem, val),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const label nElem,
    const zero
)
:
    List<T>(nElem, Zero),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const DynamicList<T, SizeMin>& lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const DynamicList<T, AnySizeMin>& lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const UList<T>& lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
template<unsigned FixedSize>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const FixedList<T, FixedSize>& lst
)
:
    capacity_(0)
{
    this->operator=(lst);
}


template<class T, int SizeMin>
template<class InputIterator>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    InputIterator begIter,
    InputIterator endIter
)
:
    List<T>(begIter, endIter),
    capacity_(List<T>::size())
{}


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
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const UIndirectList<T>& lst
)
:
    List<T>(lst),
    capacity_(List<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>::DynamicList
(
    const Xfer<List<T>>& lst
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::min_size() const
{
    return SizeMin;
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::capacity() const
{
    return capacity_;
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setCapacity
(
    const label nElem
)
{
    label nextFree = List<T>::size();
    capacity_ = nElem;

    if (nextFree > capacity_)
    {
        // Truncate addressed sizes too
        nextFree = capacity_;
    }

    // We could also enforce sizing granularity

    List<T>::setSize(capacity_);
    List<T>::size(nextFree);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::reserve
(
    const label nElem
)
{
    // Allocate more capacity if necessary
    if (nElem > capacity_)
    {
        capacity_ = max
        (
            SizeMin,
            max
            (
                nElem,
                // label(SizeInc + capacity_ * SizeMult / SizeDiv)
                label(2 * capacity_)
            )
        );

        // Adjust allocated size, leave addressed size untouched
        const label nextFree = List<T>::size();
        List<T>::setSize(capacity_);
        List<T>::size(nextFree);
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setSize
(
    const label nElem
)
{
    // Allocate more capacity if necessary
    if (nElem > capacity_)
    {
        capacity_ = max
        (
            SizeMin,
            max
            (
                nElem,
                // label(SizeInc + capacity_ * SizeMult / SizeDiv)
                label(2 * capacity_)
            )
        );

        List<T>::setSize(capacity_);
    }

    // Adjust addressed size
    List<T>::size(nElem);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::setSize
(
    const label nElem,
    const T& val
)
{
    label nextFree = List<T>::size();
    setSize(nElem);

    // Set new elements to constant value
    while (nextFree < nElem)
    {
        this->operator[](nextFree++) = val;
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize
(
    const label nElem
)
{
    this->setSize(nElem);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::resize
(
    const label nElem,
    const T& val
)
{
    this->setSize(nElem, val);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::clear()
{
    List<T>::size(0);
}


template<class T, int SizeMin>
inline void Foam::DynamicList<T, SizeMin>::clearStorage()
{
    List<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicList<T, SizeMin>::expandStorage()
{
    const label nextFree = List<T>::size();

    // Allow addressing into the entire list
    List<T>::size(capacity_);

    return nextFree;
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::shrink()
{
    const label nextFree = List<T>::size();
    if (capacity_ > nextFree)
    {
        // Use the full list when resizing
        List<T>::size(capacity_);

        // The new size
        capacity_ = nextFree;
        List<T>::setSize(capacity_);
        List<T>::size(nextFree);
    }
    return *this;
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicList<T, SizeMin>::swap
(
    DynamicList<T, AnySizeMin>& lst
)
{
    DynamicList<T, SizeMin>& cur = *this;

    // Make addressable size identical to the allocated capacity
    const label oldSize1 = cur.expandStorage();
    const label oldSize2 = lst.expandStorage();

    // Swap storage
    Foam::Swap
    (
        static_cast<UList<T>&>(cur),
        static_cast<UList<T>&>(lst)
    );

    // Match capacity to the underlying allocated list size
    cur.setCapacity(cur.size());
    lst.setCapacity(lst.size());

    // Set addressable size
    cur.setSize(oldSize2);
    lst.setSize(oldSize1);
}


template<class T, int SizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer(List<T>& lst)
{
    // Take over storage, clear addressing for lst.
    capacity_ = lst.size();
    List<T>::transfer(lst);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer
(
    DynamicList<T, AnySizeMin>& lst
)
{
    // Take over storage as-is (without shrink, without using SizeMin)
    // clear addressing and storage for old lst.
    capacity_ = lst.capacity();

    List<T>::transfer(static_cast<List<T>&>(lst));
    lst.clearStorage();  // Ensure capacity=0
}


template<class T, int SizeMin>
inline void
Foam::DynamicList<T, SizeMin>::transfer
(
    SortableList<T>& lst
)
{
    lst.shrink();  // Shrink away sort indices
    capacity_ = lst.size(); // Capacity after transfer == list size
    List<T>::transfer(lst);
}


template<class T, int SizeMin>
inline Foam::Xfer<Foam::List<T>>
Foam::DynamicList<T, SizeMin>::xfer()
{
    return xferMoveTo<List<T>>(*this);
}


template<class T, int SizeMin>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    const T& val
)
{
    const label idx = List<T>::size();
    setSize(idx + 1);

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
    setSize(idx + 1);

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
            << "Attempted appending to self" << abort(FatalError);
    }

    label idx = List<T>::size();

    setSize(idx + lst.size());

    for (const T& val : lst)
    {
        this->operator[](idx++) = val;  // copy element
    }
    return *this;
}


template<class T, int SizeMin>
template<unsigned FixedSize>
inline Foam::DynamicList<T, SizeMin>&
Foam::DynamicList<T, SizeMin>::append
(
    const FixedList<T, FixedSize>& lst
)
{
    label idx = List<T>::size();
    setSize(idx + lst.size());

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

    setSize(idx + lst.size());

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
    const UIndirectList<T>& lst
)
{
    label idx = List<T>::size();
    const label n = lst.size();

    setSize(idx + n);

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
    List<T>&& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "Attempted appending to self" << abort(FatalError);
    }

    label idx = List<T>::size();

    setSize(idx + lst.size());

    for (T& val : lst)
    {
        Foam::Swap(this->operator[](idx++), val);  // moved content
    }

    lst.clear();
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

    List<T>::size(idx);

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
        setSize(i + 1);
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
    const zero
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
template<unsigned FixedSize>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const FixedList<T, FixedSize>& lst
)
{
    const label n = lst.size();

    setSize(n);

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
        FatalErrorInFunction
            << "Attempted assignment to self" << abort(FatalError);
    }

    assignDynList(lst);
}


template<class T, int SizeMin>
template<int SizeMin2>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const DynamicList<T, SizeMin2>& lst
)
{
    assignDynList(lst);
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
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    const UIndirectList<T>& lst
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
        FatalErrorInFunction
            << "Attempted assignment to self" << abort(FatalError);
    }

    clear();
    transfer(lst);
}


template<class T, int SizeMin>
template<int SizeMin2>
inline void Foam::DynamicList<T, SizeMin>::operator=
(
    DynamicList<T, SizeMin2>&& lst
)
{
    clear();
    transfer(lst);
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


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin1, int SizeMin2>
inline void Foam::Swap(DynamicList<T, SizeMin1>& a, DynamicList<T, SizeMin2>& b)
{
    a.swap(b);
}


// ************************************************************************* //
