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

#include "autoPtr.H"
#include "refPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline constexpr Foam::PtrDynList<T, SizeMin>::PtrDynList() noexcept
:
    PtrList<T>(),
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList(const label len)
:
    PtrList<T>(len),
    capacity_(len)
{
    PtrList<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    const PtrDynList<T, SizeMin>& list
)
:
    PtrList<T>(list),
    capacity_(PtrList<T>::size())
{}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList
(
    PtrDynList<T, SizeMin>&& list
)
:
    PtrList<T>(std::move(list)),
    capacity_(list.capacity_)
{
    list.clearStorage();
}


template<class T, int SizeMin>
inline Foam::PtrDynList<T, SizeMin>::PtrDynList(UList<T*>& list)
:
    PtrList<T>(list),
    capacity_(PtrList<T>::size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::label Foam::PtrDynList<T, SizeMin>::capacity() const noexcept
{
    return capacity_;
}


template<class T, int SizeMin>
inline const T* Foam::PtrDynList<T, SizeMin>::get(const label i) const
{
    return (i >= 0 && i < PtrList<T>::size()) ? PtrList<T>::get(i) : nullptr;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::reserve(const label len)
{
    if (capacity_ < len)
    {
        // Preserve addressed size
        const label currLen = PtrList<T>::size();

        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(len, label(2*capacity_)));

        PtrList<T>::resize(capacity_);
        PtrList<T>::setAddressableSize(currLen);
    }
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::resize(const label newLen)
{
    auto& ptrs = this->ptrs_;

    const label oldLen = ptrs.size();

    if (capacity_ < newLen)
    {
        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(newLen, label(2*capacity_)));

        PtrList<T>::resize(capacity_);
    }
    else if (newLen != oldLen)
    {
        // Truncation frees old pointers
        for (label i = newLen; i < oldLen; ++i)
        {
            delete ptrs[i];
            ptrs[i] = nullptr;
        }
    }

    // Adjust addressed size
    PtrList<T>::setAddressableSize(newLen);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::clear()
{
    (this->ptrs_).free(); // free old pointers
    PtrList<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::clearStorage()
{
    PtrList<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline Foam::label Foam::PtrDynList<T, SizeMin>::expandStorage() noexcept
{
    const label currLen = PtrList<T>::size();

    // Allow addressing into the entire list
    PtrList<T>::setAddressableSize(capacity_);

    return currLen;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::shrink()
{
    const label currLen = PtrList<T>::size();
    if (currLen < capacity_)
    {
        // Adjust addressable size to trigger proper resizing
        PtrList<T>::setAddressableSize(currLen+1);

        PtrList<T>::resize(currLen);
        capacity_ = PtrList<T>::size();
    }
}


template<class T, int SizeMin>
inline Foam::label Foam::PtrDynList<T, SizeMin>::squeezeNull()
{
    const label newLen = UPtrList<T>::squeezeNull();
    resize(newLen);
    return newLen;
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::swap
(
    PtrDynList<T, AnySizeMin>& other
)
{
    // Cannot compare 'this' for different types, so use cdata()
    if (this->cdata() == other.cdata())
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UPtrList<T>::swap(other);

    // Swap capacity
    std::swap(this->capacity_, other.capacity_);
}


template<class T, int SizeMin>
template<class... Args>
inline void Foam::PtrDynList<T, SizeMin>::emplace_append(Args&&... args)
{
    this->append(new T(std::forward<Args>(args)...));
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(T* ptr)
{
    const label idx = this->size();
    resize(idx + 1);
    this->ptrs_[idx] = ptr;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(autoPtr<T>& ptr)
{
    this->append(ptr.release());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(autoPtr<T>&& ptr)
{
    this->append(ptr.release());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(std::unique_ptr<T>&& ptr)
{
    this->append(ptr.release());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(const refPtr<T>& ptr)
{
    this->append(ptr.ptr());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(const tmp<T>& ptr)
{
    this->append(ptr.ptr());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(PtrList<T>&& other)
{
    const label idx = this->size();
    const label len = other.size();

    resize(idx + len);

    for (label i=0; i < len; ++i)
    {
        set(idx + i, other.release(i));  // moves pointer
    }

    other.clear();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append
(
    PtrDynList<T, AnySizeMin>&& other
)
{
    const label idx = this->size();
    const label len = other.size();

    resize(idx + len);

    for (label i=0; i < len; ++i)
    {
        set(idx + i, other.release(i));  // moves pointer
    }

    other.clearStorage();  // Ensure capacity=0
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::remove()
{
    // Location of last element and simultaneously the new size
    const label idx = (this->size() - 1);

    if (idx < 0)
    {
        return nullptr;  // List is empty
    }

    autoPtr<T> old(this->ptrs_[idx]);
    this->ptrs_[idx] = nullptr;
    PtrList<T>::setAddressableSize(idx);

    return old;
}


template<class T, int SizeMin>
template<class... Args>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::emplace
(
    const label i,
    Args&&... args
)
{
    return this->set(i, new T(std::forward<Args>(args)...));
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    T* ptr
)
{
    if (i >= this->size())
    {
        resize(i+1);
    }

    return autoPtr<T>(UPtrList<T>::set(i, ptr));
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    autoPtr<T>& ptr
)
{
    return this->set(i, ptr.release());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    autoPtr<T>&& ptr
)
{
    return this->set(i, ptr.release());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    std::unique_ptr<T>&& ptr
)
{
    return this->set(i, ptr.release());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    const refPtr<T>& ptr
)
{
    return this->set(i, ptr.ptr());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    const tmp<T>& ptr
)
{
    return this->set(i, ptr.ptr());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::reorder(const labelUList& oldToNew)
{
    // Shrinking first is a bit annoying, but saves needing a special version.
    shrink();
    PtrList<T>::reorder(oldToNew);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrList<T>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::operator=(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrDynList<T, SizeMin>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::operator=(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrDynList<T, AnySizeMin>& list
)
{
    // Cannot compare 'this' for different types, so use cdata()
    if (this->cdata() == list.cdata())
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::operator=(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrList<T>&& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::transfer(list);
    capacity_ = PtrList<T>::size();
    list.clearStorage();
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrDynList<T, SizeMin>&& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::transfer(list);
    capacity_ = list.capacity();
    list.clearStorage();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrDynList<T, AnySizeMin>&& list
)
{
    // Cannot compare 'this' for different types, so use cdata()
    if (this->cdata() == list.cdata())
    {
        return;  // Self-assignment is a no-op
    }

    PtrList<T>::transfer(list);
    capacity_ = list.capacity();
    list.clearStorage();
}


// ************************************************************************* //
