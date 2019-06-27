/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::setAddressableSize(const label len)
{
    (this->ptrs_).setAddressableSize(len);
}


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
    setAddressableSize(0);
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
inline Foam::label Foam::PtrDynList<T, SizeMin>::capacity() const
{
    return capacity_;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::setCapacity(const label nElem)
{
    label nextFree = PtrList<T>::size();
    capacity_ = nElem;

    if (nextFree > capacity_)
    {
        // Truncate addressed sizes too
        nextFree = capacity_;
    }

    PtrList<T>::resize(capacity_);
    setAddressableSize(nextFree);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::reserve(const label nElem)
{
    if (nElem > capacity_)
    {
        // Allocate more capacity if necessary

        capacity_ = max(SizeMin, max(nElem, label(2*capacity_)));

        // Adjust allocated size, leave addressed size untouched
        const label nextFree = PtrList<T>::size();
        PtrList<T>::resize(capacity_);
        setAddressableSize(nextFree);
    }
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::resize(const label newLen)
{
    auto& ptrs = this->ptrs_;

    const label oldLen = ptrs.size();

    if (newLen > capacity_)
    {
        // Allocate more capacity if necessary
        capacity_ = max(SizeMin, max(newLen, label(2*capacity_)));

        PtrList<T>::resize(capacity_);
    }
    else if (newLen != oldLen)
    {
        // Truncation frees old pointers
        for (label i=newLen; i<oldLen; ++i)
        {
            T* ptr = ptrs[i];

            if (ptr)
            {
                delete ptr;
            }

            ptrs[i] = nullptr;
        }
    }

    setAddressableSize(newLen);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::setSize(const label newLen)
{
    this->resize(newLen);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::clear()
{
    (this->ptrs_).free(); // free old pointers
    setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::clearStorage()
{
    PtrList<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline Foam::label Foam::PtrDynList<T, SizeMin>::expandStorage()
{
    const label nextFree = PtrList<T>::size();

    // Allow addressing into the entire list
    setAddressableSize(capacity_);

    return nextFree;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::shrink()
{
    label nextFree = PtrList<T>::size();
    if (capacity_ > nextFree)
    {
        // Use the full list when resizing
        setAddressableSize(capacity_);

        // The new size
        capacity_ = nextFree;
        PtrList<T>::resize(capacity_);
        setAddressableSize(nextFree);
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
inline void Foam::PtrDynList<T, SizeMin>::append(T* ptr)
{
    const label idx = this->size();
    resize(idx + 1);
    this->ptrs_[idx] = ptr;
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(const autoPtr<T>& aptr)
{
    this->append(const_cast<autoPtr<T>&>(aptr).release());
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::append(const tmp<T>& tptr)
{
    this->append(tptr.ptr());
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
    setAddressableSize(idx);

    return old;
}


template<class T, int SizeMin>
inline const T* Foam::PtrDynList<T, SizeMin>::set(const label i) const
{
    return (i >= 0 && i < PtrList<T>::size()) ? PtrList<T>::set(i) : nullptr;
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
    const autoPtr<T>& aptr
)
{
    return this->set(i, const_cast<autoPtr<T>&>(aptr).release());
}


template<class T, int SizeMin>
inline Foam::autoPtr<T> Foam::PtrDynList<T, SizeMin>::set
(
    const label i,
    const tmp<T>& tptr
)
{
    return this->set(i, tptr.ptr());
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
    PtrList<T>::operator=(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    const PtrDynList<T, SizeMin>& list
)
{
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
    PtrList<T>::operator=(list);
    capacity_ = PtrList<T>::size();
}


template<class T, int SizeMin>
inline void Foam::PtrDynList<T, SizeMin>::operator=
(
    PtrList<T>&& list
)
{
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
    PtrList<T>::transfer(list);
    capacity_ = list.capacity();
    list.clearStorage();
}


// ************************************************************************* //
