/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::UPtrList<T>::UPtrList(const label nElem)
:
    ptrs_(nElem, reinterpret_cast<T*>(0))
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(PtrList<T>& lst)
:
    ptrs_(lst.ptrs_)
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(UPtrList<T>&& lst)
:
    ptrs_(std::move(lst.ptrs_))
{}


template<class T>
inline Foam::UPtrList<T>::UPtrList(const Xfer<UPtrList<T>>& lst)
{
    transfer(lst());
}


template<class T>
inline Foam::UPtrList<T>::UPtrList(UPtrList<T>& lst, bool reuse)
:
    ptrs_(lst.ptrs_, reuse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::UPtrList<T>::size() const
{
    return ptrs_.size();
}


template<class T>
inline bool Foam::UPtrList<T>::empty() const
{
    return ptrs_.empty();
}


template<class T>
inline void Foam::UPtrList<T>::clear()
{
    ptrs_.clear();
}


template<class T>
inline void Foam::UPtrList<T>::swap(UPtrList<T>& lst)
{
    ptrs_.swap(lst.ptrs_);
}


template<class T>
inline void Foam::UPtrList<T>::transfer(UPtrList<T>& lst)
{
    ptrs_.transfer(lst.ptrs_);
}


template<class T>
inline T& Foam::UPtrList<T>::first()
{
    return this->operator[](0);
}


template<class T>
inline const T& Foam::UPtrList<T>::first() const
{
    return this->operator[](0);
}


template<class T>
inline T& Foam::UPtrList<T>::last()
{
    return this->operator[](this->size()-1);
}


template<class T>
inline const T& Foam::UPtrList<T>::last() const
{
    return this->operator[](this->size()-1);
}


template<class T>
inline void Foam::UPtrList<T>::setSize(const label newLen)
{
    if (newLen <= 0)
    {
        clear();
        return;
    }

    const label oldLen = this->size();
    if (newLen != oldLen)
    {
        // Truncate or extend (new elements initialized to nullptr)
        ptrs_.setSize(newLen, reinterpret_cast<T*>(0));
    }
}


template<class T>
inline void Foam::UPtrList<T>::resize(const label nElem)
{
    this->setSize(nElem);
}


template<class T>
inline void Foam::UPtrList<T>::append(T* ptr)
{
    const label idx = this->size();
    ptrs_.setSize(idx + 1);
    ptrs_[idx] = ptr;
}


template<class T>
inline bool Foam::UPtrList<T>::set(const label i) const
{
    return ptrs_[i] != nullptr;
}


template<class T>
inline T* Foam::UPtrList<T>::set(const label i, T* ptr)
{
    T* old = ptrs_[i];
    ptrs_[i] = ptr;
    return old;
}


template<class T>
inline Foam::Xfer<Foam::UPtrList<T>> Foam::UPtrList<T>::xfer()
{
    return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::UPtrList<T>::operator[](const label i) const
{
    if (!ptrs_[i])
    {
        FatalErrorInFunction
            << "hanging pointer at index " << i
            << " (size " << size()
            << "), cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


template<class T>
inline T& Foam::UPtrList<T>::operator[](const label i)
{
    if (!ptrs_[i])
    {
        FatalErrorInFunction
            << "hanging pointer at index " << i
            << " (size " << size()
            << "), cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


template<class T>
inline const T* Foam::UPtrList<T>::operator()(const label i) const
{
    return ptrs_[i];
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T>
inline Foam::UPtrList<T>::iterator::iterator(T** ptr)
:
    ptr_(ptr)
{}


template<class T>
inline bool Foam::UPtrList<T>::iterator::operator==(const iterator& iter) const
{
    return ptr_ == iter.ptr_;
}


template<class T>
inline bool Foam::UPtrList<T>::iterator::operator!=(const iterator& iter) const
{
    return ptr_ != iter.ptr_;
}


template<class T>
inline T& Foam::UPtrList<T>::iterator::operator*() const
{
    return **ptr_;
}


template<class T>
inline T& Foam::UPtrList<T>::iterator::operator()() const
{
    return **ptr_;
}


template<class T>
inline typename Foam::UPtrList<T>::iterator&
Foam::UPtrList<T>::iterator::operator++()
{
    ++ptr_;
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator++(int)
{
    iterator iter(*this);
    ++ptr_;
    return iter;
}


template<class T>
inline typename Foam::UPtrList<T>::iterator&
Foam::UPtrList<T>::iterator::operator--()
{
    --ptr_;
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator--(int)
{
    iterator iter(*this);
    --ptr_;
    return iter;
}


template<class T>
inline typename Foam::UPtrList<T>::iterator&
Foam::UPtrList<T>::iterator::operator+=(label n)
{
    ptr_ += n;
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::iterator&
Foam::UPtrList<T>::iterator::operator-=(label n)
{
    ptr_ -= n;
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator+(label n) const
{
    return iterator(ptr_ + n);
}


template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::iterator::operator-(label n) const
{
    return iterator(ptr_ - n);
}


template<class T>
inline Foam::label
Foam::UPtrList<T>::iterator::operator-(const iterator& iter) const
{
    return (ptr_ - iter.ptr_);
}


template<class T>
inline T& Foam::UPtrList<T>::iterator::operator[](label n) const
{
    return *(*this + n);
}


template<class T>
inline bool Foam::UPtrList<T>::iterator::operator<(const iterator& iter) const
{
    return ptr_ < iter.ptr_;
}


template<class T>
inline bool Foam::UPtrList<T>::iterator::operator>(const iterator& iter) const
{
    return ptr_ > iter.ptr_;
}


template<class T>
inline bool Foam::UPtrList<T>::iterator::operator<=(const iterator& iter) const
{
    return ptr_ <= iter.ptr_;
}


template<class T>
inline bool Foam::UPtrList<T>::iterator::operator>=(const iterator& iter) const
{
    return ptr_ >= iter.ptr_;
}


// * * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * //

template<class T>
inline Foam::UPtrList<T>::const_iterator::const_iterator(const T* const* ptr)
:
    ptr_(ptr)
{}


template<class T>
inline Foam::UPtrList<T>::const_iterator::const_iterator(const iterator& iter)
:
    ptr_(iter.ptr_)
{}


template<class T>
inline bool Foam::UPtrList<T>::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return ptr_ == iter.ptr_;
}


template<class T>
inline bool Foam::UPtrList<T>::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return ptr_ != iter.ptr_;
}


template<class T>
inline const T& Foam::UPtrList<T>::const_iterator::operator*() const
{
    return **ptr_;
}


template<class T>
inline const T& Foam::UPtrList<T>::const_iterator::operator()() const
{
    return **ptr_;
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator&
Foam::UPtrList<T>::const_iterator::operator++()
{
    ++ptr_;
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::const_iterator::operator++(int)
{
    const_iterator iter(*this);
    ++ptr_;
    return iter;
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator&
Foam::UPtrList<T>::const_iterator::operator--()
{
    --ptr_;
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::const_iterator::operator--(int)
{
    const_iterator iter(*this);
    --ptr_;
    return iter;
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator&
Foam::UPtrList<T>::const_iterator::operator+=(label n)
{
    ptr_ += n;
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator&
Foam::UPtrList<T>::const_iterator::operator-=(label n)
{
    ptr_ -= n;
    return *this;
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::const_iterator::operator+(label n) const
{
    return const_iterator(ptr_ + n);
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::const_iterator::operator-(label n) const
{
    return const_iterator(ptr_ - n);
}


template<class T>
inline Foam::label
Foam::UPtrList<T>::const_iterator::operator-(const const_iterator& iter) const
{
    return (ptr_ - iter.ptr_);
}


template<class T>
inline const T&
Foam::UPtrList<T>::const_iterator::operator[](label n) const
{
    return *(*this + n);
}


template<class T>
inline bool Foam::UPtrList<T>::const_iterator::operator<
(
    const const_iterator& iter
) const
{
    return ptr_ < iter.ptr_;
}


template<class T>
inline bool Foam::UPtrList<T>::const_iterator::operator>
(
    const const_iterator& iter
) const
{
    return ptr_ > iter.ptr_;
}


template<class T>
inline bool Foam::UPtrList<T>::const_iterator::operator<=
(
    const const_iterator& iter
) const
{
    return ptr_ <= iter.ptr_;
}


template<class T>
inline bool Foam::UPtrList<T>::const_iterator::operator>=
(
    const const_iterator& iter
) const
{
    return ptr_ >= iter.ptr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::begin()
{
    return ptrs_.begin();
}


template<class T>
inline typename Foam::UPtrList<T>::iterator
Foam::UPtrList<T>::end()
{
    return ptrs_.end();
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::cbegin() const
{
    return ptrs_.cbegin();
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::cend() const
{
    return ptrs_.cend();
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::begin() const
{
    return ptrs_.begin();
}


template<class T>
inline typename Foam::UPtrList<T>::const_iterator
Foam::UPtrList<T>::end() const
{
    return ptrs_.end();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::UPtrList<T>::operator=(UPtrList<T>&& lst)
{
    this->clear();
    this->swap(lst);
}


// ************************************************************************* //
