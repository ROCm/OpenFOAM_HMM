/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

template<class T>
inline void Foam::List<T>::alloc()
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];
    }
}


template<class T>
inline void Foam::List<T>::reAlloc(const label s)
{
    if (this->size_ != s)
    {
        clear();
        this->size_ = s;
        alloc();
    }
}


template<class T>
template<class List2>
inline void Foam::List<T>::copyList(const List2& lst)
{
    if (this->size_)
    {
        forAll(*this, i)
        {
            this->operator[](i) = lst[i];
        }
    }
}


template<class T>
template<class List2>
inline void Foam::List<T>::allocCopyList(const List2& lst)
{
    if (this->size_)
    {
        alloc();
        copyList(lst);
    }
}


template<class T>
template<class InputIterator>
inline Foam::List<T>::List
(
    InputIterator begIter,
    InputIterator endIter,
    const label s
)
:
    UList<T>(nullptr, s)
{
    if (this->size_)
    {
        alloc();

        InputIterator iter = begIter;
        for (label i = 0; i < s; ++i)
        {
            this->operator[](i) = *iter;
            ++iter;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::List<T>::List()
{}


template<class T>
inline Foam::autoPtr<Foam::List<T>> Foam::List<T>::clone() const
{
    return autoPtr<List<T>>(new List<T>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::List<T>& Foam::List<T>::null()
{
    return NullObjectRef<List<T>>();
}


template<class T>
inline void Foam::List<T>::clear()
{
    if (this->v_)
    {
        delete[] this->v_;
        this->v_ = nullptr;
    }

    this->size_ = 0;
}


template<class T>
inline void Foam::List<T>::resize(const label newSize)
{
    this->setSize(newSize);
}


template<class T>
inline void Foam::List<T>::resize(const label newSize, const T& val)
{
    this->setSize(newSize, val);
}


template<class T>
inline T& Foam::List<T>::newElmt(const label i)
{
    const label n = this->size();

    if (i >= n)
    {
        setSize(2*n);
    }

    return UList<T>::operator[](i);
}


template<class T>
inline Foam::Xfer<Foam::List<T>> Foam::List<T>::xfer()
{
    return xferMove(*this);
}


template<class T>
inline void Foam::List<T>::append(const T& val)
{
    setSize(this->size() + 1, val);  // copy element
}


template<class T>
inline void Foam::List<T>::append(T&& val)
{
    const label idx = this->size();
    setSize(idx + 1);

    this->operator[](idx) = std::move(val);  // move assign element
}


template<class T>
inline void Foam::List<T>::append(const UList<T>& lst)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "attempted appending to self" << abort(FatalError);
    }

    label idx = this->size();
    const label n = lst.size();

    setSize(idx + n);

    for (label i=0; i<n; ++i)
    {
        this->operator[](idx++) = lst[i];  // copy element
    }
}


template<class T>
inline void Foam::List<T>::append(const UIndirectList<T>& lst)
{
    label idx = this->size();
    const label n = lst.size();

    setSize(idx + n);

    for (label i=0; i<n; ++i)
    {
        this->operator[](idx++) = lst[i];  // copy element
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::operator=(const T& val)
{
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::List<T>::operator=(const zero)
{
    UList<T>::operator=(Zero);
}


// ************************************************************************* //
