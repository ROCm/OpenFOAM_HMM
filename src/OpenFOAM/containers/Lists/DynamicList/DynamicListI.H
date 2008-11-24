/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList()
:
    List<T>(SizeInc),
    allocSize_(SizeInc)
{
    List<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const label s
)
:
    List<T>(s),
    allocSize_(s)
{
    List<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const UList<T>& lst
)
:
    List<T>(lst),
    allocSize_(lst.size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::label Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::capacity()
const
{
    return allocSize_;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setCapacity
(
    const label s
)
{
    label nextFree = List<T>::size();
    allocSize_ = s;

    if (nextFree > allocSize_)
    {
        // truncate both allocated and addressed sizes
        nextFree = allocSize_;
    }

    List<T>::setSize(allocSize_);
    List<T>::size(nextFree);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::reserve
(
    const label s
)
{
    if (s > allocSize_)
    {
        allocSize_ = max
        (
            s,
            label(SizeInc + allocSize_ * SizeMult / SizeDiv)
        );

        // adjust allocated size, leave addressed size untouched
        label nextFree = List<T>::size();
        List<T>::setSize(allocSize_);
        List<T>::size(nextFree);
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label s
)
{
    // allocate more space?
    if (s > allocSize_)
    {
        allocSize_ = max
        (
            s,
            label(SizeInc + allocSize_ * SizeMult / SizeDiv)
        );

        List<T>::setSize(allocSize_);
    }

    // adjust addressed size
    List<T>::size(s);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label s,
    const T& t
)
{
    label nextFree = List<T>::size();
    setSize(s);

    // set new elements to constant value
    while (nextFree < s)
    {
        this->operator[](nextFree++) = t;
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::clear()
{
    List<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::clearStorage()
{
    List<T>::clear();
    allocSize_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::shrink()
{
    label nextFree = List<T>::size();
    if (allocSize_ > nextFree)
    {
        // use the full list when resizing
        List<T>::size(allocSize_);
        // the new size
        allocSize_ = nextFree;
        List<T>::setSize(allocSize_);
        List<T>::size(nextFree);
    }
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer(List<T>& lst)
{
    allocSize_ = lst.size();
    List<T>::transfer(lst);   // take over storage, clear addressing for lst.
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer
(
    DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    // take over storage as-is (without shrink), clear addressing for lst.
    allocSize_ = lst.allocSize_;
    lst.allocSize_ = 0;

    List<T>::transfer(static_cast<List<T>&>(lst));
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
    const T& e
)
{
    // Work on copy free index since the size gets overwritten by setSize
    label nextFree = List<T>::size();

    reserve(nextFree+1);
    List<T>::size(nextFree+1);

    this->operator[](nextFree) = e;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
    const UList<T>& lst
)
{
    label nextFree = List<T>::size();

    if (this == &lst)
    {
        FatalErrorIn
        (
            "DynamicList<T, SizeInc, SizeMult, SizeDiv>::append(UList<T>&)"
        )   << "attempted appending to self" << abort(FatalError);
    }

    reserve(nextFree + lst.size());
    List<T>::size(nextFree + lst.size());

    forAll(lst, elemI)
    {
        this->operator[](nextFree++) = lst[elemI];
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::remove()
{
    if (List<T>::size() == 0)
    {
        FatalErrorIn
        (
            "Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::remove()"
        )   << "List is empty" << abort(FatalError);
    }

    label nextFree = List<T>::size()-1;

    const T& val = List<T>::operator[](nextFree);

    List<T>::size(nextFree);

    return val;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T& Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator()
(
    const label i
)
{
    label nextFree = List<T>::size();
    nextFree = max(nextFree, i + 1);

    reserve(nextFree);
    List<T>::size(nextFree);

    return this->operator[](i);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const T& t
)
{
    UList<T>::operator=(t);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const UList<T>& lst
)
{
    if (this == &lst)
    {
        FatalErrorIn
        (
            "DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=(UList<T>&)"
        )   << "attempted assignment to self" << abort(FatalError);
    }

    if (allocSize_ >= lst.size())
    {
        // can copy without reallocation, match sizes to avoid reallocation
        List<T>::size(lst.size());
        List<T>::operator=(lst);
    }
    else
    {
        // make the entire storage available for the copy operation
        List<T>::size(allocSize_);

        List<T>::operator=(lst);
        allocSize_ = List<T>::size();
    }
}


// ************************************************************************* //
