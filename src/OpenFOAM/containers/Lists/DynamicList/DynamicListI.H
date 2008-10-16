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

//- Construct null
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList()
:
    List<T>(SizeInc),
    allocSize_(SizeInc)
{
    List<T>::size(0);
}


//- Construct given size
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


//- Construct given size
template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
    const UList<T>& s
)
:
    List<T>(s),
    allocSize_(s.size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::label Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::allocSize()
const
{
    return allocSize_;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label s
)
{
    label nextFree = List<T>::size();
    if (s <= nextFree)
    {
        // adjust addressed size, leave allocated size untouched
        nextFree = s;
    }
    else
    {
        // adjust allocated size, leave addressed size untouched
        allocSize_ = s;
        List<T>::setSize(allocSize_);
    }
    List<T>::size(nextFree);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
    const label s,
    const T& t
)
{
    if (s <= List<T>::size())
    {
        // shrink addressable size, leave allocated size untouched
        List<T>::size() = s;
    }
    else
    {
        label nextFree = List<T>::size();
        allocSize_ = s;
        List<T>::setSize(allocSize_, t);
        List<T>::size() = nextFree;
=======
>>>>>>> olesenm:src/OpenFOAM/containers/Lists/DynamicList/DynamicListI.H
    }
    List<T>::size(nextFree);
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
    allocSize_ = List<T>::size();
    List<T>::setSize(allocSize_);
    return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer(List<T>& l)
{
    if (l.size() < List<T>::size())
    {
        FatalErrorIn
        (
            "void DynamicList<T, SizeInc, SizeMult"
            ", SizeDiv>::transfer(List<T>&)"
        )   << "Cannot replace the underlying storage of this DynamicList"
            << " of which " << List<T>::size() << " elements are used" << nl
            << "with a List of size " << l.size() << abort(FatalError);
    }
    else
    {
        allocSize_ = l.size();
        List<T>::transfer(l);   // take over storage
    }
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer
(
    DynamicList<T, SizeInc, SizeMult, SizeDiv>& l
)
{
    allocSize_ = l.allocSize();
    List<T>::transfer(l);       // take over storage. Null l.
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append(const T& e)
{
    // Work on copy free index since gets overwritten by setSize
    label nextFree = List<T>::size();

    nextFree++;

    if (nextFree > allocSize_)
    {
        allocSize_ = max
        (
            nextFree,
            label(SizeMult*allocSize_/SizeDiv + SizeInc)
        );
        List<T>::setSize(allocSize_);
    }

    List<T>::size(nextFree);

    this->operator[](nextFree - 1) = e;
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

    if (nextFree > allocSize_)
    {
        allocSize_ = max
        (
            nextFree,
            label(SizeMult*allocSize_/SizeDiv + SizeInc)
        );
        List<T>::setSize(allocSize_);
    }

    List<T>::size(nextFree);

    return this->operator[](i);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const T& t
)
{
    List<T>::operator=(t);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const List<T>& l
)
{
    List<T>::operator=(l);
    allocSize_ = List<T>::size();
}


// ************************************************************************* //
