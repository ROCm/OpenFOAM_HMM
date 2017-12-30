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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
template<class ListType>
inline void Foam::DynamicField<T, SizeMin>::assignDynField
(
    const ListType& lst
)
{
    const label newSize = lst.size();

    if (capacity_ >= newSize)
    {
        // Can copy w/o reallocating - adjust addressable size accordingly.
        Field<T>::size(lst.size());
        Field<T>::operator=(lst);
    }
    else
    {
        // Ensure list size consistency prior to copying.
        Field<T>::size(capacity_);

        Field<T>::operator=(lst);
        capacity_ = Field<T>::size();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField()
:
    Field<T>(0),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const label nElem
)
:
    Field<T>(nElem),
    capacity_(Field<T>::size())
{
    // We could also enforce sizing granularity
    Field<T>::size(0);
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& lst
)
:
    Field<T>(lst),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const Xfer<List<T>>& lst
)
:
    Field<T>(lst),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const Xfer<Field<T>>& lst
)
:
    Field<T>(lst),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& mapF,
    const labelList& mapAddressing
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


//- Construct by mapping from the given field
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
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const DynamicField<T, SizeMin>& lst
)
:
    Field<T>(lst),
    capacity_(lst.capacity())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const Xfer<DynamicField<T, SizeMin>>& lst
)
:
    Field<T>(lst),
    capacity_(Field<T>::size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::label Foam::DynamicField<T, SizeMin>::capacity()
const
{
    return capacity_;
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::setCapacity
(
    const label nElem
)
{
    label nextFree = Field<T>::size();
    capacity_ = nElem;

    if (nextFree > capacity_)
    {
        // truncate addressed sizes too
        nextFree = capacity_;
    }

    // We could also enforce sizing granularity

    Field<T>::setSize(capacity_);
    Field<T>::size(nextFree);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::reserve
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
                label(2*capacity_)
            )
        );

        // Adjust allocated size, leave addressed size untouched
        const label nextFree = Field<T>::size();
        Field<T>::setSize(capacity_);
        Field<T>::size(nextFree);
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::setSize
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
                label(2*capacity_)
            )
        );

        Field<T>::setSize(capacity_);
    }

    // Adjust addressed size
    Field<T>::size(nElem);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::setSize
(
    const label nElem,
    const T& val
)
{
    label nextFree = Field<T>::size();
    setSize(nElem);

    // Set new elements to constant value
    while (nextFree < nElem)
    {
        this->operator[](nextFree++) = val;
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::resize
(
    const label nElem
)
{
    this->setSize(nElem);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::resize
(
    const label nElem,
    const T& val
)
{
    this->setSize(nElem, val);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::clear()
{
    Field<T>::size(0);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::clearStorage()
{
    Field<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline Foam::label Foam::DynamicField<T, SizeMin>::expandStorage()
{
    const label nextFree = Field<T>::size();

    // Allow addressing into the entire list
    Field<T>::size(capacity_);

    return nextFree;
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>&
Foam::DynamicField<T, SizeMin>::shrink()
{
    label nextFree = Field<T>::size();
    if (capacity_ > nextFree)
    {
        // Use the full list when resizing
        Field<T>::size(capacity_);

        // The new size
        capacity_ = nextFree;
        Field<T>::setSize(capacity_);
        Field<T>::size(nextFree);
    }
    return *this;
}


template<class T, int SizeMin>
inline Foam::Xfer<Foam::List<T>>
Foam::DynamicField<T, SizeMin>::xfer()
{
    return xferMoveTo<List<T>>(*this);
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>&
Foam::DynamicField<T, SizeMin>::append
(
    const T& val
)
{
    const label idx = List<T>::size();
    setSize(idx + 1);

    this->operator[](idx) = val;
    return *this;
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>&
Foam::DynamicField<T, SizeMin>::append
(
    const UList<T>& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "attempted appending to self" << abort(FatalError);
    }

    label nextFree = List<T>::size();
    setSize(nextFree + lst.size());

    forAll(lst, i)
    {
        this->operator[](nextFree++) = lst[i];
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

    List<T>::size(idx);

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
        setSize(i + 1);
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
    const DynamicField<T, SizeMin>& lst
)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "Attempted assignment to self" << abort(FatalError);
    }

    assignDynField(lst);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const UList<T>& lst
)
{
    assignDynField(lst);
}


// ************************************************************************* //
