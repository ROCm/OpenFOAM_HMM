/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "CompactListList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::CompactListList<T>::CompactListList(const List<List<T> >& ll)
:
    size_(ll.size()),
    offsets_(ll.size()+1)
{
    label sumSize = 0;
    offsets_[0] = 0;
    forAll(ll, i)
    {
        sumSize += ll[i].size();
        offsets_[i+1] = sumSize;
    }

    m_.setSize(sumSize);

    label k = 0;
    forAll(ll, i)
    {
        const List<T>& lli = ll[i];

        forAll(lli, j)
        {
            m_[k++] = lli[j];
        }
    }
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const UList<label>& rowSizes
)
:
    size_(rowSizes.size()),
    offsets_(rowSizes.size()+1)
{
    label sumSize = 0;
    offsets_[0] = 0;
    forAll(rowSizes, i)
    {
        sumSize += rowSizes[i];
        offsets_[i+1] = sumSize;
    }

    m_.setSize(sumSize);
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const UList<label>& rowSizes,
    const T& t
)
:
    size_(rowSizes.size()),
    offsets_(rowSizes.size()+1)
{
    label sumSize = 0;
    offsets_[0] = 0;
    forAll(rowSizes, i)
    {
        sumSize += rowSizes[i];
        offsets_[i+1] = sumSize;
    }

    m_.setSize(sumSize, t);
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const Xfer<CompactListList<T> >& lst
)
{
    transfer(lst());
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    CompactListList<T>& lst,
    bool reUse
)
:
    size_(lst.size()),
    offsets_(lst.offsets_, reUse),
    m_(lst.m_, reUse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::CompactListList<T>::setSize(const label nRows)
{
    if (nRows == 0)
    {
        clear();
    }
    if (nRows < size())
    {
        size_ = nRows;
        offsets_.setSize(nRows+1);
        m_.setSize(offsets_[nRows]);
    }
    else if (nRows > size())
    {
        FatalErrorIn("CompactListList<T>::setSize(const label nRows)")
            << "Cannot be used to extend the list from " << offsets_.size()
            << " to " << nRows << nl
            << "    Please use one of the other setSize member functions"
            << abort(FatalError);
    }
}


template<class T>
void Foam::CompactListList<T>::setSize
(
    const label nRows,
    const label nData
)
{
    size_ = nRows;
    offsets_.setSize(nRows+1);
    m_.setSize(nData);
}


template<class T>
void Foam::CompactListList<T>::setSize
(
    const label nRows,
    const label nData,
    const T& t
)
{
    size_ = nRows;
    offsets_.setSize(nRows+1);
    m_.setSize(nData, t);
}


template<class T>
void Foam::CompactListList<T>::setSize(const UList<label>& rowSizes)
{
    size_ = rowSizes.size();
    offsets_.setSize(rowSizes.size()+1);

    label sumSize = 0;
    offsets_[0] = 0;
    forAll(rowSizes, i)
    {
        sumSize += rowSizes[i];
        offsets_[i+1] = sumSize;
    }

    m_.setSize(sumSize);
}


template<class T>
Foam::labelList Foam::CompactListList<T>::sizes() const
{
    labelList rowSizes(size());

    if (rowSizes.size() > 0)
    {
        forAll(rowSizes, i)
        {
            rowSizes[i] = offsets_[i+1] - offsets_[i];
        }
    }
    return rowSizes;
}


template<class T>
void Foam::CompactListList<T>::clear()
{
    size_ = 0;
    offsets_.clear();
    m_.clear();
}


template<class T>
void Foam::CompactListList<T>::transfer(CompactListList<T>& a)
{
    size_ = a.size_;
    offsets_.transfer(a.offsets_);
    m_.transfer(a.m_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
Foam::List<Foam::List<T> > Foam::CompactListList<T>::operator()() const
{
    List<List<T> > ll(size());

    forAll(ll, i)
    {
        ll[i] = operator[](i);
    }

    return ll;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CompactListListIO.C"

// ************************************************************************* //
