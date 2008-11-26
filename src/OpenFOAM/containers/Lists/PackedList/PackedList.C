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

#include "PackedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int nBits>
Foam::PackedList<nBits>::PackedList(const label size, const unsigned int val)
:
    List<unsigned int>(intSize(size)),
    size_(size)
{
    operator=(val);
}


template<int nBits>
Foam::PackedList<nBits>::PackedList(const PackedList<nBits>& lst)
:
    List<unsigned int>(lst),
    size_(lst.size())
{}


template<int nBits>
Foam::PackedList<nBits>::PackedList(const xfer<PackedList<nBits> >& lst)
{
    transfer(lst());
}


template<int nBits>
Foam::PackedList<nBits>::PackedList(const UList<label>& lst)
:
    List<unsigned int>(intSize(lst.size()), 0),
    size_(lst.size())
{
    forAll(lst, i)
    {
        set(i, lst[i]);
    }
}


template<int nBits>
Foam::autoPtr<Foam::PackedList<nBits> > Foam::PackedList<nBits>::clone() const
{
    return autoPtr<PackedList<nBits> >(new PackedList<nBits>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int nBits>
void Foam::PackedList<nBits>::setSize(const label size)
{
    List<unsigned int>::setSize(intSize(size));
    size_ = size;
}


template<int nBits>
void Foam::PackedList<nBits>::clear()
{
    List<unsigned int>::clear();
    size_ = 0;
}


template<int nBits>
void Foam::PackedList<nBits>::transfer(PackedList<nBits>& lst)
{
    size_ = lst.size();
    List<unsigned int>::transfer(lst);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Assignment.
template<int nBits>
void Foam::PackedList<nBits>::operator=(const PackedList<nBits>& lst)
{
    setSize(lst.size());
    List<unsigned int>::operator=(lst);
}


template<int nBits>
Foam::labelList Foam::PackedList<nBits>::operator()() const
{
    labelList elems(size());

    forAll(*this, i)
    {
        elems[i] = get(i);
    }
    return elems;
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

//template<int nBits>
//Foam::Ostream& ::Foam::operator<<(Ostream& os, const PackedList<nBits>& lst)
//{
//    os << lst();
//    return os;
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
