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

#include "PackedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int nBits>
Foam::PackedList<nBits>::PackedList(const label size, const unsigned int val)
:
    List<PackedStorage>(packedLength(size), 0u),
    size_(size)
{
    operator=(val);
}


template<int nBits>
Foam::PackedList<nBits>::PackedList(const UList<label>& lst)
:
    List<PackedStorage>(packedLength(lst.size()), 0u),
    size_(lst.size())
{
    forAll(lst, i)
    {
        set(i, lst[i]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int nBits>
Foam::labelList Foam::PackedList<nBits>::values() const
{
    labelList elems(size());

    forAll(*this, i)
    {
        elems[i] = get(i);
    }
    return elems;
}


template<int nBits>
Foam::Ostream& Foam::PackedList<nBits>::const_iterator::print(Ostream& os) const
{
    os  << "iterator<" << nBits << "> ["
        << (index_ * packing() + offset_) << "]"
        << " index:" << index_ << " offset:" << offset_
        << " value:" << unsigned(*this)
        << nl;

    return os;
}


template<int nBits>
Foam::Ostream& Foam::PackedList<nBits>::print(Ostream& os) const
{
    os  << "PackedList<" << nBits << ">"
        << " max_value:" << max_value()
        << " packing:"   << packing() << nl
        << "values: " << size() << "/" << capacity() << "( ";
    forAll(*this, i)
    {
        os << get(i) << ' ';
    }

    label packLen = packedLength(size());

    os  << ")\n"
        << "storage: " << packLen << "/" << storage().size() << "( ";

    // mask for the valid bits
    unsigned int validBits = max_value();
    for (unsigned int i = 1; i < packing(); ++i)
    {
        validBits |= (validBits << nBits);
    }

    for (label i=0; i < packLen; i++)
    {
        const PackedStorage& rawBits = storage()[i];

        // the final storage may not be full, modify validBits accordingly
        if (i+1 == packLen)
        {
            label junk = size() % packing();

            if (junk)
            {
                junk = packing() - junk;
            }

            for (label j=0; j < junk; j++)
            {
                validBits >>= nBits;
            }
        }

        for (unsigned int testBit = 0x1 << max_bits(); testBit; testBit >>= 1)
        {
            if (testBit & validBits)
            {
                if (rawBits & testBit)
                {
                    os << '1';
                }
                else
                {
                    os << '0';
                }
            }
            else
            {
                os << '.';
            }
        }
        cout << ' ';
    }
    os << ")\n";

    return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<int nBits>
void Foam::PackedList<nBits>::operator=(const PackedList<nBits>& lst)
{
    setCapacity(lst.size());
    List<PackedStorage>::operator=(lst);
}


template<int nBits>
void Foam::PackedList<nBits>::operator=(const UList<label>& lst)
{
    setCapacity(lst.size());

    forAll(lst, i)
    {
        set(i, lst[i]);
    }
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
