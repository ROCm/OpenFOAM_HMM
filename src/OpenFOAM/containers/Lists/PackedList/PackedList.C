/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "PackedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned nBits>
Foam::PackedList<nBits>::PackedList(const label size, const unsigned int val)
:
    StorageList(packedLength(size), 0u),
    size_(size)
{
    if (val)
    {
        operator=(val);
    }
}


template<unsigned nBits>
Foam::PackedList<nBits>::PackedList(const UList<label>& lst)
:
    StorageList(packedLength(lst.size()), 0u),
    size_(lst.size())
{
    forAll(lst, i)
    {
        set(i, lst[i]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


#if (UINT_MAX == 0xFFFFFFFF)
// 32-bit counting, Hamming weight method
#   define COUNT_PACKEDBITS(sum, x)                                           \
{                                                                             \
    x -= (x >> 1) & 0x55555555;                                               \
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);                           \
    sum += (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;                \
}
#elif (UINT_MAX == 0xFFFFFFFFFFFFFFFF)
// 64-bit counting, Hamming weight method
#   define COUNT_PACKEDBITS(sum, x)                                           \
{                                                                             \
    x -= (x >> 1) & 0x5555555555555555;                                       \
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);           \
    sum += (((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F) * 0x0101010101010101) >> 56;\
}
#else
// Arbitrary number of bits, Brian Kernighan's method
#   define COUNT_PACKEDBITS(sum, x)    for (; x; ++sum) { x &= x - 1; }
#endif


template<unsigned nBits>
unsigned int Foam::PackedList<nBits>::count() const
{
    register unsigned int c = 0;

    if (size_)
    {
        const label packLen = packedLength(size_);

        for (label i = 0; i < packLen; ++i)
        {
            register unsigned int bits = StorageList::operator[](i);
            COUNT_PACKEDBITS(c, bits);
        }
    }

    return c;
}


template<unsigned nBits>
bool Foam::PackedList<nBits>::trim()
{
    if (!size_)
    {
        return false;
    }

    const label oldSize = size_;

    for (label storeI = packedLength(size_) - 1; storeI >= 0; --storeI)
    {
        size_ = storeI * packing();
        unsigned int bits = StorageList::operator[](storeI);

        // found some bits
        if (bits)
        {
            while (bits)
            {
                bits >>= nBits;
                ++size_;
            }
            break;
        }
    }

    return (size_ != oldSize);
}


template<unsigned nBits>
void Foam::PackedList<nBits>::flip()
{
    if (!size_)
    {
        return;
    }

    // mask value for complete segments
    const unsigned int mask = maskLower(packing());

    const label packLen = packedLength(size_);
    for (label i=0; i < packLen; ++i)
    {
        StorageList::operator[](i) = mask & ~StorageList::operator[](i);
    }

    // mask off the final partial segment
    {
        const unsigned int off = size_ % packing();
        if (off)
        {
            const unsigned int seg = size_ / packing();

            StorageList::operator[](seg) &= maskLower(off);
        }
    }
}


template<unsigned nBits>
Foam::Xfer<Foam::labelList> Foam::PackedList<nBits>::values() const
{
    labelList elems(size_);

    forAll(*this, i)
    {
        elems[i] = get(i);
    }

    return elems.xfer();
}


template<unsigned nBits>
Foam::Ostream& Foam::PackedList<nBits>::iteratorBase::print(Ostream& os) const
{
    os  << "iterator<"  << label(nBits) << "> ["
        << this->index_ << "]"
        << " segment:"  << label(this->index_ / packing())
        << " offset:"   << label(this->index_ % packing())
        << " value:"    << this->get()
        << nl;

    return os;
}


template<unsigned nBits>
Foam::Ostream& Foam::PackedList<nBits>::print
(
    Ostream& os,
    const bool fullOutput
) const
{
    const label packLen = packedLength(size_);

    os  << "PackedList<" << nBits << ">"
        << " max_value:" << max_value()
        << " packing:"   << packing() << nl
        << " count: "     << count() << nl
        << " size/capacity: " << size_ << "/" << capacity() << nl
        << " storage/capacity: " << packLen << "/" << StorageList::size()
        << "\n(\n";

    // mask value for complete segments
    unsigned int mask = maskLower(packing());
    const label outputLen = fullOutput ? StorageList::size() : packLen;

    for (label i=0; i < outputLen; ++i)
    {
        const StorageType& rawBits = StorageList::operator[](i);

        // the final segment may not be full, modify mask accordingly
        if (i == packLen-1)
        {
            const unsigned int off = size_ % packing();

            if (off)
            {
                mask = maskLower(off);
            }
        }
        else if (i == packLen)
        {
            // no mask for unaddressed bit
            mask = 0u;
        }


        for (unsigned int testBit = (1u << max_bits()); testBit; testBit >>= 1)
        {
            if (mask & testBit)
            {
                // addressable region
                if (rawBits & testBit)
                {
                    os  << '1';
                }
                else
                {
                    os  << '-';
                }
            }
            else
            {
                if (rawBits & testBit)
                {
                    os  << '!';
                }
                else
                {
                    os  << '.';
                }
            }
        }
        os  << '\n';
    }
    os << ")\n";

    return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<unsigned nBits>
Foam::PackedList<nBits>&
Foam::PackedList<nBits>::operator=(const PackedList<nBits>& lst)
{
    StorageList::operator=(lst);
    size_ = lst.size();
    return *this;
}


template<unsigned nBits>
Foam::PackedList<nBits>&
Foam::PackedList<nBits>::operator=(const UList<label>& lst)
{
    setCapacity(lst.size());
    size_ = lst.size();

    forAll(lst, i)
    {
        set(i, lst[i]);
    }
    return *this;
}


// ************************************************************************* //
