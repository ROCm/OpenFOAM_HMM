/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "PackedBoolList.H"
#include "IOstreams.H"
#include <cctype>

// * * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * * * //

// table-lookup instead of printf("%02x", char)
//! @cond localScope
static const char hexLookup[] =
{
    "000102030405060708090a0b0c0d0e0f"
    "101112131415161718191a1b1c1d1e1f"
    "202122232425262728292a2b2c2d2e2f"
    "303132333435363738393a3b3c3d3e3f"
    "404142434445464748494a4b4c4d4e4f"
    "505152535455565758595a5b5c5d5e5f"
    "606162636465666768696a6b6c6d6e6f"
    "707172737475767778797a7b7c7d7e7f"
    "808182838485868788898a8b8c8d8e8f"
    "909192939495969798999a9b9c9d9e9f"
    "a0a1a2a3a4a5a6a7a8a9aaabacadaeaf"
    "b0b1b2b3b4b5b6b7b8b9babbbcbdbebf"
    "c0c1c2c3c4c5c6c7c8c9cacbcccdcecf"
    "d0d1d2d3d4d5d6d7d8d9dadbdcdddedf"
    "e0e1e2e3e4e5e6e7e8e9eaebecedeeef"
    "f0f1f2f3f4f5f6f7f8f9fafbfcfdfeff"
};


// number of bytes required to pack nElem
static inline unsigned int packedBytes(const Foam::label nElem)
{
    return (nElem + CHAR_BIT - 1) / CHAR_BIT;
}
//! @endcond


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PackedBoolList::PackedBoolList(Istream &is)
:
    PackedList<1>()
{
    is  >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList> Foam::PackedBoolList::used() const
{
    labelList lst(this->count());

    if (lst.size())
    {
        label nElem = 0;

        forAll(*this, elemI)
        {
            if (get(elemI))
            {
                lst[nElem++] = elemI;
            }
        }

        lst.setSize(nElem);
    }

    return lst.xfer();
}


template<class LabelListType>
Foam::label Foam::PackedBoolList::setFromIndices(const LabelListType& indices)
{
    // no better information, just guess something about the size
    reserve(indices.size());

    label cnt = 0;
    forAll(indices, elemI)
    {
        if (set(indices[elemI]))
        {
            ++cnt;
        }
    }

    return cnt;
}


template<class LabelListType>
Foam::label Foam::PackedBoolList::unsetFromIndices(const LabelListType& indices)
{
    label cnt = 0;
    forAll(indices, elemI)
    {
        if (unset(indices[elemI]))
        {
            ++cnt;
        }
    }

    return cnt;
}


Foam::label Foam::PackedBoolList::set(const UList<label>& indices)
{
    return setFromIndices(indices);
}


Foam::label Foam::PackedBoolList::set(const UIndirectList<label>& indices)
{
    return setFromIndices(indices);
}


Foam::label Foam::PackedBoolList::unset(const UList<label>& indices)
{
    return unsetFromIndices(indices);
}


Foam::label Foam::PackedBoolList::unset(const UIndirectList<label>& indices)
{
    return unsetFromIndices(indices);
}


void Foam::PackedBoolList::modulo(const PackedList<1>& lst)
{
    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    const label len = min(rhs.size(), lhs.size());

    for (label i=0; i < len; ++i)
    {
        lhs[i] &= ~rhs[i];
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const UList<bool>& lst)
{
    this->setSize(lst.size());

    forAll(*this, elemI)
    {
        set(elemI, lst[elemI]);
    }

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const UList<label>& indices)
{
    clear();
    set(indices);

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator=(const UIndirectList<label>& indices)
{
    clear();
    set(indices);

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator|=(const PackedList<1>& lst)
{
    // extend addressable area if needed
    if (this->size() < lst.size())
    {
        this->resize(lst.size());
    }

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    forAll(rhs, i)
    {
        lhs[i] |= rhs[i];
    }

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator&=(const PackedList<1>& lst)
{
    // shrink addressable area if needed
    if (this->size() > lst.size())
    {
        this->resize(lst.size());
    }

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    forAll(lhs, i)
    {
        lhs[i] &= rhs[i];
    }

    // trim to bits actually used
    this->trim();

    return *this;
}


Foam::PackedBoolList&
Foam::PackedBoolList::operator^=(const PackedList<1>& lst)
{
    // extend addressable area if needed
    if (this->size() < lst.size())
    {
        this->resize(lst.size());
    }

    // operate directly with the underlying storage
    StorageList& lhs = this->storage();
    const StorageList& rhs = lst.storage();

    forAll(rhs, i)
    {
        lhs[i] ^= rhs[i];
    }

    return *this;
}

// * * * * * * * * * * * * * *  Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& istr,
    PackedBoolList& lst
)
{
    // Takes into account that 'a' (or 'A') is 10
    static const label alphaOffset = toupper('A') - 10;
    // Takes into account that '0' is 0
    static const label zeroOffset = int('0');

    ISstream& is = dynamicCast<ISstream>(istr);

    lst.clear();
    is.fatalCheck("operator>>(Istream&, PackedBoolList&)");

    token firstTok(is);
    is.fatalCheck
    (
        "operator>>(Istream&, PackedBoolList&) : reading first token"
    );

    if (firstTok.isLabel())
    {
        const label sz = firstTok.labelToken();

        // Set list length to that read
        lst.resize(sz);

        // Read list contents as ASCII

        // Read beginning of contents
        const char delimiter = is.readBeginList("PackedBoolList");

        if (sz)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                // number of bytes when packed:
                unsigned int nBytes = packedBytes(lst.size());

                for (label i=0, storeI=0; i < sz; ++storeI)
                {
                    PackedList<1>::StorageType& stored = lst.storage()[storeI];

                    // byte-wise read
                    for
                    (
                        unsigned byte=0;
                        byte < sizeof(PackedList<1>::StorageType);
                        ++byte
                    )
                    {
                        PackedList<1>::StorageType result = 0;
                        char c = 0;

                        // Get next non-whitespace character
                        while (is.get(c) && isspace(c))
                        {}

                        for (label nibble=0; nibble < 2; ++nibble)
                        {
                            if (!isxdigit(c))
                            {
                                FatalIOErrorIn
                                (
                                    "operator>>(Istream&, PackedBoolList&) : "
                                    "reading first token",
                                    is
                                )
                                    << "Illegal hex digit: '" << c << "'"
                                    << exit(FatalIOError);
                            }

                            result <<= 4;

                            if (isdigit(c))
                            {
                                result += int(c) - zeroOffset;
                            }
                            else
                            {
                                result += toupper(c) - alphaOffset;
                            }

                            // Get character for the lower part of the byte
                            if (!nibble)
                            {
                                is.get(c);
                            }
                        }

                        stored |= result << (byte*CHAR_BIT);

                        if (!--nBytes)
                        {
                            break;
                        }
                    }
                    i += PackedList<1>::packing();
                }

                // trim possible trailing junk
                // mask off the final partial segment
                {
                    const unsigned int off = sz % PackedList<1>::packing();
                    if (off)
                    {
                        const unsigned int seg = sz / PackedList<1>::packing();

                        lst.storage()[seg] &= PackedList<1>::maskLower(off);
                    }
                }


                // skip over all trailing whitespace and zeroes
                char c = 0;
                while (is.get(c) && (isspace(c) || c == '0'))
                {}

                // put back for later readEndList() to deal with
                is.putback(c);
            }
            else
            {
                const label val = readLabel(is);

                is.fatalCheck
                (
                    "operator>>(Istream&, PackedBoolList&) : "
                    "reading the single entry"
                );

                lst = val;
            }
        }

        // Read end of contents
        is.readEndList("PackedList");
    }
    else if (firstTok.isPunctuation())
    {
        if (firstTok.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorIn
            (
                "operator>>(Istream&, PackedBoolList&)",
                is
            )
                << "incorrect first token, expected '(', found "
                << firstTok.info()
                << exit(FatalIOError);
        }

        char c = '(';

        for (label storeI=0; c && c != ')'; ++storeI)
        {
            lst.resize((storeI+1)*PackedList<1>::packing());

            PackedList<1>::StorageType& stored = lst.storage()[storeI];

            // byte-wise read
            for
            (
                unsigned byte=0;
                byte < sizeof(PackedList<1>::StorageType);
                ++byte
            )
            {
                PackedList<1>::StorageType result = 0;
                c = 0;

                // Get next non-whitespace character
                while (is.get(c) && isspace(c))
                {}

                if (!c || c == ')')
                {
                    break;
                }

                for (label nibble=0; nibble < 2; ++nibble)
                {
                    if (!isxdigit(c))
                    {
                        FatalIOErrorIn
                        (
                            "operator>>(Istream&, PackedBoolList&)",
                            is
                        )
                            << "Illegal hex digit: '" << c << "'"
                            << exit(FatalIOError);
                    }

                    result *= 16;

                    if (isdigit(c))
                    {
                        result += int(c) - zeroOffset;
                    }
                    else
                    {
                        result += toupper(c) - alphaOffset;
                    }

                    // Get character for the lower part of the byte
                    if (!nibble)
                    {
                        is.get(c);
                    }
                }

                stored |= result << (byte*CHAR_BIT);
            }
        }

        // put back for later readEndList() to deal with
        is.putback(c);

        // Read end of contents
        is.readEndList("PackedList");
    }
    else
    {
        FatalIOErrorIn
        (
            "operator>>(Istream&, PackedBoolList&)",
            is
        )
            << "incorrect first token, expected <int> or '(', found "
            << firstTok.info()
            << exit(FatalIOError);
    }

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PackedBoolList& lst
)
{
    const label sz = lst.size();

    unsigned int nBytes = packedBytes(sz);

    os  << sz << token::BEGIN_LIST;

    for (label storeI=0; nBytes; ++storeI)
    {
        PackedList<1>::StorageType stored = lst.storage()[storeI];

        for
        (
            unsigned byte=0;
            byte < sizeof(PackedList<1>::StorageType);
            ++byte
        )
        {
            os  << hexLookup[((stored & 0xFF) << 1)]
                << hexLookup[((stored & 0xFF) << 1) + 1];

            if (!--nBytes)
            {
                break;
            }
            stored >>= 8;
        }
    }

    os  << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * *  Global Operators * * * * * * * * * * * * * * //

Foam::PackedBoolList Foam::operator|
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
)
{
    PackedBoolList result(lst1);
    result |= lst2;
    return result;
}


Foam::PackedBoolList Foam::operator&
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
)
{
    PackedBoolList result(lst1);
    result &= lst2;
    return result;
}


Foam::PackedBoolList Foam::operator^
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
)
{
    PackedBoolList result(lst1);
    result ^= lst2;
    return result;
}


// ************************************************************************* //
