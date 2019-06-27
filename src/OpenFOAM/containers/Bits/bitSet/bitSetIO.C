/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "bitSet.H"
#include "IOstreams.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::bitSet::writeEntry(Ostream& os) const
{
    os  << *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::bitSet::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const bitSet& list = *this;
    const label len = list.size();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII)
    {
        if (len > 1 && list.uniform())
        {
            // Two or more entries, and all entries have identical values.
            os  << len << token::BEGIN_BLOCK << list[0] << token::END_BLOCK;
        }
        else if (!shortLen || len <= shortLen)
        {
            // Shorter list, or line-breaks suppressed
            os  << len << token::BEGIN_LIST;
            for (label i=0; i < len; ++i)
            {
                if (i) os << token::SPACE;
                os  << list[i];
            }
            os  << token::END_LIST;
        }
        else
        {
            // Longer list
            os << nl << len << nl << token::BEGIN_LIST << nl;
            for (label i=0; i < len; ++i)
            {
                os << list[i] << nl;
            }
            os << token::END_LIST << nl;
        }
    }
    else
    {
        // Contents are binary and contiguous
        os  << nl << len << nl;

        if (len)
        {
            // write(...) includes surrounding start/end delimiters
            os.write
            (
                reinterpret_cast<const char*>(list.storage().cdata()),
                list.byteSize()
            );
        }
    }

    return os;
}


void Foam::bitSet::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.writeKeyword(keyword);
    writeEntry(os);
    os  << token::END_STATEMENT << endl;
}


// * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const bitSet& bitset)
{
    return bitset.writeList(os, 40);
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<bitSet>& iproxy
)
{
    const bitSet& bitset = iproxy.t_;

    os  << "bitSet<" << bitSet::elem_per_block
        << "> size=" << bitset.size() << "/" << bitset.capacity()
        << " count=" << bitset.count()
        << nl;

    return os;
}


// ************************************************************************* //
