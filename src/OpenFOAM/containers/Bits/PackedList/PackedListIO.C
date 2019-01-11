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

#include "PackedList.H"
#include "IOstreams.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<unsigned Width>
void Foam::PackedList<Width>::writeEntry(Ostream& os) const
{
    os  << *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned Width>
Foam::Ostream& Foam::PackedList<Width>::printBits
(
    Ostream& os,
    bool debugOutput
) const
{
    os << token::BEGIN_LIST << nl;

    const label nblocks = debugOutput ? blocks_.size() : num_blocks(size());
    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        BitOps::print(os, blocks_[blocki], '.') << nl;
    }

    os << token::END_LIST << nl;

    return os;
}


template<unsigned Width>
Foam::Istream& Foam::PackedList<Width>::read(Istream& is)
{
    PackedList<Width>& list = *this;

    list.clear();
    is.fatalCheck(FUNCTION_NAME);

    token firstTok(is);
    is.fatalCheck
    (
        "PackedList::read(Istream&) : "
        "reading first token"
    );

    if (firstTok.isLabel())
    {
        const label len = firstTok.labelToken();

        // Set list length to that read
        list.resize(len);

        // Read list contents depending on data format
        if (is.format() == IOstream::ASCII)
        {
            // Read beginning of contents
            const char delimiter = is.readBeginList("PackedList");

            if (len)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (label i=0; i<len; ++i)
                    {
                        list[i] = list.readValue(is);

                        is.fatalCheck
                        (
                            "PackedList::read(Istream&) : "
                            "reading entry"
                        );
                    }
                }
                else if (delimiter == token::BEGIN_BLOCK)
                {
                    // Assign for all entries
                    list = list.readValue(is);

                    is.fatalCheck
                    (
                        "PackedList::read(Istream&) : "
                        "reading the single entry"
                    );
                }
                else
                {
                    FatalIOErrorInFunction(is)
                        << "incorrect list token, expected '(' or '{', found "
                        << firstTok.info()
                        << exit(FatalIOError);
                }
            }

            // Read end of contents
            is.readEndList("PackedList");
        }
        else
        {
            if (len)
            {
                is.read
                (
                    reinterpret_cast<char*>(list.storage().data()),
                    list.byteSize()
                );

                is.fatalCheck
                (
                    "PackedList::read(Istream&) : "
                    "reading the binary block"
                );
            }
        }
    }
    else if (firstTok.isPunctuation())
    {
        if (firstTok.pToken() == token::BEGIN_LIST)
        {
            token nextTok(is);
            is.fatalCheck(FUNCTION_NAME);

            while
            (
                !(   nextTok.isPunctuation()
                  && nextTok.pToken() == token::END_LIST
                 )
            )
            {
                is.putBack(nextTok);
                list.append(list.readValue(is));

                is  >> nextTok;
                is.fatalCheck(FUNCTION_NAME);
            }
        }
        else if (firstTok.pToken() == token::BEGIN_BLOCK)
        {
            token nextTok(is);
            is.fatalCheck(FUNCTION_NAME);

            while
            (
                !(   nextTok.isPunctuation()
                  && nextTok.pToken() == token::END_BLOCK
                 )
            )
            {
                is.putBack(nextTok);
                list.setPair(is);

                is  >> nextTok;
                is.fatalCheck(FUNCTION_NAME);
            }
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected '(', found "
                << firstTok.info()
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int>, '(' or '{', found "
            << firstTok.info()
            << exit(FatalIOError);
    }

    return is;
}


template<unsigned Width>
Foam::Ostream& Foam::PackedList<Width>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const PackedList<Width>& list = *this;
    const label len = list.size();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII)
    {
        if (len > 1 && list.uniform())
        {
            // Two or more entries, and all have identical values.
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


template<unsigned Width>
void Foam::PackedList<Width>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.writeKeyword(keyword);
    writeEntry(os);
    os  << token::END_STATEMENT << endl;
}


// * * * * * * * * * * * * * *  Friend Operators * * * * * * * * * * * * * * //

template<unsigned Width>
Foam::Istream& Foam::operator>>(Istream& is, PackedList<Width>& list)
{
    return list.read(is);
}


template<unsigned Width>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<PackedList<Width>>& iproxy
)
{
    const PackedList<Width>& list = iproxy.t_;

    os  << "PackedList<" << Width
        << "> size=" << list.size() << "/" << list.capacity()
        << " (limits: max=" << PackedList<Width>::max_value
        << ", elem_per_block=" << PackedList<Width>::elem_per_block
        << ")"
        << nl;

    return os;
}


// ************************************************************************* //
