/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "List.H"
#include "Istream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::CircularBuffer<T>::CircularBuffer(Istream& is)
{
    this->readList(is);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::CircularBuffer<T>::info(Ostream& os) const
{
    os  << "size=" << size() << '/' << capacity()
        << " begin=" << begin_
        << " end=" << end_
        /// << " one=" << this->range_one() << this->array_one()
        /// << " two=" << this->range_two() << this->array_two()
        << nl;

    return os;
}


template<class T>
Foam::Istream& Foam::CircularBuffer<T>::readList(Istream& is)
{
    // Clear list
    storage_.clear();
    begin_ = 0;
    end_ = 0;

    // More work than it should be. We avoid completely filled buffers!

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck
    (
        "CircularBuffer<T>::readList(Istream&) : "
        "reading first token"
    );

    if (tok.isCompound())
    {
        // Compound: simply transfer contents

        storage_.transfer
        (
            dynamicCast<token::Compound<List<T>>>
            (
                tok.transferCompoundToken(is)
            )
        );

        end_ = storage_.size();
        if (end_)
        {
            // Resize larger to avoid full buffer
            storage_.resize(end_ + min_size());
        }
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..), int{...} or just a plain '0'

        const label len = tok.labelToken();

        end_ = len;
        if (end_)
        {
            // Resize larger to avoid full buffer
            storage_.resize(end_ + min_size());
        }

        // Dispatch to UList reading...

        UList<T> list(storage_.data(), end_);

        is.putBack(tok);
        list.readList(is);
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        // "(...)" : read as SLList and transfer contents

        is.putBack(tok);    // Putback the opening bracket
        SLList<T> sll(is);  // Read as singly-linked list

        const label len = sll.size();

        end_ = len;
        if (end_)
        {
            // Resize larger to avoid full buffer
            storage_.resize(end_ + min_size());

            // Move assign each list element
            for (label i = 0; i < len; ++i)
            {
                storage_[i] = std::move(sll.removeHead());
            }
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }

    return is;
}


template<class T>
Foam::Ostream& Foam::CircularBuffer<T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const label len = this->size();
    const auto list1 = this->array_one();
    const auto list2 = this->array_two();

    #ifdef FULLDEBUG
    if (len != (list1.size() + list2.size()))
    {
        FatalErrorInFunction
            << "Size check failed"
            << abort(FatalError);
    }
    #endif

    if (os.format() == IOstream::BINARY && is_contiguous<T>::value)
    {
        // Binary and contiguous

        os << nl << len << nl;

        if (len)
        {
            // The TOTAL number of bytes to be written.
            // - possibly add start delimiter
            // This is much like IndirectListBase output

            os.beginRawWrite(len*sizeof(T));

            if (!list1.empty())
            {
                os.writeRaw(list1.cdata_bytes(), list1.size_bytes());
            }
            if (!list2.empty())
            {
                os.writeRaw(list2.cdata_bytes(), list2.size_bytes());
            }

            // End delimiter and/or cleanup.
            os.endRawWrite();
        }
    }
    else if
    (
        (len <= 1 || !shortLen)
     ||
        (
            (len <= shortLen)
         &&
            (
                is_contiguous<T>::value
             || Detail::ListPolicy::no_linebreak<T>::value
            )
        )
    )
    {
        // Single-line output

        // Size and start delimiter
        os << len << token::BEGIN_LIST;

        // Contents
        label i = 0;
        for (const T& val : list1)
        {
            if (i++) os << token::SPACE;
            os << val;
        }
        for (const T& val : list2)
        {
            if (i++) os << token::SPACE;
            os << val;
        }

        // End delimiter
        os << token::END_LIST;
    }
    else
    {
        // Multi-line output

        // Size and start delimiter
        os << nl << len << nl << token::BEGIN_LIST << nl;

        // Contents
        for (const T& val : list1)
        {
            os << val << nl;
        }
        for (const T& val : list2)
        {
            os << val << nl;
        }

        // End delimiter
        os << token::END_LIST << nl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
