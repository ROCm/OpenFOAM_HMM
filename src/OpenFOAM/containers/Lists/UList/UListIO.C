/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "UList.H"
#include "Ostream.H"
#include "token.H"
#include "SLList.H"
#include "contiguous.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::writeEntry(Ostream& os) const
{
    if (size())
    {
        const word tag = "List<" + word(pTraits<T>::typeName) + '>';
        if (token::compound::isCompound(tag))
        {
            os  << tag << token::SPACE;
        }
        os << *this;
    }
    else if (os.format() == IOstream::ASCII)
    {
        // Zero-sized ASCII - Write size and delimiters
        os  << 0 << token::BEGIN_LIST << token::END_LIST;
    }
    else
    {
        // Zero-sized binary - Write size only
        os  << 0;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::writeEntry(const word& keyword, Ostream& os) const
{
    if (keyword.size())
    {
        os.writeKeyword(keyword);
    }
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}


template<class T>
Foam::Ostream& Foam::UList<T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const UList<T>& list = *this;

    const label len = list.size();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !is_contiguous<T>::value)
    {
        if (len > 1 && is_contiguous<T>::value && list.uniform())
        {
            // Two or more entries, and all entries have identical values.
            os  << len << token::BEGIN_BLOCK << list[0] << token::END_BLOCK;
        }
        else if
        (
            (len <= 1 || !shortLen)
         ||
            (
                (len <= shortLen)
             &&
                (
                    Detail::ListPolicy::no_linebreak<T>::value
                 || is_contiguous<T>::value
                )
            )
        )
        {
            // Size and start delimiter
            os << len << token::BEGIN_LIST;

            // Contents
            for (label i=0; i < len; ++i)
            {
                if (i) os << token::SPACE;
                os << list[i];
            }

            // End delimiter
            os << token::END_LIST;
        }
        else
        {
            // Size and start delimiter
            os << nl << len << nl << token::BEGIN_LIST << nl;

            // Contents
            for (label i=0; i < len; ++i)
            {
                os << list[i] << nl;
            }

            // End delimiter
            os << token::END_LIST << nl;
        }
    }
    else
    {
        // Contents are binary and contiguous
        os << nl << len << nl;

        if (len)
        {
            // write(...) includes surrounding start/end delimiters
            os.write
            (
                reinterpret_cast<const char*>(list.cdata()),
                list.byteSize()
            );
        }
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::operator>>(Istream& is, UList<T>& list)
{
    // The target list length - must match with sizes read
    const label len = list.size();

    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck("operator>>(Istream&, UList<T>&) : reading first token");

    // Compound: simply transfer contents
    if (firstToken.isCompound())
    {
        List<T> elems;
        elems.transfer
        (
            dynamicCast<token::Compound<List<T>>>
            (
                firstToken.transferCompoundToken(is)
            )
        );

        const label inputLen = elems.size();

        // List lengths must match
        if (inputLen != len)
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read "
                << inputLen << " expected " << len
                << exit(FatalIOError);
        }

        for (label i = 0; i < len; ++i)
        {
            list[i] = std::move(elems[i]);
        }

        return is;
    }


    // Label: could be int(..), int{...} or just a plain '0'
    if (firstToken.isLabel())
    {
        const label inputLen = firstToken.labelToken();

        // List lengths must match
        if (inputLen != len)
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read "
                << inputLen << " expected " << len
                << exit(FatalIOError);
        }

        // Read list contents depending on data format

        if (is.format() == IOstream::ASCII || !is_contiguous<T>::value)
        {
            // Read beginning of contents
            const char delimiter = is.readBeginList("List");

            if (len)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (label i=0; i<len; ++i)
                    {
                        is >> list[i];

                        is.fatalCheck
                        (
                            "operator>>(Istream&, UList<T>&) : reading entry"
                        );
                    }
                }
                else
                {
                    // Uniform content (delimiter == token::BEGIN_BLOCK)

                    T element;
                    is >> element;

                    is.fatalCheck
                    (
                        "operator>>(Istream&, UList<T>&) : "
                        "reading the single entry"
                    );

                    for (label i=0; i<len; ++i)
                    {
                        list[i] = element;  // Copy the value
                    }
                }
            }

            // Read end of contents
            is.readEndList("List");
        }
        else if (len)
        {
            // Non-empty, binary, contiguous

            Detail::readContiguous<T>
            (
                is,
                reinterpret_cast<char*>(list.data()),
                len*sizeof(T)
            );

            is.fatalCheck
            (
                "operator>>(Istream&, UList<T>&) : reading the binary block"
            );
        }

        return is;
    }


    // "(...)" : read as SLList and transfer contents
    if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected '(', found "
                << firstToken.info()
                << exit(FatalIOError);
        }

        is.putBack(firstToken); // Putback the opening bracket

        SLList<T> sll(is);      // Read as singly-linked list

        // List lengths must match
        if (sll.size() != len)
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read "
                << sll.size() << " expected " << len
                << exit(FatalIOError);
        }

        // Move assign each list element
        for (label i = 0; i < len; ++i)
        {
            list[i] = std::move(sll.removeHead());
        }

        return is;
    }


    FatalIOErrorInFunction(is)
        << "incorrect first token, expected <int> or '(', found "
        << firstToken.info()
        << exit(FatalIOError);

    return is;
}


// ************************************************************************* //
