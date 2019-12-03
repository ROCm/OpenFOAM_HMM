/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "FixedList.H"
#include "Istream.H"
#include "Ostream.H"
#include "token.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T, unsigned N>
void Foam::FixedList<T, N>::writeEntry(Ostream& os) const
{
    const word tag = "List<" + word(pTraits<T>::typeName) + '>';
    if (token::compound::isCompound(tag))
    {
        os  << tag << token::SPACE;
    }
    os << *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T, unsigned N>
void Foam::FixedList<T, N>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    if (keyword.size())
    {
        os.writeKeyword(keyword);
    }
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}


template<class T, unsigned N>
Foam::Ostream& Foam::FixedList<T, N>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const FixedList<T, N>& list = *this;

    // Write list contents depending on data format

    // Unlike UList, no compact output form since a FixedList is generally
    // small and we desire a consistent appearance.
    // Eg, FixedList<T,2> or Pair<T> as "(-1 -1)", not as "2{-1}"

    if (os.format() == IOstream::ASCII || !is_contiguous<T>::value)
    {
        if
        (
            (N <= 1 || !shortLen)
         ||
            (
                (N <= unsigned(shortLen))
             &&
                (
                    Detail::ListPolicy::no_linebreak<T>::value
                 || is_contiguous<T>::value
                )
            )
        )
        {
            // Start delimiter
            os << token::BEGIN_LIST;

            // Contents
            for (unsigned i=0; i<N; ++i)
            {
                if (i) os << token::SPACE;
                os << list[i];
            }

            // End delimiter
            os << token::END_LIST;
        }
        else
        {
            // Start delimiter
            os << nl << token::BEGIN_LIST << nl;

            // Contents
            for (unsigned i=0; i<N; ++i)
            {
                os << list[i] << nl;
            }

            // End delimiter
            os << token::END_LIST << nl;
        }
    }
    else
    {
        // Binary, contiguous

        // write(...) includes surrounding start/end delimiters
        os.write(reinterpret_cast<const char*>(list.cdata()), N*sizeof(T));
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, unsigned N>
Foam::FixedList<T, N>::FixedList(Istream& is)
{
    operator>>(is, *this);
}


template<class T, unsigned N>
Foam::Istream& Foam::operator>>(Foam::Istream& is, FixedList<T, N>& list)
{
    is.fatalCheck(FUNCTION_NAME);

    if (is.format() == IOstream::ASCII || !is_contiguous<T>::value)
    {
        token firstToken(is);

        is.fatalCheck
        (
            "operator>>(Istream&, FixedList<T, N>&) : "
            "reading first token"
        );

        if (firstToken.isCompound())
        {
            list = dynamicCast<token::Compound<List<T>>>
            (
                firstToken.transferCompoundToken(is)
            );
        }
        else if (firstToken.isLabel())
        {
            const label len = firstToken.labelToken();

            // List lengths must match
            list.checkSize(len);
        }
        else if (!firstToken.isPunctuation())
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected <label> "
                   "or '(' or '{', found "
                << firstToken.info()
                << exit(FatalIOError);
        }
        else
        {
            // Putback the opening bracket
            is.putBack(firstToken);
        }

        // Read beginning of contents
        const char delimiter = is.readBeginList("FixedList");

        if (delimiter == token::BEGIN_LIST)
        {
            for (unsigned i=0; i<N; ++i)
            {
                is >> list[i];

                is.fatalCheck
                (
                    "operator>>(Istream&, FixedList<T, N>&) : "
                    "reading entry"
                );
            }
        }
        else
        {
            T val;
            is >> val;

            is.fatalCheck
            (
                "operator>>(Istream&, FixedList<T, N>&) : "
                "reading the single entry"
            );

            for (unsigned i=0; i<N; ++i)
            {
                list[i] = val;  // Copy the value
            }
        }

        // Read end of contents
        is.readEndList("FixedList");
    }
    else
    {
        // Binary and contiguous

        Detail::readContiguous<T>
        (
            is,
            reinterpret_cast<char*>(list.data()),
            N*sizeof(T)
        );

        is.fatalCheck
        (
            "operator>>(Istream&, FixedList<T, N>&) : "
            "reading the binary block"
        );
    }

    return is;
}


// ************************************************************************* //
