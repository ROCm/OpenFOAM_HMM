/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
#include "contiguous.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T, unsigned Size>
void Foam::FixedList<T, Size>::writeEntry(Ostream& os) const
{
    const word tag = "List<" + word(pTraits<T>::typeName) + '>';
    if (token::compound::isCompound(tag))
    {
        os  << tag << ' ';
    }

    os << *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T, unsigned Size>
void Foam::FixedList<T, Size>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.writeKeyword(keyword);
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}


template<class T, unsigned Size>
Foam::Ostream& Foam::FixedList<T, Size>::writeList
(
    Ostream& os,
    const label shortListLen
) const
{
    const FixedList<T, Size>& list = *this;

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<T>())
    {
        if (contiguous<T>() && list.uniform())
        {
            // Two or more entries, and all entries have identical values.

            os << Size << token::BEGIN_BLOCK << list[0] << token::END_BLOCK;
        }
        else if
        (
            Size <= 1 || !shortListLen
         || (Size <= unsigned(shortListLen) && contiguous<T>())
        )
        {
            // Start delimiter
            os << token::BEGIN_LIST;

            // Contents
            for (unsigned i=0; i<Size; ++i)
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
            for (unsigned i=0; i<Size; ++i)
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
        os.write(reinterpret_cast<const char*>(list.cdata()), Size*sizeof(T));
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, unsigned Size>
Foam::FixedList<T, Size>::FixedList(Istream& is)
{
    operator>>(is, *this);
}


template<class T, unsigned Size>
Foam::Istream& Foam::operator>>(Foam::Istream& is, FixedList<T, Size>& list)
{
    is.fatalCheck(FUNCTION_NAME);

    if (is.format() == IOstream::ASCII || !contiguous<T>())
    {
        token firstToken(is);

        is.fatalCheck
        (
            "operator>>(Istream&, FixedList<T, Size>&) : reading first token"
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
            for (unsigned i=0; i<Size; ++i)
            {
                is >> list[i];

                is.fatalCheck
                (
                    "operator>>(Istream&, FixedList<T, Size>&) : "
                    "reading entry"
                );
            }
        }
        else
        {
            T element;
            is >> element;

            is.fatalCheck
            (
                "operator>>(Istream&, FixedList<T, Size>&) : "
                "reading the single entry"
            );

            for (unsigned i=0; i<Size; ++i)
            {
                list[i] = element;  // Copy the value
            }
        }

        // Read end of contents
        is.readEndList("FixedList");
    }
    else
    {
        // Binary and contiguous

        is.read(reinterpret_cast<char*>(list.data()), Size*sizeof(T));

        is.fatalCheck
        (
            "operator>>(Istream&, FixedList<T, Size>&) : "
            "reading the binary block"
        );
    }

    return is;
}


template<class T, unsigned Size>
Foam::Ostream& Foam::operator<<(Ostream& os, const FixedList<T, Size>& list)
{
    return list.writeList(os, 10);
}


// ************************************************************************* //
