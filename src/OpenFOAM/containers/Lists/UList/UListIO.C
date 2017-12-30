/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
            os  << tag << ' ';
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
    os.writeKeyword(keyword);
    writeEntry(os);
    os << token::END_STATEMENT << endl;
}


template<class T>
Foam::Ostream& Foam::UList<T>::writeList
(
    Ostream& os,
    const label shortListLen
) const
{
    const UList<T>& L = *this;

    const label sz = L.size();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<T>())
    {
        // Can the contents be considered 'uniform' (ie, identical)?
        bool uniform = (sz > 1 && contiguous<T>());
        if (uniform)
        {
            for (label i=1; i < sz; ++i)
            {
                if (L[i] != L[0])
                {
                    uniform = false;
                    break;
                }
            }
        }

        if (uniform)
        {
            // Size and start delimiter
            os << sz << token::BEGIN_BLOCK;

            // Contents
            os << L[0];

            // End delimiter
            os << token::END_BLOCK;
        }
        else if
        (
            sz <= 1 || !shortListLen
         || (sz <= shortListLen && contiguous<T>())
        )
        {
            // Size and start delimiter
            os << sz << token::BEGIN_LIST;

            // Contents
            for (label i=0; i < sz; ++i)
            {
                if (i) os << token::SPACE;
                os << L[i];
            }

            // End delimiter
            os << token::END_LIST;
        }
        else
        {
            // Size and start delimiter
            os << nl << sz << nl << token::BEGIN_LIST << nl;

            // Contents
            for (label i=0; i < sz; ++i)
            {
                os << L[i] << nl;
            }

            // End delimiter
            os << token::END_LIST << nl;
        }
    }
    else
    {
        // Contents are binary and contiguous
        os << nl << sz << nl;

        if (sz)
        {
            // write(...) includes surrounding start/end delimiters
            os.write
            (
                reinterpret_cast<const char*>(L.cdata()),
                L.byteSize()
            );
        }
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const UList<T>& lst)
{
    return lst.writeList(os, 10);
}


template<class T>
Foam::Istream& Foam::operator>>(Istream& is, UList<T>& L)
{
    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck("operator>>(Istream&, UList<T>&) : reading first token");

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
        // Check list length
        const label sz = elems.size();

        if (sz != L.size())
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read " << sz
                << " expected " << L.size()
                << exit(FatalIOError);
        }
        for (label i=0; i<sz; ++i)
        {
            L[i] = elems[i];
        }
    }
    else if (firstToken.isLabel())
    {
        const label sz = firstToken.labelToken();

        // Set list length to that read
        if (sz != L.size())
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read " << sz
                << " expected " << L.size()
                << exit(FatalIOError);
        }

        // Read list contents depending on data format

        if (is.format() == IOstream::ASCII || !contiguous<T>())
        {
            // Read beginning of contents
            const char delimiter = is.readBeginList("List");

            if (sz)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (label i=0; i<sz; ++i)
                    {
                        is >> L[i];

                        is.fatalCheck
                        (
                            "operator>>(Istream&, UList<T>&) : reading entry"
                        );
                    }
                }
                else
                {
                    // uniform content (delimiter == token::BEGIN_BLOCK)

                    T element;
                    is >> element;

                    is.fatalCheck
                    (
                        "operator>>(Istream&, UList<T>&) : "
                        "reading the single entry"
                    );

                    for (label i=0; i<sz; ++i)
                    {
                        L[i] = element;
                    }
                }
            }

            // Read end of contents
            is.readEndList("List");
        }
        else
        {
            // contents are binary and contiguous

            if (sz)
            {
                is.read(reinterpret_cast<char*>(L.data()), sz*sizeof(T));

                is.fatalCheck
                (
                    "operator>>(Istream&, UList<T>&) : reading the binary block"
                );
            }
        }
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected '(', found "
                << firstToken.info()
                << exit(FatalIOError);
        }

        // Putback the opening bracket
        is.putBack(firstToken);

        // Now read as a singly-linked list
        SLList<T> sll(is);

        if (sll.size() != L.size())
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read " << sll.size()
                << " expected " << L.size()
                << exit(FatalIOError);
        }

        // Convert the singly-linked list to this list
        label i = 0;
        for
        (
            typename SLList<T>::const_iterator iter = sll.begin();
            iter != sll.end();
            ++iter, ++i
        )
        {
            L[i] = iter();
        }

    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    return is;
}


// ************************************************************************* //
