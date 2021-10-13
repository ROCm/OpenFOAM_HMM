/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
#include "token.H"
#include "SLList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::List(Istream& is)
:
    UList<T>(nullptr, 0)
{
    this->readList(is);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::List<T>::readList(Istream& is)
{
    List<T>& list = *this;

    // Anull list
    list.clear();

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck("List<T>::readList(Istream&) : reading first token");

    if (tok.isCompound())
    {
        // Compound: simply transfer contents

        list.transfer
        (
            dynamicCast<token::Compound<List<T>>>
            (
                tok.transferCompoundToken(is)
            )
        );
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..), int{...} or just a plain '0'

        const label len = tok.labelToken();

        // Resize to actual length read
        list.resize(len);

        if (is.format() == IOstream::BINARY && is_contiguous<T>::value)
        {
            // Binary and contiguous

            if (len)
            {
                Detail::readContiguous<T>
                (
                    is,
                    list.data_bytes(),
                    list.size_bytes()
                );

                is.fatalCheck
                (
                    "List<T>::readList(Istream&) : "
                    "reading the binary block"
                );
            }
        }
        else
        {
            // Begin of contents marker
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
                            "List<T>::readList(Istream&) : "
                            "reading entry"
                        );
                    }
                }
                else
                {
                    // Uniform content (delimiter == token::BEGIN_BLOCK)

                    T elem;
                    is >> elem;

                    is.fatalCheck
                    (
                        "List<T>::readList(Istream&) : "
                        "reading the single entry"
                    );

                    for (label i=0; i<len; ++i)
                    {
                        list[i] = elem;  // Copy the value
                    }
                }
            }

            // End of contents marker
            is.readEndList("List");
        }
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        // "(...)" : read as SLList and transfer contents

        is.putBack(tok);    // Putback the opening bracket
        SLList<T> sll(is);  // Read as singly-linked list

        // Reallocate and move assign list elements
        list = std::move(sll);
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


// ************************************************************************* //
