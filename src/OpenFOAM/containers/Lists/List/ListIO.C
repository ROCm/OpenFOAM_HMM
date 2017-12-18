/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "List.H"
#include "Istream.H"
#include "token.H"
#include "SLList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::List<T>::List(Istream& is)
:
    UList<T>(nullptr, 0)
{
    operator>>(is, *this);
}


template<class T>
Foam::Istream& Foam::operator>>(Istream& is, List<T>& L)
{
    // Anull list
    L.setSize(0);

    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck(FUNCTION_NAME);

    if (firstToken.isCompound())
    {
        L.transfer
        (
            dynamicCast<token::Compound<List<T>>>
            (
                firstToken.transferCompoundToken(is)
            )
        );
    }
    else if (firstToken.isLabel())
    {
        const label sz = firstToken.labelToken();

        // Set list length to that read
        L.setSize(sz);

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
                            "operator>>(Istream&, List<T>&) : reading entry"
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
                        "operator>>(Istream&, List<T>&) : "
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
            // Contents are binary and contiguous

            if (sz)
            {
                is.read(reinterpret_cast<char*>(L.data()), sz*sizeof(T));

                is.fatalCheck
                (
                    "operator>>(Istream&, List<T>&) : reading the binary block"
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

        // Read as a singly-linked list
        SLList<T> sll(is);

        // Convert the singly-linked list to this list
        L = sll;
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
