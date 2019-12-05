/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
Foam::Istream& Foam::operator>>(Istream& is, List<T>& list)
{
    // Anull list
    list.resize(0);

    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck(FUNCTION_NAME);

    // Compound: simply transfer contents
    if (firstToken.isCompound())
    {
        list.transfer
        (
            dynamicCast<token::Compound<List<T>>>
            (
                firstToken.transferCompoundToken(is)
            )
        );

        return is;
    }


    // Label: could be int(..), int{...} or just a plain '0'
    if (firstToken.isLabel())
    {
        const label len = firstToken.labelToken();

        // Resize to length read
        list.resize(len);

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
                            "operator>>(Istream&, List<T>&) : "
                            "reading entry"
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
                "operator>>(Istream&, List<T>&) : "
                "reading the binary block"
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

        // Reallocate and move assign list elements
        list = std::move(sll);

        return is;
    }


    FatalIOErrorInFunction(is)
        << "incorrect first token, expected <int> or '(', found "
        << firstToken.info()
        << exit(FatalIOError);

    return is;
}


// ************************************************************************* //
