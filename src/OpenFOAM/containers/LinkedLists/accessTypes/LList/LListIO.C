/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "LList.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::LList<LListBase, T>::LList(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Istream& Foam::operator>>(Istream& is, LList<LListBase, T>& lst)
{
    // Anull list
    lst.clear();

    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck("LList::readList : reading first token");

    if (firstToken.isLabel())
    {
        const label len = firstToken.labelToken();

        // Read beginning of contents
        const char delimiter = is.readBeginList("LList");

        if (len)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<len; ++i)
                {
                    T element;
                    is >> element;
                    lst.append(element);
                }
            }
            else
            {
                T element;
                is >> element;

                for (label i=0; i<len; ++i)
                {
                    lst.append(element);
                }
            }
        }

        // Read end of contents
        is.readEndList("LList");
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        token lastToken(is);
        is.fatalCheck(FUNCTION_NAME);

        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);

            T element;
            is >> element;
            lst.append(element);

            is >> lastToken;
            is.fatalCheck(FUNCTION_NAME);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    is.fatalCheck(FUNCTION_NAME);
    return is;
}


template<class LListBase, class T>
Foam::Ostream& Foam::LList<LListBase, T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    const label len = this->size();

    if
    (
        (len <= 1 || !shortLen)
     || (len <= shortLen)
    )
    {
        // Size and start delimiter
        os << len << token::BEGIN_LIST;

        // Contents
        bool space = false;
        for (const T& val : *this)
        {
            if (space) os << token::SPACE;
            os << val;
            space = true;
        }

        // End delimiter
        os << token::END_LIST;
    }
    else
    {
        // Size and start delimiter
        os << nl << len << nl << token::BEGIN_LIST << nl;

        // Contents
        for (const T& val : *this)
        {
            os << val << nl;
        }

        // End delimiter
        os << token::END_LIST;
    }

    os.check(FUNCTION_NAME);
    return os;
}


template<class LListBase, class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const LList<LListBase, T>& lst)
{
    return lst.writeList(os, -1);  // Always with line breaks
}


// ************************************************************************* //
