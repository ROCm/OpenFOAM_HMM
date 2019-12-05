/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "ILList.H"
#include "Istream.H"
#include "INew.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LListBase, class T>
template<class INew>
void Foam::ILList<LListBase, T>::readIstream(Istream& is, const INew& inew)
{
    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck
    (
        "ILList::readIstream : "
        "reading first token"
    );

    if (firstToken.isLabel())
    {
        const label len = firstToken.labelToken();

        // Read beginning of contents
        const char delimiter = is.readBeginList("ILList");

        if (len)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<len; ++i)
                {
                    T* p = inew(is).ptr();
                    this->append(p);

                    is.fatalCheck
                    (
                        "ILList::readIstream : "
                        "reading entry"
                    );
                }
            }
            else
            {
                T* p = inew(is).ptr();
                this->append(p);

                is.fatalCheck
                (
                    "ILList::readIstream : "
                    "reading the single entry"
                );

                for (label i=1; i<len; ++i)
                {
                    this->append(new T(*p));  // Copy construct
                }
            }
        }

        // Read end of contents
        is.readEndList("ILList");
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

            T* p = inew(is).ptr();
            this->append(p);

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
}


template<class LListBase, class T>
template<class INew>
Foam::ILList<LListBase, T>::ILList(Istream& is, const INew& inew)
{
    this->readIstream(is, inew);
}


template<class LListBase, class T>
Foam::ILList<LListBase, T>::ILList(Istream& is)
{
    this->readIstream(is, INew<T>());
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Istream& Foam::operator>>(Istream& is, ILList<LListBase, T>& list)
{
    list.clear();
    list.readIstream(is, INew<T>());

    return is;
}


// ************************************************************************* //
