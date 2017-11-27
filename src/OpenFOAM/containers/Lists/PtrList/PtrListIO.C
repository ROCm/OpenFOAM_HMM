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

#include "PtrList.H"
#include "SLList.H"
#include "Istream.H"
#include "Ostream.H"
#include "INew.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T>
template<class INew>
void Foam::PtrList<T>::read(Istream& is, const INew& inew)
{
    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck
    (
        "PtrList::read(Istream&) : "
        "reading first token"
    );


    // Label: could be int(..), int{...} or just a plain '0'
    if (firstToken.isLabel())
    {
        // Read size of list
        const label len = firstToken.labelToken();

        // Set list length to that read
        setSize(len);

        // Read beginning of contents
        const char delimiter = is.readBeginList("PtrList");

        if (len)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<len; ++i)
                {
                    T* p = inew(is).ptr();
                    set(i, p);

                    is.fatalCheck
                    (
                        "PtrList::read(Istream&) : "
                        "reading entry"
                    );
                }
            }
            else
            {
                T* p = inew(is).ptr();
                set(0, p);

                is.fatalCheck
                (
                    "PtrList::read(Istream&) : "
                    "reading the single entry"
                );

                for (label i=1; i<len; ++i)
                {
                    set(i, p->clone());
                }
            }
        }

        // Read end of contents
        is.readEndList("PtrList");

        return;
    }


    // "(...)" : read as SLList and transfer contents
    if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction
            (
                is
            )   << "incorrect first token, '(', found " << firstToken.info()
                << exit(FatalIOError);
        }

        SLList<T*> sllPtrs;

        token lastToken(is);
        while
        (
           !(
                lastToken.isPunctuation()
             && lastToken.pToken() == token::END_LIST
            )
        )
        {
            is.putBack(lastToken);

            if (is.eof())
            {
                FatalIOErrorInFunction
                (
                    is
                )   << "Premature EOF after reading " << lastToken.info()
                    << exit(FatalIOError);
            }

            sllPtrs.append(inew(is).ptr());
            is >> lastToken;
        }

        setSize(sllPtrs.size());

        // Pointers - can simply copy
        label i = 0;
        for
        (
            typename SLList<T*>::const_iterator iter = sllPtrs.cbegin();
            iter != sllPtrs.cend();
            ++iter
        )
        {
            set(i++, *iter);
        }

        return;
    }

    FatalIOErrorInFunction
    (
        is
    )   << "incorrect first token, expected <int> or '(', found "
        << firstToken.info()
        << exit(FatalIOError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
template<class INew>
Foam::PtrList<T>::PtrList(Istream& is, const INew& inew)
{
    read(is, inew);
}


template<class T>
Foam::PtrList<T>::PtrList(Istream& is)
{
    read(is, INew<T>());
}


// * * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * //

template<class T>
Foam::Istream& Foam::operator>>(Istream& is, PtrList<T>& lst)
{
    lst.clear();
    lst.read(is, INew<T>());

    return is;
}


// ************************************************************************* //
