/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "keyType.H"
#include "regExp.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::keyType Foam::keyType::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::keyType::keyType(Istream& is)
:
    word(),
    type_(option::LITERAL)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::keyType::match(const std::string& text, bool literal) const
{
    if (!literal && isPattern())
    {
        return regExp(*this).match(text);  // Match as regex
    }

    return !compare(text);  // Compare as literal
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, keyType& val)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get a word/regex"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    if (t.isWord())
    {
        val = t.wordToken();
        val.setType(keyType::LITERAL);
    }
    else if (t.isString())
    {
        // Assign from string, treat as regular expression
        val = t.stringToken();
        val.setType(keyType::REGEX);

        // Flag empty strings as an error
        if (val.empty())
        {
            FatalIOErrorInFunction(is)
                << "Empty word/expression"
                << exit(FatalIOError);
            is.setBad();
            return is;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected word or string, found "
            << t.info()
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const keyType& val)
{
    os.writeQuoted(val, val.isPattern());
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
