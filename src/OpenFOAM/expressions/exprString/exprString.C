/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "exprString.H"
#include "stringOps.H"
#include "expressionEntry.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::expressions::exprString::inplaceExpand
(
    std::string& str,
    const dictionary& dict,
    const bool stripComments
)
{
    if (stripComments)
    {
        stringOps::inplaceRemoveComments(str);
    }

    exprTools::expressionEntry::inplaceExpand(str, dict);
}


Foam::expressions::exprString
Foam::expressions::exprString::getEntry
(
    const word& key,
    const dictionary& dict,
    const bool stripComments
)
{
    exprString expr;
    expr.readEntry(key, dict, true, stripComments);  // mandatory

    return expr;
}


Foam::expressions::exprString
Foam::expressions::exprString::getOptional
(
    const word& key,
    const dictionary& dict,
    const bool stripComments
)
{
    exprString expr;
    expr.readEntry(key, dict, false, stripComments);  // optional

    return expr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::exprString::expand
(
    const dictionary& dict,
    const bool stripComments
)
{
    inplaceExpand(*this, dict, stripComments);

    #ifdef FULLDEBUG
    (void)valid();
    #endif
}


void Foam::expressions::exprString::trim()
{
    stringOps::inplaceTrim(*this);
}


bool Foam::expressions::exprString::readEntry
(
    const word& keyword,
    const dictionary& dict,
    bool mandatory,
    const bool stripComments
)
{
    const bool ok = dict.readEntry(keyword, *this, keyType::LITERAL, mandatory);

    if (ok && !empty())
    {
        this->expand(dict, stripComments);  // strip comments
    }
    else
    {
        clear();
    }

    return ok;
}


bool Foam::expressions::exprString::writeEntry
(
    const word& keyword,
    Ostream& os,
    bool writeEmpty
) const
{
    const bool ok = (writeEmpty || !empty());

    if (ok)
    {
        if (!keyword.empty())
        {
            os.writeKeyword(keyword);
        }

        // Write as regular or verbatim string

        token tok(*this);
        if (!empty())
        {
            tok.setType(token::tokenType::VERBATIM);
        }

        os.write(tok);
        os.endEntry();
    }

    return ok;
}


// ************************************************************************* //
