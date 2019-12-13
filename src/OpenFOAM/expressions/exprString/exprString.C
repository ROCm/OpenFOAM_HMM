/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
Foam::expressions::exprString::getExpression
(
    const word& name,
    const dictionary& dict,
    const bool stripComments
)
{
    string orig(dict.get<string>(name));

    // No validation
    expressions::exprString expr;
    expr.assign(std::move(orig));

    inplaceExpand(expr, dict, stripComments);

    return expr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::expressions::exprString&
Foam::expressions::exprString::expand
(
    const dictionary& dict,
    const bool stripComments
)
{
    inplaceExpand(*this, dict, stripComments);

    #ifdef FULLDEBUG
    (void)valid();
    #endif

    return *this;
}


// ************************************************************************* //
