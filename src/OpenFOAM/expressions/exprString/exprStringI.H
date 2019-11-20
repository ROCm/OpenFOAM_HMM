/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Original code Copyright (C) 2012-2018 Bernhard Gschaider
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

#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::expressions::exprString::exprString
(
    const std::string& s,
    bool doValidate
)
:
    string(s)
{
    #ifdef FULLDEBUG
    if (doValidate)
    {
        (void)valid();
    }
    #endif
}


inline Foam::expressions::exprString::exprString
(
    std::string&& s,
    bool doValidate
)
:
    string(std::move(s))
{
    #ifdef FULLDEBUG
    if (doValidate)
    {
        (void)valid();
    }
    #endif
}


inline Foam::expressions::exprString::exprString
(
    const char* s,
    bool doValidate
)
:
    string(s)
{
    #ifdef FULLDEBUG
    if (doValidate)
    {
        (void)valid();
    }
    #endif
}


inline Foam::expressions::exprString::exprString
(
    const std::string& str,
    const dictionary& dict,
    const bool removeComments
)
:
    string(str)
{
    expand(dict, removeComments);
}


inline Foam::expressions::exprString::exprString
(
    std::string&& str,
    const dictionary& dict,
    const bool removeComments
)
:
    string(std::move(str))
{
    expand(dict, removeComments);
}


inline Foam::expressions::exprString::exprString
(
    Istream& is,
    const dictionary& dict,
    const bool removeComments
)
:
    string(is)
{
    expand(dict, removeComments);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::expressions::exprString::valid() const
{
    const bool ok = (std::string::npos == find('$'));

    #ifdef FULLDEBUG
    if (!ok)
    {
        FatalErrorInFunction
            << "Unexpanded '$' in " << *this << nl
            << exit(FatalError);
    }
    #endif

    return ok;
}


inline Foam::expressions::exprString
Foam::expressions::exprString::toExpr(const std::string& str)
{
    exprString expr;

    expr.string::operator=(str);

    return expr;
}


inline Foam::expressions::exprString
Foam::expressions::exprString::toExpr(std::string&& str)
{
    exprString expr;

    expr.string::operator=(std::move(str));

    return expr;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::expressions::exprString&
Foam::expressions::exprString::operator=(const char* str)
{
    string::operator=(str);

    #ifdef FULLDEBUG
    (void)valid();
    #endif

    return *this;
}


inline Foam::expressions::exprString&
Foam::expressions::exprString::operator=(const std::string& str)
{
    string::operator=(str);

    #ifdef FULLDEBUG
    (void)valid();
    #endif

    return *this;
}


inline Foam::expressions::exprString&
Foam::expressions::exprString::operator=(std::string&& str)
{
    string::operator=(std::move(str));

    #ifdef FULLDEBUG
    (void)valid();
    #endif

    return *this;
}


// ************************************************************************* //
