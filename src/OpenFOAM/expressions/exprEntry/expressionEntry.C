/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Original code Copyright (C) 2014-2018 Bernhard Gschaider
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

#include "expressionEntry.H"
#include "stringOps.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace exprTools
{

defineTypeName(expressionEntry);
defineRunTimeSelectionTable(expressionEntry, empty);

// Various types can be used directly without any changes

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    expressionEntry,
    empty,
    direct
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    expressionEntry,
    empty,
    label
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    expressionEntry,
    empty,
    scalar
);

addNamedToRunTimeSelectionTable
(
    expressionEntry,
    expressionEntry,
    empty,
    word
);

} // End namespace exprTools
} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{
// Same code as in stringOps.C

// Acceptable values for $variable names.
//
// Similar to word::valid(), except we don't have the benefit of a parser
// to filter out other unacceptable entries for us.
//
// Does not currently accept '/' in a variable name.
// We would like "$file/$name" to expand as two variables.
static inline bool validVariableChar(char c)
{
    return
    (
        std::isalnum(c)
     || c == '.'
     || c == ':'
     || c == '_'
    );
}


// For input string of "$variable with other" return the length of
// the variable.
//
// Intentionally will not capture ':+', ':-' alterations. Use ${ .. } for that
static inline std::string::size_type findVariableLen
(
    const std::string& s,
    std::string::size_type pos,
    const char sigil = '$'
)
{
    std::string::size_type len = 0;

    if (pos < s.length())
    {
        if (s[pos] == sigil)
        {
            // Skip leading '$' in the count!
            ++pos;
        }

        for
        (
            auto iter = s.cbegin() + pos;
            iter != s.cend() && validVariableChar(*iter);
            ++iter
        )
        {
            ++len;
        }
    }

    return len;
}

} // End anonymous namespace


namespace Foam
{

inline static const entry* getVariableOrDie
(
    const word& name,
    const dictionary& dict
)
{
    const entry* eptr = dict.findScoped(name, keyType::LITERAL_RECURSIVE);

    if (!eptr)
    {
        FatalIOErrorInFunction(dict)
            << "No dictionary entry " << name << nl
            << exit(FatalIOError);
    }

    if (eptr->isDict())
    {
        FatalIOErrorInFunction(dict)
            << "Found dictionary " << name << " instead of entry" << nl
            << exit(FatalIOError);
    }

    return eptr;
}


} // End namespace Foam


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::exprTools::expressionEntry>
Foam::exprTools::expressionEntry::New
(
    const word& name
)
{
    auto cstrIter = emptyConstructorTablePtr_->cfind(name);

    if (!cstrIter.found())
    {
        FatalErrorInLookup
        (
            "expressionEntry",
            name,
            *emptyConstructorTablePtr_
        )  << exit(FatalError);
    }

    return autoPtr<expressionEntry>(cstrIter()());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::expressions::exprString
Foam::exprTools::expressionEntry::expand
(
    const std::string& orig,
    const dictionary& dict
)
{
    // This is much like stringOps::inplaceExpand
    constexpr const char sigil = '$';

    // Copy to a exprString, without any validation (using assign)
    expressions::exprString s;
    s.assign(orig);

    std::string::size_type varBeg = 0;

    while
    (
        (varBeg = s.find(sigil, varBeg)) != std::string::npos
     // && varBeg < s.size()-1
    )
    {
        // No handling of escape characters

        if (varBeg == s.size()-1)
        {
            // Die if we ended with a '$'
            FatalErrorInFunction
                << "'" << sigil << "' found at end of " << s
                << "(originally " << orig << ')' << nl
                << exit(FatalError);
        }

        std::string::size_type varEnd = varBeg;
        std::string::size_type delim = 0;

        word castTo, varName;

        if (s[varBeg+1] == '[')
        {
            // An expression pattern with $[...]

            varEnd = s.find(']', varBeg);
            delim = 1;

            if (varEnd == std::string::npos)
            {
                FatalErrorInFunction
                    << "No correct terminating ']' found in " << s
                    << " (originally " << orig << ")" << nl
                    << exit(FatalError);
            }

            // Look for embedded (type) cast

            const auto lparen = varBeg+2;
            if (lparen < s.size() && s[lparen] == '(')
            {
                const auto rparen = s.find(')', lparen);

                if (rparen > varEnd)
                {
                    // Handles both "$[( ...]" and "$[( ...])" cases

                    auto& err = FatalErrorInFunction;

                    if (rparen == std::string::npos)
                    {
                        err << "No closing ')' found in ";
                    }
                    else
                    {
                        err << "Closing ')' found outside of";
                    }

                    err << " substring "
                        << s.substr(varBeg, varEnd-varBeg)
                        << " (" << orig << ')' << nl
                        << exit(FatalError);
                }

                castTo.assign(s.substr(lparen+1, rparen - lparen - 1));
                varName.assign(s.substr(rparen+1, varEnd - rparen - 1));
            }
            else
            {
                varName.assign
                (
                    s.substr(varBeg + 1 + delim, varEnd - varBeg - 2*delim)
                );
            }

            stringOps::inplaceTrim(varName);
        }
        else
        {
            if (s[varBeg+1] == '{')
            {
                varEnd = s.find('}', varBeg);
                delim = 1;
            }
            else
            {
                // Handling regular $var construct
                varEnd += findVariableLen(s, varBeg, sigil);
            }

            if (varEnd == std::string::npos)
            {
                // Likely parsed '${...' without closing '}' - abort
                break;
            }
            else if (varEnd == varBeg)
            {
                // Parsed '${}' or $badChar  - skip over or die?
                FatalErrorInFunction
                    << "No valid character after the $ in " << s
                    << "(originally " << orig << ")" << endl
                    << exit(FatalError);
            }
            else
            {
                // Assign - assumed to be validated with findVariableLen()
                varName.assign
                (
                    s.substr(varBeg + 1 + delim, varEnd - varBeg - 2*delim)
                );
            }
        }


        // Length of original text to replace (incl. decorators)
        const auto replaceLen = (varEnd - varBeg + 1);

        const entry* eptr = getVariableOrDie(varName, dict);

        std::string varValue;

        if (castTo.empty())
        {
            // Serialized with spaces
            varValue = eptr->stream().toString();
        }
        else
        {
            varValue = expressionEntry::New(castTo)->toExpr(*eptr);
        }

        s.std::string::replace(varBeg, replaceLen, varValue);
        varBeg += varValue.size();
    }

    return s;
}


Foam::expressions::exprString
Foam::exprTools::expressionEntry::getExpression
(
    const word& name,
    const dictionary& dict,
    const bool removeComments
)
{
    string str(dict.get<string>(name));

    if (removeComments)
    {
        stringOps::inplaceRemoveComments(str);
    }

    return expand(str, dict);
}


// ************************************************************************* //
