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

#include "evalEntry.H"
#include "dictionary.H"
#include "OTstream.H"
#include "stringOps.H"
#include "fieldExprDriver.H"
#include "addToMemberFunctionSelectionTable.H"
#include <cctype>

#undef  DetailInfo
#define DetailInfo  if (::Foam::infoDetailLevel > 0) InfoErr


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        evalEntry,
        execute,
        primitiveEntryIstream,
        eval
    );

} // End namespace functionEntry
} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{
    // This is akin to a SafeIOWarning, which does not yet exist
    inline void safeIOWarning
    (
        const Foam::IOstream& is,
        const std::string& msg
    )
    {
        std::cerr
            << "--> FOAM Warning :\n"
            << "    Reading \"" << is.name() << "\" at line "
            << is.lineNumber() << '\n'
            << "    " << msg << std::endl;
    }

} // End anonymous namespace


namespace Foam
{

// Slurp a string until a closing '}' is found.
// Track balanced bracket/brace pairs, with max stack depth of 60.
static bool slurpUntilBalancedBrace(ISstream& is, std::string& str)
{
    constexpr const unsigned bufLen = 1024;
    static char buf[bufLen];

    is.fatalCheck(FUNCTION_NAME);

    unsigned nChar = 0;
    unsigned depth = 1; // Initial '{' already seen by caller
    char c;

    str.clear();
    while (is.get(c))
    {
        buf[nChar++] = c;

        if (c == token::BEGIN_BLOCK)
        {
            ++depth;
        }
        else if (c == token::END_BLOCK)
        {
            --depth;
            if (!depth)
            {
                // Closing '}' character - do not include in output
                --nChar;
                str.append(buf, nChar);
                return true;
            }
        }
        else if (c == '/')
        {
            // Strip C/C++ comments from expressions
            // Note: could also peek instead of get/putback

            if (!is.get(c))
            {
                break;  // Premature end of stream
            }
            else if (c == '/')
            {
                --nChar;  // Remove initial '/' from buffer

                // C++ comment: discard through newline
                (void) is.getLine(nullptr, '\n');
            }
            else if (c == '*')
            {
                --nChar;  // Remove initial '/' from buffer

                // C-style comment: discard through to "*/" ending
                if (!is.seekCommentEnd_Cstyle())
                {
                    break;  // Premature end of stream
                }
            }
            else
            {
                // Reanalyze the char
                is.putback(c);
            }
        }

        if (nChar == bufLen)
        {
            str.append(buf, nChar);  // Flush full buffer
            nChar = 0;
        }
    }


    // Abnormal exit of the loop

    str.append(buf, nChar);  // Finalize pending content

    safeIOWarning(is, "Premature end while reading expression - missing '}'?");

    is.fatalCheck(FUNCTION_NAME);
    return false;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tokenList Foam::functionEntries::evalEntry::evaluate
(
    const dictionary& parentDict,
    const string& inputExpr,
    label fieldWidth,
    const Istream& is
)
{
    // Field width for the result
    if (fieldWidth < 1)
    {
        FatalIOErrorInFunction(is)
            << "Invalid field width: " << fieldWidth << nl << endl
            << exit(FatalIOError);
    }

    #ifdef FULLDEBUG
    DetailInfo
        << "input: " << inputExpr << endl;
    #endif

    // Expand with env=true, empty=true, subDict=false
    // with comments stripped.
    // Special handling of $[...] syntax enabled.

    string s;

    // Passed '${{ expr }}' by accident, or on purpuse
    if
    (
        inputExpr[0] == token::DOLLAR
     && inputExpr[1] == token::BEGIN_BLOCK
     && inputExpr[2] == token::BEGIN_BLOCK
     && inputExpr[inputExpr.length()-1] == token::END_BLOCK
     && inputExpr[inputExpr.length()-2] == token::END_BLOCK
    )
    {
        s.assign(inputExpr, 3, inputExpr.length()-5);
    }
    else
    {
        s.assign(inputExpr);
    }

    expressions::exprString::inplaceExpand(s, parentDict, true);
    stringOps::inplaceTrim(s);

    // An extraneous trailing ';' is a common input error.
    // - trim if it does not influence the result

    const auto trailing = s.find(';');
    if (std::string::npos != trailing)
    {
        bool ignore = true;
        for (size_t other = trailing; ignore && other < s.length(); ++other)
        {
            ignore = s[other] == ';' || std::isspace(s[other]);
        }

        if (ignore)
        {
            // Can trim trailing without semantical change
            s.erase(trailing);
            stringOps::inplaceTrim(s);
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Invalid input (after trailing ';') for #eval" << nl
                << s << endl
                << exit(FatalIOError);
        }
    }

    #ifdef FULLDEBUG
    DetailInfo
        << "expanded: " << s << endl;
    #endif

    if (s.empty())
    {
        InfoErr
            << "Empty #eval - line "
            << is.lineNumber() << " in file " <<  parentDict.name() << nl;

        return tokenList();
    }

    expressions::exprResult result;
    {
        expressions::fieldExprDriver driver(fieldWidth);
        driver.parse(s);
        result = std::move(driver.result());
    }

    if (!result.hasValue() || !result.size())
    {
        InfoErr
            << "Failed #eval - line "
            << is.lineNumber() << " in file " <<  parentDict.name() << nl;

        return tokenList();
    }

    OTstream toks;
    if (result.size() <= 1)
    {
        result.writeValue(toks);
    }
    else
    {
        result.writeField(toks);
    }

    return std::move(toks);
}


Foam::tokenList Foam::functionEntries::evalEntry::evaluate
(
    const dictionary& parentDict,
    Istream& is
)
{
    #ifdef FULLDEBUG
    DetailInfo
        << "Using #eval - line "
        << is.lineNumber() << " in file " <<  parentDict.name() << nl;
    #endif

    token tok(is);
    label fieldWidth(1);  // Field width for the result
    if (tok.isLabel())
    {
        // - #eval INT "expr"
        // - #eval INT { expr }
        // - #eval INT #{ expr #}
        fieldWidth = max(1, tok.labelToken());
        is >> tok;
    }

    string str;  // The string to evaluate
    if (tok.isString())
    {
        // - #eval "expr"
        // - #eval #{ expr #}
        // - #eval ${{ expr }} - wierd but handled
        str = tok.stringToken();
    }
    else if (tok.isPunctuation(token::BEGIN_BLOCK))
    {
        // - #eval { expr }
        slurpUntilBalancedBrace(dynamic_cast<ISstream&>(is), str);
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Invalid input for #eval."
               " Expecting a string or block to evaluate, but found" << nl
            << tok.info() << endl
            << exit(FatalIOError);
    }

    tokenList toks
    (
        evalEntry::evaluate(parentDict, str, fieldWidth, is)
    );

    return toks;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::evalEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    tokenList toks(evaluate(parentDict, is));

    entry.append(std::move(toks), true);  // Lazy resizing

    return true;
}


bool Foam::functionEntries::evalEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    const string& inputExpr,
    label fieldWidth,
    Istream& is
)
{
    tokenList toks(evaluate(parentDict, inputExpr, fieldWidth, is));

    entry.append(std::move(toks), true);  // Lazy resizing

    return true;
}


// ************************************************************************* //
