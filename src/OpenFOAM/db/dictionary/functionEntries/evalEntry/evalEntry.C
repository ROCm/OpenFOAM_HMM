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

#include "evalEntry.H"
#include "dictionary.H"
#include "OTstream.H"
#include "stringOps.H"
#include "fieldExprDriver.H"
#include "addToMemberFunctionSelectionTable.H"

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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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

    // String to evaluate
    string s;

    token tok(is);

    if (!tok.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get string to evaluate"
            << exit(FatalIOError);

        return tokenList();
    }

    if (tok.isString())
    {
        s = tok.stringToken();
    }
    else if (tok == token::BEGIN_BLOCK)
    {
        dynamic_cast<ISstream&>(is).getLine(s, token::END_BLOCK);
    }
    else
    {
        is.putBack(tok);

        FatalIOErrorInFunction(is)
            << "Invalid input for #eval" << nl
            << exit(FatalIOError);
    }

    #ifdef FULLDEBUG
    DetailInfo
        << "input: " << s << endl;
    #endif

    // Expand with env=true, empty=true, subDict=false
    // with comments stripped.
    // Special handling of $[...] syntax enabled.
    expressions::exprString::inplaceExpand(s, parentDict, true);
    stringOps::inplaceTrim(s);

    // An extraneous trailing ';' is a common input error, catch it now.
    // May need to relax in the future, trim or something else

    if (std::string::npos != s.find(';'))
    {
        FatalIOErrorInFunction(is)
            << "Invalid input for #eval" << nl
            << s << endl
            << exit(FatalIOError);
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
        expressions::fieldExprDriver driver(1);
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

    // Could average/reduce to a single value, but probably not needed
    //// result.testIfSingleValue(false);  // No parallel check

    OTstream toks;
    result.writeValue(toks);

    return std::move(toks);
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


// ************************************************************************* //
