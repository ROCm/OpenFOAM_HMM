/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider <bgschaid@hfd-research.com>
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

#include "exprDriver.H"
#include "expressionEntry.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

    defineTypeNameAndDebug(exprDriver, 0);

} // End namespace expressions
} // End namespace Foam



// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
#if 0
static string getEntryString
(
    const dictionary& dict,
    const string& key
)
{
    const entry* eptr = dict.findEntry(key, keyType::REGEX_RECURSIVE);

    if (!eptr)
    {
        FatalErrorInFunction
            << "Entry " << key << " not found in "
            << dict.name() << nl
            << exit(FatalError);
    }
    else if (eptr->isDict())
    {
        FatalErrorInFunction
            << "Entry " << key << " found in "
            << dict.name() << " but is a dictionary" << nl
            << exit(FatalError);
    }

    return exprTools::expressionEntry::evaluate(*eptr);
}
#endif
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::exprDriver::exprDriver
(
    bool cacheReadFields,
    bool searchInMemory,
    bool searchFiles,
    const dictionary& dict
)
:
    dict_(dict),
    result_(),
    variableStrings_(),
    variables_(),
    stashedTokenId_(0),

    // Controls
    debugScanner_(dict.lookupOrDefault("debugScanner", false)),
    debugParser_(dict.lookupOrDefault("debugParser", false)),
    allowShadowing_
    (
        dict.lookupOrDefault("allowShadowing", false)
    ),
    prevIterIsOldTime_
    (
        dict.lookupOrDefault("prevIterIsOldTime", false)
    ),
    cacheReadFields_(cacheReadFields),
    searchInMemory_(searchInMemory || cacheReadFields),
    searchFiles_(searchFiles)
{}


Foam::expressions::exprDriver::exprDriver
(
    const exprDriver& rhs
)
:
    dict_(rhs.dict_),
    result_(rhs.result_),
    variableStrings_(rhs.variableStrings_),
    variables_(rhs.variables_),
    stashedTokenId_(0),

    debugScanner_(rhs.debugScanner_),
    debugParser_(rhs.debugParser_),
    allowShadowing_(rhs.allowShadowing_),
    prevIterIsOldTime_(rhs.prevIterIsOldTime_),

    cacheReadFields_(rhs.cacheReadFields_),
    searchInMemory_(rhs.searchInMemory_),
    searchFiles_(rhs.searchFiles_)
{}


Foam::expressions::exprDriver::exprDriver
(
    const dictionary& dict
)
:
    exprDriver
    (
        dict.lookupOrDefault("cacheReadFields", false),
        dict.lookupOrDefault("searchInMemory", true),
        dict.lookupOrDefault("searchFiles", false),
        dict
    )
{
    readDict(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::exprDriver::readDict
(
    const dictionary& dict
)
{
    dict.readIfPresent("debugBaseDriver", debug);

    // Regular variables
    variableStrings_ = readVariableStrings(dict);

    // Other tables?
    // readTable("timelines", dict, lines_);
    // readTable("lookuptables", dict, lookup_);
    // readTable("lookuptables2D", dict, lookup2D_);

    return true;
}


void Foam::expressions::exprDriver::clearResult()
{
    result_.clear();
}


bool Foam::expressions::exprDriver::update()
{
    return true;
}


void Foam::expressions::exprDriver::updateSpecialVariables(bool force)
{}


void Foam::expressions::exprDriver::clearVariables()
{
    variables_.clear();
    addVariables(variableStrings_, false);
}


void Foam::expressions::exprDriver::evaluateVariable
(
    const word& varName,
    const expressions::exprString& expr
)
{
    parse(expr);
    result_.testIfSingleValue();

    DebugInfo
        << "Evaluating: " << expr << " -> " << varName << endl
        << result_;

    // Overwrite with a copy
    variables_.set(varName, exprResult(result_));
}


void Foam::expressions::exprDriver::evaluateVariableRemote
(
    string remote,
    const word& varName,
    const expressions::exprString& expr
)
{
    NotImplemented;
}


Foam::expressions::exprResult
Foam::expressions::exprDriver::getRemoteResult
(
    const exprDriver& other
) const
{
    // With warnings (noWarn = false)
    return other.result().getUniform(this->size(), false);
}


void Foam::expressions::exprDriver::addVariables
(
    const expressions::exprString& expr,
    bool clear
)
{
    if (clear)
    {
        clearVariables();
    }

    // Allow inline list of semicolon-separated variables
    const auto varExpressions =
        stringOps::split<expressions::exprString>(expr, ';');

    for (const auto& subMatch : varExpressions)
    {
        string varExpr(stringOps::trim(subMatch.str()));
        if (varExpr.empty())
        {
            continue;
        }

        // Split on '=' for lhsExpr = rhsExpr
        //
        // varName = rhsExpr
        // varName{where} = rhsExpr

        const auto eqPos = varExpr.find('=');

        if (eqPos == std::string::npos)
        {
            FatalIOErrorInFunction(dict_)
                << "No '=' found in expression " << varExpr << nl << nl
                << exit(FatalIOError);
        }

        // The RHS
        expressions::exprString rhsExpr
        (
            expressions::exprString::toExpr
            (
                stringOps::trim(varExpr.substr(eqPos+1))
            )
        );

        // The LHS
        varExpr.resize(eqPos);
        stringOps::inplaceTrim(varExpr);

        // Check for varName{where}
        const auto lbrace = varExpr.find('{');

        if (lbrace != std::string::npos)
        {
            const auto rbrace = varExpr.find('}');

            if (rbrace == std::string::npos || rbrace < lbrace)
            {
                FatalErrorInFunction
                    // << "Context: " << driverContext_ << nl
                    << "No closing '}' found in " << varExpr << nl
                    << exit(FatalError);
            }
            else if (lbrace+1 == rbrace)
            {
                FatalErrorInFunction
                    // << "Context: " << driverContext_ << nl
                    << "Empty '{}' location in " << varExpr << nl
                    << exit(FatalError);
            }

            const word varName(word::validate(varExpr.substr(0, lbrace)));

            const expressions::exprString remoteExpr
            (
                expressions::exprString::toExpr
                (
                    varExpr.substr(lbrace+1, rbrace-lbrace-1)
                )
            );

            // Fails if derived class does not implement!

            evaluateVariableRemote(remoteExpr, varName, rhsExpr);
        }
        else
        {
            const word varName(word::validate(varExpr));

            evaluateVariable(varName, rhsExpr);
        }
    }
}


void Foam::expressions::exprDriver::addVariables
(
    const UList<expressions::exprString>& list,
    bool clear
)
{
    if (clear)
    {
        clearVariables();
    }

    for (const auto& expr : list)
    {
        addVariables(expr, false); // No clear (already done)
    }
}


void Foam::expressions::exprDriver::setDebugging
(
    bool scannerDebug,
    bool parserDebug
)
{
    debugScanner_ = scannerDebug;
    debugParser_ = parserDebug;
}


void Foam::expressions::exprDriver::setDebugging
(
    const exprDriver& rhs
)
{
    debugScanner_ = rhs.debugScanner_;
    debugParser_ = rhs.debugParser_;
}


void Foam::expressions::exprDriver::setSearchBehaviour
(
    bool cacheReadFields,
    bool searchInMemory,
    bool searchFiles
)
{
    searchInMemory_ = searchInMemory_ || cacheReadFields_;

    #ifdef FULLDEBUG
    Info<< "Searching "
        << " registry:" << searchInMemory_
        << " disk:" << searchFiles_
        << " cache-read:" << cacheReadFields_ << nl;
    #endif
}


void Foam::expressions::exprDriver::setSearchBehaviour
(
    const exprDriver& rhs
)
{
    setSearchBehaviour
    (
        rhs.cacheReadFields_,
        rhs.searchInMemory_,
        rhs.searchFiles_
    );
}


// ************************************************************************* //
