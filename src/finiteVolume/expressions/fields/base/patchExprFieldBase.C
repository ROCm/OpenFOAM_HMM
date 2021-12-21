/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 Bernhard Gschaider
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

#include "patchExprFieldBase.H"
#include "facePointPatch.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "pointMesh.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::expressions::patchExprFieldBase::readExpressions
(
    const dictionary& dict,
    enum expectedTypes expectedType,
    bool wantPointData
)
{
    if (debug_)
    {
        Info<< "Expression BC with " << dict << nl;
    }

    valueExpr_.clear();
    gradExpr_.clear();
    fracExpr_.clear();

    string exprValue, exprGrad, exprFrac;
    bool evalValue = false, evalGrad = false, evalFrac = false;

    if (expectedTypes::VALUE_TYPE == expectedType)
    {
        // Mandatory
        evalValue = dict.readEntry("valueExpr", exprValue, keyType::LITERAL);
    }
    else if (expectedTypes::GRADIENT_TYPE == expectedType)
    {
        // Mandatory
        evalGrad = dict.readEntry("gradientExpr", exprGrad, keyType::LITERAL);
    }
    else
    {
        // MIXED_TYPE
        evalValue =
            dict.readIfPresent("valueExpr", exprValue, keyType::LITERAL);

        evalGrad =
            dict.readIfPresent("gradientExpr", exprGrad, keyType::LITERAL);

        if (!evalValue && !evalGrad)
        {
            FatalIOErrorInFunction(dict)
                << "Entries 'valueExpr' and 'gradientExpr' "
                   "(mixed-conditon) not found in dictionary "
                << dict.name() << nl
                << exit(FatalIOError);
        }

        if (debug_)
        {
            if (!evalValue)
            {
                Info<< "Mixed with no valueExpr" << nl;
            }
            if (!evalGrad)
            {
                Info<< "Mixed with no gradientExpr" << nl;
            }
        }
    }


    // When both value/gradient specified (ie, mixed) expect a fraction
    // - if missing, defer treatment to inherited BC

    if (evalValue && evalGrad && dict.readIfPresent("fractionExpr", exprFrac))
    {
        stringOps::inplaceTrim(exprFrac);

        if (exprFrac == "0" || exprFrac == "1")
        {
            // Special cases, handled with more efficiency
            fracExpr_ = exprFrac;
        }
        else if (!exprFrac.empty())
        {
            evalFrac = true;
            if (wantPointData)
            {
                exprFrac = "toPoint(" + exprFrac + ")";
            }
        }
    }


    // Expansions

    if (evalValue)
    {
        valueExpr_ = expressions::exprString(exprValue, dict);
    }
    if (evalGrad)
    {
        gradExpr_ = expressions::exprString(exprGrad, dict);
    }
    if (evalFrac)
    {
        fracExpr_ = expressions::exprString(exprFrac, dict);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::patchExprFieldBase::patchExprFieldBase()
:
    debug_(false),
    evalOnConstruct_(false),
    valueExpr_(),
    gradExpr_(),
    fracExpr_()
{}


Foam::expressions::patchExprFieldBase::patchExprFieldBase
(
    const dictionary& dict,
    enum expectedTypes expectedType,
    bool wantPointData
)
:
    debug_(dict.getOrDefault("debug", false)),
    evalOnConstruct_(dict.getOrDefault("evalOnConstruct", false)),
    valueExpr_(),
    gradExpr_(),
    fracExpr_()
{
    readExpressions(dict, expectedType, wantPointData);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::patchExprFieldBase::write(Ostream& os) const
{
    os.writeEntryIfDifferent<bool>("evalOnConstruct", false, evalOnConstruct_);

    // Do not emit debug_ value

    // Write expression, but not empty ones
    valueExpr_.writeEntry("valueExpr", os, false);
    gradExpr_.writeEntry("gradientExpr", os, false);
    fracExpr_.writeEntry("fractionExpr", os, false);
}


// ************************************************************************* //
