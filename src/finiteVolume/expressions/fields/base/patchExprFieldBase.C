/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Original code Copyright (C) 2011-2018 Bernhard Gschaider
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

#include "patchExprFieldBase.H"
#include "facePointPatch.H"
#include "fvMesh.H"
#include "fvPatch.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::fvPatch&
Foam::expressions::patchExprFieldBase::getFvPatch(const facePointPatch& pp)
{
    const polyMesh& pmesh = pp.boundaryMesh().mesh().mesh();

    const fvMesh* meshptr = isA<fvMesh>(pmesh);

    if (!meshptr)
    {
        FatalErrorInFunction
            << "Point patch not attached to a base fvMesh, "
            << "cannot use patch expressions" << nl << endl
            << exit(FatalError);
    }

    return meshptr->boundary()[pp.index()];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::patchExprFieldBase::patchExprFieldBase()
:
    patchExprFieldBase(false)
{}


Foam::expressions::patchExprFieldBase::patchExprFieldBase
(
    bool allowGradient
)
:
    debug_(false),
    allowGradient_(allowGradient),
    evalOnConstruct_(false),
    valueExpr_(),
    gradExpr_(),
    fracExpr_()
{}


Foam::expressions::patchExprFieldBase::patchExprFieldBase
(
    const dictionary& dict,
    bool allowGradient,
    bool isPointVal
)
:
    debug_(dict.lookupOrDefault("debug", false)),
    allowGradient_(allowGradient),
    evalOnConstruct_(dict.lookupOrDefault<bool>("evalOnConstruct", false)),
    valueExpr_(),
    gradExpr_(),
    fracExpr_()
{
    if (debug_)
    {
        Info<< "Expression BC with " << dict << nl;
    }

    string expr;

    if (dict.readIfPresent("valueExpr", expr))
    {
        valueExpr_ = expressions::exprString(expr, dict);
    }
    else
    {
        // No value expression - same as Zero
        if (debug_)
        {
            Info<< "No valueExpr" << nl;
        }
    }

    if (allowGradient)
    {
        if (dict.readIfPresent("gradientExpr", expr))
        {
            gradExpr_ = expressions::exprString(expr, dict);
        }
        else
        {
            // No gradient expression - same as Zero

            if (debug_)
            {
                Info<< "No gradientExpr" << nl;
            }
        }

        if (dict.readIfPresent("fractionExpr", expr))
        {
            if (!expr.empty() && expr != "0")
            {
                if (isPointVal)
                {
                    expr = "toPoint(" + expr + ")";
                }

                fracExpr_ = expressions::exprString(expr, dict);
            }
        }
        else
        {
            // No fraction expression - same as 1 (one)
            // Mixed BC may elect to simply ignore gradient expression
        }
    }
}


Foam::expressions::patchExprFieldBase::patchExprFieldBase
(
    const patchExprFieldBase& rhs
)
:
    debug_(rhs.debug_),
    allowGradient_(rhs.allowGradient_),
    evalOnConstruct_(rhs.evalOnConstruct_),
    valueExpr_(rhs.valueExpr_),
    gradExpr_(rhs.gradExpr_),
    fracExpr_(rhs.fracExpr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::patchExprFieldBase::write(Ostream& os) const
{
    os.writeEntryIfDifferent<bool>("evalOnConstruct", false, evalOnConstruct_);

    // Do not emit debug_ value

    if (!valueExpr_.empty())
    {
        os.writeEntry("valueExpr", valueExpr_);
    }
    if (!gradExpr_.empty())
    {
        os.writeEntry("gradientExpr", gradExpr_);
    }
    if (!fracExpr_.empty())
    {
        os.writeEntry("fractionExpr", fracExpr_);
    }
}


// ************************************************************************* //
