/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "pimpleControl.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::pimpleControl::read()
{
    solutionControl::read(false);

    // Read solution controls
    const dictionary& pimpleDict = dict();
    nOuterCorr_ = pimpleDict.lookupOrDefault<label>("nOuterCorrectors", 1);
    nCorr_ = pimpleDict.lookupOrDefault<label>("nCorrectors", 1);
    turbOnFinalIterOnly_ =
        pimpleDict.lookupOrDefault<Switch>("turbOnFinalIterOnly", true);
}


bool Foam::pimpleControl::criteriaSatisfied()
{
    if ((corr_ == 0) || residualControl_.empty() || finalIter())
    {
        return false;
    }

    bool firstIter = corr_ == 1;

    bool achieved = true;
    const dictionary& solverDict = mesh_.solverPerformanceDict();

    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        label fieldI = applyToField(variableName);
        if (fieldI != -1)
        {
            const List<lduMatrix::solverPerformance> sp(iter().stream());
            const scalar residual = sp.last().initialResidual();

            if (firstIter)
            {
                residualControl_[fieldI].initialResidual =
                    sp.first().initialResidual();
            }

            bool absCheck = residual < residualControl_[fieldI].absTol;

            bool relCheck = false;

            scalar relative = 0.0;
            if (!firstIter)
            {
                scalar iniRes =
                    residualControl_[fieldI].initialResidual
                  + ROOTVSMALL;

                relative = residual/iniRes;

                relCheck = relative < residualControl_[fieldI].relTol;
            }

            achieved = achieved && (absCheck || relCheck);

            if (debug)
            {
                Info<< algorithmName_ << "loop statistics:" << endl;

                Info<< "    " << variableName << " iter " << corr_
                    << ": ini res = "
                    << residualControl_[fieldI].initialResidual
                    << ", abs tol = " << residual
                    << " (" << residualControl_[fieldI].absTol << ")"
                    << ", rel tol = " << relative
                    << " (" << residualControl_[fieldI].relTol << ")"
                    << endl;
            }
        }
    }

    return achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimpleControl::pimpleControl(fvMesh& mesh)
:
    solutionControl(mesh, "PIMPLE"),
    nOuterCorr_(0),
    nCorr_(0),
    corr_(0),
    turbOnFinalIterOnly_(true)
{
    read();

    if (nOuterCorr_ > 1)
    {
        Info<< nl;
        if (!residualControl_.empty())
        {
            Info<< algorithmName_ << ": max iterations = " << nOuterCorr_
                << endl;
            forAll(residualControl_, i)
            {
                Info<< "    field " << residualControl_[i].name << token::TAB
                    << ": relTol " << residualControl_[i].relTol
                    << ", tolerance " << residualControl_[i].absTol
                    << nl;
            }
            Info<< endl;
        }
        else
        {
            Info<< algorithmName_ << ": no residual control data found. " << nl
                << "Calculations will employ " << nOuterCorr_
                << " corrector loops" << nl << endl;
        }
    }
    else
    {
        Info<< nl << algorithmName_ << ": Operating solver in PISO mode" << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleControl::~pimpleControl()
{}


// ************************************************************************* //
