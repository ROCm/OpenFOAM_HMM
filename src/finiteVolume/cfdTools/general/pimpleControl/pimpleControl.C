/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimpleControl, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::pimpleControl::applyToField(const word& fieldName) const
{
    forAll(residualControl_, i)
    {
        if (residualControl_[i].name.match(fieldName))
        {
            return i;
        }
    }

    return -1;
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
                Info<< "PIMPLE loop statistics:" << endl;

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
    mesh_(mesh),
    nCorr_(0),
    corr_(0),
    residualControl_()
{
    read();

    if (residualControl_.size() > 0)
    {
        Info<< "PIMPLE: max iterations = " << nCorr_ << endl;
        forAll(residualControl_, i)
        {
            Info<< "    field " << residualControl_[i].name << token::TAB
                << ": relTol " << residualControl_[i].relTol
                << ", absTol " << residualControl_[i].absTol
                << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "No PIMPLE residual control data found. "
            << "Calculations will employ a fixed number of corrector loops"
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pimpleControl::~pimpleControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pimpleControl::read()
{
    const dictionary& dict = mesh_.solutionDict().subDict("PIMPLE");
    nCorr_ = dict.lookupOrDefault<label>("nOuterCorrectors", 1);
    const dictionary residualDict(dict.subOrEmptyDict("residualControl"));
    DynamicList<fieldData> data(residualDict.toc().size());
    wordHashSet fieldNames;
    forAllConstIter(dictionary, residualDict, iter)
    {
        if (fieldNames.insert(iter().keyword()))
        {
            fieldData fd;
            fd.name = iter().keyword().c_str();
            if (iter().isDict())
            {
                const dictionary& fieldDict(iter().dict());
                fd.relTol = readScalar(fieldDict.lookup("relTol"));
                fd.absTol = readScalar(fieldDict.lookup("absTol"));
                fd.initialResidual = 0.0;
                data.append(fd);
            }
            else
            {
                FatalErrorIn("bool Foam::pimpleControl::read()")
                    << "Residual data for " << iter().keyword()
                    << " must be specified as a dictionary";
            }
        }
    }

    residualControl_.transfer(data);
}


// ************************************************************************* //
