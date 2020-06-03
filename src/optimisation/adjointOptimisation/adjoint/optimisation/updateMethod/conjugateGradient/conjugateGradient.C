/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "conjugateGradient.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conjugateGradient, 0);
    addToRunTimeSelectionTable
    (
        updateMethod,
        conjugateGradient,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::conjugateGradient::allocateFields()
{
    // Set active design variables, if necessary
    if (activeDesignVars_.empty())
    {
        activeDesignVars_ = identity(objectiveDerivatives_.size());
    }

    // Allocate old fields
    dxOld_ = scalarField(activeDesignVars_.size(), Zero);
    sOld_ = scalarField(activeDesignVars_.size(), Zero);
}


void Foam::conjugateGradient::readFromDict()
{
    if (optMethodIODict_.headerOk())
    {
        optMethodIODict_.readEntry("dxOld", dxOld_);
        optMethodIODict_.readEntry("sOld", sOld_);
        optMethodIODict_.readEntry("counter", counter_);
        optMethodIODict_.readEntry("eta", eta_);

        label nDVs = optMethodIODict_.get<label>("nDVs");
        correction_ = scalarField(nDVs, Zero);

        if (activeDesignVars_.empty())
        {
            activeDesignVars_ = identity(nDVs);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conjugateGradient::conjugateGradient
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    updateMethod(mesh, dict),

    activeDesignVars_(0),
    dxOld_(0),
    sOld_(0),
    counter_(0),
    betaType_
    (
        coeffsDict().getOrDefault<word>("betaType", "FletcherReeves")
    )
{
    if
    (
        !coeffsDict().readIfPresent("activeDesignVariables", activeDesignVars_)
    )
    {
        // If not, all available design variables will be used.
        // Number is not know at the moment
        Info<< "\t Did not find explicit definition of active design variables. "
            << "Treating all available ones as active " << endl;
    }

    // Check if beta type is valid
    if
    (
       !(betaType_ == "FletcherReeves")
    && !(betaType_ == "PolakRibiere")
    && !(betaType_ == "PolakRibiereRestarted")
    )
    {
        FatalErrorInFunction
           << "Invalid betaType " << betaType_ << ". Valid options are "
           << "FletcherReeves, PolakRibiere, PolakRibiereRestarted"
           << nl << nl
           << exit(FatalError);
    }

    // Read old dx and s, if present
    readFromDict();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conjugateGradient::computeCorrection()
{
    if (counter_ == 0)
    {
        allocateFields();

        Info<< "Using steepest descent for the first iteration" << endl;
        correction_ = -eta_*objectiveDerivatives_;

        dxOld_.map(-objectiveDerivatives_, activeDesignVars_);
        sOld_ = dxOld_;
    }
    else
    {
        scalarField dx = scalarField(activeDesignVars_.size(), Zero);
        dx.map(-objectiveDerivatives_, activeDesignVars_);

        scalar beta(Zero);
        if (betaType_ == "FletcherReeves")
        {
            beta = globalSum(dx*dx)/globalSum(dxOld_ * dxOld_);
        }
        else if (betaType_ == "PolakRibiere")
        {
            beta = globalSum(dx*(dx - dxOld_))/globalSum(dxOld_ * dxOld_);
        }
        else if (betaType_ == "PolakRibiereRestarted")
        {
            beta =
            max
            (
                scalar(0),
                globalSum(dx*(dx - dxOld_))/globalSum(dxOld_ * dxOld_)
            );
            if (beta == scalar(0))
            {
                Info<< "Computed negative beta. Resetting to zero" << endl;
            }
        }

        scalarField s(dx + beta*sOld_);

        correction_ = Zero;
        forAll(activeDesignVars_, varI)
        {
            correction_[activeDesignVars_[varI]] = eta_*s[varI];
        }

        // Store fields for the next iteration
        dxOld_ = dx;
        sOld_ = s;
    }

    ++counter_;
}


void Foam::conjugateGradient::updateOldCorrection
(
    const scalarField& oldCorrection
)
{
    sOld_.map(oldCorrection, activeDesignVars_);
    sOld_ /= eta_;
    correction_ = oldCorrection;
}


void Foam::conjugateGradient::write()
{
    optMethodIODict_.add<scalarField>("dxOld", dxOld_, true);
    optMethodIODict_.add<scalarField>("sOld", sOld_, true);
    optMethodIODict_.add<label>("counter", counter_, true);
    optMethodIODict_.add<label>("nDVs", objectiveDerivatives_.size(), true);

    updateMethod::write();
}


// ************************************************************************* //
