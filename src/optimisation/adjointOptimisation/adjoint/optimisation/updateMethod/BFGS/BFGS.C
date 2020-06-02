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

#include "BFGS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BFGS, 0);
    addToRunTimeSelectionTable
    (
        updateMethod,
        BFGS,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::BFGS::allocateMatrices()
{
    // Set active design variables, if necessary
    if (activeDesignVars_.empty())
    {
        activeDesignVars_ = identity(objectiveDerivatives_.size());
    }

    // Set previous HessianInv to be a diagonal matrix
    SquareMatrix<scalar> temp(activeDesignVars_.size(), I);

    // Allocate correct size and content to HessianInv matrices
    // has a max. capability of approximately 34000 variables.
    HessianInvOld_ = temp;
    HessianInv_ = temp;
}


void Foam::BFGS::updateHessian()
{
    // Vectors needed to construct the inverse HessianInv matrix
    scalarField y(activeDesignVars_.size(), Zero);
    scalarField s(activeDesignVars_.size(), Zero);
    y.map(objectiveDerivatives_ - derivativesOld_, activeDesignVars_);
    s.map(correctionOld_, activeDesignVars_);

    scalar ys = globalSum(s*y);

    if (counter_ == 1 && scaleFirstHessian_)
    {
        scalar scaleFactor = ys/globalSum(y*y);
        Info<< "Scaling Hessian with factor " << scaleFactor << endl;
        forAll(activeDesignVars_, varI)
        {
            HessianInvOld_[varI][varI] *= scaleFactor;
        }
    }

    // Check curvature condition
    if (ys < curvatureThreshold_)
    {
        WarningInFunction
            << " y*s is below threshold! y*s=" << ys << endl;
    }

    DebugInfo
        << "Hessian curvature index " << ys << endl;

    // Construct the inverse HessianInv
    HessianInv_ =
        HessianInvOld_
      + (ys + globalSum(leftMult(y, HessianInvOld_)*y))/sqr(ys)*outerProd(s, s)
      - (scalar(1)/ys)*
        (
           outerProd(rightMult(HessianInvOld_, y), s)
         + outerProd(s, leftMult(y, HessianInvOld_))
        );
}


void Foam::BFGS::update()
{
    // In the first few iterations, use steepest descent but update the Hessian
    // matrix
    if (counter_ < nSteepestDescent_)
    {
        Info<< "Using steepest descent to update design variables" << endl;
        correction_ = -eta_*objectiveDerivatives_;
    }
    // else use BFGS formula to update the design variables
    else
    {
        // Compute correction for active design variables
        scalarField activeDerivs(activeDesignVars_.size(), Zero);
        activeDerivs.map(objectiveDerivatives_, activeDesignVars_);
        scalarField activeCorrection
        (
            -etaHessian_*rightMult(HessianInv_, activeDerivs)
        );

        // Transfer correction to the global list
        correction_ = Zero;
        forAll(activeDesignVars_, varI)
        {
            correction_[activeDesignVars_[varI]] = activeCorrection[varI];
        }
    }

    // Store fields for the next iteration
    derivativesOld_ = objectiveDerivatives_;
    correctionOld_ = correction_;
    HessianInvOld_ = HessianInv_;
}


void Foam::BFGS::readFromDict()
{
    if (optMethodIODict_.headerOk())
    {
        optMethodIODict_.readEntry("HessianInvOld", HessianInvOld_);
        optMethodIODict_.readEntry("derivativesOld", derivativesOld_);
        optMethodIODict_.readEntry("correctionOld", correctionOld_);
        optMethodIODict_.readEntry("counter", counter_);
        optMethodIODict_.readEntry("eta", eta_);

        const label n(HessianInvOld_.n());
        HessianInv_ = SquareMatrix<scalar>(n, Zero);
        correction_ = scalarField(correctionOld_.size(), Zero);

        if (activeDesignVars_.empty())
        {
            activeDesignVars_ = identity(n);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BFGS::BFGS
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    updateMethod(mesh, dict),
    etaHessian_
    (
        coeffsDict().getOrDefault<scalar>("etaHessian", 1)
    ),
    nSteepestDescent_
    (
        coeffsDict().getOrDefault<label>("nSteepestDescent", 1)
    ),
    activeDesignVars_(0),
    scaleFirstHessian_
    (
        coeffsDict().getOrDefault<bool>("scaleFirstHessian", false)
    ),
    curvatureThreshold_
    (
        coeffsDict().getOrDefault<scalar>("curvatureThreshold", 1e-10)
    ),
    // Construct null matrix since we dont know the dimension yet
    HessianInv_(),
    HessianInvOld_(),
    derivativesOld_(0),
    correctionOld_(0),
    counter_(Zero)
{
    if
    (
        !coeffsDict().readIfPresent("activeDesignVariables", activeDesignVars_)
    )
    {
        // If not, all available design variables will be used. Number is not
        // know at the moment
        Info<< "\t Did not find explicit definition of active design variables."
            << " Treating all available ones as active" << endl;
    }

    // Read old hessian, correction and derivatives, if present
    readFromDict();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BFGS::computeCorrection()
{
    if (counter_ == 0)
    {
        allocateMatrices();
    }
    else
    {
        updateHessian();
    }

    update();
    ++counter_;
}


void Foam::BFGS::updateOldCorrection(const scalarField& oldCorrection)
{
    updateMethod::updateOldCorrection(oldCorrection);
    correctionOld_ = oldCorrection;
}


void Foam::BFGS::write()
{
    optMethodIODict_.add<SquareMatrix<scalar>>
    (
        "HessianInvOld",
        HessianInvOld_,
        true
    );
    optMethodIODict_.add<scalarField>("derivativesOld", derivativesOld_, true);
    optMethodIODict_.add<scalarField>("correctionOld", correctionOld_, true);
    optMethodIODict_.add<label>("counter", counter_, true);

    updateMethod::write();
}


// ************************************************************************* //
