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

#include "DBFGS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DBFGS, 0);
    addToRunTimeSelectionTable
    (
        updateMethod,
        DBFGS,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::DBFGS::allocateMatrices()
{
    // Set active design variables, if necessary
    if (activeDesignVars_.empty())
    {
        activeDesignVars_ = identity(objectiveDerivatives_.size());
    }

    // Set previous Hessian to be a diagonal matrix
    SquareMatrix<scalar> temp(activeDesignVars_.size(), I);

    // Allocate correct size and content to Hessian matrices
    // has a max. capability of approximately 34000 variables.
    HessianOld_ = temp;
    Hessian_ = temp;
}


void Foam::DBFGS::updateHessian()
{
    // Vectors needed to construct the inverse Hessian matrix
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
            HessianOld_[varI][varI] /= scaleFactor;
        }
    }

    scalar sBs = globalSum(leftMult(s, HessianOld_)*s);

    // Check curvature condition and apply dampening is necessary
    scalar theta(1);
    if (ys < gamma_*sBs)
    {
        WarningInFunction
            << " y*s is below threshold. Using damped form" << endl;
        theta = (scalar(1)-gamma_)*sBs/(sBs - ys);
    }

    DebugInfo
        << "Hessian curvature index " << ys << endl;

    scalarField r(theta*y + (scalar(1)-theta)*rightMult(HessianOld_, s));

    // Construct the inverse Hessian
    Hessian_ =
        HessianOld_
      - outerProd(rightMult(HessianOld_, s), leftMult(s/sBs, HessianOld_))
      + outerProd(r, r/globalSum(s*r));
}


void Foam::DBFGS::update()
{
    SquareMatrix<scalar> HessianInv = inv(Hessian_);

    // In the first few iterations, use steepest descent but update the Hessian
    // matrix
    if (counter_ < nSteepestDescent_)
    {
        Info<< "Using steepest descent to update design variables" << endl;
        correction_ = -eta_*objectiveDerivatives_;
    }
    // Else use DBFGS formula to update design variables
    else
    {
        // Compute correction for active design variables
        scalarField activeDerivs(activeDesignVars_.size(), Zero);
        activeDerivs.map(objectiveDerivatives_, activeDesignVars_);
        scalarField activeCorrection
        (
            -etaHessian_*rightMult(HessianInv, activeDerivs)
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
    HessianOld_ = Hessian_;
}


void Foam::DBFGS::readFromDict()
{
    if (optMethodIODict_.headerOk())
    {
        optMethodIODict_.readEntry("HessianOld",  HessianOld_);
        optMethodIODict_.readEntry("derivativesOld", derivativesOld_);
        optMethodIODict_.readEntry("correctionOld", correctionOld_);
        optMethodIODict_.readEntry("counter", counter_);
        optMethodIODict_.readEntry("eta", eta_);

        label n = HessianOld_.n();
        Hessian_ = SquareMatrix<scalar>(n, Zero);
        correction_ = scalarField(correctionOld_.size(), Zero);

        if (activeDesignVars_.empty())
        {
            activeDesignVars_ = identity(n);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DBFGS::DBFGS
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    updateMethod(mesh, dict),

    // Construct null matrix since we dont know the dimension yet
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
    Hessian_(),
    HessianOld_(),
    derivativesOld_(0),
    correctionOld_(0),
    counter_(0),
    gamma_(coeffsDict().getOrDefault<scalar>("gamma", 0.2))

{
    if
    (
        !coeffsDict().readIfPresent("activeDesignVariables", activeDesignVars_)
    )
    {
        // If not, all available design variables will be used. Number is not
        // know at the moment
        Info<< "\t Did not find explicit definition of active design variables."
               " Treating all available ones as active " << endl;
    }

    // read old hessian, correction and derivatives, if present
    readFromDict();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DBFGS::computeCorrection()
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


void Foam::DBFGS::updateOldCorrection(const scalarField& oldCorrection)
{
    updateMethod::updateOldCorrection(oldCorrection);
    correctionOld_ = oldCorrection;
}


void Foam::DBFGS::write()
{
    optMethodIODict_.add<SquareMatrix<scalar>>("HessianOld", HessianOld_, true);
    optMethodIODict_.add<scalarField>("derivativesOld", derivativesOld_, true);
    optMethodIODict_.add<scalarField>("correctionOld", correctionOld_, true);
    optMethodIODict_.add<label>("counter", counter_, true);

    updateMethod::write();
}


// ************************************************************************* //
