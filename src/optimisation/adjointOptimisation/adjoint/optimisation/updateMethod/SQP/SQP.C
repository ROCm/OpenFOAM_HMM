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

#include "SQP.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SQP, 1);
    addToRunTimeSelectionTable
    (
        updateMethod,
        SQP,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        constrainedOptimisationMethod,
        SQP,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::SQP::allocateMatrices()
{
    // Set active design variables, if necessary
    if (activeDesignVars_.empty())
    {
        activeDesignVars_ = identity(objectiveDerivatives_.size());
    }

    // Set previous Hessian to be a diagonal matrix
    SquareMatrix<scalar> temp(activeDesignVars_.size(), I);

    // Allocate correct size and content to Hessian matrices
    // Has a max. capability of approximately 34000 variables.
    HessianOld_ = temp;
    Hessian_ = temp;

    // Set size of Lagrange multipliers
    lamdas_.setSize(constraintDerivatives_.size());
    lamdas_ = Zero;

    // Set corerction size
    correction_.setSize(objectiveDerivatives_.size());
    correction_ = Zero;
}


void Foam::SQP::updateHessian()
{
    // Vectors needed to construct the (inverse) Hessian matrix
    scalarField y(activeDesignVars_.size(), Zero);
    scalarField s(activeDesignVars_.size(), Zero);
    scalarField LagrangianDerivativesOld = objectiveDerivativesOld_;
    forAll(constraintDerivatives_, cI)
    {
        LagrangianDerivatives_ -= lamdas_[cI] * constraintDerivatives_[cI];
        LagrangianDerivativesOld -= lamdas_[cI] * constraintDerivativesOld_[cI];
    }
    y.map(LagrangianDerivatives_ - LagrangianDerivativesOld, activeDesignVars_);
    s.map(correctionOld_, activeDesignVars_);

    scalar ys = globalSum(s*y);
    if (counter_ == 1 && scaleFirstHessian_)
    {
        if (ys > scalar(0))
        {
            scalar scaleFactor = ys/globalSum(y*y);
            Info<< "Scaling Hessian with factor " << scaleFactor << endl;
            forAll(activeDesignVars_, varI)
            {
                HessianOld_[varI][varI] /= scaleFactor;
            }
        }
        else
        {
            WarningInFunction
                << " y*s is negative. Skipping the scaling of the first Hessian"
                << endl;
        }
    }
    scalar sBs = globalSum(leftMult(s, HessianOld_)*s);

    // Check curvature condition
    scalar theta(1);
    if (ys < dumpingThreshold_*sBs)
    {
        WarningInFunction
            << " y*s is below threshold. Using damped form" << endl;
        theta = (1 - dumpingThreshold_)*sBs/(sBs - ys);
    }
    scalarField r(theta*y + (scalar(1) - theta)*rightMult(HessianOld_, s));
    DebugInfo
        << "Unmodified Hessian curvature index " << ys << endl;
    DebugInfo
        << "Modified Hessian curvature index " << globalSum(r*s) << endl;

    // Update the Hessian
    Hessian_ =
        HessianOld_
      - outerProd(rightMult(HessianOld_, s), leftMult(s/sBs, HessianOld_))
      + outerProd(r, r/globalSum(s*r));
}


void Foam::SQP::computeLagrangeMultipliersAndCorrect()
{
    SquareMatrix<scalar> HessianInv = inv(Hessian_);  //also denoted below as W
    if (debug > 1)
    {
        Info<< "Hessian " << Hessian_ << endl;
        Info<< "HessianInv " << HessianInv << endl;
        label n = Hessian_.n();
        SquareMatrix<scalar> test(n, Zero);
        for (label k = 0; k < n; k++)
        {
            for (label l = 0; l < n; l++)
            {
                scalar elem(Zero);
                for (label i = 0; i < n; i++)
                {
                    elem += Hessian_[k][i] * HessianInv[i][l];
                }
                test[k][l]=elem;
            }
        }
        Info<< "Validation " << test << endl;
    }

    // Compute new Lagrange multipliers
    label nc = constraintDerivatives_.size();
    scalarField activeDerivs(activeDesignVars_.size(), Zero);

    // activeDerivs.map(objectiveDerivatives_, activeDesignVars_);
    activeDerivs.map(LagrangianDerivatives_, activeDesignVars_);
    scalarField WgradL = rightMult(HessianInv, activeDerivs);

    scalarField lamdaRHS(nc, Zero);
    forAll(lamdaRHS, cI)
    {
        scalarField activeConsDerivs(activeDesignVars_.size(), Zero);
        activeConsDerivs.map(constraintDerivatives_[cI], activeDesignVars_);
        lamdaRHS[cI] = globalSum(activeConsDerivs * WgradL) - cValues_[cI];
        if (debug > 1)
        {
            Info<< "lamdaRHS total|deriv part|constraint part "
                << lamdaRHS[cI] << " " << globalSum(activeConsDerivs * WgradL)
                << " " << cValues_[cI] << endl;
        }
    }

    // lhs for the lamda system
    SquareMatrix<scalar> AWA(nc, Zero);
    PtrList<scalarField> WA(nc);
    for (label j = 0; j < nc; j++)
    {
        scalarField gradcJ(activeDesignVars_.size(), Zero);
        gradcJ.map(constraintDerivatives_[j], activeDesignVars_);
        WA.set(j, new scalarField(rightMult(HessianInv, gradcJ)));
        for (label i = 0; i < nc; i++)
        {
            scalarField gradcI(activeDesignVars_.size(), Zero);
            gradcI.map(constraintDerivatives_[i], activeDesignVars_);
            AWA[i][j] = globalSum(gradcI * WA[j]);
        }
    }
    SquareMatrix<scalar> invAWA = inv(AWA);
    scalarField deltaLamda = rightMult(invAWA, lamdaRHS);
    if (debug > 1)
    {
        Info<< "AWA " << AWA << endl;
        Info<< "AWAInv " << invAWA << endl;
        Info<< "lamda update " << deltaLamda << endl;
    }
    lamdas_ += deltaLamda;

    // Compute design variables correction
    scalarField activeCorrection(-WgradL);
    forAll(WA, cI)
    {
        activeCorrection += WA[cI]*deltaLamda[cI];
    }
    activeCorrection *= etaHessian_;
    // Transfer correction to the global list
    correction_ = Zero;
    forAll(activeDesignVars_, varI)
    {
        correction_[activeDesignVars_[varI]] = activeCorrection[varI];
    }
    if (counter_ == 0)
    {
        correction_ *= eta_;
    }
}


void Foam::SQP::storeOldFields()
{
    objectiveDerivativesOld_ = objectiveDerivatives_;
    if (constraintDerivativesOld_.empty())
    {
        constraintDerivativesOld_.setSize(constraintDerivatives_.size());
    }
    forAll(constraintDerivativesOld_, cI)
    {
        constraintDerivativesOld_[cI] = constraintDerivatives_[cI];
    }
    correctionOld_ = correction_;
    HessianOld_ = Hessian_;
}


void Foam::SQP::readFromDict()
{
    if (optMethodIODict_.headerOk())
    {
        optMethodIODict_.readEntry("Hessian", Hessian_);
        optMethodIODict_.readEntry("HessianOld", HessianOld_);
        optMethodIODict_.readEntry
        (
            "objectiveDerivativesOld",
            objectiveDerivativesOld_
        );
        optMethodIODict_.readEntry
        (
            "constraintDerivativesOld",
            constraintDerivativesOld_
        );
        optMethodIODict_.readEntry("correctionOld", correctionOld_);
        optMethodIODict_.readEntry("lamdas", lamdas_);
        optMethodIODict_.readEntry("counter", counter_);
        optMethodIODict_.readEntry("eta", eta_);

        correction_ = scalarField(correctionOld_.size(), Zero);

        if (activeDesignVars_.empty())
        {
            activeDesignVars_ = identity(correction_.size());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SQP::SQP(const fvMesh& mesh, const dictionary& dict)
:
    constrainedOptimisationMethod(mesh, dict),

    etaHessian_
    (
        coeffsDict().getOrDefault<scalar>("etaHessian", 1)
    ),
    activeDesignVars_(0),
    scaleFirstHessian_
    (
        coeffsDict().getOrDefault<bool>("scaleFirstHessian", false)
    ),
    dumpingThreshold_
    (
        coeffsDict().getOrDefault<scalar>("dumpingThreshold", 0.2)
    ),
    LagrangianDerivatives_(0),
    Hessian_(),  // construct null matrix since we dont know the dimension yet
    HessianOld_(),
    objectiveDerivativesOld_(0),
    constraintDerivativesOld_(0),
    correctionOld_(0),
    lamdas_(0),
    counter_(0),
    objFunctionFolder_
    (
        mesh_.time().globalPath()/"optimisation"/"objective"/
        mesh_.time().timeName()
    ),
    meritFunctionFile_(nullptr),
    mu_(Zero),
    delta_
    (
        coeffsDict().getOrDefault<scalar>("delta", 0.1)
    )
{
    if
    (
        !coeffsDict().readIfPresent("activeDesignVariables", activeDesignVars_)
    )
    {
        // If not, all available design variables will be used. Number is not
        // know at the moment
        Info<< "\t Did not find explicit definition of active design "
            << "variables. Treating all available ones as active " << endl;
    }

    // Create folder to merit function
    if (Pstream::master())
    {
        mkDir(objFunctionFolder_);
    }

    // Read old hessian, correction and derivatives, if present
    readFromDict();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SQP::computeCorrection()
{
    // Allocate correct sizes in first update
    if (counter_ == 0)
    {
        allocateMatrices();
    }

    // The first iteration uses a unitary Hessian. No need to update
    LagrangianDerivatives_ = objectiveDerivatives_;
    if (counter_ != 0)
    {
        updateHessian();
    }

    // Update lamdas and desing vars
    computeLagrangeMultipliersAndCorrect();

    // Store fields for the next iteration and write them to file
    storeOldFields();

    counter_++;
}


Foam::scalar Foam::SQP::computeMeritFunction()
{
    // If condition is not met, update mu value
    if (mu_ < max(mag(lamdas_)) + delta_)
    {
        mu_ = max(mag(lamdas_)) + 2*delta_;
        if (debug > 1)
        {
            Info<< "Updated mu value to " << mu_ << endl;
        }
    }
    scalar L = objectiveValue_ + mu_*sum(mag(cValues_));

    return L;
}


Foam::scalar Foam::SQP::meritFunctionDirectionalDerivative()
{
    scalar deriv =
        globalSum(objectiveDerivatives_*correction_)
      - mu_*sum(mag(cValues_));

    return deriv;
}


void Foam::SQP::updateOldCorrection(const scalarField& oldCorrection)
{
    updateMethod::updateOldCorrection(oldCorrection);
    correctionOld_ = oldCorrection;
}


void Foam::SQP::write()
{
    // Write updateMethod dictionary
    optMethodIODict_.add<SquareMatrix<scalar>>("Hessian", Hessian_, true);
    optMethodIODict_.add<SquareMatrix<scalar>>("HessianOld", HessianOld_, true);
    optMethodIODict_.
        add<scalarField>
        (
            "objectiveDerivativesOld", objectiveDerivativesOld_, true
        );
    optMethodIODict_.
        add<List<scalarField>>
        (
            "constraintDerivativesOld", constraintDerivativesOld_, true
        );
    optMethodIODict_.add<scalarField>("correctionOld", correctionOld_, true);
    optMethodIODict_.add<scalarField>("lamdas", lamdas_, true);
    optMethodIODict_.add<label>("counter", counter_, true);

    updateMethod::write();

    // Write merit function
    scalar constraintPart = sum(mag(cValues_));
    scalar merit = objectiveValue_ + mu_*constraintPart;
    if (Pstream::master())
    {
        unsigned int width = IOstream::defaultPrecision() + 6;
        unsigned int constraintsSize = lamdas_.size();
        constraintsSize = constraintsSize*(width + 1) + 2;

        // Open file and write header
        if (!meritFunctionFile_)
        {
            meritFunctionFile_.reset
            (
                new OFstream(objFunctionFolder_/word("meritFunction"))
            );

            meritFunctionFile_()
                << setw(1) << "#" << " "
                << setw(width) << "merit" << " "
                << setw(width) << "J" << " "
                << setw(constraintsSize) << "lamdas" << " "
                << setw(constraintsSize) << "constraints" << " "
                << setw(width) << "mu" << " "
                << setw(width) << "constraintContr" << endl;

        }

        meritFunctionFile_()
            << setw(1) << mesh_.time().value() -1 << " "
            << setw(width) << merit << " "
            << setw(width) << objectiveValue_ << " "
            << setw(1) << "(";

        forAll(lamdas_, cI)
        {
            meritFunctionFile_()
                << setw(width) << lamdas_[cI] << setw(1) << " ";
        }
        meritFunctionFile_() << setw(3) << ")(";
        forAll(cValues_, cI)
        {
            meritFunctionFile_()
                << setw(width) << cValues_[cI] << setw(1) << " ";
        }
        meritFunctionFile_() << setw(2) << ") ";
        meritFunctionFile_() << setw(width) << mu_ << " ";
        meritFunctionFile_() << setw(width) << constraintPart << endl;
    }
}


// ************************************************************************* //
