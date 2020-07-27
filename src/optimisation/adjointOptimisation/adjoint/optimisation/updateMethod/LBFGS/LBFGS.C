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

#include "LBFGS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LBFGS, 0);
    addToRunTimeSelectionTable
    (
        updateMethod,
        LBFGS,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::LBFGS::allocateMatrices()
{
    // Set active design variables, if necessary
    if (activeDesignVars_.empty())
    {
        activeDesignVars_ = identity(objectiveDerivatives_.size());
    }

    // Allocate vectors
    label nVars(activeDesignVars_.size());
    for (label i = 0; i < nPrevSteps_; i++)
    {
        y_.set(i, new scalarField(nVars, Zero));
        s_.set(i, new scalarField(nVars, Zero));
    }
}


void Foam::LBFGS::pivotFields(PtrList<scalarField>& list, const scalarField& f)
{
    if (counter_ > nPrevSteps_)
    {
        // Reorder list by moving pointers down the line
        labelList newOrder(nPrevSteps_, -1);
        newOrder[0] = nPrevSteps_ - 1;
        for (label i = 1; i < nPrevSteps_; ++i)
        {
            newOrder[i] = i - 1;
        }
        list.reorder(newOrder);

        // Fill in last element with the provided field
        list[nPrevSteps_ - 1] = f;
    }
    else
    {
        list[counter_ - 1] = f;
    }
}


void Foam::LBFGS::updateVectors()
{
    // Update list of y. Can only be done here since objectiveDerivatives_
    // was not known at the end of the previous loop
    scalarField yRecent
        (objectiveDerivatives_ - derivativesOld_, activeDesignVars_);
    pivotFields(y_, yRecent);
    // Update list of s.
    // correction_ holds the previous correction
    scalarField sActive(correctionOld_, activeDesignVars_);
    pivotFields(s_, sActive);

    DebugInfo
        << "y fields" << nl << y_ << endl;
    DebugInfo
        << "s fields" << nl << s_ << endl;
}


void Foam::LBFGS::steepestDescentUpdate()
{
    Info<< "Using steepest descent to update design variables" << endl;
    correction_ = -eta_*objectiveDerivatives_;
}


void Foam::LBFGS::LBFGSUpdate()
{
    // L-BFGS two loop recursion
    //~~~~~~~~~~~~~~~~~~~~~~~~~~
    label nSteps(min(counter_, nPrevSteps_));
    label nLast(nSteps - 1);
    scalarField q(objectiveDerivatives_, activeDesignVars_);
    scalarField a(nSteps, Zero);
    scalarField r(nSteps, Zero);
    for (label i = nLast; i > -1; --i)
    {
        r[i] = 1./globalSum(y_[i]*s_[i]);
        a[i] = r[i]*globalSum(s_[i]*q);
        q -= a[i]*y_[i];
    }

    scalar gamma =
        globalSum(y_[nLast]*s_[nLast])/globalSum(y_[nLast]*y_[nLast]);
    q *= gamma;

    scalarField b(activeDesignVars_.size(), Zero);
    for (label i = 0; i < nSteps; ++i)
    {
        b = r[i]*globalSum(y_[i]*q);
        q += s_[i]*(a[i] -b);
    }

    // Update correction
    forAll(activeDesignVars_, varI)
    {
        correction_[activeDesignVars_[varI]] = -etaHessian_*q[varI];
    }
}


void Foam::LBFGS::update()
{
    // In the first few iterations, use steepest descent but update the Hessian
    // matrix
    if (counter_ < nSteepestDescent_)
    {
        steepestDescentUpdate();
    }
    // else use LBFGS formula to update the design variables
    else
    {
        LBFGSUpdate();
    }

    // Store fields for the next iteration
    derivativesOld_ = objectiveDerivatives_;
    correctionOld_ = correction_;
}


void Foam::LBFGS::readFromDict()
{
    if (optMethodIODict_.headerOk())
    {
        optMethodIODict_.readEntry("y", y_);
        optMethodIODict_.readEntry("s", s_);
        optMethodIODict_.readEntry("derivativesOld", derivativesOld_);
        optMethodIODict_.readEntry("counter", counter_);
        optMethodIODict_.readEntry("eta", eta_);
        optMethodIODict_.readEntry("correctionOld", correctionOld_);

        correction_ = scalarField(correctionOld_.size(), Zero);

        if (activeDesignVars_.empty())
        {
            activeDesignVars_ = identity(derivativesOld_.size());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LBFGS::LBFGS
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
    nPrevSteps_
    (
        coeffsDict().getOrDefault<label>("nPrevSteps", 10)
    ),
    y_(nPrevSteps_),
    s_(nPrevSteps_),
    derivativesOld_(0),
    counter_(Zero)
{
    if
    (
        !coeffsDict().readIfPresent("activeDesignVariables", activeDesignVars_)
    )
    {
        // If not, all available design variables will be used. Number is not
        // know at the moment
        Info<< "\t Did not find explicit definition of active design variables. "
            << "Treating all available ones as active " << endl;
    }

    // Read old Hessian, correction and derivatives, if present
    readFromDict();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LBFGS::computeCorrection()
{
    if (counter_ == 0)
    {
        allocateMatrices();
    }
    else
    {
        updateVectors();
    }

    update();
    ++counter_;
}


void Foam::LBFGS::updateOldCorrection(const scalarField& oldCorrection)
{
    updateMethod::updateOldCorrection(oldCorrection);
    correctionOld_ = oldCorrection;
}


void Foam::LBFGS::write()
{
    optMethodIODict_.add<PtrList<scalarField>>("y", y_, true);
    optMethodIODict_.add<PtrList<scalarField>>("s", s_, true);
    optMethodIODict_.add<scalarField>("derivativesOld", derivativesOld_, true);
    optMethodIODict_.add<scalarField>("correctionOld", correctionOld_, true);
    optMethodIODict_.add<label>("counter", counter_, true);

    updateMethod::write();
}


// ************************************************************************* //
