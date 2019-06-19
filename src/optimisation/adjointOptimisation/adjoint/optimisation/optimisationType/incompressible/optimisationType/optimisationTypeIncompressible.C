/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "optimisationTypeIncompressible.H"
#include "constrainedOptimisationMethod.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(optimisationType, 0);
defineRunTimeSelectionTable(optimisationType, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

optimisationType::optimisationType
(
    fvMesh& mesh,
    const dictionary& dict,
    PtrList<adjointSolverManager>& adjointSolverManagers
)
:
    mesh_(mesh),
    dict_(dict),
    adjointSolvManagers_(adjointSolverManagers),
    updateMethod_
    (
        updateMethod::New(mesh_, dict_.subDict("updateMethod"))
    ),
    sourcePtr_(nullptr),
    lineSearch_(lineSearch::New(dict_.subDict("updateMethod"), mesh.time()))
{
    // Figure out number of adjoint solvers corresponding to constraints.
    // Looks in all operating poitns
    label nConstraints(0);
    for (const adjointSolverManager& adjManagerI : adjointSolvManagers_)
    {
        nConstraints += adjManagerI.nConstraints();
    }

    // Sanity checks for combinations of number of constraints and
    // optimisation methods
    if
    (
        nConstraints
     && !isA<constrainedOptimisationMethod>(updateMethod_())
    )
    {
        // Has constraints but is not a constraint optimisation method
        auto cstTablePtr
        (
            constrainedOptimisationMethod::dictionaryConstructorTablePtr_
        );
        FatalErrorInFunction
            << "Found " << nConstraints << " adjoint solvers corresponding to "
            << "constraints but the optimisation method used "
            << "(" << updateMethod_().type() << ") "
            << "is not a constrainedOptimisationMethod. " << nl
            << "Available constrainedOptimisationMethods are :" << nl
            << cstTablePtr->sortedToc()
            << exit(FatalError);
    }
    else if
    (
        !nConstraints
     && isA<constrainedOptimisationMethod>(updateMethod_())
    )
    {
        // Does not have constraints but is a constrained optimisation method
        WarningInFunction
            << "Did not find any adjoint solvers corresponding to "
            << "constraints but the optimisation method used "
            << "(" << updateMethod_().type() << ") "
            << "is a constrainedOptimisationMethod. " << nl << nl
            << "This can cause some constraintOptimisationMethods to misbehave."
            << nl << nl
            << "Either the isConstraint bool is not set in one of the adjoint "
            << "solvers or you should consider using an updateMethod "
            << "that is not a constrainedOptimisationMethod"
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<optimisationType> optimisationType::New
(
    fvMesh& mesh,
    const dictionary& dict,
    PtrList<adjointSolverManager>& adjointSolverManagers
)
{
    const word modelType(dict.subDict("optimisationType").get<word>("type"));

    Info<< "optimisationType type : " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown optimisationType type " << modelType << nl << nl
            << "Valid optimisationType types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<optimisationType>
    (
        cstrIter()(mesh, dict, adjointSolverManagers)
    );
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

tmp<scalarField> optimisationType::computeDirection()
{
    // Sum contributions
    scalarField objectiveSens;
    PtrList<scalarField> constraintSens;
    scalar objectiveValue(Zero);
    scalarField constraintValues;

    for (adjointSolverManager& adjSolvManager : adjointSolvManagers_)
    {
        const scalar opWeight = adjSolvManager.operatingPointWeight();

        // Allocate objective sens size if necessary
        tmp<scalarField> tadjointSolverManagerSens =
            adjSolvManager.aggregateSensitivities();

        if (objectiveSens.empty())
        {
            objectiveSens.setSize(tadjointSolverManagerSens().size(), Zero);
        }

        objectiveSens += opWeight*tadjointSolverManagerSens();
        objectiveValue += opWeight*adjSolvManager.objectiveValue();

        // Allocate constraint sens size if necessary
        PtrList<scalarField> adjointSolverManagerConstSens =
            adjSolvManager.constraintSensitivities();

        tmp<scalarField> cValues = adjSolvManager.constraintValues();

        if (constraintSens.empty())
        {
            constraintSens.setSize(adjointSolverManagerConstSens.size());
            forAll(constraintSens, cI)
            {
                constraintSens.set
                (
                    cI,
                    new scalarField
                    (
                        adjointSolverManagerConstSens[cI].size(),
                        Zero
                    )
                );
                constraintValues.setSize(cValues().size());
                constraintValues = Zero;
            }
        }

        forAll(constraintSens, cI)
        {
            constraintSens[cI] += opWeight*adjointSolverManagerConstSens[cI];
        }
        constraintValues += opWeight*cValues();
    }

    // Based on the sensitivities, return design variables correction
    updateMethod_->setObjectiveDeriv(objectiveSens);
    updateMethod_->setConstraintDeriv(constraintSens);
    updateMethod_->setObjectiveValue(objectiveValue);
    updateMethod_->setConstraintValues(constraintValues);
    tmp<scalarField> tcorrection
    (
        new scalarField(objectiveSens.size(), Zero)
    );
    scalarField& correction = tcorrection.ref();
    correction = updateMethod_->returnCorrection();

    return tcorrection;
}


scalar optimisationType::computeMeritFunction()
{
    // Compute new objective and constraint values and update the ones
    // in updateMethod
    scalar objectiveValue(Zero);
    scalarField constraintValues;

    for (adjointSolverManager& adjSolvManager : adjointSolvManagers_)
    {
        const scalar opWeight = adjSolvManager.operatingPointWeight();

        objectiveValue += opWeight*adjSolvManager.objectiveValue();
        tmp<scalarField> cValues = adjSolvManager.constraintValues();

        if (constraintValues.empty())
        {
            constraintValues.setSize(cValues().size(), Zero);
        }
        constraintValues += opWeight*cValues();
    }
    updateMethod_->setObjectiveValue(objectiveValue);
    updateMethod_->setConstraintValues(constraintValues);

    return updateMethod_->computeMeritFunction();
}


scalar optimisationType::meritFunctionDirectionalDerivative()
{
    return updateMethod_->meritFunctionDirectionalDerivative();
}


void optimisationType::updateOldCorrection(const scalarField& oldCorrection)
{
    updateMethod_->updateOldCorrection(oldCorrection);
}


void optimisationType::write()
{
    updateMethod_->write();
}


const autoPtr<volScalarField>& optimisationType::sourcePtr()
{
    return sourcePtr_;
}


autoPtr<lineSearch>& optimisationType::getLineSearch()
{
    return lineSearch_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
