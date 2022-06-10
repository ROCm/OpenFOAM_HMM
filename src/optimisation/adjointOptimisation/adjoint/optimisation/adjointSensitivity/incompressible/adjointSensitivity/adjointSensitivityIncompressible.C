/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
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

#include "runTimeSelectionTables.H"
#include "adjointSensitivityIncompressible.H"
#include "boundaryAdjointContribution.H"
#include "incompressibleAdjointSolver.H"
#include "wallFvPatch.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointSensitivity, 0);
defineRunTimeSelectionTable(adjointSensitivity, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointSensitivity::adjointSensitivity
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleAdjointSolver& adjointSolver
)
:
    sensitivity(mesh, dict),
    derivatives_(0),
    adjointSolver_(adjointSolver),
    primalVars_(adjointSolver.getPrimalVars()),
    adjointVars_(adjointSolver.getAdjointVars()),
    objectiveManager_(adjointSolver.getObjectiveManager())
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<adjointSensitivity> adjointSensitivity::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleAdjointSolver& adjointSolver
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "adjointSensitivity type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "adjointSensitivity",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<adjointSensitivity>
    (
        ctorPtr(mesh, dict, adjointSolver)
    );
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

const scalarField& adjointSensitivity::calculateSensitivities()
{
    assembleSensitivities();
    write(type());
    return derivatives_;
}


const scalarField& adjointSensitivity::getSensitivities() const
{
    return derivatives_;
}


void adjointSensitivity::clearSensitivities()
{
    derivatives_ = scalar(0);
    if (fieldSensPtr_)
    {
        fieldSensPtr_().primitiveFieldRef() = scalar(0);
    }
}


void adjointSensitivity::write(const word& baseName)
{
    sensitivity::write(baseName);
}


tmp<volTensorField> adjointSensitivity::computeGradDxDbMultiplier()
{
    return adjointSolver_.computeGradDxDbMultiplier();
}


tmp<volVectorField> adjointSensitivity::adjointMeshMovementSource()
{
    tmp<volTensorField> tgradDxDbMult = computeGradDxDbMultiplier();
    volTensorField& gradDxDbMult = tgradDxDbMult.ref();

    tmp<volVectorField> tadjointMeshMovementSource
    (
        new volVectorField
        (
            IOobject
            (
               "adjointMeshMovementSource",
               mesh_.time().timeName(),
               mesh_,
               IOobject::NO_READ,
               IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(gradDxDbMult.dimensions()/dimLength, Zero)
        )
    );

    volVectorField& source = tadjointMeshMovementSource.ref();

    source -= fvc::div(gradDxDbMult.T());

    // Terms from fvOptions
    fv::options::New(this->mesh_).postProcessSens
    (
        source.primitiveFieldRef(), adjointVars_.solverName()
    );

    return (tadjointMeshMovementSource);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
