/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "optMeshMovement.H"
#include "cellQuality.H"
#include "createZeroField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(optMeshMovement, 0);
    defineRunTimeSelectionTable(optMeshMovement, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::optMeshMovement::getMaxAllowedDisplacement() const
{
    if (!maxAllowedDisplacement_)
    {
        FatalErrorInFunction
            << "maxAllowedDisplacement requested but not set" << nl
            << exit(FatalError);
    }

    return maxAllowedDisplacement_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optMeshMovement::optMeshMovement
(
    fvMesh& mesh,
    const dictionary& dict,
    const labelList& patchIDs
)
:
    maxAllowedDisplacement_(nullptr),
    mesh_(mesh),
    dict_(dict),
    correction_(0),
    patchIDs_(patchIDs),
    pointsInit_(mesh.points()),
    displMethodPtr_(displacementMethod::New(mesh_, patchIDs_)),
    writeMeshQualityMetrics_
    (
        dict.getOrDefault("writeMeshQualityMetrics", false)
    )
{
    // Set maxAllowedDisplacement if provided
    if (dict.found("maxAllowedDisplacement"))
    {
        maxAllowedDisplacement_.reset
        (
            new scalar(dict.get<scalar>("maxAllowedDisplacement"))
        );
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::optMeshMovement> Foam::optMeshMovement::New
(
    fvMesh& mesh,
    const dictionary& dict,
    const labelList& patchIDs
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "optMeshMovement type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "type",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<optMeshMovement>(ctorPtr(mesh, dict, patchIDs));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::optMeshMovement::setCorrection(const scalarField& correction)
{
    correction_ = correction;
}


void Foam::optMeshMovement::moveMesh()
{
    // Move mesh
    displMethodPtr_->update();

    // Check mesh quality
    mesh_.checkMesh(true);

    // If needed, plot mesh quality metrics
    writeMeshQualityMetrics();
}


Foam::autoPtr<Foam::displacementMethod>&
Foam::optMeshMovement::returnDisplacementMethod()
{
    return displMethodPtr_;
}


const Foam::labelList& Foam::optMeshMovement::getPatchIDs()
{
    return patchIDs_;
}


void Foam::optMeshMovement::writeMeshQualityMetrics()
{
    if (writeMeshQualityMetrics_)
    {
        cellQuality cellQualityEngine(mesh_);
        tmp<scalarField> cellNonOrtho(cellQualityEngine.nonOrthogonality());
        tmp<scalarField> cellSkewness(cellQualityEngine.skewness());
        Info<< "Average, Max cell non - orthogonality "
            << gAverage(cellNonOrtho())
            << " " << gMax(cellNonOrtho()) << endl;
        Info<< "Average, Max cell skewness " << gAverage(cellSkewness())
            << " " << gMax(cellSkewness()) << endl;
        autoPtr<volScalarField> nonOrthoPtr
        (
           createZeroFieldPtr<scalar>(mesh_, "nonOrtho", dimless)
        );
        autoPtr<volScalarField> skewnessPtr
        (
           createZeroFieldPtr<scalar>(mesh_, "skewness", dimless)
        );
        nonOrthoPtr().primitiveFieldRef() = cellNonOrtho();
        skewnessPtr().primitiveFieldRef() = cellSkewness();
        nonOrthoPtr().write();
        skewnessPtr().write();
    }
}


void Foam::optMeshMovement::storeDesignVariables()
{
    pointsInit_ = mesh_.points();
}


void Foam::optMeshMovement::resetDesignVariables()
{
    Info<< "optMeshMovement:: resetting mesh points" << endl;
    mesh_.movePoints(pointsInit_);
}


bool Foam::optMeshMovement::maxAllowedDisplacementSet() const
{
    return maxAllowedDisplacement_.valid();
}


Foam::labelList Foam::optMeshMovement::getActiveDesignVariables() const
{
    NotImplemented;
    return labelList(0);
}


// ************************************************************************* //
