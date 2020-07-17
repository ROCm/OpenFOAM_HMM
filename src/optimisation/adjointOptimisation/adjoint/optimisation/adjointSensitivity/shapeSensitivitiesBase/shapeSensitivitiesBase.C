/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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

#include "shapeSensitivitiesBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shapeSensitivitiesBase, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::shapeSensitivitiesBase::writeFaceBasedSens() const
{
    // Wall face sensitivity projected to normal
    if (wallFaceSensNormalPtr_)
    {
        constructAndWriteSensitivityField<scalar>
        (
            wallFaceSensNormalPtr_,
            "faceSensNormal" + surfaceFieldSuffix_
        );
    }

    if (writeAllSurfaceFiles_)
    {
        // Wall face sensitivity vectors
        if (wallFaceSensVecPtr_)
        {
            constructAndWriteSensitivityField<vector>
            (
                wallFaceSensVecPtr_,
                "faceSensVec" + surfaceFieldSuffix_
            );
        }

        // Normal sens as vectors
        if (wallFaceSensNormalVecPtr_)
        {
            constructAndWriteSensitivityField<vector>
            (
                wallFaceSensNormalVecPtr_,
                "faceSensNormalVec" + surfaceFieldSuffix_
            );
        }
    }
}


void Foam::shapeSensitivitiesBase::writePointBasedSens() const
{
    // Wall point sensitivity projected to normal
    if (wallPointSensNormalPtr_)
    {
        constructAndWriteSensitivtyPointField<scalar>
        (
            wallPointSensNormalPtr_,
            "pointSensNormal" + surfaceFieldSuffix_
        );
    }

    // Write point-based sensitivities, if present
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (writeAllSurfaceFiles_)
    {
        // Wall point sensitivity vectors
        if (wallPointSensVecPtr_)
        {
            constructAndWriteSensitivtyPointField<vector>
            (
                wallPointSensVecPtr_,
                "pointSensVec" + surfaceFieldSuffix_
            );
        }

        // Normal point as vectors
        if (wallPointSensNormalVecPtr_)
        {
            constructAndWriteSensitivtyPointField<vector>
            (
                wallPointSensNormalVecPtr_,
                "pointSensNormalVec" + surfaceFieldSuffix_
            );
        }
    }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shapeSensitivitiesBase::shapeSensitivitiesBase
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    meshShape_(mesh),
    surfaceFieldSuffix_(word::null),
    writeAllSurfaceFiles_
    (
        dict.getOrDefault<bool>
        (
            "writeAllSurfaceFiles",
            false
        )
    ),
    sensitivityPatchIDs_
    (
        mesh.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches", keyType::REGEX_RECURSIVE)
        )
    ),
    wallFaceSensVecPtr_(nullptr),
    wallFaceSensNormalPtr_(nullptr),
    wallFaceSensNormalVecPtr_(nullptr),

    wallPointSensVecPtr_(nullptr),
    wallPointSensNormalPtr_(nullptr),
    wallPointSensNormalVecPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelHashSet&
Foam::shapeSensitivitiesBase::sensitivityPatchIDs() const
{
    return sensitivityPatchIDs_;
}


void Foam::shapeSensitivitiesBase::setSensitivityPatchIDs
(
    const labelHashSet& sensPatchIDs
)
{
    sensitivityPatchIDs_ = sensPatchIDs;
}


void Foam::shapeSensitivitiesBase::clearSensitivities()
{
    // Face-based boundary sens
    if (wallFaceSensVecPtr_)
    {
        wallFaceSensVecPtr_() = vector::zero;
    }
    if (wallFaceSensNormalVecPtr_)
    {
        wallFaceSensNormalVecPtr_() = vector::zero;
    }
    if (wallFaceSensNormalPtr_)
    {
        wallFaceSensNormalPtr_() = scalar(0);
    }

    // Point-based boundary sens
    if (wallPointSensVecPtr_)
    {
        for (vectorField& patchSens : wallPointSensVecPtr_())
        {
            patchSens = vector::zero;
        }
    }
    if (wallPointSensNormalVecPtr_)
    {
        for (vectorField& patchSens : wallPointSensNormalVecPtr_())
        {
            patchSens = vector::zero;
        }
    }
    if (wallPointSensNormalPtr_)
    {
        for (scalarField& patchSens : wallPointSensNormalPtr_())
        {
            patchSens = scalar(0);
        }
    }
}


void Foam::shapeSensitivitiesBase::write()
{
    writeFaceBasedSens();
    writePointBasedSens();
}


void Foam::shapeSensitivitiesBase::setSuffix(const word& suffix)
{
    surfaceFieldSuffix_ = suffix;
}


Foam::tmp<Foam::volVectorField>
Foam::shapeSensitivitiesBase::getWallFaceSensVec()
{
    if (wallFaceSensVecPtr_)
    {
        return
            constructVolSensitivtyField<vector>
            (
                wallFaceSensVecPtr_,
                "faceSensVec" + surfaceFieldSuffix_
            );
    }
    else
    {
        WarningInFunction
            << " no faceSensVec boundary field. Returning zero" << endl;

        return
            tmp<volVectorField>
            (
                createZeroFieldPtr<vector>
                (
                    meshShape_,
                    "faceSensVec" + surfaceFieldSuffix_,
                    dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::shapeSensitivitiesBase::getWallFaceSensNormal()
{
    if (wallFaceSensNormalPtr_)
    {
        return
            constructVolSensitivtyField<scalar>
            (
                wallFaceSensNormalPtr_,
                "faceSensNormal" + surfaceFieldSuffix_
            );
    }
    else
    {
        WarningInFunction
            << " no wallFaceSensNormal boundary field. Returning zero" << endl;

        return
            tmp<volScalarField>
            (
                createZeroFieldPtr<scalar>
                (
                    meshShape_,
                    "faceSensNormal" + surfaceFieldSuffix_, dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::volVectorField>
Foam::shapeSensitivitiesBase::getWallFaceSensNormalVec()
{
    if (wallFaceSensNormalVecPtr_)
    {
        return
            constructVolSensitivtyField<vector>
            (
                wallFaceSensNormalVecPtr_,
                "faceSensNormalVec" + surfaceFieldSuffix_
            );
    }
    else
    {
        WarningInFunction
            << " no wallFaceSensNormalVec boundary field. Returning zero"
            << endl;

        return
            tmp<volVectorField>
            (
                createZeroFieldPtr<vector>
                (
                    meshShape_,
                    "faceSensNormalVec" + surfaceFieldSuffix_,
                    dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::pointVectorField>
Foam::shapeSensitivitiesBase::getWallPointSensVec()
{
    tmp<volVectorField> tWallFaceSensVec = getWallFaceSensVec();
    volPointInterpolation volPointInter(meshShape_);

    return (volPointInter.interpolate(tWallFaceSensVec));
}


Foam::tmp<Foam::pointScalarField>
Foam::shapeSensitivitiesBase::getWallPointSensNormal()
{
    tmp<volScalarField> tWallFaceSensNormal = getWallFaceSensNormal();
    volPointInterpolation volPointInter(meshShape_);

    return (volPointInter.interpolate(tWallFaceSensNormal));
}


Foam::tmp<Foam::pointVectorField>
Foam::shapeSensitivitiesBase::getWallPointSensNormalVec()
{
    tmp<volVectorField> tWallFaceSensNormalVec = getWallFaceSensNormalVec();
    volPointInterpolation volPointInter(meshShape_);

    return (volPointInter.interpolate(tWallFaceSensNormalVec));
}


const Foam::boundaryVectorField&
Foam::shapeSensitivitiesBase::getWallFaceSensVecBoundary() const
{
    return wallFaceSensVecPtr_();
}


const Foam::boundaryScalarField&
Foam::shapeSensitivitiesBase::getWallFaceSensNormalBoundary() const
{
    return wallFaceSensNormalPtr_();
}


const Foam::boundaryVectorField&
Foam::shapeSensitivitiesBase::getWallFaceSensNormalVecBoundary() const
{
    return wallFaceSensNormalVecPtr_();
}


// ************************************************************************* //
