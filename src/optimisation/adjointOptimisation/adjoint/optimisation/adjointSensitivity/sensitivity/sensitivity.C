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

#include "runTimeSelectionTables.H"
#include "sensitivity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sensitivity, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sensitivity::writeFaceBasedSens() const
{
    const word suffix(adjointSolverName_ + surfaceFieldSuffix_);

    // Wall face sensitivity projected to normal
    if (wallFaceSensNormalPtr_.valid())
    {
        constructAndWriteSensitivityField<scalar>
        (
            wallFaceSensNormalPtr_,
            "faceSensNormal" + suffix
        );
    }

    if (writeAllSurfaceFiles_)
    {
        // Wall face sensitivity vectors
        if (wallFaceSensVecPtr_.valid())
        {
            constructAndWriteSensitivityField<vector>
            (
                wallFaceSensVecPtr_,
                "faceSensVec" + suffix
            );
        }

        // Normal sens as vectors
        if (wallFaceSensNormalVecPtr_.valid())
        {
            constructAndWriteSensitivityField<vector>
            (
                wallFaceSensNormalVecPtr_,
                "faceSensNormalVec" + suffix
            );
        }
    }
}


void Foam::sensitivity::writePointBasedSens() const
{
    const word suffix(adjointSolverName_ + surfaceFieldSuffix_);

    // Wall point sensitivity projected to normal
    if (wallPointSensNormalPtr_.valid())
    {
        constructAndWriteSensitivtyPointField<scalar>
        (
            wallPointSensNormalPtr_,
            "pointSensNormal" + suffix
        );
    }

    // Write point-based sensitivities, if present
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (writeAllSurfaceFiles_)
    {
        // Wall point sensitivity vectors
        if (wallPointSensVecPtr_.valid())
        {
            constructAndWriteSensitivtyPointField<vector>
            (
                wallPointSensVecPtr_,
                "pointSensVec" + suffix
            );
        }

        // Normal point as vectors
        if (wallPointSensNormalVecPtr_.valid())
        {
            constructAndWriteSensitivtyPointField<vector>
            (
                wallPointSensNormalVecPtr_,
                "pointSensNormalVec" + suffix
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sensitivity::sensitivity
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName
)
:
    mesh_(mesh),
    dict_(dict),
    sensitivityPatchIDs_(0),
    adjointSolverName_(adjointSolverName),
    surfaceFieldSuffix_(word::null),
    writeAllSurfaceFiles_
    (
        dict.lookupOrDefault<bool>
        (
            "writeAllSurfaceFiles",
            false
        )
    ),

    wallFaceSensVecPtr_(nullptr),
    wallFaceSensNormalPtr_(nullptr),
    wallFaceSensNormalVecPtr_(nullptr),

    wallPointSensVecPtr_(nullptr),
    wallPointSensNormalPtr_(nullptr),
    wallPointSensNormalVecPtr_(nullptr),

    fieldSensPtr_(nullptr)
{
    labelHashSet patches
    (
        mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches"))
    );

    if (patches.empty())
    {
        WarningInFunction
            << "There is no patch on which to compute sensitivities. "
            << "Check optimisationDict" << nl
            << endl;
    }
    sensitivityPatchIDs_ = patches.toc();
};


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

const Foam::dictionary& Foam::sensitivity::dict() const
{
    return dict_;
}


bool Foam::sensitivity::readDict(const dictionary& dict)
{
    dict_ = dict;

    return true;
}


const Foam::labelList& Foam::sensitivity::sensitivityPatchIDs() const
{
    return sensitivityPatchIDs_;
}


void Foam::sensitivity::setSensitivityPatchIDs(const labelList& sensPatchIDs)
{
    sensitivityPatchIDs_ = sensPatchIDs;
}


void Foam::sensitivity::computeDerivativesSize()
{
    // Does nothing
}


void Foam::sensitivity::write(const word& baseName)
{
    writeFaceBasedSens();

    writePointBasedSens();

    if (fieldSensPtr_.valid())
    {
        fieldSensPtr_().write();
    }
}


Foam::tmp<Foam::volVectorField> Foam::sensitivity::getWallFaceSensVec()
{
    if (wallFaceSensVecPtr_.valid())
    {
        return
            constructVolSensitivtyField<vector>
            (
                wallFaceSensVecPtr_,
                "faceSensVec" + adjointSolverName_
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
                    mesh_,
                    "faceSensVec" + adjointSolverName_,
                    dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::volScalarField> Foam::sensitivity::getWallFaceSensNormal()
{
    if (wallFaceSensNormalPtr_.valid())
    {
        return
            constructVolSensitivtyField<scalar>
            (
                wallFaceSensNormalPtr_,
                "faceSensNormal" + adjointSolverName_
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
                    mesh_,
                    "faceSensNormal" + adjointSolverName_, dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::volVectorField> Foam::sensitivity::getWallFaceSensNormalVec()
{
    if (wallFaceSensNormalVecPtr_.valid())
    {
        return
            constructVolSensitivtyField<vector>
            (
                wallFaceSensNormalVecPtr_,
                "faceSensNormalVec" + adjointSolverName_
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
                    mesh_,
                    "faceSensNormalVec" + adjointSolverName_,
                    dimless
                ).ptr()
            );
    }
}


Foam::tmp<Foam::pointVectorField> Foam::sensitivity::getWallPointSensVec()
{
    tmp<volVectorField> tWallFaceSensVec = getWallFaceSensVec();
    volPointInterpolation volPointInter(mesh_);

    return (volPointInter.interpolate(tWallFaceSensVec));
}


Foam::tmp<Foam::pointScalarField> Foam::sensitivity::getWallPointSensNormal()
{
    tmp<volScalarField> tWallFaceSensNormal = getWallFaceSensNormal();
    volPointInterpolation volPointInter(mesh_);

    return (volPointInter.interpolate(tWallFaceSensNormal));
}


Foam::tmp<Foam::pointVectorField>
Foam::sensitivity::getWallPointSensNormalVec()
{
    tmp<volVectorField> tWallFaceSensNormalVec = getWallFaceSensNormalVec();
    volPointInterpolation volPointInter(mesh_);

    return (volPointInter.interpolate(tWallFaceSensNormalVec));
}


const Foam::boundaryVectorField&
Foam::sensitivity::getWallFaceSensVecBoundary() const
{
    return wallFaceSensVecPtr_();
}


const Foam::boundaryScalarField&
Foam::sensitivity::getWallFaceSensNormalBoundary() const
{
    return wallFaceSensNormalPtr_();
}


const Foam::boundaryVectorField&
Foam::sensitivity::getWallFaceSensNormalVecBoundary() const
{
    return wallFaceSensNormalVecPtr_();
}


// ************************************************************************* //
