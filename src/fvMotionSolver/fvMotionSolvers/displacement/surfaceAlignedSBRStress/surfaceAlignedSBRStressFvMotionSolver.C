/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "surfaceAlignedSBRStressFvMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "pointIndexHit.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "unitConversion.H"
#include "motionDiffusivity.H"
#include "fvcSmooth.H"
#include "pointMVCWeight.H"
#include "dimensionedSymmTensor.H"
#include "quaternion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceAlignedSBRStressFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        surfaceAlignedSBRStressFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceAlignedSBRStressFvMotionSolver::
surfaceAlignedSBRStressFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementSBRStressFvMotionSolver(mesh, dict),
    surfaceNames_(coeffDict().lookup("surfaces")),
    surfaceMesh_(surfaceNames_.size()),
    cellRot_
    (
        IOobject
        (
            "cellRot",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvMesh_,
        dimensionedVector("zero", dimless, vector::zero)
    ),
    maxAng_(coeffDict().lookupOrDefault<scalar>("maxAng", 80.0)),
    minAng_(coeffDict().lookupOrDefault<scalar>("minAng", 20.0)),
    accFactor_(coeffDict().lookupOrDefault<scalar>("accFactor", 1.0)),
    smoothFactor_(coeffDict().lookupOrDefault<scalar>("smoothFactor", 0.9)),
    nNonOrthogonalCorr_(readLabel(coeffDict().lookup("nNonOrthogonalCorr"))),
    pointDisplacement_(pointDisplacement()),
    sigmaD_
    (
        IOobject
        (
            "sigmaD",
            fvMesh_.time().timeName(),
            fvMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvMesh_,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    minSigmaDiff_(coeffDict().lookupOrDefault<scalar>("minSigmaDiff", 1e-4))
{
    forAll (surfaceNames_, i)
    {
        surfaceMesh_.set
        (
            i,
            new triSurfaceMesh
            (
                IOobject
                (
                    surfaceNames_[i],
                    mesh.time().constant(),
                    "triSurface",
                    mesh.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceAlignedSBRStressFvMotionSolver::
~surfaceAlignedSBRStressFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::surfaceAlignedSBRStressFvMotionSolver::calculateCellRot()
{
    cellRot_.internalField() = vector::zero;
    pointDisplacement_.internalField() = vector::zero;

    // Find intersections
    pointField start(fvMesh_.nInternalFaces());
    pointField end(start.size());

    for (label faceI = 0; faceI < fvMesh_.nInternalFaces(); faceI++)
    {
        start[faceI] = fvMesh_.cellCentres()[fvMesh_.faceOwner()[faceI]];
        end[faceI] = fvMesh_.cellCentres()[fvMesh_.faceNeighbour()[faceI]];
    }

    DynamicList<label> hitCells;

    forAll (surfaceMesh_, surfI)
    {
        List<pointIndexHit> hit(start.size());
        surfaceMesh_[surfI].findLineAny(start, end, hit);

        labelField pointsCount(fvMesh_.nPoints(), 1);

        const vectorField& nf = surfaceMesh_[surfI].faceNormals();

        const vectorField& SfMesh = fvMesh_.faceAreas();

        const vectorField nSfMesh(SfMesh/mag(SfMesh));

        DynamicList<label> cellsHit;

        forAll(hit, faceI)
        {
            if (hit[faceI].hit())
            {
                label rotCellId(-1);
                const vector hitPoint = hit[faceI].hitPoint();

                if (fvMesh_.isInternalFace(faceI))
                {
                    const vector cCOne =
                        fvMesh_.cellCentres()[fvMesh_.faceOwner()[faceI]];

                    const vector cCTwo =
                        fvMesh_.cellCentres()[fvMesh_.faceNeighbour()[faceI]];

                    if (mag(cCOne - hitPoint) < mag(cCTwo - hitPoint))
                    {
                        rotCellId = fvMesh_.faceOwner()[faceI];
                    }
                    else
                    {
                        rotCellId = fvMesh_.faceNeighbour()[faceI];
                    }
                }
                else
                {
                    label patchI = fvMesh_.boundaryMesh().whichPatch(faceI);
                    if
                    (
                        isA<processorPolyPatch>(fvMesh_.boundaryMesh()[patchI])
                    )
                    {
                        const point& ownCc =
                            fvMesh_.cellCentres()[fvMesh_.faceOwner()[faceI]];

                        const vector cCentreOne = ownCc - hitPoint;

                        const vector nbrCc =
                            refCast<const processorPolyPatch>
                            (
                                fvMesh_.boundaryMesh()[patchI]
                            ).neighbFaceCellCentres()[faceI];

                        const vector cCentreTwo = nbrCc - hitPoint;

                        if (cCentreOne < cCentreTwo)
                        {
                            rotCellId = fvMesh_.faceOwner()[faceI];
                        }
                    }
                }

// Note: faces on boundaries that get hit are noy included as the
// pointDisplacement on boundaries is usually zero for this solver.

                // Search for closest direction on faces on mesh
                // and surface normal.
                if (rotCellId != -1)
                {
                    const labelList& cFaces = fvMesh_.cells()[rotCellId];

                    scalar cosMax(-GREAT);
                    label faceId(-1);
                    forAll(cFaces, k)
                    {
                        scalar tmp =
                            mag(nf[hit[faceI].index()] & nSfMesh[cFaces[k]]);

                        if (tmp > cosMax)
                        {
                            cosMax = tmp;
                            faceId = cFaces[k];
                        }
                    }
                    if
                    (
                        faceId != -1
                    &&
                        (
                            ::cos(degToRad(minAng_)) > cosMax
                            || cosMax > ::cos(degToRad(maxAng_))

                        )
                    )
                    {
                        cellRot_[rotCellId] =
                            nSfMesh[faceId] ^ nf[hit[faceI].index()];

                        const scalar magRot = mag(cellRot_[rotCellId]);

                        if (magRot > 0)
                        {
                            const scalar theta = ::asin(magRot);
                            quaternion q(cellRot_[rotCellId]/magRot, theta);
                            const tensor R = q.R();
                            const labelList& cPoints =
                                fvMesh_.cellPoints(rotCellId);

                            forAll (cPoints, j)
                            {
                                const label pointId = cPoints[j];

                                pointsCount[pointId]++;

                                const vector pointPos =
                                    fvMesh_.points()[pointId];

                                pointDisplacement_[pointId] +=
                                    (R & (pointPos - hitPoint))
                                  - (pointPos - hitPoint);
                            }
                        }
                    }
                }
            }
        }

        forAll (pointDisplacement_.internalField(), iPoint)
        {
            vector& point = pointDisplacement_.internalField()[iPoint];
            point /= pointsCount[iPoint];
        }


    }
}


void Foam::surfaceAlignedSBRStressFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the mtionSolver accordingly
    this->movePoints(fvMesh_.points());

    volVectorField& cellDisp = cellDisplacement();

    diffusivity().correct();

    // Calculate rotations on surface intersection
    calculateCellRot();

    tmp<volVectorField> tUd
    (
        new volVectorField
        (
            IOobject
            (
                "Ud",
                fvMesh_.time().timeName(),
                fvMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvMesh_,
            dimensionedVector("zero", dimLength, vector::zero),
            cellMotionBoundaryTypes<vector>
            (
                pointDisplacement().boundaryField()
            )
        )
    );

    volVectorField& Ud = tUd();

    const vectorList& C = fvMesh_.C();
    forAll (Ud, i)
    {
        pointMVCWeight pointInter(fvMesh_, C[i], i);
        Ud[i] = pointInter.interpolate(pointDisplacement_);
    }

    volScalarField Udx(Ud.component(vector::X));
    fvc::smooth(Udx, smoothFactor_);

    volScalarField Udy(Ud.component(vector::Y));
    fvc::smooth(Udy, smoothFactor_);

    volScalarField Udz(Ud.component(vector::Z));
    fvc::smooth(Udz, smoothFactor_);

    Ud.replace(vector::X, Udx);
    Ud.replace(vector::Y, Udy);
    Ud.replace(vector::Z, Udz);

    const volTensorField gradD("gradD", fvc::grad(Ud));

    tmp<volScalarField> tmu
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                fvMesh_.time().timeName(),
                fvMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvMesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    volScalarField& mu =  tmu();

    const scalarList& V = fvMesh_.V();
    mu.internalField() = (1.0/V);

    const volScalarField lambda(-mu);

    const volSymmTensorField newSigmaD
    (
        mu*twoSymm(gradD) + (lambda*I)*tr(gradD)
    );

    const volSymmTensorField magNewSigmaD(sigmaD_ + accFactor_*newSigmaD);

    const scalar diffSigmaD =
        gSum(mag(sigmaD_.oldTime().internalField()))
     -  gSum(mag(magNewSigmaD.internalField()));

    if (mag(diffSigmaD) > minSigmaDiff_)
    {
        sigmaD_ = magNewSigmaD;
    }

    const surfaceScalarField Df(diffusivity().operator()());
    pointDisplacement_.boundaryField().updateCoeffs();

    const volTensorField gradCd(fvc::grad(cellDisp));

    for (int nonOrth=0; nonOrth<=nNonOrthogonalCorr_; nonOrth++)
    {
        fvVectorMatrix DEqn
        (
            fvm::laplacian
            (
                2*Df*fvc::interpolate(mu),
                cellDisp,
                "laplacian(diffusivity,cellDisplacement)"
            )
          + fvc::div
            (
                Df*
                (
                    fvc::interpolate(mu)
                  * (fvMesh_.Sf() & fvc::interpolate(gradCd.T() - gradCd))
                  - fvc::interpolate(lambda)*fvMesh_.Sf()
                  * fvc::interpolate(tr(gradCd))
                )
            )
          ==
            fvc::div(sigmaD_)
        );

        DEqn.solve();
    }
}

// ************************************************************************* //
