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

#include "sensitivitySurfacePointsIncompressible.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivitySurfacePoints, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivitySurfacePoints,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void sensitivitySurfacePoints::read()
{
    includeSurfaceArea_ =
        dict().getOrDefault<bool>("includeSurfaceArea", false);
    includePressureTerm_ =
        dict().getOrDefault<bool>("includePressure", true);
    includeGradStressTerm_ =
        dict().getOrDefault<bool>("includeGradStressTerm", true);
    includeTransposeStresses_ =
        dict().getOrDefault<bool>("includeTransposeStresses", true);
    useSnGradInTranposeStresses_ =
        dict().getOrDefault<bool>("useSnGradInTranposeStresses", false);
    includeDivTerm_ =
        dict().getOrDefault<bool>("includeDivTerm", false);
    includeDistance_ =
        dict().getOrDefault<bool>
        (
            "includeDistance",
            adjointVars_.adjointTurbulence().ref().includeDistance()
        );
    includeMeshMovement_ =
        dict().getOrDefault<bool>("includeMeshMovement", true);
    includeObjective_ =
        dict().getOrDefault<bool>("includeObjectiveContribution", true);

    // Allocate new solvers if necessary
    if (includeDistance_ && !eikonalSolver_)
    {
        eikonalSolver_.reset
        (
            new adjointEikonalSolver
            (
                mesh_,
                dict(),
                primalVars_.RASModelVariables(),
                adjointVars_,
                sensitivityPatchIDs_
            )
        );
    }

    if (includeMeshMovement_ && !meshMovementSolver_)
    {
        meshMovementSolver_.reset
        (
            new adjointMeshMovementSolver
            (
                mesh_,
                dict(),
                *this,
                sensitivityPatchIDs_,
                eikonalSolver_
            )
        );
    }
}


void sensitivitySurfacePoints::finaliseFaceMultiplier()
{
    // Solve extra equations if necessary
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    autoPtr<boundaryVectorField> distanceSensPtr(nullptr);
    if (includeDistance_)
    {
        eikonalSolver_->solve();
        distanceSensPtr.reset(createZeroBoundaryPtr<vector>(mesh_));
        const boundaryVectorField& sens =
            eikonalSolver_->distanceSensitivities();
        for (const label patchI : sensitivityPatchIDs_)
        {
            distanceSensPtr()[patchI] = sens[patchI];
        }
    }

    autoPtr<boundaryVectorField> meshMovementSensPtr(nullptr);
    if (includeMeshMovement_)
    {
        meshMovementSolver_->solve();
        meshMovementSensPtr.reset(createZeroBoundaryPtr<vector>(mesh_));
        const boundaryVectorField& sens =
            meshMovementSolver_->meshMovementSensitivities();
        for (const label patchI : sensitivityPatchIDs_)
        {
            meshMovementSensPtr()[patchI] = sens[patchI];
        }
    }

    // Add to other terms multiplying dxFace/dxPoints
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const scalarField& magSf = patch.magSf();
        // Distance related terms
        if (includeDistance_)
        {
            wallFaceSens_()[patchI] += distanceSensPtr()[patchI];
        }

        // Mesh movement related terms
        if (includeMeshMovement_)
        {
            wallFaceSens_()[patchI] += meshMovementSensPtr()[patchI];
        }

        // Add local face area
        //~~~~~~~~~~~~~~~~~~~~
        // Sensitivities DO include locale surface area, to get
        // the correct weighting from the contributions of various faces.
        // Normalized at the end.
        // dSfdbMult already includes the local area. No need to re-multiply
        wallFaceSens_()[patchI] *= magSf;
        dnfdbMult_()[patchI] *= magSf;
    }
}


void sensitivitySurfacePoints::finalisePointSensitivities()
{
    // Geometric (or "direct") sensitivities are better computed directly on
    // the points. Compute them and add the ones that depend on dxFace/dxPoint
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        vectorField nf(patch.nf());

        // Point sens result for patch
        vectorField& pointPatchSens = wallPointSensVecPtr_()[patchI];

        // Face sens for patch
        const vectorField& facePatchSens = wallFaceSens_()[patchI];

        // Geometry variances
        const vectorField& dSfdbMultPatch = dSfdbMult_()[patchI];
        const vectorField& dnfdbMultPatch = dnfdbMult_()[patchI];

        // Correspondance of local point addressing to global point addressing
        const labelList& meshPoints = patch.patch().meshPoints();

        // List with mesh faces. Global addressing
        const faceList& faces = mesh_.faces();

        // Each local patch point belongs to these local patch faces
        // (local numbering)
        const labelListList& patchPointFaces = patch.patch().pointFaces();

        // Index of first face in patch
        const label patchStartIndex = patch.start();

        // Geometry differentiation engine
        deltaBoundary dBoundary(mesh_);

        // Loop over patch points.
        // Collect contributions from each boundary face this point belongs to
        forAll(meshPoints, ppI)
        {
            const labelList& pointFaces = patchPointFaces[ppI];
            forAll(pointFaces, pfI)
            {
                label localFaceIndex = pointFaces[pfI];
                label globalFaceIndex = patchStartIndex + localFaceIndex;
                const face& faceI = faces[globalFaceIndex];

                // Point coordinates. All indices in global numbering
                pointField p(faceI.points(mesh_.points()));
                tensorField p_d(faceI.size(), Zero);
                forAll(faceI, facePointI)
                {
                    if (faceI[facePointI] == meshPoints[ppI])
                    {
                        p_d[facePointI] = tensor::I;
                    }
                }
                tensorField deltaNormals =
                    dBoundary.makeFaceCentresAndAreas_d(p, p_d);

                // Element [0] is the variation in the face center
                // (dxFace/dxPoint)
                const tensor& deltaCf = deltaNormals[0];
                pointPatchSens[ppI] += facePatchSens[localFaceIndex] & deltaCf;

                // Term multiplying d(Sf)/d(point displacement) and
                // d(nf)/d(point displacement)
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (includeObjective_)
                {
                    // Element [1] is the variation in the (dimensional) normal
                    const tensor& deltaSf = deltaNormals[1];
                    pointPatchSens[ppI] +=
                        dSfdbMultPatch[localFaceIndex] & deltaSf;

                    // Element [2] is the variation in the unit normal
                    const tensor& deltaNf = deltaNormals[2];
                    pointPatchSens[ppI] +=
                        dnfdbMultPatch[localFaceIndex] & deltaNf;
                }
            }
        }
    }
}


void sensitivitySurfacePoints::constructGlobalPointNormalsAndAreas
(
    vectorField& pointNormals,
    scalarField& pointMagSf
)
{
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const scalarField& magSf = patch.magSf();
        vectorField nf(patch.nf());

        // Correspondance of local point addressing to global point addressing
        const labelList& meshPoints = patch.patch().meshPoints();

        // Each local patch point belongs to these local patch faces
        // (local numbering)
        const labelListList& patchPointFaces = patch.patch().pointFaces();

        // Loop over patch points
        forAll(meshPoints, ppI)
        {
            const labelList& pointFaces = patchPointFaces[ppI];
            forAll(pointFaces, pfI)
            {
                const label localFaceIndex = pointFaces[pfI];

                // Accumulate information for point normals
                pointNormals[meshPoints[ppI]] += nf[localFaceIndex];
                pointMagSf[meshPoints[ppI]] += magSf[localFaceIndex];
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pointNormals,
        plusEqOp<vector>(),
        vector::zero
    );
    syncTools::syncPointList
    (
        mesh_,
        pointMagSf,
        plusEqOp<scalar>(),
        scalar(0)
    );
}


void sensitivitySurfacePoints::setSuffixName()
{
    word suffix(dict().getOrDefault<word>("suffix", word::null));
    // Determine suffix for fields holding the sens
    if (includeMeshMovement_)
    {
        shapeSensitivitiesBase::setSuffix
        (
            adjointVars_.solverName() + "ESI" + suffix
        );
    }
    else
    {
        shapeSensitivitiesBase::setSuffix
        (
            adjointVars_.solverName() + "SI" + suffix
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivitySurfacePoints::sensitivitySurfacePoints
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager
)
:
    adjointSensitivity
    (
        mesh,
        dict,
        primalVars,
        adjointVars,
        objectiveManager
    ),
    shapeSensitivitiesBase(mesh, dict),
    includeSurfaceArea_(false),
    includePressureTerm_(false),
    includeGradStressTerm_(false),
    includeTransposeStresses_(false),
    useSnGradInTranposeStresses_(false),
    includeDivTerm_(false),
    includeDistance_(false),
    includeMeshMovement_(false),
    includeObjective_(false),
    eikonalSolver_(nullptr),
    meshMovementSolver_(nullptr),
    wallFaceSens_(createZeroBoundaryPtr<vector>(mesh_)),
    dSfdbMult_(createZeroBoundaryPtr<vector>(mesh_)),
    dnfdbMult_(createZeroBoundaryPtr<vector>(mesh_))

{
    read();

    // Allocate boundary field pointer
    wallPointSensVecPtr_.reset(createZeroBoundaryPointFieldPtr<vector>(mesh_));
    wallPointSensNormalPtr_.reset
    (
        createZeroBoundaryPointFieldPtr<scalar>(mesh_)
    );
    wallPointSensNormalVecPtr_.reset
    (
        createZeroBoundaryPointFieldPtr<vector>(mesh_)
    );

    // Allocate appropriate space for sensitivities
    label nTotalPoints(0);
    for (const label patchI : sensitivityPatchIDs_)
    {
        label nPoints = mesh_.boundaryMesh()[patchI].nPoints();
        nTotalPoints += returnReduce(nPoints, sumOp<label>());
    }

    // Derivatives for all (x,y,z) components of the displacement are kept
    derivatives_ = scalarField(3*nTotalPoints, Zero);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sensitivitySurfacePoints::readDict(const dictionary& dict)
{
    if (sensitivity::readDict(dict))
    {
        if (eikonalSolver_)
        {
            eikonalSolver_().readDict(dict);
        }

        if (meshMovementSolver_)
        {
            meshMovementSolver_().readDict(dict);
        }

        return true;
    }

    return false;
}


void sensitivitySurfacePoints::accumulateIntegrand(const scalar dt)
{
    // Grab references
    const volScalarField& p = primalVars_.p();
    const volVectorField& U = primalVars_.U();

    const volScalarField& pa = adjointVars_.pa();
    const volVectorField& Ua = adjointVars_.Ua();
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars_.adjointTurbulence();

    DebugInfo
        << "    Calculating auxilary quantities " << endl;

    // Fields needed to calculate adjoint sensitivities
    volScalarField nuEff(adjointTurbulence->nuEff());
    volTensorField gradUa(fvc::grad(Ua));
    volTensorField gradU(fvc::grad(U));

    // Explicitly correct the boundary gradient to get rid of the
    // tangential component
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<wallFvPatch>(patch))
        {
            tmp<vectorField> tnf = mesh_.boundary()[patchI].nf();
            const vectorField& nf = tnf();
            gradU.boundaryFieldRef()[patchI] =
                nf*U.boundaryField()[patchI].snGrad();
        }
    }

    // Auxiliary terms
    volVectorField gradp(fvc::grad(p));
    volTensorField stress(nuEff*(gradU + T(gradU)));
    autoPtr<volVectorField> stressXPtr
    (
        createZeroFieldPtr<vector>(mesh_, "stressX", stress.dimensions())
    );
    autoPtr<volVectorField> stressYPtr
    (
        createZeroFieldPtr<vector>(mesh_, "stressY", stress.dimensions())
    );
    autoPtr<volVectorField> stressZPtr
    (
        createZeroFieldPtr<vector>(mesh_, "stressZ", stress.dimensions())
    );

    stressXPtr().replace(0, stress.component(0));
    stressXPtr().replace(1, stress.component(1));
    stressXPtr().replace(2, stress.component(2));

    stressYPtr().replace(0, stress.component(3));
    stressYPtr().replace(1, stress.component(4));
    stressYPtr().replace(2, stress.component(5));

    stressZPtr().replace(0, stress.component(6));
    stressZPtr().replace(1, stress.component(7));
    stressZPtr().replace(2, stress.component(8));

    volTensorField gradStressX(fvc::grad(stressXPtr()));
    volTensorField gradStressY(fvc::grad(stressYPtr()));
    volTensorField gradStressZ(fvc::grad(stressZPtr()));

    // Solve extra equations if necessary
    if (includeDistance_)
    {
        eikonalSolver_->accumulateIntegrand(dt);
    }

    autoPtr<boundaryVectorField> meshMovementSensPtr(nullptr);
    if (includeMeshMovement_)
    {
        meshMovementSolver_->accumulateIntegrand(dt);
    }

    // Terms from the adjoint turbulence model
    const boundaryVectorField& adjointTMsensitivities =
        adjointTurbulence->wallShapeSensitivities();

    // Objective references
    PtrList<objective>& functions(objectiveManager_.getObjectiveFunctions());

    DebugInfo
        << "    Calculating adjoint sensitivity. " << endl;

    // The face-based part of the sensitivities, i.e. terms that multiply
    // dxFace/dxPoint.
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();

        // Adjoint stress term
        // vectorField stressTerm
        //     (
        //        -(nf & DUa.boundaryField()[patchI])
        //        *nuEff.boundaryField()[patchI]
        //       & gradU.boundaryField()[patchI].T();
        //     )

        vectorField stressTerm
        (
          - (
                Ua.boundaryField()[patchI].snGrad()
              & U.boundaryField()[patchI].snGrad()
            )
          * nuEff.boundaryField()[patchI]
          * nf
        );

        vectorField gradStressTerm(patch.size(), Zero);
        if (includeGradStressTerm_)
        {
            // Terms corresponding to contributions from converting delta to
            // thetas are added through the corresponding adjoint boundary
            // conditions instead of grabing contributions from the objective
            // function. Useful to have a unified formulation for low- and
            // high-re meshes
            const fvPatchVectorField& Uab = Ua.boundaryField()[patchI];
            gradStressTerm = (-((Uab & nf)*gradp.boundaryField()[patchI]));
            gradStressTerm +=
            (
                Uab.component(0)*gradStressX.boundaryField()[patchI]
              + Uab.component(1)*gradStressY.boundaryField()[patchI]
              + Uab.component(2)*gradStressZ.boundaryField()[patchI]
            ) & nf;
        }

        if (includeTransposeStresses_)
        {
            vectorField gradUaNf
                (
                    useSnGradInTranposeStresses_ ?
                    (Ua.boundaryField()[patchI].snGrad() & nf)*nf :
                    (gradUa.boundaryField()[patchI] & nf)
                );
            stressTerm -=
                nuEff.boundaryField()[patchI]
               *(gradUaNf & U.boundaryField()[patchI].snGrad())
               *nf;
        }

        if (includeDivTerm_)
        {
            stressTerm +=
                scalar(1./3.)*nuEff.boundaryField()[patchI]
              * (
                    ((Ua.boundaryField()[patchI].snGrad() &nf)*nf)
                  & U.boundaryField()[patchI].snGrad()
                )
               *nf;
        }

        // Adjoint pressure terms
        vectorField pressureTerm(patch.size(), Zero);
        if (includePressureTerm_)
        {
            pressureTerm =
            (
                (nf * pa.boundaryField()[patchI])
              & U.boundaryField()[patchI].snGrad()
            )
           *nf;
        }

        vectorField dxdbMultiplierTot(patch.size(), Zero);
        if (includeObjective_)
        {
            // Term from objectives multiplying dxdb
            forAll(functions, funcI)
            {
                const scalar wei = functions[funcI].weight();
                // dt added in wallFaceSens_
                dxdbMultiplierTot +=
                    wei*functions[funcI].dxdbDirectMultiplier(patchI);

                // Fill in multipliers of d(Sf)/db and d(nf)/db
                dSfdbMult_()[patchI] +=
                    wei*dt*functions[funcI].dSdbMultiplier(patchI);
                dnfdbMult_()[patchI] +=
                    wei*dt*functions[funcI].dndbMultiplier(patchI);
            }
        }

        // Fill in dxFace/dxPoint multiplier.
        // Missing geometric contributions which are directly computed on the
        // points
        wallFaceSens_()[patchI] +=
        (
            stressTerm
          + gradStressTerm
          + pressureTerm
          + adjointTMsensitivities[patchI]
          + dxdbMultiplierTot
        )*dt;
    }
}


void sensitivitySurfacePoints::assembleSensitivities()
{
    // Add remaining parts to term multiplying dxFace/dxPoints
    // Solves for post-processing PDEs
    finaliseFaceMultiplier();

    // Geometric (or "direct") sensitivities are better computed directly on
    // the points. Compute them and add the ones that depend on dxFace/dxPoint
    finalisePointSensitivities();

    // polyPatch::pointNormals will give the wrong result for points
    // belonging to multiple patches or patch-processorPatch intersections.
    // Keeping a mesh-wide field to allow easy reduction using syncTools.
    // A bit expensive? Better way?
    vectorField pointNormals(mesh_.nPoints(), Zero);
    scalarField pointMagSf(mesh_.nPoints(), Zero);
    constructGlobalPointNormalsAndAreas(pointNormals, pointMagSf);

    // Do parallel communications to avoid wrong values at processor boundaries
    // Global field for accumulation
    vectorField pointSensGlobal(mesh_.nPoints(), Zero);
    for (const label patchI : sensitivityPatchIDs_)
    {
        const labelList& meshPoints = mesh_.boundaryMesh()[patchI].meshPoints();
        forAll(meshPoints, ppI)
        {
            const label globaPointI = meshPoints[ppI];
            pointSensGlobal[globaPointI] +=
                wallPointSensVecPtr_()[patchI][ppI];
        }
    }

    // Accumulate dJ/dx_i
    syncTools::syncPointList
    (
        mesh_,
        pointSensGlobal,
        plusEqOp<vector>(),
        vector::zero
    );

    // Transfer back to local fields
    for (const label patchI : sensitivityPatchIDs_)
    {
        const labelList& meshPoints =
            mesh_.boundaryMesh()[patchI].meshPoints();
        wallPointSensVecPtr_()[patchI].map(pointSensGlobal, meshPoints);
    }

    // Compute normal sens and append to return field
    label nPassedDVs(0);
    for (const label patchI : sensitivityPatchIDs_)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        List<scalarField> procPatchSens(Pstream::nProcs());
        //if (patch.size()>0)
        {
            const labelList& meshPoints = patch.meshPoints();

            // Avoid storing unit point normals in the global list since we
            // might divide multiple times with the number of faces belonging
            // to the point. Instead do the division locally, per patch use
            vectorField patchPointNormals(pointNormals, meshPoints);
            patchPointNormals /= mag(patchPointNormals) + VSMALL;
            if (!includeSurfaceArea_)
            {
                wallPointSensVecPtr_()[patchI] /=
                    scalarField(pointMagSf, meshPoints);
            }
            wallPointSensNormalPtr_()[patchI] =
                wallPointSensVecPtr_()[patchI]
              & patchPointNormals;
            wallPointSensNormalVecPtr_()[patchI] =
                wallPointSensNormalPtr_()[patchI]
               *patchPointNormals;

            // 1. Gather sens from all processors for this patch and communicate
            // them back. Potentially large memory overhead but the rest of the
            // code structure assumes that all procs know all sensitivity
            // derivatives
            //
            // 2. Transfer vectorial sensitivities to scalarField.
            // Needed since the normal point vector is wrongly computed at patch
            // boundaries and cannot be used to reconstruct a vectorial movement
            // from just its normal component
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            procPatchSens[Pstream::myProcNo()].setSize
            (
                3*wallPointSensNormalVecPtr_()[patchI].size()
            );
            scalarField& patchScalarSens = procPatchSens[Pstream::myProcNo()];
            forAll(wallPointSensNormalVecPtr_()[patchI], ptI)
            {
                patchScalarSens[3*ptI] =
                    wallPointSensNormalVecPtr_()[patchI][ptI].x();
                patchScalarSens[3*ptI + 1] =
                    wallPointSensNormalVecPtr_()[patchI][ptI].y();
                patchScalarSens[3*ptI + 2] =
                    wallPointSensNormalVecPtr_()[patchI][ptI].z();
            }
            Pstream::gatherList(procPatchSens);
            Pstream::scatterList(procPatchSens);

            forAll(procPatchSens, procI)
            {
                const scalarField& procSens = procPatchSens[procI];
                forAll(procSens, dvI)
                {
                    derivatives_[nPassedDVs + dvI] = procSens[dvI];
                }
                nPassedDVs += procSens.size();
            }
        }
    }
}


void sensitivitySurfacePoints::clearSensitivities()
{
    // Reset terms in post-processing PDEs
    if (includeDistance_)
    {
        eikonalSolver_->reset();
    }
    if (includeMeshMovement_)
    {
        meshMovementSolver_->reset();
    }

    // Reset local fields to zero
    wallFaceSens_() = vector::zero;
    dSfdbMult_() = vector::zero;
    dnfdbMult_() = vector::zero;

    // Reset sensitivity fields
    adjointSensitivity::clearSensitivities();
    shapeSensitivitiesBase::clearSensitivities();
}


void sensitivitySurfacePoints::write(const word& baseName)
{
    setSuffixName();
    adjointSensitivity::write();
    shapeSensitivitiesBase::write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace incompressible

// ************************************************************************* //
