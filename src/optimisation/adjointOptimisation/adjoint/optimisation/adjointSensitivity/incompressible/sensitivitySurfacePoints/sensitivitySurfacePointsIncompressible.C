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

// * * * * * * * * * * * Private  Member Functions  * * * * * * * * * * * * * //

void sensitivitySurfacePoints::read()
{
    includeSurfaceArea_ =
        dict().lookupOrDefault<bool>("includeSurfaceArea", false);
    includePressureTerm_ =
        dict().lookupOrDefault<bool>("includePressure", true);
    includeGradStressTerm_ =
        dict().lookupOrDefault<bool>("includeGradStressTerm", true);
    includeTransposeStresses_ =
        dict().lookupOrDefault<bool>("includeTransposeStresses", true);
    includeDivTerm_ =
        dict().lookupOrDefault<bool>("includeDivTerm", false);
    includeDistance_ =
        dict().lookupOrDefault<bool>
        (
            "includeDistance",
            adjointVars_.adjointTurbulence().ref().includeDistance()
        );
    includeMeshMovement_ =
        dict().lookupOrDefault<bool>("includeMeshMovement", true);
    includeObjective_ =
        dict().lookupOrDefault<bool>("includeObjectiveContribution", true);

    // Allocate new solvers if necessary
    if (includeDistance_ && eikonalSolver_.empty())
    {
        eikonalSolver_.reset
        (
            new adjointEikonalSolver
            (
                mesh_,
                dict(),
                primalVars_.RASModelVariables(),
                adjointVars_.adjointTurbulence(),
                sensitivityPatchIDs_
            )
        );
    }

    if (includeMeshMovement_ && meshMovementSolver_.empty())
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivitySurfacePoints::sensitivitySurfacePoints
(
    const fvMesh& mesh,
    const dictionary& dict,
    incompressibleVars& primalVars,
    incompressibleAdjointVars& adjointVars,
    objectiveManager& objectiveManager,
    fv::optionAdjointList& fvOptionsAdjoint
)
:
    adjointSensitivity
    (
        mesh,
        dict,
        primalVars,
        adjointVars,
        objectiveManager,
        fvOptionsAdjoint
    ),
    derivatives_(0),
    includeSurfaceArea_(false),
    includePressureTerm_(false),
    includeGradStressTerm_(false),
    includeTransposeStresses_(false),
    includeDivTerm_(false),
    includeDistance_(false),
    includeMeshMovement_(false),
    includeObjective_(false),
    eikonalSolver_(nullptr),
    meshMovementSolver_(nullptr)
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
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sensitivitySurfacePoints::readDict(const dictionary& dict)
{
    if (sensitivity::readDict(dict))
    {
        if (eikonalSolver_.valid())
        {
            eikonalSolver_().readDict(dict);
        }

        if (meshMovementSolver_.valid())
        {
            meshMovementSolver_().readDict(dict);
        }

        return true;
    }

    return false;
}


const scalarField& sensitivitySurfacePoints::calculateSensitivities()
{
    // Grab references
    const volScalarField& p = primalVars_.p();
    const volVectorField& U = primalVars_.U();

    const volScalarField& pa = adjointVars_.pa();
    const volVectorField& Ua = adjointVars_.Ua();
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars_.adjointTurbulence();

    // Restore to zero
    derivatives_ = Zero;
    forAll(mesh_.boundary(), patchI)
    {
        wallPointSensVecPtr_()[patchI] = vector::zero;
        wallPointSensNormalPtr_()[patchI] = Zero;
        wallPointSensNormalVecPtr_()[patchI] = vector::zero;
    }

    Info<< "    Calculating auxilary quantities " << endl;

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

    // solve extra equations if necessary
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

    // Terms from the adjoint turbulence model
    const boundaryVectorField& adjointTMsensitivities =
        adjointTurbulence->wallShapeSensitivities();

    // Objective references
    PtrList<objective>& functions(objectiveManager_.getObjectiveFunctions());

    Info<< "    Calculating adjoint sensitivity. " << endl;

    // The face-based part of the sensitivities, i.e. terms that multiply
    // dxFace/dxPoint. Sensitivities DO include locale surface area, to get
    // the correct weighting from the contributions of various faces.
    // Normalized at the end.
    autoPtr<boundaryVectorField> wallFaceSens
    (
        createZeroBoundaryPtr<vector>(mesh_)
    );

    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();
        const scalarField& magSf = patch.magSf();

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

        vectorField gradStressTerm(patch.size(), vector::zero);
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
            stressTerm -=
                nuEff.boundaryField()[patchI]
               *(
                    // Note: in case of laminar or low-Re flows,
                    // includes a spurious tangential gradUa component
                    // (gradUa.boundaryField()[patchI] & nf)
                    ((Ua.boundaryField()[patchI].snGrad() &nf)*nf)
                  & U.boundaryField()[patchI].snGrad()
                )
              * nf;
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
        vectorField pressureTerm(patch.size(), vector::zero);
        if (includePressureTerm_)
        {
            pressureTerm =
            (
                (nf * pa.boundaryField()[patchI])
              & U.boundaryField()[patchI].snGrad()
            )
           *nf;
        }

        // Distance related terms
        vectorField distanceTerm(pressureTerm.size(), vector::zero);
        if (includeDistance_)
        {
            distanceTerm = distanceSensPtr()[patchI];
        }

        // Mesh movement related terms
        vectorField meshMovementTerm(pressureTerm.size(), vector::zero);
        if (includeMeshMovement_)
        {
            meshMovementTerm = meshMovementSensPtr()[patchI];
        }


        vectorField dxdbMultiplierTot
        (
            mesh_.boundary()[patchI].size(), vector::zero
        );
        if (includeObjective_)
        {
            // Term from objectives multiplying dxdb
            forAll(functions, funcI)
            {
                dxdbMultiplierTot +=
                    functions[funcI].weight()
                  * functions[funcI].dxdbDirectMultiplier(patchI);
            }
        }

        // Fill in dxFace/dxPoint multiplier.
        // Missing geometric contributions which are directly computed on the
        // points
        wallFaceSens()[patchI] =
            stressTerm
          + gradStressTerm
          + pressureTerm
          + distanceTerm
          + meshMovementTerm
          + adjointTMsensitivities[patchI]
          + dxdbMultiplierTot;
        wallFaceSens()[patchI] *= magSf;
    }

    // polyPatch::pointNormals will give the wrong result for points
    // belonging to multiple patches or patch-processorPatch intersections.
    // Keeping a mesh-wide field to allow easy reduction using syncTools.
    // A bit expensive? Better way?
    vectorField pointNormals(mesh_.nPoints(), vector::zero);
    scalarField pointMagSf(mesh_.nPoints(), Zero);

    // Geometric (or "direct") sensitivities are better computed directly on
    // the points. Compute them and add the ones that depend on dxFace/dxPoint
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const scalarField& magSf = patch.magSf();
        vectorField nf(patch.nf());

        // Point sens result for patch
        vectorField& pointPatchSens = wallPointSensVecPtr_()[patchI];

        // Face sens for patch
        const vectorField& facePatchSens = wallFaceSens()[patchI];

        vectorField dSdbMultiplierTot(patch.size(), vector::zero);
        vectorField dndbMultiplierTot(patch.size(), vector::zero);
        forAll(functions, funcI)
        {
            dSdbMultiplierTot +=
                functions[funcI].weight() //includes surface by itself
               *functions[funcI].dSdbMultiplier(patchI);
            dndbMultiplierTot +=
                functions[funcI].weight()
               *functions[funcI].dndbMultiplier(patchI)
               *magSf;
        }

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
                tensorField p_d(faceI.size(), tensor::zero);
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
                        dSdbMultiplierTot[localFaceIndex] & deltaSf;

                    // Element [2] is the variation in the unit normal
                    const tensor& deltaNf = deltaNormals[2];
                    pointPatchSens[ppI] +=
                        dndbMultiplierTot[localFaceIndex] & deltaNf;
                }

                // Accumulate information for point normals
                pointNormals[meshPoints[ppI]] += nf[localFaceIndex];
                pointMagSf[meshPoints[ppI]] += magSf[localFaceIndex];
            }
        }
    }

    // Do parallel communications to avoid wrong values at processor boundaries
    // - global field for accumulation
    vectorField pointSensGlobal(mesh_.nPoints(), vector::zero);
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

    // Accumulate dJ/dx_i, pointNormals and pointFaces number
    syncTools::syncPointList
    (
        mesh_,
        pointSensGlobal,
        plusEqOp<vector>(),
        vector::zero
    );
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

            // avoid storing unit point normals in the global list since we
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

    // Write sens fields
    write(type());

    return (derivatives_);
}


void sensitivitySurfacePoints::write(const word& baseName)
{
    //determine suffix for fields holding the sens
    if (includeMeshMovement_)
    {
        surfaceFieldSuffix_ = "ESI";
    }
    else
    {
        surfaceFieldSuffix_ = "SI";
    }
    adjointSensitivity::write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace incompressible

// ************************************************************************* //
