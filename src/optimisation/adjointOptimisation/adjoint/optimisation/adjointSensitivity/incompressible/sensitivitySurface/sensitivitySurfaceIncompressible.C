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

#include "sensitivitySurfaceIncompressible.H"
#include "PrimitivePatchInterpolation.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivitySurface, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivitySurface,
    dictionary
);

// * * * * * * * * * * * Private  Member Functions  * * * * * * * * * * * * * //

void sensitivitySurface::read()
{
    includeSurfaceArea_ =
        dict().lookupOrDefault<bool>("includeSurfaceArea", true);
    includePressureTerm_ =
        dict().lookupOrDefault<bool>("includePressure", true);
    includeGradStressTerm_ =
        dict().lookupOrDefault<bool>("includeGradStressTerm", true);
    includeTransposeStresses_ =
        dict().lookupOrDefault<bool>("includeTransposeStresses", true);
    includeDivTerm_ = dict().lookupOrDefault<bool>("includeDivTerm", false);
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
    writeGeometricInfo_ =
        dict().lookupOrDefault<bool>("writeGeometricInfo", false);

    // Allocate new solvers if necessary
    if (includeDistance_ && eikonalSolver_.empty())
    {
        eikonalSolver_.reset
        (
            new adjointEikonalSolver
            (
                mesh_,
                dict_,
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
                dict_,
                *this,
                sensitivityPatchIDs_,
                eikonalSolver_
            )
        );
    }
}


void sensitivitySurface::addGeometricSens()
{
    if (includeObjective_)
    {
        // Grab objective refs
        PtrList<objective>& functions
            (objectiveManager_.getObjectiveFunctions());
        // Compute sens for all points in parameterized patches.
        // Interfacing points will be accumulated later
        autoPtr<pointBoundaryVectorField> pointSensdSdb
        (
            createZeroBoundaryPointFieldPtr<vector>(mesh_)
        );
        autoPtr<pointBoundaryVectorField> pointSensdndb
        (
            createZeroBoundaryPointFieldPtr<vector>(mesh_)
        );
        // Geometric (or "direct") sensitivities are better computed directly
        // on the points
        forAll(sensitivityPatchIDs_, pI)
        {
            const label patchI = sensitivityPatchIDs_[pI];
            const fvPatch& patch = mesh_.boundary()[patchI];
            vectorField nf(patch.nf());

            // point sens result for patch
            vectorField& patchdSdb = pointSensdSdb()[patchI];
            vectorField& patchdndb = pointSensdndb()[patchI];

            vectorField dSdbMultiplierTot(patch.size(), vector::zero);
            vectorField dndbMultiplierTot(patch.size(), vector::zero);
            forAll(functions, funcI)
            {
                dSdbMultiplierTot +=
                    functions[funcI].weight()
                   *functions[funcI].dSdbMultiplier(patchI);
                dndbMultiplierTot +=
                    functions[funcI].weight()
                   *functions[funcI].dndbMultiplier(patchI);
            }
            // Correspondance of local point addressing to global point
            // addressing
            const labelList& meshPoints = patch.patch().meshPoints();
            //  List with mesh faces. Global addressing
            const faceList& faces = mesh_.faces();
            //  Each local patch point belongs to these local patch faces
            //  (local numbering)
            const labelListList& patchPointFaces = patch.patch().pointFaces();
            //  index of first face in patch
            const label patchStartIndex = patch.start();
            //  geometry differentiation engine
            deltaBoundary dBoundary(mesh_);
            //  Loop over patch points.
            //  Collect contributions from each boundary face this point
            //  belongs to
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

                    // Element [1] is the variation in the (dimensional) normal
                    const tensor& deltaSf = deltaNormals[1];
                    patchdSdb[ppI] +=
                        dSdbMultiplierTot[localFaceIndex] & deltaSf;

                    // Element [2] is the variation in the unit normal
                    const tensor& deltaNf = deltaNormals[2];
                    patchdndb[ppI]   +=
                        dndbMultiplierTot[localFaceIndex] & deltaNf;
                }
            }
        }
        // Do parallel communications to avoid wrong values at processor
        // boundaries
        vectorField dSdbGlobal(mesh_.nPoints(), vector::zero);
        vectorField dndbGlobal(mesh_.nPoints(), vector::zero);
        forAll(sensitivityPatchIDs_, pI)
        {
            const label patchI = sensitivityPatchIDs_[pI];
            const labelList& meshPoints =
                mesh_.boundaryMesh()[patchI].meshPoints();
            forAll(meshPoints, ppI)
            {
                const label globaPointI = meshPoints[ppI];
                dSdbGlobal[globaPointI] += pointSensdSdb()[patchI][ppI];
                dndbGlobal[globaPointI] += pointSensdndb()[patchI][ppI];
            }
        }
        // Accumulate over processors
        syncTools::syncPointList
        (
            mesh_, dSdbGlobal, plusEqOp<vector>(), vector::zero
        );
        syncTools::syncPointList
        (
            mesh_, dndbGlobal, plusEqOp<vector>(), vector::zero
        );
        // Transfer back to local fields and map to faces
        forAll(sensitivityPatchIDs_, pI)
        {
            const label patchI = sensitivityPatchIDs_[pI];
            const fvPatch& patch = mesh_.boundary()[patchI];
            const labelList& meshPoints = patch.patch().meshPoints();
            const scalarField& magSf = patch.magSf();
            pointSensdSdb()[patchI].map(dSdbGlobal, meshPoints);
            pointSensdndb()[patchI].map(dndbGlobal, meshPoints);
            // Map dSf/dx and dnf/dx term from points to faces.  Ideally, all
            // sensitivities should be computed at points rather than faces.
            PrimitivePatchInterpolation<polyPatch> patchInter(patch.patch());
            vectorField dSdbFace
            (
                patchInter.pointToFaceInterpolate(pointSensdSdb()[patchI])
            );
            // dSdb already contains the face area. Divide with it to make it
            // compatible with the rest of the terms
            dSdbFace /= magSf;

            tmp<vectorField> tdndbFace =
                patchInter.pointToFaceInterpolate(pointSensdndb()[patchI]);
            const vectorField& dndbFace = tdndbFace();

            // Add to sensitivity fields
            wallFaceSensVecPtr_()[patchI] += dSdbFace + dndbFace;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivitySurface::sensitivitySurface
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
    writeGeometricInfo_(false),
    eikonalSolver_(nullptr),
    meshMovementSolver_(nullptr),

    nfOnPatchPtr_(nullptr),
    SfOnPatchPtr_(nullptr),
    CfOnPatchPtr_(nullptr)
{
    read();

    // Allocate boundary field pointer
    wallFaceSensVecPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    wallFaceSensNormalPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    wallFaceSensNormalVecPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));

    // Allocate fields to contain geometric info
    if (writeGeometricInfo_)
    {
        nfOnPatchPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "nfOnPatch",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                vector::zero
            )
        );

        SfOnPatchPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "SfOnPatch",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                vector::zero
            )
        );

        CfOnPatchPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "CfOnPatch",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                vector::zero
            )
        );
    }

    // Allocate appropriate space for the sensitivity field
    computeDerivativesSize();
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool sensitivitySurface::readDict(const dictionary& dict)
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


void sensitivitySurface::computeDerivativesSize()
{
    label nFaces(0);
    forAll(sensitivityPatchIDs_, pI)
    {
        const label patchI = sensitivityPatchIDs_[pI];
        nFaces += mesh_.boundary()[patchI].size();
    }
    derivatives_.setSize(nFaces);
}


const scalarField& sensitivitySurface::calculateSensitivities()
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

    // Update geometric fields for use by external users
    if (writeGeometricInfo_)
    {
        forAll(sensitivityPatchIDs_, pI)
        {
            const label patchI = sensitivityPatchIDs_[pI];
            const fvPatch& patch = mesh_.boundary()[patchI];
            tmp<vectorField> tnf = patch.nf();
            const vectorField& nf = tnf();
            const vectorField& Sf = patch.Sf();
            const vectorField& Cf = patch.Cf();

            nfOnPatchPtr_().boundaryFieldRef()[patchI] = nf;
            SfOnPatchPtr_().boundaryFieldRef()[patchI] = Sf;
            CfOnPatchPtr_().boundaryFieldRef()[patchI] = Cf;
        }
    }

    Info<< "    Calculating auxilary quantities " << endl;
    // Fields needed to calculate adjoint sensitivities
    const autoPtr<incompressible::RASModelVariables>&
       turbVars = primalVars_.RASModelVariables();
    const singlePhaseTransportModel& lamTransp = primalVars_.laminarTransport();
    volScalarField nuEff(lamTransp.nu() + turbVars->nutRef());
    volTensorField gradUa(fvc::grad(Ua));
    volTensorField gradU(fvc::grad(U));

    // Explicitly correct the boundary gradient to get rid of the tangential
    // component
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

    Info<< "    Calculating adjoint sensitivity. " << endl;

    // Sensitivities do not include locale surface area by default.
    // Part of the sensitivities that multiplies dxFace/db
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();

        // Includes spurious tangential gradU part. Deprecated
        /*
        vectorField stressAndPressureTerm =
              (
                - (
                       Ua.boundaryField()[patchI].snGrad()
                    + (gradUa.boundaryField()[patchI] & nf)
                  ) * nuEff.boundaryField()[patchI]
                + pa.boundaryField()[patchI] *nf
              ) & gradU.boundaryField()[patchI].T();
        */

        // Adjoint stress term
        vectorField stressTerm
        (
          - (
                Ua.boundaryField()[patchI].snGrad()
              & U.boundaryField()[patchI].snGrad()
            )
          * nuEff.boundaryField()[patchI]
          * nf
        );


        if (includeTransposeStresses_)
        {
            stressTerm -=
                nuEff.boundaryField()[patchI]
              * (
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
              * nf;
        }

        vectorField gradStressTerm(patch.size(), vector::zero);
        if (includeGradStressTerm_)
        {
            // Terms corresponding to contributions from converting delta to
            // thetas are added through the corresponding adjoint boundary
            // conditions instead of grabing contributions from the objective
            // function.  Useful to have a unified formulation for low- and
            // high-re meshes
            const fvPatchVectorField& Uab = Ua.boundaryField()[patchI];
            gradStressTerm = - ((Uab & nf)*gradp.boundaryField()[patchI]);
            gradStressTerm +=
            (
                Uab.component(0) * gradStressX.boundaryField()[patchI]
              + Uab.component(1) * gradStressY.boundaryField()[patchI]
              + Uab.component(2) * gradStressZ.boundaryField()[patchI]
            ) & nf;
        }

        // Adjoint pressure terms
        vectorField pressureTerm(patch.size(), vector::zero);
        if (includePressureTerm_)
        {
            pressureTerm =
            (
                (nf*pa.boundaryField()[patchI])
                & U.boundaryField()[patchI].snGrad()
            )* nf;
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

        PtrList<objective>& functions
            (objectiveManager_.getObjectiveFunctions());

        // Term from objectives including x directly (e.g. moments)
        vectorField dxdbMultiplierTot(pressureTerm.size(), vector::zero);
        if (includeObjective_)
        {
            forAll(functions, funcI)
            {
                dxdbMultiplierTot +=
                    functions[funcI].weight()
                  * (
                        functions[funcI].dxdbDirectMultiplier(patchI)
                    );
            }
        }

        // Fill in sensitivity fields
        wallFaceSensVecPtr_()[patchI] =
            stressTerm
          + gradStressTerm
          + pressureTerm
          + distanceTerm
          + meshMovementTerm
          + adjointTMsensitivities[patchI]
          + dxdbMultiplierTot;
    }

    // Add the sensitivity part corresponding to changes of the normal vector
    // Computed at points and mapped to faces
    addGeometricSens();

    // Project to normal face vector
    label nPassedFaces(0);
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf(patch.nf());
        const vectorField& nf = tnf();
        const scalarField& magSf = patch.magSf();

        if (includeSurfaceArea_)
        {
            wallFaceSensVecPtr_()[patchI] *= magSf;
        }

        wallFaceSensNormalPtr_()[patchI] = wallFaceSensVecPtr_()[patchI] & nf;
        wallFaceSensNormalVecPtr_()[patchI] =
            wallFaceSensNormalPtr_()[patchI] * nf;

        forAll(patch, fI)
        {
            derivatives_[nPassedFaces + fI]
                = wallFaceSensNormalPtr_()[patchI][fI];
        }
        nPassedFaces += patch.size();
    }

    // Write sens fields
    write(type());

    return (derivatives_);
}


autoPtr<adjointEikonalSolver>& sensitivitySurface::getAdjointEikonalSolver()
{
    return eikonalSolver_;
}


void sensitivitySurface::write(const word& baseName)
{
    // Determine suffix for fields holding the sens
    if (includeMeshMovement_)
    {
        surfaceFieldSuffix_ = word("ESI");
    }
    else
    {
        surfaceFieldSuffix_ = word("SI");
    }
    adjointSensitivity::write();

    if (writeGeometricInfo_)
    {
        nfOnPatchPtr_().write();
        SfOnPatchPtr_().write();
        CfOnPatchPtr_().write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace incompressible

// ************************************************************************* //
