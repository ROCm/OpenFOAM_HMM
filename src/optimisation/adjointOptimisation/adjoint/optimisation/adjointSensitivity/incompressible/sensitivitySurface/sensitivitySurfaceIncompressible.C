/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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

#include "sensitivitySurfaceIncompressible.H"
#include "PrimitivePatchInterpolation.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"
#include "faMatrices.H"
#include "famSup.H"
#include "famLaplacian.H"
#include "volSurfaceMapping.H"
#include "fixedValueFaPatchFields.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivitySurface, 1);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivitySurface,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

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
        for (const label patchI : sensitivityPatchIDs_)
        {
            const fvPatch& patch = mesh_.boundary()[patchI];
            const vectorField nf(patch.nf());

            // point sens result for patch
            vectorField& patchdSdb = pointSensdSdb()[patchI];
            vectorField& patchdndb = pointSensdndb()[patchI];

            vectorField dSdbMultiplierTot(patch.size(), Zero);
            vectorField dndbMultiplierTot(patch.size(), Zero);
            for (auto& fun : functions)
            {
                dSdbMultiplierTot += fun.weight()*fun.dSdbMultiplier(patchI);
                dndbMultiplierTot += fun.weight()*fun.dndbMultiplier(patchI);
            }
            // Correspondence of local point addressing to global point
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
                for (label localFaceIndex : pointFaces)
                {
                    label globalFaceIndex = patchStartIndex + localFaceIndex;
                    const face& faceI = faces[globalFaceIndex];
                    // Point coordinates. All indices in global numbering
                    const pointField p(faceI.points(mesh_.points()));
                    tensorField p_d(faceI.size(), Zero);
                    forAll(faceI, facePointI)
                    {
                        if (faceI[facePointI] == meshPoints[ppI])
                        {
                            p_d[facePointI] = tensor::I;
                        }
                    }
                    const tensorField deltaNormals
                    (
                        dBoundary.makeFaceCentresAndAreas_d(p, p_d)
                    );

                    // Element [1] is the variation in the (dimensional) normal
                    const tensor& deltaSf = deltaNormals[1];
                    patchdSdb[ppI] +=
                        dSdbMultiplierTot[localFaceIndex] & deltaSf;

                    // Element [2] is the variation in the unit normal
                    const tensor& deltaNf = deltaNormals[2];
                    patchdndb[ppI] +=
                        dndbMultiplierTot[localFaceIndex] & deltaNf;
                }
            }
        }
        // Do parallel communications to avoid wrong values at processor
        // boundaries
        vectorField dSdbGlobal(mesh_.nPoints(), Zero);
        vectorField dndbGlobal(mesh_.nPoints(), Zero);
        for (const label patchI : sensitivityPatchIDs_)
        {
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
        for (const label patchI : sensitivityPatchIDs_)
        {
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


void sensitivitySurface::setSuffixName()
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


void sensitivitySurface::smoothSensitivities()
{
    // Read in parameters
    const label iters(dict().getOrDefault<label>("iters", 500));
    const scalar tolerance(dict().getOrDefault<scalar>("tolerance", 1.e-06));
    autoPtr<faMesh> aMeshPtr(nullptr);

    IOobject faceLabels
    (
        "faceLabels",
        mesh_.time().findInstance
        (
            mesh_.dbDir()/faMesh::meshSubDir,
            "faceLabels",
            IOobject::READ_IF_PRESENT
        ),
        faMesh::meshSubDir,
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    // If the faMesh already exists, read it
    if (faceLabels.typeHeaderOk<labelIOList>(false))
    {
        Info<< "Reading the already constructed faMesh" << endl;
        aMeshPtr.reset(new faMesh(mesh_));
    }
    else
    {
        // Dictionary used to construct the faMesh
        dictionary faMeshDefinition;

        IOobject faMeshDefinitionDict
        (
            "faMeshDefinition",
            mesh_.time().caseSystem(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        // If the faMeshDefinitionDict exists, use it to construct the mesh
        if (faMeshDefinitionDict.typeHeaderOk<IOdictionary>(false))
        {
            Info<< "Reading faMeshDefinition from system " << endl;
            faMeshDefinition = IOdictionary(faMeshDefinitionDict);
        }
        // Otherwise, faMesh is generated from all patches on which we compute
        // sensitivities
        else
        {
            Info<< "Constructing faMeshDefinition from sensitivity patches"
                << endl;
            wordList polyMeshPatches(sensitivityPatchIDs_.size());
            label i(0);
            for (const label patchID : sensitivityPatchIDs_)
            {
                polyMeshPatches[i++] = mesh_.boundary()[patchID].name();
            }
            faMeshDefinition.add<wordList>("polyMeshPatches", polyMeshPatches);
            (void)faMeshDefinition.subDictOrAdd("boundary");
        }

        // Construct faMesh
        aMeshPtr.reset(new faMesh(mesh_, faMeshDefinition));
    }
    faMesh& aMesh = aMeshPtr.ref();

    // Physical radius of the smoothing, provided either directly or computed
    // based on the average 'length' of boundary faces
    const scalar Rphysical
        (dict().getOrDefault<scalar>("radius", computeRadius(aMesh)));
    DebugInfo
        << "Physical radius of the sensitivity smoothing "
        << Rphysical << nl << endl;

    // Radius used as the diffusivity in the Helmholtz filter, computed as a
    // function of the physical radius
    const dimensionedScalar RpdeSqr
    (
        "RpdeSqr", dimArea, sqr(Rphysical/(2.*::sqrt(3.)))
    );

    dimensionedScalar one("1", dimless, 1.);

    // Mapping engine
    volSurfaceMapping vsm(aMesh);

    // Source term in faMatrix needs to be an areaField
    areaScalarField sens
    (
        IOobject
        (
            "sens",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFaPatchField<scalar>::typeName
    );

    // Copy sensitivities to area field
    sens.primitiveFieldRef() =
        vsm.mapToSurface<scalar>(wallFaceSensNormalPtr_());

    // Initialisation of the smoothed sensitivities field based on the original
    // sensitivities
    areaScalarField smoothedSens("smoothedSens", sens);
    for (label iter = 0; iter < iters; ++iter)
    {
        Info<< "Sensitivity smoothing iteration " << iter << endl;

        faScalarMatrix smoothEqn
        (
            fam::Sp(one, smoothedSens)
          - fam::laplacian(RpdeSqr, smoothedSens)
         ==
            sens
        );

        smoothEqn.relax();

        const scalar residual(mag(smoothEqn.solve().initialResidual()));

        DebugInfo
            << "Max smoothSens " << gMax(mag(smoothedSens)()) << endl;

        // Print execution time
        mesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance)
        {
            Info<< "\n***Reached smoothing equation convergence limit, "
                   "iteration " << iter << "***\n\n";
            break;
        }
    }

    // Field used to write the smoothed sensitivity field to file
    volScalarField volSmoothedSens
    (
        IOobject
        (
            "smoothedSurfaceSens" + surfaceFieldSuffix_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    // Transfer result back to volField and write
    vsm.mapToVolume(smoothedSens, volSmoothedSens.boundaryFieldRef());
    volSmoothedSens.write();
}


scalar sensitivitySurface::computeRadius(const faMesh& aMesh)
{
    scalar averageArea(gAverage(aMesh.S().field()));
    const Vector<label>& geometricD = mesh_.geometricD();
    const boundBox& bounds = mesh_.bounds();
    forAll(geometricD, iDir)
    {
        if (geometricD[iDir] == -1)
        {
            averageArea /= bounds.span()[iDir];
        }
    }
    scalar mult = dict().getOrDefault<scalar>("meanRadiusMultiplier", 10);

    return mult*pow(averageArea, scalar(1)/scalar(mesh_.nGeometricD() - 1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivitySurface::sensitivitySurface
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
    writeGeometricInfo_(false),
    smoothSensitivities_(false),
    eikonalSolver_(nullptr),
    meshMovementSolver_(nullptr),

    nfOnPatchPtr_(nullptr),
    SfOnPatchPtr_(nullptr),
    CfOnPatchPtr_(nullptr)
{
    read();
    setSuffixName();

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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sensitivitySurface::read()
{
    includeSurfaceArea_ =
        dict().getOrDefault<bool>("includeSurfaceArea", true);
    includePressureTerm_ =
        dict().getOrDefault<bool>("includePressure", true);
    includeGradStressTerm_ =
        dict().getOrDefault<bool>("includeGradStressTerm", true);
    includeTransposeStresses_ =
        dict().getOrDefault<bool>("includeTransposeStresses", true);
    useSnGradInTranposeStresses_ =
        dict().getOrDefault<bool>("useSnGradInTranposeStresses", false);
    includeDivTerm_ = dict().getOrDefault<bool>("includeDivTerm", false);
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
    writeGeometricInfo_ =
        dict().getOrDefault<bool>("writeGeometricInfo", false);
    smoothSensitivities_ =
        dict().getOrDefault<bool>("smoothSensitivities", false);

    // Allocate new solvers if necessary
    if (includeDistance_ && !eikonalSolver_)
    {
        eikonalSolver_.reset
        (
            new adjointEikonalSolver
            (
                mesh_,
                dict_,
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
                dict_,
                *this,
                sensitivityPatchIDs_,
                eikonalSolver_
            )
        );
    }
}


bool sensitivitySurface::readDict(const dictionary& dict)
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


void sensitivitySurface::computeDerivativesSize()
{
    label nFaces(0);
    for (const label patchI : sensitivityPatchIDs_)
    {
        nFaces += mesh_.boundary()[patchI].size();
    }
    derivatives_.setSize(nFaces);
}


void sensitivitySurface::accumulateIntegrand(const scalar dt)
{
    // Grab references
    const volScalarField& p = primalVars_.p();
    const volVectorField& U = primalVars_.U();

    const volScalarField& pa = adjointVars_.pa();
    const volVectorField& Ua = adjointVars_.Ua();
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars_.adjointTurbulence();

    Info<< "    Calculating auxiliary quantities " << endl;
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

    // Accumulate source for additional post-processing PDEs, if necessary
    if (includeDistance_)
    {
        eikonalSolver_->accumulateIntegrand(dt);
    }

    if (includeMeshMovement_)
    {
        meshMovementSolver_->accumulateIntegrand(dt);
    }

    // Terms from the adjoint turbulence model
    const boundaryVectorField& adjointTMsensitivities =
        adjointTurbulence->wallShapeSensitivities();

    DebugInfo
        << "    Calculating adjoint sensitivity. " << endl;

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
              * nf;
        }

        vectorField gradStressTerm(patch.size(), Zero);
        if (includeGradStressTerm_)
        {
            // Terms corresponding to contributions from converting delta to
            // thetas are added through the corresponding adjoint boundary
            // conditions instead of grabbing contributions from the objective
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
        vectorField pressureTerm(patch.size(), Zero);
        if (includePressureTerm_)
        {
            pressureTerm =
            (
                (nf*pa.boundaryField()[patchI])
                & U.boundaryField()[patchI].snGrad()
            )* nf;
        }

        PtrList<objective>& functions
            (objectiveManager_.getObjectiveFunctions());

        // Term from objectives including x directly (e.g. moments)
        vectorField dxdbMultiplierTot(pressureTerm.size(), Zero);
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
        wallFaceSensVecPtr_()[patchI] +=
        (
            stressTerm
          + gradStressTerm
          + pressureTerm
          + adjointTMsensitivities[patchI]
          + dxdbMultiplierTot
        )*dt;
    }

    // Add the sensitivity part corresponding to changes of the normal vector
    // Computed at points and mapped to faces
    addGeometricSens();
}


void sensitivitySurface::assembleSensitivities()
{
    // Update geometric fields for use by external users
    if (writeGeometricInfo_)
    {
        for (const label patchI : sensitivityPatchIDs_)
        {
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

    // Solve extra equations if necessary
    // Solved using accumulated sources over time
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


    // Project to normal face vector
    label nPassedFaces(0);
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> tnf(patch.nf());
        const vectorField& nf = tnf();

        // Distance related terms
        if (includeDistance_)
        {
            wallFaceSensVecPtr_()[patchI] += distanceSensPtr()[patchI];
        }

        // Mesh movement related terms
        if (includeMeshMovement_)
        {
            wallFaceSensVecPtr_()[patchI] += meshMovementSensPtr()[patchI];
        }

        if (includeSurfaceArea_)
        {
            wallFaceSensVecPtr_()[patchI] *= patch.magSf();
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

    // Smooth sensitivities if needed
    if (smoothSensitivities_)
    {
        smoothSensitivities();
    }
}


void sensitivitySurface::clearSensitivities()
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
    // Reset sensitivity fields
    adjointSensitivity::clearSensitivities();
    shapeSensitivitiesBase::clearSensitivities();
}


autoPtr<adjointEikonalSolver>& sensitivitySurface::getAdjointEikonalSolver()
{
    return eikonalSolver_;
}


void sensitivitySurface::write(const word& baseName)
{
    setSuffixName();
    adjointSensitivity::write();
    shapeSensitivitiesBase::write();

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
