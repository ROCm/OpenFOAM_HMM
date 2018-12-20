/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "solarLoad.H"
#include "surfaceFields.H"
#include "vectorList.H"
#include "addToRunTimeSelectionTable.H"
#include "boundaryRadiationProperties.H"
#include "gravityMeshObject.H"
#include "cyclicAMIPolyPatch.H"
#include "mappedPatchBase.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(solarLoad, 0);
        addToRadiationRunTimeSelectionTables(solarLoad);
    }
}

const Foam::word Foam::radiation::solarLoad::viewFactorWalls = "viewFactorWall";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::radiation::solarLoad::updateHitFaces()
{
    if (hitFaces_.empty())
    {
         hitFaces_.reset(new faceShading(mesh_, solarCalc_.direction()));
         return true;
    }
    else
    {
        switch (solarCalc_.sunDirectionModel())
        {
            case solarCalculator::mSunDirConstant:
            {
                return false;
                break;
            }
            case solarCalculator::mSunDirTracking:
            {
                label updateIndex = label
                (
                    mesh_.time().value()/solarCalc_.sunTrackingUpdateInterval()
                );

                if (updateIndex > updateTimeIndex_)
                {
                    Info << "Updating Sun position..." << endl;
                    updateTimeIndex_ = updateIndex;
                    solarCalc_.correctSunDirection();
                    hitFaces_->direction() = solarCalc_.direction();
                    hitFaces_->correct();
                    return true;
                    break;
                }
            }
        }
    }

    return false;
}


void Foam::radiation::solarLoad::updateAbsorptivity
(
    const labelHashSet& includePatches
)
{
    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(mesh_);

    for (const label patchID : includePatches)
    {
        absorptivity_[patchID].setSize(nBands_);
        for (label bandI = 0; bandI < nBands_; bandI++)
        {
            absorptivity_[patchID][bandI] =
                boundaryRadiation.absorptivity(patchID, bandI);
        }
    }
}


void Foam::radiation::solarLoad::updateDirectHitRadiation
(
    const labelList& hitFacesId,
    const labelHashSet& includeMappedPatchBasePatches
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const scalarField& V = mesh_.V();
    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    forAll(hitFacesId, i)
    {
        const label faceI = hitFacesId[i];
        label patchID = patches.whichPatch(faceI);
        const polyPatch& pp = patches[patchID];
        const label localFaceI = faceI - pp.start();
        const vector qPrim = solarCalc_.directSolarRad()*solarCalc_.direction();

        if (includeMappedPatchBasePatches[patchID])
        {
            const vectorField n = pp.faceNormals();

            for (label bandI = 0; bandI < nBands_; bandI++)
            {
                qrBf[patchID][localFaceI] +=
                    (qPrim & n[localFaceI])
                  * spectralDistribution_[bandI]
                  * absorptivity_[patchID][bandI]()[localFaceI];
            }
        }
        else
        {
            const vectorField& sf = mesh_.Sf().boundaryField()[patchID];
            const label cellI = pp.faceCells()[localFaceI];

            for (label bandI = 0; bandI < nBands_; bandI++)
            {
                Ru_[cellI] +=
                    (qPrim & sf[localFaceI])
                  * spectralDistribution_[bandI]
                  * absorptivity_[patchID][bandI]()[localFaceI]
                  / V[cellI];
            }
        }
    }
}


void Foam::radiation::solarLoad::updateSkyDiffusiveRadiation
(
    const labelHashSet& includePatches,
    const labelHashSet& includeMappedPatchBasePatches
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const scalarField& V = mesh_.V();
    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    switch(solarCalc_.sunLoadModel())
    {
        case solarCalculator::mSunLoadFairWeatherConditions:
        case solarCalculator::mSunLoadTheoreticalMaximum:
        {
            for (const label patchID : includePatches)
            {
                const polyPatch& pp = patches[patchID];
                const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

                const vectorField n = pp.faceNormals();
                const labelList& cellIds = pp.faceCells();

                forAll(n, faceI)
                {
                    const scalar cosEpsilon(verticalDir_ & -n[faceI]);

                    scalar Ed(0.0);
                    scalar Er(0.0);
                    const scalar cosTheta(solarCalc_.direction() & -n[faceI]);

                    {
                        // Above the horizon
                        if (cosEpsilon == 0.0)
                        {
                            // Vertical walls
                            scalar Y(0);

                            if (cosTheta > -0.2)
                            {
                                Y = 0.55+0.437*cosTheta + 0.313*sqr(cosTheta);
                            }
                            else
                            {
                                Y = 0.45;
                            }

                            Ed = solarCalc_.C()*Y*solarCalc_.directSolarRad();
                        }
                        else
                        {
                            //Other than vertical walls
                            Ed =
                                solarCalc_.C()
                              * solarCalc_.directSolarRad()
                              * (1.0 + cosEpsilon)/2.0;
                        }

                        // Ground reflected
                        Er =
                            solarCalc_.directSolarRad()
                          * (solarCalc_.C() + Foam::sin(solarCalc_.beta()))
                          * solarCalc_.groundReflectivity()
                          * (1.0 - cosEpsilon)/2.0;
                    }

                    const label cellI = cellIds[faceI];
                    if (includeMappedPatchBasePatches[patchID])
                    {
                        for (label bandI = 0; bandI < nBands_; bandI++)
                        {
                            qrBf[patchID][faceI] +=
                                (Ed + Er)
                              * spectralDistribution_[bandI]
                              * absorptivity_[patchID][bandI]()[faceI];
                        }
                    }
                    else
                    {
                        for (label bandI = 0; bandI < nBands_; bandI++)
                        {
                            Ru_[cellI] +=
                                (Ed + Er)
                              * spectralDistribution_[bandI]
                              * absorptivity_[patchID][bandI]()[faceI]
                              * sf[faceI]/V[cellI];
                        }
                    }
                }
            }
        }
        break;

        case solarCalculator::mSunLoadConstant:
        {
            for (const label patchID : includePatches)
            {
                const polyPatch& pp = patches[patchID];
                const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

                const labelList& cellIds = pp.faceCells();
                forAll(pp, faceI)
                {
                    const label cellI = cellIds[faceI];
                    if (includeMappedPatchBasePatches[patchID])
                    {
                        for (label bandI = 0; bandI < nBands_; bandI++)
                        {
                            qrBf[patchID][faceI] +=
                                solarCalc_.diffuseSolarRad()
                              * spectralDistribution_[bandI]
                              * absorptivity_[patchID][bandI]()[faceI];
                        }
                    }
                    else
                    {
                        for (label bandI = 0; bandI < nBands_; bandI++)
                        {
                            Ru_[cellI] +=
                                (
                                    spectralDistribution_[bandI]
                                  * absorptivity_[patchID][bandI]()[faceI]
                                  * solarCalc_.diffuseSolarRad()
                                )*sf[faceI]/V[cellI];
                        }
                    }
                }
            }
            break;
        }
    }
}


void Foam::radiation::solarLoad::initialise(const dictionary& coeffs)
{
    if (coeffs.readIfPresent("gridUp", verticalDir_))
    {
         verticalDir_.normalise();
    }
    else
    {
        const uniformDimensionedVectorField& g =
            meshObjects::gravity::New(mesh_.time());
        verticalDir_ = (-g/mag(g)).value();
    }

    includePatches_ = mesh_.boundaryMesh().indices(viewFactorWalls);

    coeffs.readEntry("useVFbeamToDiffuse", useVFbeamToDiffuse_);

    coeffs.readEntry("spectralDistribution", spectralDistribution_);

    spectralDistribution_ =
        spectralDistribution_/sum(spectralDistribution_);

    nBands_ = spectralDistribution_.size();

    if (useVFbeamToDiffuse_)
    {
        map_.reset
        (
            new IOmapDistribute
            (
                IOobject
                (
                    "mapDist",
                    mesh_.facesInstance(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }

    coeffs.readIfPresent("solidCoupled", solidCoupled_);
    coeffs.readIfPresent("wallCoupled", wallCoupled_);
    coeffs.readIfPresent("updateAbsorptivity", updateAbsorptivity_);
}


void Foam::radiation::solarLoad::calculateQdiff
(
    const labelHashSet& includePatches,
    const labelHashSet& includeMappedPatchBasePatches
)
{
    scalarListIOList FmyProc
    (
        IOobject
        (
            "F",
            mesh_.facesInstance(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );


    if (finalAgglom_.size() > 0 && coarseMesh_.empty())
    {
        coarseMesh_.reset
        (
            new singleCellFvMesh
            (
                IOobject
                (
                    "coarse:" + mesh_.name(),
                    mesh_.polyMesh::instance(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                finalAgglom_
            )
        );
    }

    label nLocalVFCoarseFaces = 0;
    forAll(includePatches_, i)
    {
        const label patchI = includePatches_[i];
        nLocalVFCoarseFaces += coarseMesh_->boundaryMesh()[patchI].size();
    }

    label totalFVNCoarseFaces = nLocalVFCoarseFaces;
    reduce(totalFVNCoarseFaces, sumOp<label>());

    // Calculate weighted absorptivity on coarse patches
    List<scalar> localCoarseRave(nLocalVFCoarseFaces);
    List<scalar> localCoarsePartialArea(nLocalVFCoarseFaces);
    List<vector> localCoarseNorm(nLocalVFCoarseFaces);

    scalarField compactCoarseRave(map_->constructSize(), 0.0);
    scalarField compactCoarsePartialArea(map_->constructSize(), 0.0);
    vectorList compactCoarseNorm(map_->constructSize(), Zero);

    const boundaryRadiationProperties& boundaryRadiation =
        boundaryRadiationProperties::New(mesh_);

    coarseToFine_.setSize(includePatches_.size());

    const labelList& hitFacesId = hitFaces_->rayStartFaces();

    label startI = 0;
    label compactI = 0;
    forAll(includePatches_, i)
    {
        const label patchID = includePatches_[i];
        const polyPatch& pp = mesh_.boundaryMesh()[patchID];

        const polyPatch& cpp = coarseMesh_->boundaryMesh()[patchID];

        const labelList& agglom = finalAgglom_[patchID];
        //if (pp.size() > 0)
        if (agglom.size() > 0)
        {
            label nAgglom = max(agglom) + 1;
            coarseToFine_[i] = invertOneToMany(nAgglom, agglom);
        }

        scalarField r(pp.size(), 0.0);
        for (label bandI = 0; bandI < nBands_; bandI++)
        {
            const tmp<scalarField> tr =
                spectralDistribution_[bandI]
               *boundaryRadiation.reflectivity(patchID, bandI);
            r += tr();
        }

        scalarList Rave(cpp.size(), 0.0);
        scalarList area(cpp.size(), 0.0);

        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

        const labelList& coarsePatchFace =
            coarseMesh_->patchFaceMap()[patchID];

        forAll(cpp, coarseI)
        {
            const label coarseFaceID = coarsePatchFace[coarseI];
            const labelList& fineFaces = coarseToFine_[i][coarseFaceID];

            UIndirectList<scalar> fineSf(sf, fineFaces);
            scalar fineArea = sum(fineSf());

            scalar fullArea = 0.0;
            forAll(fineFaces, j)
            {
                label faceI = fineFaces[j];
                label globalFaceI = faceI + pp.start();

                if (hitFacesId.found(globalFaceI))
                {
                    fullArea += sf[faceI];
                }
                Rave[coarseI] += (r[faceI]*sf[faceI])/fineArea;
            }
            localCoarsePartialArea[compactI++] = fullArea/fineArea;
        }

        SubList<scalar>
        (
            localCoarseRave,
            Rave.size(),
            startI
        ) = Rave;


        const vectorList coarseNSf = cpp.faceNormals();
        SubList<vector>
        (
            localCoarseNorm,
            cpp.size(),
            startI
        ) = coarseNSf;
        startI += cpp.size();
    }


    SubList<scalar>(compactCoarsePartialArea, nLocalVFCoarseFaces) =
        localCoarsePartialArea;

    SubList<scalar>(compactCoarseRave, nLocalVFCoarseFaces) =
        localCoarseRave;

    SubList<vector>(compactCoarseNorm, nLocalVFCoarseFaces) =
        localCoarseNorm;

    map_->distribute(compactCoarsePartialArea);
    map_->distribute(compactCoarseRave);
    map_->distribute(compactCoarseNorm);


    // Calculate coarse hitFaces and Sun direct hit heat fluxes
    scalarList localqDiffusive(nLocalVFCoarseFaces, 0.0);

    label locaFaceI = 0;
    forAll(includePatches_, i)
    {
        const label patchID = includePatches_[i];
        const polyPatch& pp = coarseMesh_->boundaryMesh()[patchID];
        const polyPatch& ppf = mesh_.boundaryMesh()[patchID];

        const labelList& coarsePatchFace = coarseMesh_->patchFaceMap()[patchID];
        const scalarField& sf = mesh_.magSf().boundaryField()[patchID];

        scalarField a(ppf.size(), 0.0);
        for (label bandI = 0; bandI < nBands_; bandI++)
        {
            const tmp<scalarField> ta =
                spectralDistribution_[bandI]
              * absorptivity_[patchID][bandI]();
            a += ta();
        }

        forAll(pp, coarseI)
        {
            const label coarseFaceID = coarsePatchFace[coarseI];
            const labelList& fineFaces = coarseToFine_[i][coarseFaceID];
            UIndirectList<scalar> fineSf(sf, fineFaces);
            scalar fineArea = sum(fineSf());

            scalar aAve = 0.0;
            forAll(fineFaces, j)
            {
                label faceI = fineFaces[j];
                aAve += (a[faceI]*sf[faceI])/fineArea;
            }

            const scalarList& vf = FmyProc[locaFaceI];
            const labelList& compactFaces = visibleFaceFaces_[locaFaceI];

            forAll(compactFaces, j)
            {
                label compactI = compactFaces[j];

                localqDiffusive[locaFaceI] +=
                    compactCoarsePartialArea[compactI]
                  * aAve
                  * (solarCalc_.directSolarRad()*solarCalc_.direction())
                  & compactCoarseNorm[compactI]
                  * vf[j]
                  * compactCoarseRave[compactI];

            }
            locaFaceI++;
        }
    }

    volScalarField::Boundary& qsBf = qsecondRad_.boundaryFieldRef();

    // Fill qsecondRad_
    label compactId = 0;
    forAll(includePatches_, i)
    {
        const label patchID = includePatches_[i];
        const polyPatch& pp = coarseMesh_->boundaryMesh()[patchID];

        if (pp.size() > 0)
        {
            scalarField& qrp = qsBf[patchID];

            const labelList& coarsePatchFace =
                coarseMesh_->patchFaceMap()[patchID];

            forAll(pp, coarseI)
            {
                const label coarseFaceID = coarsePatchFace[coarseI];
                const labelList& fineFaces = coarseToFine_[i][coarseFaceID];
                forAll(fineFaces, k)
                {
                    label faceI = fineFaces[k];
                    qrp[faceI] = localqDiffusive[compactId];
                }
                compactId ++;
            }
        }
    }

    const scalarField& V = mesh_.V();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    for (const label patchID : includePatches)
    {
        const scalarField& qSecond = qsecondRad_.boundaryField()[patchID];
        if (includeMappedPatchBasePatches[patchID])
        {
            qrBf[patchID] += qSecond;
        }
        else
        {
            const polyPatch& pp = patches[patchID];
            const labelList& cellIds = pp.faceCells();
            const scalarField& sf = mesh_.magSf().boundaryField()[patchID];
            forAll(pp, faceI)
            {
                const label cellI = cellIds[faceI];
                Ru_[cellI] += qSecond[faceI]*sf[faceI]/V[cellI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::solarLoad::solarLoad(const volScalarField& T)
:
    radiationModel(typeName, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    coarseMesh_(),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qsecondRad_
    (
        IOobject
        (
            "qsecondRad",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    hitFaces_(),
    Ru_
    (
        IOobject
        (
            "Ru",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    ),
    solarCalc_(this->subDict(typeName + "Coeffs"), mesh_),
    verticalDir_(Zero),
    useVFbeamToDiffuse_(false),
    includePatches_(mesh_.boundary().size(), -1),
    coarseToFine_(),
    nBands_(1),
    spectralDistribution_(nBands_),
    map_(),
    visibleFaceFaces_
    (
        IOobject
        (
            "visibleFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    solidCoupled_(true),
    absorptivity_(mesh_.boundaryMesh().size()),
    updateAbsorptivity_(false),
    firstIter_(true),
    updateTimeIndex_(0)
{
    initialise(coeffs_);
}


Foam::radiation::solarLoad::solarLoad
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    coarseMesh_(),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qsecondRad_
    (
        IOobject
        (
            "qsecondRad",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    hitFaces_(),
    Ru_
    (
        IOobject
        (
            "Ru",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    ),
    solarCalc_(coeffs_, mesh_),
    verticalDir_(Zero),
    useVFbeamToDiffuse_(false),
    includePatches_(mesh_.boundary().size(), -1),
    coarseToFine_(),
    nBands_(1),
    spectralDistribution_(nBands_),
    map_(),
    visibleFaceFaces_
    (
        IOobject
        (
            "visibleFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    solidCoupled_(true),
    wallCoupled_(false),
    absorptivity_(mesh_.boundaryMesh().size()),
    updateAbsorptivity_(false),
    firstIter_(true),
    updateTimeIndex_(0)
{
    initialise(coeffs_);
}


Foam::radiation::solarLoad::solarLoad
(
    const dictionary& dict,
    const volScalarField& T,
    const word radWallFieldName
)
:
    radiationModel("none", T),
    finalAgglom_
    (
        IOobject
        (
            "finalAgglom",
            mesh_.facesInstance(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    coarseMesh_(),
    qr_
    (
        IOobject
        (
            radWallFieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    qsecondRad_
    (
        IOobject
        (
            "qsecondRad",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    hitFaces_(),
    Ru_
    (
        IOobject
        (
            "Ru",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    ),
    solarCalc_(dict, mesh_),
    verticalDir_(Zero),
    useVFbeamToDiffuse_(false),
    includePatches_(mesh_.boundary().size(), -1),
    coarseToFine_(),
    nBands_(1),
    spectralDistribution_(nBands_),
    map_(),
    visibleFaceFaces_
    (
        IOobject
        (
            "visibleFaceFaces",
            mesh_.facesInstance(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    ),
    solidCoupled_(true),
    wallCoupled_(false),
    absorptivity_(mesh_.boundaryMesh().size()),
    updateAbsorptivity_(false),
    firstIter_(true)
{
    initialise(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::solarLoad::~solarLoad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label  Foam::radiation::solarLoad::nBands() const
{
    return nBands_;
}


bool Foam::radiation::solarLoad::read()
{
    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::solarLoad::calculate()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    labelHashSet includePatches;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
        {
            includePatches.insert(patchI);
        }
    }

    labelHashSet includeMappedPatchBasePatches;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if
        (
            (isA<mappedPatchBase>(pp) && solidCoupled_)
         || (isA<wallPolyPatch>(pp) && wallCoupled_)
        )
        {
            includeMappedPatchBasePatches.insert(patchI);
        }
    }

    if (updateAbsorptivity_ || firstIter_)
    {
        updateAbsorptivity(includePatches);
    }

    bool facesChanged = updateHitFaces();
    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    if (facesChanged)
    {
        // Reset Ru and qr
        Ru_ = dimensionedScalar("Ru", dimMass/dimLength/pow3(dimTime), 0.0);
        qrBf = 0.0;

        // Add direct hit radiation
        const labelList& hitFacesId = hitFaces_->rayStartFaces();
        updateDirectHitRadiation(hitFacesId, includeMappedPatchBasePatches);

        // Add sky diffusive radiation
        updateSkyDiffusiveRadiation
        (
            includePatches,
            includeMappedPatchBasePatches
        );

        // Add indirect diffusive radiation
        if (useVFbeamToDiffuse_)
        {
            calculateQdiff(includePatches, includeMappedPatchBasePatches);
        }

        firstIter_ = false;
    }

    if (debug)
    {
        if (mesh_.time().writeTime())
        {
            Ru_.write();
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiation::solarLoad::Rp() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "Rp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar
        (
            dimMass/pow3(dimTime)/dimLength/pow4(dimTemperature),
            Zero
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::solarLoad::Ru() const
{
    return Ru_;
}


// ************************************************************************* //
