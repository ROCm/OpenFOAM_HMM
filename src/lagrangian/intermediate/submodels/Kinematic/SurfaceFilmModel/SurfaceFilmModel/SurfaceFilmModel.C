/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "SurfaceFilmModel.H"
#include "mathematicalConstants.H"
#include "surfaceFilmRegionModel.H"
#include "liquidFilmBase.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    g_(owner.g()),
    ejectedParcelType_(0),
    injectionOffset_(1.1),
    minDiameter_(0),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    UFilmPatch_(0),
    rhoFilmPatch_(0),
    deltaFilmPatch_(0),
    nParcelsTransferred_(0),
    nParcelsInjected_(0),
    totalMassTransferred_(0)
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    g_(owner.g()),
    ejectedParcelType_
    (
        this->coeffDict().template getOrDefault<label>("ejectedParcelType", -1)
    ),
    injectionOffset_
    (
        this->coeffDict().template getOrDefault<scalar>("injectionOffset", 1.1)
    ),
    minDiameter_
    (
        this->coeffDict().template getOrDefault<scalar>("minDiameter", -1)
    ),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    UFilmPatch_(0),
    rhoFilmPatch_(0),
    deltaFilmPatch_(owner.mesh().boundary().size()),
    nParcelsTransferred_(0),
    nParcelsInjected_(0),
    totalMassTransferred_()
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const SurfaceFilmModel<CloudType>& sfm
)
:
    CloudSubModelBase<CloudType>(sfm),
    g_(sfm.g_),
    ejectedParcelType_(sfm.ejectedParcelType_),
    injectionOffset_(sfm.injectionOffset_),
    minDiameter_(sfm.minDiameter_),
    massParcelPatch_(sfm.massParcelPatch_),
    diameterParcelPatch_(sfm.diameterParcelPatch_),
    UFilmPatch_(sfm.UFilmPatch_),
    rhoFilmPatch_(sfm.rhoFilmPatch_),
    deltaFilmPatch_(sfm.deltaFilmPatch_),
    nParcelsTransferred_(sfm.nParcelsTransferred_),
    nParcelsInjected_(sfm.nParcelsInjected_),
    totalMassTransferred_(sfm.totalMassTransferred_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::~SurfaceFilmModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class CloudTrackType>
void Foam::SurfaceFilmModel<CloudType>::injectParticles
(
    const label primaryPatchi,
    const labelList& injectorCellsPatch,
    CloudTrackType& cloud
)
{
    const fvMesh& mesh = this->owner().mesh();
    const vectorField& Cf = mesh.C().boundaryField()[primaryPatchi];
    const vectorField& Sf = mesh.Sf().boundaryField()[primaryPatchi];
    const scalarField& magSf =
        mesh.magSf().boundaryField()[primaryPatchi];

    forAll(injectorCellsPatch, j)
    {
        if (diameterParcelPatch_[j] > 0)
        {
            const label celli = injectorCellsPatch[j];

            const scalar offset =
                max
                (
                    diameterParcelPatch_[j],
                    deltaFilmPatch_[primaryPatchi][j]
                );
            const point pos = Cf[j] - injectionOffset_*offset*Sf[j]/magSf[j];

            // Create a new parcel
            parcelType* pPtr =
                new parcelType(this->owner().pMesh(), pos, celli);

            // Check/set new parcel thermo properties
            cloud.setParcelThermoProperties(*pPtr, 0.0);

            setParcelProperties(*pPtr, j);

            if (pPtr->nParticle() > 0.001)
            {
                // Check new parcel properties
                cloud.checkParcelProperties(*pPtr, 0.0, false);

                // Add the new parcel to the cloud
                cloud.addParticle(pPtr);

                nParcelsInjected_++;
            }
            else
            {
                // TODO: cache mass and re-distribute?
                delete pPtr;
            }
        }
    }
}

template<class CloudType>
template<class TrackCloudType>
void Foam::SurfaceFilmModel<CloudType>::inject(TrackCloudType& cloud)
{
    if (!this->active())
    {
        return;
    }

    const fvMesh& mesh = this->owner().mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Retrieve the film model from the owner database
    const regionModels::surfaceFilmModels::surfaceFilmRegionModel* filmModel =
        mesh.time().objectRegistry::template findObject
        <regionModels::surfaceFilmModels::surfaceFilmRegionModel>
        (
            "surfaceFilmProperties"
        );

    // Check the singleLayer type of films
    if (filmModel && filmModel->active())
    {

        const labelList& filmPatches = filmModel->intCoupledPatchIDs();
        const labelList& primaryPatches = filmModel->primaryPatchIDs();

        forAll(filmPatches, i)
        {
            const label filmPatchi = filmPatches[i];
            const label primaryPatchi = primaryPatches[i];

            const labelList& injectorCellsPatch = pbm[primaryPatchi].faceCells();

            cacheFilmFields(filmPatchi, primaryPatchi, *filmModel);

            injectParticles(primaryPatchi, injectorCellsPatch, cloud);
        }
    }

    // Check finite area films
    wordList names =
        mesh.time().objectRegistry::template
            sortedNames<regionModels::regionFaModel>();

    forAll (names, i)
    {
        const regionModels::regionFaModel* regionFa =
            mesh.time().objectRegistry::template cfindObject
            <
                regionModels::regionFaModel
            >(names[i]);

        // Check that it is a type areaFilm
        if (regionFa && isA<areaFilm>(*regionFa) && regionFa->active())
        {
            areaFilm& film =
                const_cast<areaFilm&>(refCast<const areaFilm>(*regionFa));

            const label patchId = regionFa->patchID();

            const labelList& injectorCellsPatch = pbm[patchId].faceCells();

            cacheFilmFields(patchId, film);

            injectParticles(patchId, injectorCellsPatch, cloud);

            forAll(injectorCellsPatch, facei)
            {
                if (diameterParcelPatch_[facei] > 0)
                {
                    film.addSources
                    (
                        patchId,
                        facei,
                       -massParcelPatch_[facei],    // mass
                        Zero,                       // tangential momentum
                        Zero,                       // impingement
                        Zero                        // energy
                    );
                }
            }
        }
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const label primaryPatchi,
    const regionModels::surfaceFilmModels::surfaceFilmRegionModel& filmModel
)
{
    massParcelPatch_ = filmModel.cloudMassTrans().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, massParcelPatch_);

    diameterParcelPatch_ =
        filmModel.cloudDiameterTrans().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, diameterParcelPatch_, maxEqOp<scalar>());

    UFilmPatch_ = filmModel.Us().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, UFilmPatch_);

    rhoFilmPatch_ = filmModel.rho().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, rhoFilmPatch_);

    deltaFilmPatch_[primaryPatchi] =
        filmModel.delta().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, deltaFilmPatch_[primaryPatchi]);
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const regionModels::areaSurfaceFilmModels::liquidFilmBase& filmModel
)
{
    const volSurfaceMapping& map = filmModel.region().vsm();

    massParcelPatch_.setSize(filmModel.Uf().size(), Zero);

    const scalarField& massParcelPatch =
        filmModel.cloudMassTrans().boundaryField()[filmPatchi];

    map.mapToField(massParcelPatch, massParcelPatch_);

    const scalarField&  diameterParcelPatch =
        filmModel.cloudDiameterTrans().boundaryField()[filmPatchi];

    diameterParcelPatch_.setSize(filmModel.Uf().size(), Zero);
    map.mapToField(diameterParcelPatch, diameterParcelPatch_);

    UFilmPatch_.setSize(filmModel.Uf().size(), vector::zero);
    map.mapToField(filmModel.Uf(), UFilmPatch_);

    rhoFilmPatch_.setSize(UFilmPatch_.size(), Zero);
    map.mapToField(filmModel.rho(), rhoFilmPatch_);

    deltaFilmPatch_[filmPatchi].setSize(UFilmPatch_.size(), Zero);
    map.mapToField(filmModel.h(), deltaFilmPatch_[filmPatchi]);
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFacei
) const
{
    // Set parcel properties
    scalar vol = mathematical::pi/6.0*pow3(diameterParcelPatch_[filmFacei]);
    p.d() = diameterParcelPatch_[filmFacei];
    p.U() = UFilmPatch_[filmFacei];
    p.rho() = rhoFilmPatch_[filmFacei];

    p.nParticle() = massParcelPatch_[filmFacei]/p.rho()/vol;

    if (minDiameter_ != -1)
    {
        if (p.d() < minDiameter_)
        {
            p.nParticle() = 0;
        }
    }

    if (ejectedParcelType_ >= 0)
    {
        p.typeId() = ejectedParcelType_;
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::info(Ostream& os)
{
    label nTrans0 =
        this->template getModelProperty<label>("nParcelsTransferred");

    label nInject0 =
        this->template getModelProperty<label>("nParcelsInjected");

    scalar massTransferred0 =
        this->template getModelProperty<scalar>("massTransferred");

    label nTransTotal =
        nTrans0 + returnReduce(nParcelsTransferred_, sumOp<label>());

    label nInjectTotal =
        nInject0 + returnReduce(nParcelsInjected_, sumOp<label>());

    scalar massTransferredTotal =
        massTransferred0 + returnReduce(totalMassTransferred_, sumOp<scalar>());


    os  << "    Surface film:" << nl
        << "      - parcels absorbed            = " << nTransTotal << nl
        << "      - mass absorbed               = " << massTransferredTotal << nl
        << "      - parcels ejected             = " << nInjectTotal << endl;

    if (this->writeTime())
    {
        this->setModelProperty("nParcelsTransferred", nTransTotal);
        this->setModelProperty("nParcelsInjected", nInjectTotal);
        this->setModelProperty("massTransferred", massTransferredTotal);

        nParcelsTransferred_ = 0;
        nParcelsInjected_ = 0;
        totalMassTransferred_ = 0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceFilmModelNew.C"

// ************************************************************************* //
