/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
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
#include "surfaceFilmModel.H"
#include "directMappedWallPolyPatch.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    g_(dimensionedVector("zero", dimAcceleration, vector::zero)),
    coeffDict_(dictionary::null),
    active_(false),
    injectorCellsPatch_(0),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    UFilmPatch_(0),
    rhoFilmPatch_(0),
    nParcelsTransferred_(0),
    nParcelsInjected_(0)
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const dictionary& dict,
    CloudType& owner,
    const dimensionedVector& g,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    g_(g),
    coeffDict_(dict.subDict(type + "Coeffs")),
    active_(true),
    injectorCellsPatch_(0),
    massParcelPatch_(0),
    diameterParcelPatch_(0),
    UFilmPatch_(0),
    rhoFilmPatch_(0),
    nParcelsTransferred_(0),
    nParcelsInjected_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::~SurfaceFilmModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class TrackData>
void Foam::SurfaceFilmModel<CloudType>::inject(TrackData& td)
{
    if (!active_)
    {
        return;
    }

    // Retrieve the film model from the owner database
    const surfaceFilmModels::surfaceFilmModel& filmModel =
        this->owner().db().objectRegistry::lookupObject
        <surfaceFilmModels::surfaceFilmModel>
        (
            "surfaceFilmProperties"
        );

    const labelList& filmPatches = filmModel.filmBottomPatchIDs();
    const labelList& primaryPatches = filmModel.primaryPatchIDs();

    forAll(filmPatches, i)
    {
        const label primaryPatchI = primaryPatches[i];
        const directMappedWallPolyPatch& wpp =
            refCast<const directMappedWallPolyPatch>
            (
                 this->owner().mesh().boundaryMesh()[primaryPatchI]
            );

        injectorCellsPatch_ = wpp.faceCells();

        const label filmPatchI = filmPatches[i];
        const mapDistribute& distMap = wpp.map();
        cacheFilmFields(filmPatchI, distMap, filmModel);

        forAll(injectorCellsPatch_, j)
        {
            if (diameterParcelPatch_[j] > 0)
            {
                const label cellI = injectorCellsPatch_[j];
                const point& pos = this->owner().mesh().C()[cellI];

                // Create a new parcel
                typename CloudType::parcelType* pPtr =
                    new typename CloudType::parcelType(td.cloud(), pos, cellI);
                setParcelProperties(*pPtr, j);

                // Check new parcel properties
//                td.cloud().checkParcelProperties(*pPtr, 0.0, true);
                td.cloud().checkParcelProperties(*pPtr, 0.0, false);

                // Add the new parcel to the cloud
                td.cloud().addParticle(pPtr);

                nParcelsInjected_++;
            }
        }
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::cacheFilmFields
(
    const label filmPatchI,
    const mapDistribute& distMap,
    const surfaceFilmModels::surfaceFilmModel& filmModel
)
{
    massParcelPatch_ = filmModel.massForPrimary().boundaryField()[filmPatchI];
    distMap.distribute(massParcelPatch_);

    diameterParcelPatch_ =
        filmModel.diametersForPrimary().boundaryField()[filmPatchI];
    distMap.distribute(diameterParcelPatch_);

    UFilmPatch_ = filmModel.U().boundaryField()[filmPatchI];
    distMap.distribute(UFilmPatch_);

    rhoFilmPatch_ = filmModel.rho().boundaryField()[filmPatchI];
    distMap.distribute(rhoFilmPatch_);
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFaceI
) const
{
    // Set parcel properties
    scalar vol = mathematical::pi/6.0*pow3(diameterParcelPatch_[filmFaceI]);
    p.d() = diameterParcelPatch_[filmFaceI];
    p.U() = UFilmPatch_[filmFaceI];
    p.rho() = rhoFilmPatch_[filmFaceI];

    p.nParticle() = massParcelPatch_[filmFaceI]/p.rho()/vol;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceFilmModelNew.C"

// ************************************************************************* //
