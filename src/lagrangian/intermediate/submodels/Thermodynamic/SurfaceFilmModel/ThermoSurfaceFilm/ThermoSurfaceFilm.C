/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "ThermoSurfaceFilm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const dictionary& dict,
    CloudType& owner
)
:
    KinematicSurfaceFilm<CloudType>(dict, owner, typeName,  false),
    thermo_
    (
        owner.db().objectRegistry::template lookupObject<SLGThermo>("SLGThermo")
    ),
    TFilmPatch_(0),
    CpFilmPatch_(0)
{}


template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const ThermoSurfaceFilm<CloudType>& sfm
)
:
    KinematicSurfaceFilm<CloudType>(sfm, false),
    thermo_(sfm.thermo_),
    TFilmPatch_(sfm.TFilmPatch_),
    CpFilmPatch_(sfm.CpFilmPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ThermoSurfaceFilm<CloudType>::transferParcel
(
    parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    const label patchi = pp.index();

    this->initFilmModels();

    bool bInteraction(false);

    // Check the singleLayer film models
    if (this->filmModel_)
    {
        if (this->filmModel_->isRegionPatch(patchi))
        {
            const label facei = pp.whichFace(p.face());

            switch (this->interactionType_)
            {
                case KinematicSurfaceFilm<CloudType>::itBounce:
                {
                    this->bounceInteraction(p, pp, facei, keepParticle);

                    break;
                }
                case KinematicSurfaceFilm<CloudType>::itAbsorb:
                {
                    const scalar m = p.nParticle()*p.mass();

                    this->absorbInteraction //<regionFilm>
                        (*(this->filmModel_), p, pp, facei, m, keepParticle);

                    break;
                }
                case KinematicSurfaceFilm<CloudType>::itSplashBai:
                {
                    // Local pressure
                    const scalar pc = thermo_.thermo().p()[p.cell()];
                    const liquidProperties& liq = thermo_.liquids().properties()[0];
                    const scalar sigma = liq.sigma(pc, p.T());
                    const scalar mu = liq.mu(pc, p.T());

                    bool dry = this->deltaFilmPatch_[patchi][facei] < this->deltaWet_;

                    if (dry)
                    {
                        this->drySplashInteraction //<CloudType, regionFilm>
                            (*(this->filmModel_), sigma, mu, p, pp, facei, keepParticle);
                    }
                    else
                    {
                        this->wetSplashInteraction //<regionFilm>
                            (*(this->filmModel_), sigma, mu, p, pp, facei, keepParticle);
                    }

                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown interaction type enumeration"
                        << abort(FatalError);
                }
            }

            // Transfer parcel/parcel interactions complete
            bInteraction = true;
        }
    }

    for (areaFilm& film : this->areaFilms_)
    {
        if (patchi == film.patchID())
        {
            const label facei = pp.whichFace(p.face());

            switch (this->interactionType_)
            {
                // It only supports absorp model
                case KinematicSurfaceFilm<CloudType>::itAbsorb:
                {
                    const scalar m = p.nParticle()*p.mass();

                    this->absorbInteraction //<areaFilm>
                    (
                        film, p, pp, facei, m, keepParticle
                    );
                    break;
                }
                case KinematicSurfaceFilm<CloudType>::itBounce:
                {
                    this->bounceInteraction(p, pp, facei, keepParticle);

                    break;
                }
                case KinematicSurfaceFilm<CloudType>::itSplashBai:
                {
                    // Local pressure
                    const scalar pc = thermo_.thermo().p()[p.cell()];
                    const liquidProperties& liq = thermo_.liquids().properties()[0];
                    const scalar sigma = liq.sigma(pc, p.T());
                    const scalar mu = liq.mu(pc, p.T());

                    bool dry = film.h()[facei] < this->deltaWet_;

                    if (dry)
                    {
                        this->drySplashInteraction //<areaFilm>
                            (film, sigma, mu, p, pp, facei, keepParticle);
                    }
                    else
                    {
                        this->wetSplashInteraction //<areaFilm>
                            (film, sigma, mu, p, pp, facei, keepParticle);
                    }

                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown interaction type enumeration"
                        << abort(FatalError);
                }
            }
            // Transfer parcel/parcel interactions complete
            bInteraction = true;
        }
    }

    // Parcel not interacting with film
    return bInteraction;
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const label primaryPatchi,
    const regionModels::surfaceFilmModels::surfaceFilmRegionModel& filmModel
)
{
    KinematicSurfaceFilm<CloudType>::cacheFilmFields
    (
        filmPatchi,
        primaryPatchi,
        filmModel
    );

    TFilmPatch_ = filmModel.Ts().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, TFilmPatch_);

    CpFilmPatch_ = filmModel.Cp().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, CpFilmPatch_);
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const areaFilm& filmModel
)
{
    KinematicSurfaceFilm<CloudType>::cacheFilmFields
    (
        filmPatchi,
        filmModel
    );
    const volSurfaceMapping& map = filmModel.region().vsm();

    TFilmPatch_.setSize(filmModel.Tf().size(), Zero);
    map.mapToField(filmModel.Tf(), TFilmPatch_);
    CpFilmPatch_.setSize(filmModel.Tf().size(), Zero);
    map.mapToField(filmModel.Cp(), CpFilmPatch_);
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFacei
) const
{
    KinematicSurfaceFilm<CloudType>::setParcelProperties(p, filmFacei);

    // Set parcel properties
    p.T() = TFilmPatch_[filmFacei];
    p.Cp() = CpFilmPatch_[filmFacei];
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::info(Ostream& os)
{
    KinematicSurfaceFilm<CloudType>::info(os);
}


// ************************************************************************* //
