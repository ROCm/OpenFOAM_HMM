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

#include "ThermoSurfaceFilm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const dictionary& dict,
    CloudType& owner,
    const dimensionedVector& g
)
:
    SurfaceFilmModel<CloudType>(dict, owner, g, typeName),
    TFilmPatch_(0),
    cpFilmPatch_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::~ThermoSurfaceFilm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ThermoSurfaceFilm<CloudType>::active() const
{
    return true;
}


template<class CloudType>
bool Foam::ThermoSurfaceFilm<CloudType>::transferParcel
(
    const parcelType& p,
    const label patchI
)
{
    // Retrieve the film model from the owner database
    surfaceFilmModels::surfaceFilmModel& filmModel =
        const_cast<surfaceFilmModels::surfaceFilmModel&>
        (
            this->owner().db().objectRegistry::
                lookupObject<surfaceFilmModels::surfaceFilmModel>
                (
                    "surfaceFilmProperties"
                )
        );

    if (filmModel.isFilmPatch(patchI))
    {
        const polyPatch& pp = this->owner().mesh().boundaryMesh()[patchI];

        const label faceI = pp.whichFace(p.face());

        // Patch face normal
        const vector& nf = pp.faceNormals()[faceI];

        // Relative parcel velocity
        const vector Urel =
            p.U() - this->owner().U().boundaryField()[patchI][faceI];

        // Parcel mass
        const scalar m = p.nParticle()*p.mass();

        // Add the particle properties as sources to the film model
        filmModel.addSources
        (
            patchI,
            faceI,
            m,                              // mass
            m*(Urel - nf*(Urel & nf)),      // tangential momentum
            m*mag(Urel & nf),               // impingement pressure
            m*p.hs()                        // energy
        );

        if (debug)
        {
            Info<< "ThermoSurfaceFilm<CloudType>::transferParcel:" << nl
                << "    Effective increase in film height = "
                << p.nParticle()*p.volume()/mag(pp.faceAreas()[faceI]) << endl;
        }

        this->nParcelsTransferred()++;

        // Flag to remove parcel p from owner cloud
        return true;
    }
    else
    {
        return false;
    }
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::cacheFilmFields
(
    const label filmPatchI,
    const mapDistribute& distMap,
    const surfaceFilmModels::surfaceFilmModel& filmModel
)
{
    SurfaceFilmModel<CloudType>::cacheFilmFields
    (
        filmPatchI,
        distMap,
        filmModel
    );

    TFilmPatch_ = filmModel.Ts().boundaryField()[filmPatchI];
    distMap.distribute(TFilmPatch_);

    cpFilmPatch_ = filmModel.cp().boundaryField()[filmPatchI];
    distMap.distribute(cpFilmPatch_);
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFaceI
) const
{
    SurfaceFilmModel<CloudType>::setParcelProperties(p, filmFaceI);

    // Set parcel properties
    p.T() = TFilmPatch_[filmFaceI];
    p.cp() = cpFilmPatch_[filmFaceI];
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::info(Ostream& os) const
{
    os  << "    Parcels transferred to film     = "
        << returnReduce(this->nParcelsTransferred(), sumOp<label>()) << nl
        << "    Number of film parcels added    = "
        << returnReduce(this->nParcelsInjected(), sumOp<label>()) << nl;
}


// ************************************************************************* //
