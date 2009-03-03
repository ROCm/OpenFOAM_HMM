/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "EntrainmentInjection.H"
#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::EntrainmentInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if (sum(nParcelsInjected_) < nParcelsPerInjector_*positions_.size())
    {
        return positions_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::EntrainmentInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if (sum(nParcelsInjected_) < nParcelsPerInjector_*positions_.size())
    {
        return this->volumeTotal_/nParcelsPerInjector_;
    }
    else
    {
        return 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EntrainmentInjection<CloudType>::EntrainmentInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    c_(readScalar(this->coeffDict().lookup("c"))),
    rhoc_
    (
        owner.db().objectRegistry::lookupObject<volScalarField>
        (
            this->coeffDict().lookup("rhoc")
        )
    ),
    rhol_
    (
        owner.db().objectRegistry::lookupObject<volScalarField>
        (
            this->coeffDict().lookup("rhol")
        )
    ),
    Uc_
    (
        owner.db().objectRegistry::lookupObject<volVectorField>
        (
            this->coeffDict().lookup("Uc")
        )
    ),
    positionsFile_(this->coeffDict().lookup("positionsFile")),
    positions_
    (
        IOobject
        (
            positionsFile_,
            owner.db().time().constant(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    injectorCells_(positions_.size()),
    nParcelsPerInjector_
    (
        readLabel(this->coeffDict().lookup("parcelsPerInjector"))
    ),
    nParcelsInjected_(positions_.size(), 0),
    U0_(this->coeffDict().lookup("U0")),
    diameters_(positions_.size()),
    parcelPDF_
    (
        pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    )
{
    // Construct parcel diameters - one per injector cell
    forAll(diameters_, i)
    {
        diameters_[i] = parcelPDF_->sample();
    }

    // Determine total volume of particles to inject
    this->volumeTotal_ =
         nParcelsPerInjector_
        *sum(pow3(diameters_))
        *mathematicalConstant::pi/6.0;

    // Set/cahce the injector cells
    forAll(positions_, i)
    {
        this->findCellAtPosition
        (
            injectorCells_[i],
            positions_[i]
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EntrainmentInjection<CloudType>::~EntrainmentInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::EntrainmentInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::EntrainmentInjection<CloudType>::timeEnd() const
{
    return GREAT;
}


template<class CloudType>
void Foam::EntrainmentInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const scalar,
    vector& position,
    label& cellOwner
)
{
    position = positions_[parcelI];
    cellOwner = injectorCells_[parcelI];
}


template<class CloudType>
Foam::vector Foam::EntrainmentInjection<CloudType>::velocity
(
    const label,
    const scalar
)
{
    return U0_;
}


template<class CloudType>
Foam::scalar Foam::EntrainmentInjection<CloudType>::d0
(
    const label parcelI,
    const scalar
) const
{
    return diameters_[parcelI];
}


template<class CloudType>
bool Foam::EntrainmentInjection<CloudType>::validInjection(const label parcelI)
{

    const label cellI = injectorCells_[parcelI];
    if
    (
         nParcelsInjected_[parcelI] < nParcelsPerInjector_
      && rhol_[cellI] < c_*0.5*rhoc_[cellI]*magSqr(Uc_[cellI])
    )
    {
        nParcelsInjected_[parcelI]++;
        return true;
    }

    return false;
}


// ************************************************************************* //
