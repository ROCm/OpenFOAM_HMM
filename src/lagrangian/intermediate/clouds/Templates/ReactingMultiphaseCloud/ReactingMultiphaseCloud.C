/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "ReactingMultiphaseCloud.H"

#include "DevolatilisationModel.H"
#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::preEvolve()
{
    ReactingCloud<ParcelType>::preEvolve();
}


template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::evolveCloud()
{
    const volScalarField& T = this->thermo().thermo().T();
    const volScalarField cp = this->thermo().thermo().Cp();
    const volScalarField& p = this->thermo().thermo().p();

    autoPtr<interpolation<scalar> > rhoInterp = interpolation<scalar>::New
    (
        this->solution().interpolationSchemes(),
        this->rho()
    );

    autoPtr<interpolation<vector> > UInterp = interpolation<vector>::New
    (
        this->solution().interpolationSchemes(),
        this->U()
    );

    autoPtr<interpolation<scalar> > muInterp = interpolation<scalar>::New
    (
        this->solution().interpolationSchemes(),
        this->mu()
    );

    autoPtr<interpolation<scalar> > TInterp = interpolation<scalar>::New
    (
        this->solution().interpolationSchemes(),
        T
    );

    autoPtr<interpolation<scalar> > cpInterp = interpolation<scalar>::New
    (
        this->solution().interpolationSchemes(),
        cp
    );

    autoPtr<interpolation<scalar> > pInterp = interpolation<scalar>::New
    (
        this->solution().interpolationSchemes(),
        p
    );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterp(),
        UInterp(),
        muInterp(),
        TInterp(),
        cpInterp(),
        pInterp(),
        this->g().value()
    );

    label preInjectionSize = this->size();

    this->surfaceFilm().inject(td);

    // Update the cellOccupancy if the size of the cloud has changed
    // during the injection.
    if (preInjectionSize != this->size())
    {
        this->updateCellOccupancy();

        preInjectionSize = this->size();
    }

    this->injection().inject(td);

    if (this->solution().coupled())
    {
        resetSourceTerms();
    }

    // Assume that motion will update the cellOccupancy as necessary
    // before it is required.
    motion(td);
}


template<class ParcelType>
void  Foam::ReactingMultiphaseCloud<ParcelType>::motion
(
    typename ParcelType::trackData& td
)
{
    ReactingCloud<ParcelType>::motion(td);
}


template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::postEvolve()
{
    ReactingCloud<ParcelType>::postEvolve();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseCloud<ParcelType>::ReactingMultiphaseCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    ReactingCloud<ParcelType>(cloudName, rho, U, g, thermo, false),
    reactingMultiphaseCloud(),
    constProps_(this->particleProperties()),
    devolatilisationModel_
    (
        DevolatilisationModel<ReactingMultiphaseCloud<ParcelType> >::New
        (
            this->subModelProperties(),
            *this
        )
    ),
    surfaceReactionModel_
    (
        SurfaceReactionModel<ReactingMultiphaseCloud<ParcelType> >::New
        (
            this->subModelProperties(),
            *this
        )
    ),
    dMassDevolatilisation_(0.0)
{
    if (readFields)
    {
        ParcelType::readFields(*this);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseCloud<ParcelType>::~ReactingMultiphaseCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    ReactingCloud<ParcelType>::checkParcelProperties
    (
        parcel,
        lagrangianDt,
        fullyDescribed
    );

    label idGas = this->composition().idGas();
    label idLiquid = this->composition().idLiquid();
    label idSolid = this->composition().idSolid();

    if (!fullyDescribed)
    {
        parcel.YGas() = this->composition().Y0(idGas);
        parcel.YLiquid() = this->composition().Y0(idLiquid);
        parcel.YSolid() = this->composition().Y0(idSolid);
    }
    else
    {
        this->checkSuppliedComposition
        (
            parcel.YGas(),
            this->composition().Y0(idGas),
            "YGas"
        );
        this->checkSuppliedComposition
        (
            parcel.YLiquid(),
            this->composition().Y0(idLiquid),
            "YLiquid"
        );
        this->checkSuppliedComposition
        (
            parcel.YSolid(),
            this->composition().Y0(idSolid),
            "YSolid"
        );
    }
}


template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::resetSourceTerms()
{
    ReactingCloud<ParcelType>::resetSourceTerms();
}


template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::evolve()
{
    if (this->solution().active())
    {
        preEvolve();

        evolveCloud();

        postEvolve();

        info();
        Info<< endl;
    }
}


template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::info() const
{
    ReactingCloud<ParcelType>::info();
    Info<< "    Mass transfer devolatilisation  = "
        << returnReduce(dMassDevolatilisation_, sumOp<scalar>()) << nl;
    Info<< "    Mass transfer surface reaction  = "
        << returnReduce(dMassSurfaceReaction_, sumOp<scalar>()) << nl;
}


template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::addToMassDevolatilisation
(
    const scalar dMass
)
{
    dMassDevolatilisation_ += dMass;
}


template<class ParcelType>
void Foam::ReactingMultiphaseCloud<ParcelType>::addToMassSurfaceReaction
(
    const scalar dMass
)
{
    dMassSurfaceReaction_ += dMass;
}


// ************************************************************************* //
