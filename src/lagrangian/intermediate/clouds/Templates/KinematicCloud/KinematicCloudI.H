/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::word& Foam::KinematicCloud<ParcelType>::cloudType() const
{
    return cloudType_;
}


template<class ParcelType>
inline const Foam::label Foam::KinematicCloud<ParcelType>::parcelTypeId() const
{
    return parcelTypeId_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::KinematicCloud<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::Time& Foam::KinematicCloud<ParcelType>::runTime() const
{
    return runTime_;
}


template<class ParcelType>
inline const Foam::volPointInterpolation&
Foam::KinematicCloud<ParcelType>::vpi() const
{
    return vpi_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const Foam::Switch Foam::KinematicCloud<ParcelType>::coupled() const
{
    return coupled_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::KinematicCloud<ParcelType>::rho() const
{
    return rho_;
}


template<class ParcelType>
inline const Foam::volVectorField& Foam::KinematicCloud<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::KinematicCloud<ParcelType>::mu() const
{
    return mu_;
}


template<class ParcelType>
inline const Foam::dimensionedVector&
Foam::KinematicCloud<ParcelType>::g() const
{
    return g_;
}


template<class ParcelType>
inline const Foam::DispersionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::dispersion() const
{
    return dispersionModel_;
}


template<class ParcelType>
inline Foam::DispersionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::dispersion()
{
    return dispersionModel_();
}


template<class ParcelType>
inline const Foam::DragModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::drag() const
{
    return dragModel_;
}


template<class ParcelType>
inline const Foam::InjectionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::injection() const
{
    return injectionModel_;
}


template<class ParcelType>
inline const Foam::WallInteractionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::wallInteraction() const
{
    return wallInteractionModel_;
}


template<class ParcelType>
inline const Foam::label Foam::KinematicCloud<ParcelType>::nInjections() const
{
    return nInjections_;
}


template<class ParcelType>
inline const Foam::label Foam::KinematicCloud<ParcelType>::nParcelsAddedTotal() const
{
    return nParcelsAddedTotal_;
}


template<class ParcelType>
inline const Foam::scalar Foam::KinematicCloud<ParcelType>::time0() const
{
    return time0_;
}


template<class ParcelType>
inline const Foam::scalar Foam::KinematicCloud<ParcelType>::massTotal() const
{
    return massTotal_;
}


template<class ParcelType>
inline const Foam::scalar Foam::KinematicCloud<ParcelType>::massInjected() const
{
    return massInjected_;
}


template<class ParcelType>
inline const Foam::scalar Foam::KinematicCloud<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
         const ParcelType& p = iter();
         sysMass += p.mass()*p.nParticle();
    }

    return sysMass;
}


template<class ParcelType>
inline Foam::Random& Foam::KinematicCloud<ParcelType>::rndGen()
{
    return rndGen_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<ParcelType>::UTrans()
{
    return UTrans_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<ParcelType>::UCoeff()
{
    return UCoeff_;
}


template<class ParcelType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::KinematicCloud<ParcelType>::SU(volVectorField& U) const
{
    Info<< "UTrans min/max = "
        << min(UTrans_) << ", " << max(UTrans_) << endl;
    Info<< "UCoeff min/max = "
        << min(UCoeff_) << ", " << max(UCoeff_) << endl;

    return UTrans_/(mesh_.V()*runTime_.deltaT())
         - fvm::Sp(UCoeff_/mesh_.V(), U)
         + UCoeff_/mesh_.V()*U;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<ParcelType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->cloudName() + "Theta",
                runTime_.timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    );

    scalarField& theta = ttheta().internalField();
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        theta[cellI] += p.nParticle()*p.volume();
    }

    theta /= mesh().cellVolumes();

    return ttheta;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<ParcelType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->cloudName() + "Alpha",
                runTime_.timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
        )
    );

    scalarField& alpha = talpha().internalField();
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        alpha[cellI] += p.nParticle()*p.mass();
    }

    alpha /= (mesh().cellVolumes()*rho_);

    return talpha;
}


// ************************************************************************* //
