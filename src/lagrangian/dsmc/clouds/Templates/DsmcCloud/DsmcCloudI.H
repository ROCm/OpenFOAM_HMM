/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
inline const Foam::word& Foam::DsmcCloud<ParcelType>::cloudType() const
{
    return cloudType_;
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::DsmcCloud<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::DsmcCloud<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const Foam::CollisionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::collision() const
{
    return collisionModel_;
}


template<class ParcelType>
inline const Foam::InjectionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::injection() const
{
    return injectionModel_;
}


template<class ParcelType>
inline Foam::InjectionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::injection()
{
    return injectionModel_();
}


template<class ParcelType>
inline const Foam::WallInteractionModel<Foam::DsmcCloud<ParcelType> >&
Foam::DsmcCloud<ParcelType>::wallInteraction() const
{
    return wallInteractionModel_;
}


template<class ParcelType>
inline const Foam::vectorIntegrationScheme&
Foam::DsmcCloud<ParcelType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class ParcelType>
inline Foam::scalar Foam::DsmcCloud<ParcelType>::massInSystem() const
{
    scalar sysMass = 0.0;

    return sysMass;
}


template<class ParcelType>
inline Foam::Random& Foam::DsmcCloud<ParcelType>::rndGen()
{
    return rndGen_;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::rhoN() const
{
    tmp<volScalarField> trhoN
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "rhoN",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, -3, 0, 0, 0), 0.0)
        )
    );

    return trhoN;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::rhoM() const
{
    tmp<volScalarField> trhoM
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "rhoM",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0), 0.0)
        )
    );

    return trhoM;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::rhoNP() const
{
    tmp<volScalarField> trhoNP
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "rhoNP",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, -3, 0, 0, 0), 0.0)
        )
    );

    return trhoNP;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volVectorField>
Foam::DsmcCloud<ParcelType>::U() const
{
    tmp<volScalarField> tU
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "U",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector
            (
                "zero",
                dimensionSet(0, 1, -1, 0, 0),
                vector::zero
            )
        )
    );

    return tU;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::DsmcCloud<ParcelType>::T() const
{
    tmp<volScalarField> tT
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "T",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0), 0.0)
        )
    );

    return tT;
}

// template<class ParcelType>
// inline const Foam::tmp<Foam::volScalarField>
// Foam::DsmcCloud<ParcelType>::alpha() const
// {
//     tmp<volScalarField> talpha
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 this->name() + "Alpha",
//                 this->db().time().timeName(),
//                 this->db(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE,
//                 false
//             ),
//             mesh_,
//             dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
//         )
//     );

//     scalarField& alpha = talpha().internalField();
//     forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
//     {
//         const ParcelType& p = iter();
//         const label cellI = p.cell();

//         alpha[cellI] += p.nParticle()*p.mass();
//     }

//     alpha /= (mesh().cellVolumes()*rho_);

//     return talpha;
// }


// ************************************************************************* //
