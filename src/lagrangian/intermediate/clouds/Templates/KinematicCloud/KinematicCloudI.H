/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

#include "fvmSup.H"

// * * * * * * * * * * * cloudSolution Member Functions  * * * * * * * * * * //

template<class CloudType>
inline const Foam::fvMesh&
Foam::KinematicCloud<CloudType>::cloudSolution::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicCloud<CloudType>::cloudSolution::dict() const
{
    return dict_;
}


template<class CloudType>
inline const Foam::Switch
Foam::KinematicCloud<CloudType>::cloudSolution::active() const
{
    return active_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicCloud<CloudType>::cloudSolution::sourceTermDict() const
{
    return dict_.subDict("sourceTerms");
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicCloud<CloudType>::cloudSolution::interpolationSchemes() const
{
    return dict_.subDict("interpolationSchemes");
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicCloud<CloudType>::cloudSolution::integrationSchemes() const
{
    return dict_.subDict("integrationSchemes");
}


template<class CloudType>
inline const Foam::Switch
Foam::KinematicCloud<CloudType>::cloudSolution::transient() const
{
    return transient_;
}


template<class CloudType>
inline const Foam::Switch
Foam::KinematicCloud<CloudType>::cloudSolution::steadyState() const
{
    return !transient_;
}


template<class CloudType>
inline Foam::label
Foam::KinematicCloud<CloudType>::cloudSolution::calcFrequency() const
{
    return calcFrequency_;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicCloud<CloudType>::cloudSolution::maxCo() const
{
    return maxCo_;
}


template<class CloudType>
inline Foam::label Foam::KinematicCloud<CloudType>::cloudSolution::iter() const
{
    return iter_;
}


template<class CloudType>
inline Foam::label Foam::KinematicCloud<CloudType>::cloudSolution::nextIter()
{
    return ++iter_;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicCloud<CloudType>::cloudSolution::deltaT() const
{
    return deltaT_;
}


template<class CloudType>
inline const Foam::Switch
Foam::KinematicCloud<CloudType>::cloudSolution::coupled() const
{
    return coupled_;
}


template<class CloudType>
inline const Foam::Switch
Foam::KinematicCloud<CloudType>::cloudSolution::cellValueSourceCorrection()
const
{
    return cellValueSourceCorrection_;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicCloud<CloudType>::cloudSolution::maxTrackTime() const
{
    return maxTrackTime_;
}


template<class CloudType>
inline const Foam::Switch
Foam::KinematicCloud<CloudType>::cloudSolution::resetSourcesOnStartup() const
{
    return resetSourcesOnStartup_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::KinematicCloud<CloudType>&
Foam::KinematicCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline bool Foam::KinematicCloud<CloudType>::hasWallImpactDistance() const
{
    return !collision().controlsWallInteraction();
}


template<class CloudType>
inline const Foam::fvMesh& Foam::KinematicCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const typename Foam::KinematicCloud<CloudType>::cloudSolution&
Foam::KinematicCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline typename Foam::KinematicCloud<CloudType>::cloudSolution&
Foam::KinematicCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::KinematicCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::KinematicCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicCloud<CloudType>::mu() const
{
    return mu_;
}


template<class CloudType>
inline const Foam::dimensionedVector& Foam::KinematicCloud<CloudType>::g() const
{
    return g_;
}


template<class CloudType>
//inline const typename CloudType::parcelType::forceType&
inline const typename Foam::KinematicCloud<CloudType>::forceType&
Foam::KinematicCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline const Foam::CollisionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::collision() const
{
    return collisionModel_();
}


template<class CloudType>
inline Foam::CollisionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::collision()
{
    return collisionModel_();
}


template<class CloudType>
inline const Foam::DispersionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::dispersion() const
{
    return dispersionModel_;
}


template<class CloudType>
inline Foam::DispersionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::dispersion()
{
    return dispersionModel_();
}


template<class CloudType>
inline const Foam::InjectionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::injection() const
{
    return injectionModel_;
}


template<class CloudType>
inline const Foam::PatchInteractionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::patchInteraction() const
{
    return patchInteractionModel_;
}


template<class CloudType>
inline Foam::PatchInteractionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::patchInteraction()
{
    return patchInteractionModel_();
}


template<class CloudType>
inline Foam::InjectionModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::injection()
{
    return injectionModel_();
}


template<class CloudType>
inline Foam::PostProcessingModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::postProcessing()
{
    return postProcessingModel_();
}


template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::surfaceFilm() const
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::KinematicCloud<CloudType> >&
Foam::KinematicCloud<CloudType>::surfaceFilm()
{
    return surfaceFilmModel_();
}


template<class CloudType>
inline const Foam::vectorIntegrationScheme&
Foam::KinematicCloud<CloudType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
         const parcelType& p = iter();
         sysMass += p.mass()*p.nParticle();
    }

    return sysMass;
}


template<class CloudType>
inline Foam::vector
Foam::KinematicCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearMomentum += p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        linearKineticEnergy += 0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicCloud<CloudType>::rotationalKineticEnergyOfSystem() const
{
    scalar rotationalKineticEnergy = 0.0;

    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();

        rotationalKineticEnergy +=
            0.5*p.momentOfInertia()*(p.omega() & p.omega());
    }

    return rotationalKineticEnergy;
}


template<class CloudType>
inline Foam::cachedRandom& Foam::KinematicCloud<CloudType>::rndGen()
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*> >&
Foam::KinematicCloud<CloudType>::cellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UTrans()
{
    return UTrans_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UTrans() const
{
    return UTrans_();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UCoeff()
{
    return UCoeff_();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<CloudType>::UCoeff() const
{
    return UCoeff_();
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::KinematicCloud<CloudType>::SU(volVectorField& U) const
{
    if (debug)
    {
        Info<< "UTrans min/max = " << min(UTrans()).value() << ", "
            << max(UTrans()).value() << nl
            << "UCoeff min/max = " << min(UCoeff()).value() << ", "
            << max(UCoeff()).value() << endl;
    }

    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            return
                UTrans()/(mesh_.V()*this->db().time().deltaT())
              + fvm::SuSp(-UCoeff()/mesh_.V(), U)
              + UCoeff()/mesh_.V()*U;
        }
        else
        {
            tmp<fvVectorMatrix> tfvm(new fvVectorMatrix(U, dimForce));
            fvVectorMatrix& fvm = tfvm();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>(new fvVectorMatrix(U, dimForce));
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::theta() const
{
    tmp<volScalarField> ttheta
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Theta",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& theta = ttheta().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        theta[cellI] += p.nParticle()*p.volume();
    }

    theta /= mesh().V();

    return ttheta;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::alpha() const
{
    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "Alpha",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    scalarField& alpha = talpha().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        alpha[cellI] += p.nParticle()*p.mass();
    }

    alpha /= (mesh().V()*rho_);

    return talpha;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<CloudType>::rhoEff() const
{
    tmp<volScalarField> trhoEff
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "RhoEff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    scalarField& rhoEff = trhoEff().internalField();
    forAllConstIter(typename KinematicCloud<CloudType>, *this, iter)
    {
        const parcelType& p = iter();
        const label cellI = p.cell();

        rhoEff[cellI] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh().V();

    return trhoEff;
}


// ************************************************************************* //
