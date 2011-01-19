/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

template<class ParcelType>
inline const Foam::fvMesh&
Foam::KinematicCloud<ParcelType>::cloudSolution::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicCloud<ParcelType>::cloudSolution::dict() const
{
    return dict_;
}


template<class ParcelType>
inline const Foam::Switch
Foam::KinematicCloud<ParcelType>::cloudSolution::active() const
{
    return active_;
}


template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicCloud<ParcelType>::cloudSolution::sourceTermDict() const
{
    return dict_.subDict("sourceTerms");
}


template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicCloud<ParcelType>::cloudSolution::interpolationSchemes() const
{
    return dict_.subDict("interpolationSchemes");
}


template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicCloud<ParcelType>::cloudSolution::integrationSchemes() const
{
    return dict_.subDict("integrationSchemes");
}


template<class ParcelType>
inline const Foam::Switch
Foam::KinematicCloud<ParcelType>::cloudSolution::transient() const
{
    return transient_;
}


template<class ParcelType>
inline const Foam::Switch
Foam::KinematicCloud<ParcelType>::cloudSolution::steadyState() const
{
    return !transient_;
}


template<class ParcelType>
inline Foam::label
Foam::KinematicCloud<ParcelType>::cloudSolution::calcFrequency() const
{
    return calcFrequency_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicCloud<ParcelType>::cloudSolution::maxCo() const
{
    return maxCo_;
}


template<class ParcelType>
inline Foam::label Foam::KinematicCloud<ParcelType>::cloudSolution::iter() const
{
    return iter_;
}


template<class ParcelType>
inline Foam::label Foam::KinematicCloud<ParcelType>::cloudSolution::nextIter()
{
    return ++iter_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicCloud<ParcelType>::cloudSolution::deltaT() const
{
    return deltaT_;
}


template<class ParcelType>
inline const Foam::Switch
Foam::KinematicCloud<ParcelType>::cloudSolution::coupled() const
{
    return coupled_;
}


template<class ParcelType>
inline const Foam::Switch
Foam::KinematicCloud<ParcelType>::cloudSolution::cellValueSourceCorrection()
const
{
    return cellValueSourceCorrection_;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicCloud<ParcelType>::cloudSolution::maxTrackTime() const
{
    return maxTrackTime_;
}


template<class ParcelType>
inline const Foam::Switch
Foam::KinematicCloud<ParcelType>::cloudSolution::resetSourcesOnStartup() const
{
    return resetSourcesOnStartup_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::KinematicCloud<ParcelType>&
Foam::KinematicCloud<ParcelType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class ParcelType>
inline bool Foam::KinematicCloud<ParcelType>::hasWallImpactDistance() const
{
    return !collision().controlsWallInteraction();
}


template<class ParcelType>
inline const Foam::fvMesh& Foam::KinematicCloud<ParcelType>::mesh() const
{
    return mesh_;
}


template<class ParcelType>
inline const Foam::IOdictionary&
Foam::KinematicCloud<ParcelType>::particleProperties() const
{
    return particleProperties_;
}


template<class ParcelType>
inline const typename Foam::KinematicCloud<ParcelType>::cloudSolution&
Foam::KinematicCloud<ParcelType>::solution() const
{
    return solution_;
}


template<class ParcelType>
inline typename Foam::KinematicCloud<ParcelType>::cloudSolution&
Foam::KinematicCloud<ParcelType>::solution()
{
    return solution_;
}


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::KinematicCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicCloud<ParcelType>::subModelProperties() const
{
    return subModelProperties_;
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
inline const Foam::particleForces&
Foam::KinematicCloud<ParcelType>::forces() const
{
    return forces_;
}


template<class ParcelType>
inline const Foam::CollisionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::collision() const
{
    return collisionModel_();
}


template<class ParcelType>
inline Foam::CollisionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::collision()
{
    return collisionModel_();
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
inline const Foam::PatchInteractionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::patchInteraction() const
{
    return patchInteractionModel_;
}


template<class ParcelType>
inline Foam::InjectionModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::injection()
{
    return injectionModel_();
}


template<class ParcelType>
inline Foam::PostProcessingModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::postProcessing()
{
    return postProcessingModel_();
}


template<class ParcelType>
inline const Foam::SurfaceFilmModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::surfaceFilm() const
{
    return surfaceFilmModel_();
}


template<class ParcelType>
inline Foam::SurfaceFilmModel<Foam::KinematicCloud<ParcelType> >&
Foam::KinematicCloud<ParcelType>::surfaceFilm()
{
    return surfaceFilmModel_();
}


template<class ParcelType>
inline const Foam::vectorIntegrationScheme&
Foam::KinematicCloud<ParcelType>::UIntegrator() const
{
    return UIntegrator_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicCloud<ParcelType>::massInSystem() const
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
inline Foam::vector
Foam::KinematicCloud<ParcelType>::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        linearMomentum += p.mass()*p.U();
    }

    return linearMomentum;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicCloud<ParcelType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        linearKineticEnergy += 0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicCloud<ParcelType>::
rotationalKineticEnergyOfSystem() const
{
    scalar rotationalKineticEnergy = 0.0;

    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        rotationalKineticEnergy +=
            0.5*p.momentOfInertia()*(p.omega() & p.omega());
    }

    return rotationalKineticEnergy;
}


template<class ParcelType>
inline Foam::cachedRandom& Foam::KinematicCloud<ParcelType>::rndGen()
{
    return rndGen_;
}


template<class ParcelType>
inline Foam::List<Foam::DynamicList<ParcelType*> >&
Foam::KinematicCloud<ParcelType>::cellOccupancy()
{
    if (cellOccupancyPtr_.empty())
    {
        buildCellOccupancy();
    }

    return cellOccupancyPtr_();
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<ParcelType>::UTrans()
{
    return UTrans_();
}


template<class ParcelType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicCloud<ParcelType>::UTrans() const
{
    return UTrans_();
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<ParcelType>::UCoeff()
{
    return UCoeff_();
}


template<class ParcelType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicCloud<ParcelType>::UCoeff() const
{
    return UCoeff_();
}


template<class ParcelType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::KinematicCloud<ParcelType>::SU(volVectorField& U) const
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
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        theta[cellI] += p.nParticle()*p.volume();
    }

    theta /= mesh().V();

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
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        alpha[cellI] += p.nParticle()*p.mass();
    }

    alpha /= (mesh().V()*rho_);

    return talpha;
}


template<class ParcelType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicCloud<ParcelType>::rhoEff() const
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
    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        const label cellI = p.cell();

        rhoEff[cellI] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh().V();

    return trhoEff;
}


// ************************************************************************* //
