/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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

#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::ThermoCloud<ParcelType>&
Foam::ThermoCloud<ParcelType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::ThermoCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::SLGThermo& Foam::ThermoCloud<ParcelType>::thermo() const
{
    return thermo_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::ThermoCloud<ParcelType>::T() const
{
    return T_;
}


template<class ParcelType>
inline const Foam::volScalarField& Foam::ThermoCloud<ParcelType>::p() const
{
    return p_;
}


template<class ParcelType>
inline const Foam::HeatTransferModel<Foam::ThermoCloud<ParcelType> >&
Foam::ThermoCloud<ParcelType>::heatTransfer() const
{
    return heatTransferModel_;
}


template<class ParcelType>
inline const Foam::scalarIntegrationScheme&
Foam::ThermoCloud<ParcelType>::TIntegrator() const
{
    return TIntegrator_;
}


template<class ParcelType>
inline bool Foam::ThermoCloud<ParcelType>::radiation() const
{
    return radiation_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<ParcelType>::hsTrans()
{
    return hsTrans_();
}


template<class ParcelType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<ParcelType>::hsTrans() const
{
    return hsTrans_();
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<ParcelType>::hsCoeff()
{
    return hsCoeff_();
}


template<class ParcelType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ThermoCloud<ParcelType>::hsCoeff() const
{
    return hsCoeff_();
}


template<class ParcelType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ThermoCloud<ParcelType>::Sh(volScalarField& hs) const
{
    if (debug)
    {
        Info<< "hsTrans min/max = " << min(hsTrans()).value() << ", "
            << max(hsTrans()).value() << nl
            << "hsCoeff min/max = " << min(hsCoeff()).value() << ", "
            << max(hsCoeff()).value() << endl;
    }

    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("hs"))
        {
            const volScalarField Cp(thermo_.thermo().Cp());

            return
                hsTrans()/(this->mesh().V()*this->db().time().deltaT())
              + fvm::SuSp(-hsCoeff()/(Cp*this->mesh().V()), hs)
              + hsCoeff()/(Cp*this->mesh().V())*hs;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(hs, dimEnergy/dimTime));
            fvScalarMatrix& fvm = tfvm();

            fvm.source() = -hsTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(hs, dimEnergy/dimTime));
}


template<class ParcelType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<ParcelType>::Ep() const
{
    tmp<volScalarField> tEp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "radiation::Ep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    // Need to check if coupled as field is created on-the-fly
    if (radiation_ && this->solution().coupled())
    {
        scalarField& Ep = tEp().internalField();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();

        forAllConstIter(typename ThermoCloud<ParcelType>, *this, iter)
        {
            const ParcelType& p = iter();
            const label cellI = p.cell();
            Ep[cellI] += p.nParticle()*p.areaP()*pow4(p.T());
        }

        Ep *= epsilon*physicoChemical::sigma.value()/V;
    }

    return tEp;
}


template<class ParcelType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<ParcelType>::ap() const
{
    tmp<volScalarField> tap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "radiation::ap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    // Need to check if coupled as field is created on-the-fly
    if (radiation_ && this->solution().coupled())
    {
        scalarField& ap = tap().internalField();
        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();

        forAllConstIter(typename ThermoCloud<ParcelType>, *this, iter)
        {
            const ParcelType& p = iter();
            const label cellI = p.cell();
            ap[cellI] += p.nParticle()*p.areaP();
        }

        ap *= epsilon/V;
    }

    return tap;
}


template<class ParcelType>
inline Foam::tmp<Foam::volScalarField>
Foam::ThermoCloud<ParcelType>::sigmap() const
{
    tmp<volScalarField> tsigmap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "radiation::sigmap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

    // Need to check if coupled as field is created on-the-fly
    if (radiation_ && this->solution().coupled())
    {
        scalarField& sigmap = tsigmap().internalField();

        const scalarField& V = this->mesh().V();
        const scalar epsilon = constProps_.epsilon0();
        const scalar f = constProps_.f0();

        forAllConstIter(typename ThermoCloud<ParcelType>, *this, iter)
        {
            const ParcelType& p = iter();
            const label cellI = p.cell();
            sigmap[cellI] += p.nParticle()*p.areaP();
        }

        sigmap *= (1.0 - f)*(1.0 - epsilon)/V;
    }

    return tsigmap;
}


// ************************************************************************* //
