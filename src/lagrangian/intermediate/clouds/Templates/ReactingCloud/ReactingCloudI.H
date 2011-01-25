/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::ReactingCloud<ParcelType>&
Foam::ReactingCloud<ParcelType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::ReactingCloud<ParcelType>::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::CompositionModel<Foam::ReactingCloud<ParcelType> >&
Foam::ReactingCloud<ParcelType>::composition() const
{
    return compositionModel_;
}


template<class ParcelType>
inline const Foam::PhaseChangeModel<Foam::ReactingCloud<ParcelType> >&
Foam::ReactingCloud<ParcelType>::phaseChange() const
{
    return phaseChangeModel_;
}


template<class ParcelType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::ReactingCloud<ParcelType>::rhoTrans(const label i)
{
    return rhoTrans_[i];
}


template<class ParcelType>
inline
const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<ParcelType>::rhoTrans() const
{
    return rhoTrans_;
}


template<class ParcelType>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >&
Foam::ReactingCloud<ParcelType>::rhoTrans()
{
    return rhoTrans_;
}


template<class ParcelType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingCloud<ParcelType>::SYi
(
    const label i,
    volScalarField& Yi
) const
{
    if (this->solution().coupled())
    {
        if (this->solution().semiImplicit("Yi"))
        {
            tmp<volScalarField> trhoTrans
            (
                new volScalarField
                (
                    IOobject
                    (
                        this->name() + "rhoTrans",
                        this->db().time().timeName(),
                        this->db(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    this->mesh(),
                    dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
                )
            );

            volScalarField& sourceField = trhoTrans();

            sourceField.internalField() =
                rhoTrans_[i]/(this->db().time().deltaT()*this->mesh().V());

            const dimensionedScalar YiSMALL("YiSMALL", dimless, SMALL);

            return
                fvm::Sp(neg(sourceField)*sourceField/(Yi + YiSMALL), Yi)
              + pos(sourceField)*sourceField;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Yi, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm();

            fvm.source() = -rhoTrans_[i]/this->db().time().deltaT();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(Yi, dimMass/dimTime));
}


template<class ParcelType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<ParcelType>::Srho(const label i) const
{
    tmp<DimensionedField<scalar, volMesh> > tRhoi
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                "zero",
                rhoTrans_[0].dimensions()/dimTime/dimVolume,
                0.0
            )
        )
    );

    if (this->solution().coupled())
    {
        scalarField& rhoi = tRhoi();
        rhoi = rhoTrans_[i]/(this->db().time().deltaT()*this->mesh().V());
    }

    return tRhoi;
}


template<class ParcelType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::ReactingCloud<ParcelType>::Srho() const
{
    tmp<DimensionedField<scalar, volMesh> > trhoTrans
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                "zero",
                rhoTrans_[0].dimensions()/dimTime/dimVolume,
                0.0
            )
        )
    );

    if (this->solution().coupled())
    {
        scalarField& sourceField = trhoTrans();
        forAll(rhoTrans_, i)
        {
            sourceField += rhoTrans_[i];
        }

        sourceField /= this->db().time().deltaT()*this->mesh().V();
    }

    return trhoTrans;
}


template<class ParcelType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::ReactingCloud<ParcelType>::Srho(volScalarField& rho) const
{
    if (this->solution().coupled())
    {
        tmp<volScalarField> trhoTrans
        (
            new volScalarField
            (
                IOobject
                (
                    this->name() + "rhoTrans",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
            )
        );

        scalarField& sourceField = trhoTrans();

        if (this->solution().semiImplicit("rho"))
        {

            forAll(rhoTrans_, i)
            {
                sourceField += rhoTrans_[i];
            }
            sourceField /= this->db().time().deltaT()*this->mesh().V();

            return fvm::SuSp(trhoTrans()/rho, rho);
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm();

            forAll(rhoTrans_, i)
            {
                sourceField += rhoTrans_[i];
            }

            fvm.source() = -trhoTrans()/this->db().time().deltaT();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>(new fvScalarMatrix(rho, dimMass/dimTime));
}


// ************************************************************************* //
