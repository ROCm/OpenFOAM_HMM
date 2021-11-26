/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "SemiImplicitSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "Constant.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::Enum
<
    typename Foam::fv::SemiImplicitSource<Type>::volumeModeType
>
Foam::fv::SemiImplicitSource<Type>::volumeModeTypeNames_
({
    { volumeModeType::vmAbsolute, "absolute" },
    { volumeModeType::vmSpecific, "specific" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::SemiImplicitSource<Type>::setFieldData(const dictionary& dict)
{
    label count = dict.size();

    fieldNames_.resize(count);
    Su_.resize(fieldNames_.size());
    Sp_.resize(fieldNames_.size());

    fv::option::resetApplied();

    count = 0;
    for (const entry& dEntry : dict)
    {
        fieldNames_[count] = dEntry.keyword();

        if (!dEntry.isDict())
        {
            Tuple2<Type, scalar> injectionRate;
            dEntry.readEntry(injectionRate);

            Su_.set
            (
                count,
                new Function1Types::Constant<Type>
                (
                    "Su",
                    injectionRate.first()
                )
            );
            Sp_.set
            (
                count,
                new Function1Types::Constant<scalar>
                (
                    "Sp",
                    injectionRate.second()
                )
            );
        }
        else
        {
            const dictionary& Sdict = dEntry.dict();
            Su_.set(count, Function1<Type>::New("Su", Sdict, &mesh_));
            Sp_.set(count, Function1<scalar>::New("Sp", Sdict, &mesh_));
        }

        ++count;
    }

    // Set volume normalisation
    if (volumeMode_ == vmAbsolute)
    {
        VDash_ = V_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::SemiImplicitSource<Type>::SemiImplicitSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    volumeMode_(vmAbsolute),
    VDash_(1.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "SemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    typename GeometricField<Type, fvPatchField, volMesh>::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Su",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<Type>(eqn.dimensions()/dimVolume, Zero),
        false
    );

    const scalar tmVal = mesh_.time().timeOutputValue();

    UIndirectList<Type>(Su, cells_) = Su_[fieldi].value(tmVal)/VDash_;

    volScalarField::Internal Sp
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensioned<scalar>(Su.dimensions()/psi.dimensions(), Zero),
        false
    );

    UIndirectList<scalar>(Sp, cells_) = Sp_[fieldi].value(tmVal)/VDash_;

    eqn += Su + fvm::SuSp(Sp, psi);
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "SemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    return this->addSup(eqn, fieldi);
}


template<class Type>
bool Foam::fv::SemiImplicitSource<Type>::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        volumeMode_ = volumeModeTypeNames_.get("volumeMode", coeffs_);
        setFieldData(coeffs_.subDict("injectionRateSuSp"));

        return true;
    }

    return false;
}


// ************************************************************************* //
