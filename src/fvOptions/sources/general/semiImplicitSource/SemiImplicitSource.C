/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
#include "Tuple2.H"

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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::SemiImplicitSource<Type>::setFieldInjectionRates
(
    const dictionary& dict
)
{
    label count = dict.size();

    fieldNames_.resize_nocopy(count);

    Su_.resize(count);
    Sp_.resize(count);

    fv::option::resetApplied();

    count = 0;
    for (const entry& dEntry : dict)
    {
        const word& fieldName = dEntry.keyword();

        if (dEntry.isDict())
        {
            const dictionary& subdict = dEntry.dict();
            Su_.set(count, Function1<Type>::New("Su", subdict, &mesh_));
            Sp_.set(count, Function1<scalar>::New("Sp", subdict, &mesh_));
        }
        else
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

        fieldNames_[count] = fieldName;
        ++count;
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
    VDash_(1)
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
    return this->addSup(volScalarField::null(), eqn, fieldi);
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

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    const word& fieldName = fieldNames_[fieldi];

    const scalar tmVal = mesh_.time().timeOutputValue();

    const dimensionSet SuDims(eqn.dimensions()/dimVolume);
    const dimensionSet SpDims(SuDims/psi.dimensions());


    // Explicit source
    {
        const dimensioned<Type> SuValue
        (
            "Su",
            SuDims,
            Su_[fieldi].value(tmVal)/VDash_
        );

        if (mag(SuValue.value()) <= ROOTVSMALL)
        {
            // No-op
        }
        else if (this->useSubMesh())
        {
            auto tsu = DimensionedField<Type, volMesh>::New
            (
                name_ + fieldName + "Su",
                mesh_,
                dimensioned<Type>(SuDims, Zero)
            );
            UIndirectList<Type>(tsu.ref(), cells_) = SuValue.value();

            eqn += tsu;
        }
        else
        {
            eqn += SuValue;
        }
    }


    // Implicit source
    {
        const dimensioned<scalar> SpValue
        (
            "Sp",
            SpDims,
            Sp_[fieldi].value(tmVal)/VDash_
        );

        if (mag(SpValue.value()) <= ROOTVSMALL)
        {
            // No-op
        }
        else if (this->useSubMesh())
        {
            auto tsp = DimensionedField<scalar, volMesh>::New
            (
                name_ + fieldName + "Sp",
                mesh_,
                dimensioned<scalar>(SpDims, Zero)
            );
            UIndirectList<scalar>(tsp.ref(), cells_) = SpValue.value();

            eqn += fvm::SuSp(tsp, psi);
        }
        else
        {
            eqn += fvm::SuSp(SpValue, psi);
        }
    }
}


template<class Type>
bool Foam::fv::SemiImplicitSource<Type>::read(const dictionary& dict)
{
    VDash_ = 1;

    if (fv::cellSetOption::read(dict))
    {
        volumeMode_ = volumeModeTypeNames_.get("volumeMode", coeffs_);

        // Set volume normalisation
        if (volumeMode_ == vmAbsolute)
        {
            VDash_ = V_;
        }

        {
            setFieldInjectionRates
            (
                coeffs_.subDict("injectionRateSuSp", keyType::LITERAL)
            );
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
