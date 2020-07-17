/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "freestreamFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::freestreamFvPatchField<Type>::freestreamFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(p, iF),
    freestreamBCPtr_()
{}


template<class Type>
Foam::freestreamFvPatchField<Type>::freestreamFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchField<Type>(p, iF),
    freestreamBCPtr_()
{
    this->patchType() = dict.getOrDefault<word>("patchType", word::null);

    this->phiName_ = dict.getOrDefault<word>("phi", "phi");

    if (dict.found("freestreamValue"))
    {
        freestreamValue() =
            Field<Type>("freestreamValue", dict, p.size());

        if (dict.found("value"))
        {
            fvPatchField<Type>::operator=
            (
                Field<Type>("value", dict, p.size())
            );
        }
        else
        {
            fvPatchField<Type>::operator=(freestreamValue());
        }
    }
    else
    {
        // Freestream value provided by another patch
        freestreamBCPtr_ =
            fvPatchField<Type>::New(p, iF, dict.subDict("freestreamBC"));

        // Force user to supply an initial value
        // - we do not know if the supplied BC has all dependencies available
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
}


template<class Type>
Foam::freestreamFvPatchField<Type>::freestreamFvPatchField
(
    const freestreamFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchField<Type>(ptf, p, iF, mapper),
    freestreamBCPtr_()
{
    if (ptf.freestreamBCPtr_)
    {
        freestreamBCPtr_ =
            fvPatchField<Type>::New(ptf.freestreamBCPtr_(), p, iF, mapper);
    }
}


template<class Type>
Foam::freestreamFvPatchField<Type>::freestreamFvPatchField
(
    const freestreamFvPatchField<Type>& ptf
)
:
    inletOutletFvPatchField<Type>(ptf),
    freestreamBCPtr_()
{
    if (ptf.freestreamBCPtr_)
    {
        freestreamBCPtr_ = ptf.freestreamBCPtr_->clone();
    }
}


template<class Type>
Foam::freestreamFvPatchField<Type>::freestreamFvPatchField
(
    const freestreamFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(ptf, iF),
    freestreamBCPtr_()
{
    if (ptf.freestreamBCPtr_)
    {
        freestreamBCPtr_ = ptf.freestreamBCPtr_->clone();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::freestreamFvPatchField<Type>::autoMap(const fvPatchFieldMapper& m)
{
    inletOutletFvPatchField<Type>::autoMap(m);
    if (freestreamBCPtr_)
    {
        freestreamBCPtr_->autoMap(m);
    }
}


template<class Type>
void Foam::freestreamFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    inletOutletFvPatchField<Type>::rmap(ptf, addr);

    const auto& fsptf = refCast<const freestreamFvPatchField<Type>>(ptf);

    if (fsptf.freestreamBCPtr_)
    {
        freestreamBCPtr_->rmap(fsptf.freestreamBCPtr_(), addr);
    }
}


template<class Type>
void Foam::freestreamFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (freestreamBCPtr_)
    {
        freestreamBCPtr_->evaluate();
        freestreamValue() = freestreamBCPtr_();
    }

    inletOutletFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::freestreamFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);

    if (freestreamBCPtr_)
    {
        os.beginBlock("freestreamBC");
        freestreamBCPtr_->write(os);
        os.endBlock();
    }
    else
    {
        freestreamValue().writeEntry("freestreamValue", os);
    }
    this->writeEntry("value", os);
}


// ************************************************************************* //
