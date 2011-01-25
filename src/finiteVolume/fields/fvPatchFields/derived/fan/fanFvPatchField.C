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

#include "fanFvPatchField.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicFvPatchField<Type>(p, iF),
    f_(0),
    jump_(this->size(), 0.0)
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpCyclicFvPatchField<Type>(ptf, p, iF, mapper),
    f_(ptf.f_),
    jump_(ptf.jump_, mapper)
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicFvPatchField<Type>(p, iF),
    f_(),
    jump_(this->size(), 0.0)
{
    {
        Istream& is = dict.lookup("f");
        is.format(IOstream::ASCII);
        is >> f_;

        // Check that f_ table is same on both sides.?
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::blocking);
    }
}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    jumpCyclicFvPatchField<Type>(ptf),
    f_(ptf.f_),
    jump_(ptf.jump_)
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicFvPatchField<Type>(ptf, iF),
    f_(ptf.f_),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fanFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpCyclicFvPatchField<Type>::autoMap(m);
    jump_.autoMap(m);
}


template<class Type>
void Foam::fanFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    jumpCyclicFvPatchField<Type>::rmap(ptf, addr);

    const fanFvPatchField<Type>& tiptf =
        refCast<const fanFvPatchField<Type> >(ptf);
    jump_.rmap(tiptf.jump_, addr);
}


template<class Type>
void Foam::fanFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("patchType") << "cyclic" << token::END_STATEMENT << nl;

    IOstream::streamFormat fmt0 = os.format(IOstream::ASCII);
    os.writeKeyword("f") << f_ << token::END_STATEMENT << nl;
    os.format(fmt0);

    this->writeEntry("value", os);
}


// ************************************************************************* //
