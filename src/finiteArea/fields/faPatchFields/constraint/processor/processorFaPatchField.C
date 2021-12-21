/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "processorFaPatchField.H"
#include "processorFaPatch.H"
#include "IPstream.H"
#include "OPstream.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFaPatch>(p))
{}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    coupledFaPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFaPatch>(p))
{}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFaPatch>(p))
{
    if (!isType<processorFaPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    coupledFaPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFaPatch>(p, dict))
{
    if (!isType<processorFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFaPatchField<Type>(ptf),
    procPatch_(refCast<const processorFaPatch>(ptf.patch()))
{}


template<class Type>
Foam::processorFaPatchField<Type>::processorFaPatchField
(
    const processorFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFaPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFaPatchField<Type>::~processorFaPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::processorFaPatchField<Type>::patchNeighbourField() const
{
    return *this;
}


template<class Type>
void Foam::processorFaPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        procPatch_.send(commsType, this->patchInternalField()());
    }
}


template<class Type>
void Foam::processorFaPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        procPatch_.receive<Type>(commsType, *this);

        if (doTransform())
        {
            transform(*this, procPatch_.forwardT(), *this);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorFaPatchField<Type>::snGrad() const
{
    return this->patch().deltaCoeffs()*(*this - this->patchInternalField());
}


template<class Type>
void Foam::processorFaPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    procPatch_.send
    (
        commsType,
        this->patch().patchInternalField(psiInternal)()
    );
}


template<class Type>
void Foam::processorFaPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    solveScalarField pnf
    (
        procPatch_.receive<solveScalar>(commsType, this->size())()
    );

    // Transform according to the transformation tensor
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result

    const labelUList& edgeFaces = this->patch().edgeFaces();

    if (add)
    {
        forAll(edgeFaces, elemI)
        {
            result[edgeFaces[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(edgeFaces, elemI)
        {
            result[edgeFaces[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


template<class Type>
void Foam::processorFaPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    procPatch_.send
    (
        commsType,
        this->patch().patchInternalField(psiInternal)()
    );
}


template<class Type>
void Foam::processorFaPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>&,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    Field<Type> pnf
    (
        procPatch_.receive<Type>(commsType, this->size())()
    );

    // Multiply the field by coefficients and add into the result

    const labelUList& edgeFaces = this->patch().edgeFaces();

    if (add)
    {
        forAll(edgeFaces, elemI)
        {
            result[edgeFaces[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(edgeFaces, elemI)
        {
            result[edgeFaces[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


// ************************************************************************* //
