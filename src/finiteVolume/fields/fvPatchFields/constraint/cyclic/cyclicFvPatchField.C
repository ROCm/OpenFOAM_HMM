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

#include "cyclicFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isType<cyclicFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicFvPatchField<Type>::cyclicFvPatchField\n"
            "(\n"
            "    const cyclicFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{
    if (!isType<cyclicFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicFvPatchField<Type>::cyclicFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    Pout<< "Construct from dictionary for " << p.name() << endl
        << "Underlying cyclic:" << cyclicPatch_.name()
        << " with parallel:" << cyclicPatch_.parallel() << endl;

    this->coupledFvPatchField<Type>::evaluate(Pstream::blocking);
}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Referred patch functionality
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<class Type>
void cyclicFvPatchField<Type>::snGrad
(
    Field<Type>& exchangeBuf,
    const Field<Type>& subFld,
    const coupledFvPatch& referringPatch,
    const label size,
    const label start
) const
{
    if (subFld.size() != size)
    {
        FatalErrorIn("cyclicFvPatchField<Type>::snGrad(..)")
            << "Call with correct slice size." << abort(FatalError);
    }

    // Slice delta coeffs
    SubField<scalar> subDc(referringPatch.deltaCoeffs(), size, start);

    // Get internal field
    tmp<Field<Type> > patchFld(new Field<Type>(size));
    this->patchInternalField(patchFld(), referringPatch, size, start);

    exchangeBuf = subDc * (subFld - patchFld);
}


template<class Type>
void cyclicFvPatchField<Type>::initEvaluate
(
    Field<Type>& exchangeBuf,
    const coupledFvPatch& referringPatch,
    const label size,
    const label start
) const
{
    //? What about updateCoeffs? What if patch holds face-wise additional
    // information? Where is it stored? Who updates it?
//    if (!this->updated())
//    {
//        this->updateCoeffs();
//    }

    if (exchangeBuf.size() != size)
    {
        FatalErrorIn("cyclicFvPatchField<Type>::initEvaluate(..)")
            << "Call with correct slice size." << abort(FatalError);
    }

    // Get sender side. Equivalent to (1-w)*patchNeighbourField() in non-remote
    // version (so includes the transformation!).

    //Pout<< "initEvaluate name:" << cyclicPatch_.name()
    //    << " size:" << size << " start:" << start
    //    << " referringPatch.weights():" << referringPatch.weights() << endl;

    const SubField<scalar> subWeights
    (
        referringPatch.weights(),
        size,
        start
    );

    tmp<Field<Type> > patchFld(new Field<Type>(size));
    this->patchInternalField(patchFld(), referringPatch, size, start);

    //Pout<< "initEvaluate name:" << cyclicPatch_.name()
    //    << " patchFld:" << patchFld()
    //    << " subWeights:" << subWeights << endl;

    if (doTransform())
    {
        //Pout<< "initEvaluate name:" << cyclicPatch_.name()
        //    << " reverseT:" << reverseT() << endl;
        tmp<Field<Type> > tfld =
            (1.0-subWeights)
          * transform(reverseT(), patchFld);

        forAll(tfld(), i)
        {
            exchangeBuf[i] = tfld()[i];
        }
        //Pout<< "initEvaluate name:" << cyclicPatch_.name()
        //    << " exchangeBuf:" << exchangeBuf << endl;
    }
    else
    {
        //Pout<< "initEvaluate name:" << cyclicPatch_.name()
        //    << " no transform" << endl;
        tmp<Field<Type> > tfld = (1.0-subWeights)*patchFld;

        forAll(tfld(), i)
        {
            exchangeBuf[i] = tfld()[i];
        }
        //Pout<< "initEvaluate name:" << cyclicPatch_.name()
        //    << " exchangeBuf:" << exchangeBuf << endl;
    }
}


template<class Type>
void cyclicFvPatchField<Type>::evaluate
(
    Field<Type>& exchangeBuf,
    const coupledFvPatch& referringPatch,
    const label size,
    const label start
) const
{
//    if (!this->updated())
//    {
//        this->updateCoeffs();
//    }

    if (exchangeBuf.size() != size)
    {
        FatalErrorIn("cyclicFvPatchField<Type>::evaluate(..)")
            << "Call with correct slice size." << abort(FatalError);
    }

    const SubField<scalar> subWeights
    (
        referringPatch.weights(),
        size,
        start
    );

    tmp<Field<Type> > patchFld(new Field<Type>(size));
    this->patchInternalField(patchFld(), referringPatch, size, start);

    exchangeBuf += subWeights * patchFld;

    //?? fvPatchField<Type>::evaluate();
}


template<class Type>
void cyclicFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction,
    const coupledFvPatch& referringPatch,
    const label size,
    const label start,
    scalarField& exchangeBuf
) const
{
    const unallocLabelList& faceCells = referringPatch.faceCells();

    label facei = start;

    forAll(exchangeBuf, elemI)
    {
        exchangeBuf[elemI] = psiInternal[faceCells[facei++]];
    }
}


template<class Type>
void cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const coupledFvPatch& referringPatch,
    const label size,
    const label start,
    scalarField& exchangeBuf
) const
{
    // Transform according to the transformation tensor
    transformCoupleField(exchangeBuf, cmpt);

    // Multiply the field by coefficients and add into the result

    const unallocLabelList& faceCells = referringPatch.faceCells();

    label facei = start;

    forAll(exchangeBuf, elemI)
    {
        result[faceCells[facei]] -= coeffs[facei]*exchangeBuf[elemI];
        facei++;
    }
}


// Local patch functionality
// ~~~~~~~~~~~~~~~~~~~~~~~~~

template<class Type>
tmp<Field<Type> > cyclicFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->internalField();
    const unallocLabelList& nbrFaceCells =
        cyclicPatch().cyclicPatch().neighbPatch().faceCells();

    tmp<Field<Type> > tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf();


    if (doTransform())
    {
        forAll(pnf, facei)
        {
            pnf[facei] = transform
            (
                forwardT(), iField[nbrFaceCells[facei]]
            );
        }
    }
    else
    {
        forAll(pnf, facei)
        {
            pnf[facei] = iField[nbrFaceCells[facei]];
        }
    }

    return tpnf;
}


template<class Type>
void cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    scalarField pnf(this->size());

    const unallocLabelList& nbrFaceCells =
        cyclicPatch().cyclicPatch().neighbPatch().faceCells();

    forAll(pnf, facei)
    {
        pnf[facei] = psiInternal[nbrFaceCells[facei]];
    }

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& faceCells = cyclicPatch_.faceCells();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
