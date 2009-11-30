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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "demandDrivenData.H"
#include "transformField.H"
#include "diagTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const coupledFvPatchField<Type>& processorFvPatchField<Type>::patchField
(
    const label patchID
) const
{
    //const GeometricField<Type, fvPatchField, volMesh>& field =
    //this->db().objectRegistry::
    //lookupObject<GeometricField<Type, fvPatchField, volMesh> >
    //(
    //    this->dimensionedInternalField().name()
    //);
    const GeometricField<Type, fvPatchField, volMesh>& field =
        static_cast
        <
            const GeometricField<Type, fvPatchField, volMesh>&
        >(this->dimensionedInternalField());

    return refCast<const coupledFvPatchField<Type> >
    (
        field.boundaryField()[patchID]
    );
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p))
{}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p))
{}


// Construct by mapping given processorFvPatchField<Type>
template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isType<processorFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
            "(\n"
            "    const processorFvPatchField<Type>& ptf,\n"
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
processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isType<processorFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "processorFvPatchField<Type>::processorFvPatchField\n"
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
}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{}


template<class Type>
processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
processorFvPatchField<Type>::~processorFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > processorFvPatchField<Type>::patchNeighbourField() const
{
    return *this;
}


template<class Type>
void processorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        const processorPolyPatch& pp = procPatch_.procPolyPatch();

        // Get reference to proc patch built-in buffer
        List<Type>& sendBuf = procPatch_.setSendBuf<Type>(this->size());

        forAll(pp.patchIDs(), i)
        {
            SubField<Type> subSendFld(pp.subSlice(sendBuf, i));
            Field<Type>& subFld = static_cast<Field<Type>&>
            (
                static_cast<UList<Type>&>(subSendFld)
            );

            label patchI = pp.patchIDs()[i];

            //Pout<< "initEvaluate on "
            //    << this->dimensionedInternalField().name()
            //    << " patch:" << pp.name()
            //    << " subSize:" << (pp.starts()[i+1]-pp.starts()[i])
            //    << " subStart:" << pp.starts()[i]
            //    << " subPatch:" << patchI << endl;

            if (patchI == -1)
            {
                // Assign internal field
                this->patchInternalField
                (
                    subFld,
                    procPatch_,
                    subSendFld.size(),
                    pp.starts()[i]
                );
            }
            else
            {
                // Assign evaluation of referred patch
                patchField(patchI).initEvaluate
                (
                    subFld,
                    procPatch_,
                    subSendFld.size(),
                    pp.starts()[i]
                );
            }
        }
        procPatch_.compressedBufferSend<Type>(commsType);
    }
}


template<class Type>
void processorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        procPatch_.compressedReceive<Type>(commsType, *this);

        const processorPolyPatch& pp = procPatch_.procPolyPatch();

        forAll(pp.patchIDs(), i)
        {
            label patchI = pp.patchIDs()[i];

            //Pout<< "evaluate on " << this->dimensionedInternalField().name()
            //    << " patch:" << pp.name()
            //    << " subSize:" << (pp.starts()[i+1]-pp.starts()[i])
            //    << " subStart:" << pp.starts()[i]
            //    << " subPatch:" << patchI << endl;

            if (patchI == -1)
            {
                // No evaluation needed.
            }
            else
            {
                SubField<Type> subRecvFld(pp.subSlice(*this, i));

                Field<Type>& subFld = static_cast<Field<Type>&>
                (
                    static_cast<UList<Type>&>(subRecvFld)
                );

                patchField(patchI).evaluate
                (
                    subFld,
                    procPatch_,
                    subRecvFld.size(),
                    pp.starts()[i]
                );
            }
        }
    }
}


template<class Type>
tmp<Field<Type> > processorFvPatchField<Type>::snGrad() const
{
    tmp<Field<Type> > tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf();

    const processorPolyPatch& pp = procPatch_.procPolyPatch();

    forAll(pp.patchIDs(), i)
    {
        label patchI = pp.patchIDs()[i];
        label subStart = pp.starts()[i];
        label subSize = pp.starts()[i+1] - pp.starts()[i];

        const SubField<Type> subThis(pp.subSlice(*this, i));

        SubField<Type> subPnf(pp.subSlice(pnf, i));
        Field<Type>& subFld = static_cast<Field<Type>&>
        (
            static_cast<UList<Type>&>(subPnf)
        );

        //Pout<< "snGrad on " << this->dimensionedInternalField().name()
        //    << " patch:" << pp.name()
        //    << " subSize:" << (pp.starts()[i+1]-pp.starts()[i])
        //    << " subStart:" << pp.starts()[i]
        //    << " subPatch:" << patchI << endl;


        if (patchI == -1)
        {
            // Slice delta coeffs
            const SubField<scalar> subDc
            (
                pp.subSlice(procPatch_.deltaCoeffs(), i)
            );
            tmp<Field<Type> > subInt(new Field<Type>(subSize));
            this->patchInternalField(subInt(), procPatch_, subSize, subStart);

            subFld = (subDc*(subThis-subInt))();
        }
        else
        {
            patchField(patchI).snGrad
            (
                subFld,
                subThis,
                procPatch_,
                subSize,
                subStart
            );
        }
    }

    return tpnf;
}


template<class Type>
void processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Get reference to proc patch built-in buffer
    List<scalar>& sendFld = procPatch_.setSendBuf<scalar>(this->size());

    const processorPolyPatch& pp = procPatch_.procPolyPatch();

    forAll(pp.patchIDs(), i)
    {
        label subStart = pp.starts()[i];
        label subSize = pp.starts()[i+1] - pp.starts()[i];
        SubField<scalar> subSendFld(sendFld, subSize, subStart);
        Field<scalar>& subFld = static_cast<Field<scalar>&>
        (
            static_cast<UList<scalar>&>(subSendFld)
        );

        label patchI = pp.patchIDs()[i];

        //Pout<< "initInterfaceMatrixUpdate on "
        //    << this->dimensionedInternalField().name()
        //    << " patch:" << pp.name()
        //    << " subSize:" << (pp.starts()[i+1]-pp.starts()[i])
        //    << " subStart:" << pp.starts()[i]
        //    << " subPatch:" << patchI << endl;

        if (patchI == -1)
        {
            const unallocLabelList& faceCells = pp.faceCells();

            label facei = subStart;

            forAll(subFld, i)
            {
                subFld[i] = psiInternal[faceCells[facei++]];
            }
        }
        else
        {
            patchField(patchI).initInterfaceMatrixUpdate
            (
                psiInternal,
                result,
                m,
                coeffs,
                cmpt,
                procPatch_,
                subSize,
                subStart,
                subFld
            );
        }
    }

    procPatch_.compressedBufferSend<scalar>(commsType);
}


template<class Type>
void processorFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix& m,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    const List<scalar>& recvFld = procPatch_.compressedBufferReceive<scalar>
    (
        commsType,
        this->size()
    );

    const processorPolyPatch& pp = procPatch_.procPolyPatch();

    forAll(pp.patchIDs(), i)
    {
        label subStart = pp.starts()[i];
        label subSize = pp.starts()[i+1] - pp.starts()[i];

        SubField<scalar> subRecvFld(recvFld, subSize, subStart);
        Field<scalar>& subFld = static_cast<Field<scalar>&>
        (
            static_cast<UList<scalar>&>(subRecvFld)
        );

        label patchI = pp.patchIDs()[i];

        //Pout<< "updateInterfaceMatrix on "
        //    << this->dimensionedInternalField().name()
        //    << " patch:" << pp.name()
        //    << " subSize:" << (pp.starts()[i+1]-pp.starts()[i])
        //    << " subStart:" << pp.starts()[i]
        //    << " subPatch:" << patchI << endl;

        if (patchI == -1)
        {
            const unallocLabelList& faceCells = pp.faceCells();

            label facei = subStart;

            forAll(subFld, elemI)
            {
                result[faceCells[facei]] -= coeffs[facei]*subFld[elemI];

                facei++;
            }
        }
        else
        {
            patchField(patchI).updateInterfaceMatrix
            (
                psiInternal,
                result,
                m,
                coeffs,
                cmpt,
                procPatch_,
                subSize,
                subStart,
                subFld
            );
        }
    }
}


//template<class Type>
//void Foam::processorFvPatchField<Type>::transformCoupleField
//(
//    scalarField& f,
//    const direction cmpt
//) const
//{
//    if (pTraits<Type>::rank > 0)
//    {
//        const processorPolyPatch& pp = procPatch_.procPolyPatch();
//
//        forAll(pp.patchIDs(), i)
//        {
//            label patchI = pp.patchIDs()[i];
//
//            if (patchI == -1)
//            {
//                // ? anything needs to be transformed?
//            }
//            else
//            {
//                const coupledPolyPatch& cpp =
//                    refCast<const coupledPolyPatch>(pp.boundaryMesh()[patchI]);
//
//                if (!cpp.parallel())
//                {
//                    const tensor& T = cpp.forwardT();
//
//                    SubField<scalar> subFld(pp.subSlice(f, i));
//                    const scalarField& fld =
//                        static_cast<const scalarField&>(subFld);
//
//                    const_cast<scalarField&>(fld) *=
//                        pow(diag(T).component(cmpt), rank());
//                }
//            }
//        }
//    }
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
