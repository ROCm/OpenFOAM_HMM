/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(refCast<const oversetFvPatch>(p)),
    zeroGradientFvPatchField<Type>(p, iF),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    LduInterfaceField<Type>(refCast<const oversetFvPatch>(p)),
    zeroGradientFvPatchField<Type>(ptf, p, iF, mapper),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{
    if (!isA<oversetFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    LduInterfaceField<Type>(refCast<const oversetFvPatch>(p)),
    zeroGradientFvPatchField<Type>(p, iF, dict),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{
    if (!isA<oversetFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value") && this->coupled())
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf
)
:
    LduInterfaceField<Type>(ptf.oversetPatch_),
    zeroGradientFvPatchField<Type>(ptf),
    oversetPatch_(ptf.oversetPatch_)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(ptf.oversetPatch_),
    zeroGradientFvPatchField<Type>(ptf, iF),
    oversetPatch_(ptf.oversetPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::oversetFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (oversetPatch_.master())
    {
        // Trigger interpolation
        const fvMesh& mesh = this->internalField().mesh();
        const dictionary& fvSchemes = mesh.schemesDict();
        const word& fldName = this->internalField().name();

        if (&mesh.lduAddr() != &mesh.fvMesh::lduAddr())
        {
            // Running extended addressing. Interpolate always.
            if (debug)
            {
                Info<< "Interpolating solved-for field " << fldName << endl;
            }

            // Interpolate without bc update (since would trigger infinite
            // recursion back to oversetFvPatchField<Type>::evaluate)
            // The only problem is bcs that use the new cell values
            // (e.g. zeroGradient, processor). These need to appear -after-
            // the 'overset' bc.
            mesh.interpolate
            (
                const_cast<Field<Type>&>
                (
                    this->primitiveField()
                )
            );
        }
        else if
        (
           !fvSchemes.found("oversetInterpolation")
        || !fvSchemes.found("oversetInterpolationRequired")
        )
        {
            IOWarningInFunction(fvSchemes)
                << "Missing required dictionary entries"
                << " 'oversetInterpolation' and 'oversetInterpolationRequired'"
                << ". Skipping overset interpolation for field "
                << fldName << endl;
        }
        else
        {
            const dictionary& intDict = fvSchemes.subDict
            (
                "oversetInterpolationRequired"
            );
            if (intDict.found(fldName))
            {
                if (debug)
                {
                    Info<< "Interpolating field " << fldName << endl;
                }

                // Interpolate without bc update (since would trigger infinite
                // recursion back to oversetFvPatchField<Type>::evaluate)
                // The only problem is bcs that use the new cell values
                // (e.g. zeroGradient, processor). These need to appear -after-
                // the 'overset' bc.
                mesh.interpolate
                (
                    const_cast<Field<Type>&>
                    (
                        this->primitiveField()
                    )
                );
            }
            else if (debug)
            {
                Info<< "Skipping overset interpolation for field "
                    << fldName << endl;
            }
        }
    }
}


template<class Type>
void Foam::oversetFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    zeroGradientFvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::oversetFvPatchField<Type>::initInterfaceMatrixUpdate
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{}


template<class Type>
void Foam::oversetFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // Add remote values
    if (oversetPatch_.master())
    {
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();

        // Try to find out if the solve routine comes from the mesh
        // TBD. This should be cleaner.
        if (&mesh.lduAddr() == &mesh.fvMesh::lduAddr())
        {
            //Pout<< "** not running extended addressing..." << endl;
            return;
        }

        const labelListList& stencil = oversetPatch_.stencil();

        if (stencil.size() != psiInternal.size())
        {
            FatalErrorInFunction << "psiInternal:" << psiInternal.size()
                << " stencil:" << stencil.size() << exit(FatalError);
        }

        const mapDistribute& map = oversetPatch_.cellInterpolationMap();
        const List<scalarList>& wghts =
            oversetPatch_.cellInterpolationWeights();
        const labelList& cellIDs = oversetPatch_.interpolationCells();
        const scalarList& factor = oversetPatch_.cellInterpolationWeight();
        const scalarField& normalisation = oversetPatch_.normalisation();


        // Since we're inside initEvaluate/evaluate there might be processor
        // comms underway. Change the tag we use.
        scalarField work(psiInternal);
        map.mapDistributeBase::distribute(work, UPstream::msgType()+1);

        forAll(cellIDs, i)
        {
            label celli = cellIDs[i];

            const scalarList& w = wghts[celli];
            const labelList& nbrs = stencil[celli];
            scalar f = factor[celli];
            const scalar norm = normalisation[celli];

            // Add the non-local matrix contribution to psi. Note that the
            // matrix coefficients are -weights

            if (add)
            {
                f = -1.0*f;
            }

            forAll(nbrs, nbrI)
            {
                label slotI = nbrs[nbrI];
                if (slotI >= psiInternal.size())
                {
                    result[celli] += f*norm*w[nbrI]*work[slotI];
                }
            }
        }
    }
}


template<class Type>
void Foam::oversetFvPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>&,
    const bool add,
    const Field<Type>&,
    const scalarField&,
    const Pstream::commsTypes commsType
) const
{
    NotImplemented;
}


template<class Type>
void Foam::oversetFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField&,
    const Pstream::commsTypes
) const
{
    NotImplemented;
}


template<class Type>
void Foam::oversetFvPatchField<Type>::write(Ostream& os) const
{
    zeroGradientFvPatchField<Type>::write(os);
    // Make sure to write the value for ease of postprocessing.
    this->writeEntry("value", os);
}


// ************************************************************************* //
