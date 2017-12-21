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
#include "cellCellStencil.H"
#include "dynamicOversetFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    semiImplicitOversetFvPatchField<Type>(p, iF)
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
    semiImplicitOversetFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    semiImplicitOversetFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf
)
:
    semiImplicitOversetFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::oversetFvPatchField<Type>::oversetFvPatchField
(
    const oversetFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    semiImplicitOversetFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::oversetFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (this->oversetPatch_.master())
    {
        // Trigger interpolation
        const fvMesh& mesh = this->internalField().mesh();
        const dictionary& fvSchemes = mesh.schemesDict();
        const word& fldName = this->internalField().name();

        if (&mesh.lduAddr() != &mesh.fvMesh::lduAddr())
        {
            // Running extended addressing. Flux correction already done
            // in the linear solver (through the initUpdateInterfaceMatrix
            // method below)
            if (debug)
            {
                Info<< "Skipping overset interpolation for solved-for field "
                    << fldName << endl;
            }
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

    const oversetFvPatch& ovp = this->oversetPatch_;

    if (ovp.master())
    {
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();

        // Try to find out if the solve routine comes from the mesh
        // TBD. This should be cleaner.
        if (&mesh.lduAddr() == &mesh.fvMesh::lduAddr())
        {
            return;
        }

        const labelListList& stencil = ovp.stencil();

        if (stencil.size() != psiInternal.size())
        {
            FatalErrorInFunction << "psiInternal:" << psiInternal.size()
                << " stencil:" << stencil.size() << exit(FatalError);
        }

        const mapDistribute& map = ovp.cellInterpolationMap();
        const List<scalarList>& wghts = ovp.cellInterpolationWeights();
        const labelList& cellIDs = ovp.interpolationCells();
        const scalarList& factor = ovp.cellInterpolationWeight();
        const scalarField& normalisation = ovp.normalisation();


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


// ************************************************************************* //
