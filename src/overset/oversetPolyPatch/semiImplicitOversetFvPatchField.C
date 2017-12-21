/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
Foam::semiImplicitOversetFvPatchField<Type>::semiImplicitOversetFvPatchField
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
Foam::semiImplicitOversetFvPatchField<Type>::semiImplicitOversetFvPatchField
(
    const semiImplicitOversetFvPatchField<Type>& ptf,
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
Foam::semiImplicitOversetFvPatchField<Type>::semiImplicitOversetFvPatchField
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
Foam::semiImplicitOversetFvPatchField<Type>::semiImplicitOversetFvPatchField
(
    const semiImplicitOversetFvPatchField<Type>& ptf
)
:
    LduInterfaceField<Type>(ptf.oversetPatch_),
    zeroGradientFvPatchField<Type>(ptf),
    oversetPatch_(ptf.oversetPatch_)
{}


template<class Type>
Foam::semiImplicitOversetFvPatchField<Type>::semiImplicitOversetFvPatchField
(
    const semiImplicitOversetFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(ptf.oversetPatch_),
    zeroGradientFvPatchField<Type>(ptf, iF),
    oversetPatch_(ptf.oversetPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::semiImplicitOversetFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << " field " <<  this->internalField().name()
            << " patch " << this->patch().name() << endl;
    }


    if (!this->updated())
    {
        this->updateCoeffs();
    }

    zeroGradientFvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::semiImplicitOversetFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << " field " <<  this->internalField().name()
            << " patch " << this->patch().name() << endl;
    }

    if (this->oversetPatch().master())
    {
        // Trigger interpolation
        const fvMesh& mesh = this->internalField().mesh();
        //const dictionary& fvSchemes = mesh.schemesDict();
        const word& fldName = this->internalField().name();

        if (&mesh.lduAddr() != &mesh.fvMesh::lduAddr())
        {
            // Running extended addressing so called from within fvMatrix::solve
            FatalErrorInFunction
                << "Running extended addressing is not allowed when solving "
                << fldName
                << " Please choose a dynamicFvMesh without matrix adaptation"
                << exit(FatalError);
        }
        else
        {
            if (debug)
            {
                Info<< "Interpolating field " << fldName << endl;
            }

            // Interpolate without bc update (since would trigger infinite
            // recursion back to
            // semiImplicitOversetFvPatchField<Type>::evaluate)
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
    }
}


template<class Type>
void Foam::semiImplicitOversetFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << " field " <<  this->internalField().name()
            << " patch " << this->patch().name() << endl;
    }

    const oversetFvPatch& ovp = this->oversetPatch();

    // Set all interpolated cells
    if (ovp.master())
    {
        const labelListList& stencil = ovp.stencil();

        if (stencil.size() != psiInternal.size())
        {
            FatalErrorInFunction << "psiInternal:" << psiInternal.size()
                << " stencil:" << stencil.size() << exit(FatalError);
        }

        const mapDistribute& map = ovp.cellInterpolationMap();
        const List<scalarList>& wghts = ovp.cellInterpolationWeights();
        const labelList& cellIDs = ovp.interpolationCells();
        //const scalarList& factor = ovp.cellInterpolationWeight();

        // Since we're inside initEvaluate/evaluate there might be processor
        // comms underway. Change the tag we use.
        scalarField work(psiInternal);
        map.mapDistributeBase::distribute(work, UPstream::msgType()+1);

        forAll(cellIDs, i)
        {
            label celli = cellIDs[i];

            const scalarList& w = wghts[celli];
            const labelList& nbrs = stencil[celli];

            //scalar f = factor[celli];

            scalar s(0.0);
            forAll(nbrs, nbrI)
            {
                s += w[nbrI]*work[nbrs[nbrI]];
            }

            //Pout<< "cell:" << celli << " interpolated value:" << s << endl;
            result[celli] = s;  //(1.0-f)*result[celli] + f*s;
        }
    }
}


template<class Type>
void Foam::semiImplicitOversetFvPatchField<Type>::write(Ostream& os) const
{
    zeroGradientFvPatchField<Type>::write(os);
    // Make sure to write the value for ease of postprocessing.
    this->writeEntry("value", os);
}


// ************************************************************************* //
