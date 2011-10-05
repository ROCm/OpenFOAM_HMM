/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "rmap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::rmap
(
    fvPatchField<Type>& destBC,
    const labelList& reverseAddressing,
    const fvPatchField<Type>& sourceBC
)
{
    // Assign value
    destBC.Field<Type>::rmap(sourceBC, reverseAddressing);

    // Assign other properties
    if (isA<mixedFvPatchField<Type> >(destBC))
    {
        mixedFvPatchField<Type>& mp =
            refCast<mixedFvPatchField<Type> >(destBC);

        if (isA<mixedFvPatchField<Type> >(sourceBC))
        {
            const mixedFvPatchField<Type>& Tp =
                refCast<const mixedFvPatchField<Type> >(sourceBC);

            mp.refValue().rmap(Tp.refValue(), reverseAddressing);
            mp.refGrad().rmap(Tp.refGrad(), reverseAddressing);
            mp.valueFraction().rmap(Tp.valueFraction(), reverseAddressing);
        }
        else if (isA<fixedGradientFvPatchField<Type> >(sourceBC))
        {
            const fixedGradientFvPatchField<Type>& Tp =
                refCast<const fixedGradientFvPatchField<Type> >
                (
                    sourceBC
                );
            // Make pure fixedGradient
            mp.refValue().rmap(Tp, reverseAddressing); // unused
            mp.refGrad().rmap(Tp.gradient(), reverseAddressing);
            mp.valueFraction().rmap
            (
                Field<Type>(reverseAddressing.size(), 0.0),
                reverseAddressing
            );
        }
        else if (isA<zeroGradientFvPatchField<Type> >(sourceBC))
        {
            // Make pure fixedGradient with gradient = 0
            mp.refValue().rmap(sourceBC, reverseAddressing); // unused
            mp.refGrad().rmap
            (
                Field<Type>(reverseAddressing.size(), 0.0),
                reverseAddressing
            );
            mp.valueFraction().rmap
            (
                Field<Type>(reverseAddressing.size(), 0.0),
                reverseAddressing
            );
        }
        else if (isA<fixedValueFvPatchField<Type> >(sourceBC))
        {
            // Make pure fixedValue
            mp.refValue().rmap(sourceBC, reverseAddressing);
            mp.refGrad().rmap
            (
                Field<Type>(reverseAddressing.size(), 0.0),
                reverseAddressing
            ); // unused
            mp.valueFraction().rmap
            (
                Field<Type>(reverseAddressing.size(), 1.0),
                reverseAddressing
            );
        }
        else if (isA<calculatedFvPatchField<Type> >(sourceBC))
        {
            // Make pure fixedValue
            mp.refValue().rmap(sourceBC, reverseAddressing);
            mp.refGrad().rmap
            (
                Field<Type>(reverseAddressing.size(), 0.0),
                reverseAddressing
            ); // unused
            mp.valueFraction().rmap
            (
                Field<Type>(reverseAddressing.size(), 1.0),
                reverseAddressing
            );
        }
        else
        {
            FatalErrorIn("rmap(..)")
                << "Don't know how to map source bc "
                << sourceBC.type()
                << " into a mixed boundary condition at "
                << destBC.patch().name()
                << exit(FatalError);
        }
    }
    else if (isA<fixedGradientFvPatchField<Type> >(destBC))
    {
        fixedGradientFvPatchField<Type>& mp =
            refCast<fixedGradientFvPatchField<Type> >(destBC);

        if (isA<fixedGradientFvPatchField<Type> >(sourceBC))
        {
            const fixedGradientFvPatchField<Type>& Tp =
                refCast<const fixedGradientFvPatchField<Type> >
                (
                    sourceBC
                );
            mp.gradient().rmap(Tp.gradient(), reverseAddressing);
        }
        else if (isA<mixedFvPatchField<Type> >(sourceBC))
        {
            const mixedFvPatchField<Type>& Tp =
                refCast<const mixedFvPatchField<Type> >(sourceBC);
            mp.gradient().rmap(Tp.snGrad(), reverseAddressing);
        }
        else if (isA<zeroGradientFvPatchField<Type> >(sourceBC))
        {
            mp.gradient().rmap
            (
                Field<Type>(reverseAddressing.size(), 0.0),
                reverseAddressing
            );
        }
        else
        {
            FatalErrorIn("rmap(..)")
                << "Don't know how to map source bc "
                << sourceBC.type()
                << " into a fixedGradient boundary condition at "
                << destBC.patch().name()
                << exit(FatalError);
        }
    }
}


template<class Type>
void Foam::rmap
(
    GeometricField<Type, fvPatchField, volMesh>& dest,
    const GeometricField<Type, fvPatchField, volMesh>& source,
    const labelList& faceProcAddressing,
    const labelList& cellProcAddressing,
    const labelList& boundaryProcAddressing
)
{
    if (dest.dimensions() != source.dimensions())
    {
        FatalErrorIn("rmap(..)")
            << "Different dimensions for = for fields " << dest.name()
            << " and " << source.name() << endl
            << "     dimensions : " << dest.dimensions()
            << " = " << source.dimensions() << endl
            << exit(FatalError);
    }

    // Copy internal field
    dest.internalField().rmap(source.internalField(), cellProcAddressing);

    // Copy boundary properties as mixed
    forAll(source.boundaryField(), patchI)
    {
        label curBPatch = boundaryProcAddressing[patchI];

        if (curBPatch == -1)
        {
            // Unknown patch. Do not change any values.
        }
        else
        {
            // Get addressing slice for this patch
            const labelList::subList cp =
                source.mesh().boundary()[patchI].patchSlice
                (
                    faceProcAddressing
                );

            const label curPatchStart =
                dest.mesh().boundaryMesh()[curBPatch].start();

            labelList reverseAddressing(cp.size());

            forAll(cp, faceI)
            {
                // Subtract one to take into account offsets for
                // face direction.
                if (cp[faceI] <= 0)
                {
                    FatalErrorIn("rmap(..)")
                        << "Problem:"
                        << " patch:" << source.mesh().boundary()[patchI].name()
                        << " field:" << source.name()
                        << " local face:" << faceI
                        << " mapped to:" << cp[faceI] << exit(FatalError);
                }

                reverseAddressing[faceI] = cp[faceI] - 1 - curPatchStart;
            }

            // Map curBPatch from source patch. Is like rmap but also
            // copies non-value properties from alike patchFields.
            rmap
            (
                dest.boundaryField()[curBPatch],
                reverseAddressing,
                source.boundaryField()[patchI]
            );
        }
    }

    // Copy timeIndex
    dest.timeIndex() = source.timeIndex();
}


template<class Type>
void Foam::rmap
(
    fvMatrix<Type>& dest,
    const fvMatrix<Type>& source,
    const labelList& faceProcAddressing,
    const labelList& cellProcAddressing,
    const labelList& boundaryProcAddressing
)
{
    dest.source().rmap(source.source(), cellProcAddressing);

    FieldField<Field, Type>& sourceInternal =
        const_cast<fvMatrix<Type>&>(source).internalCoeffs();
    FieldField<Field, Type>& sourceBoundary =
        const_cast<fvMatrix<Type>&>(source).boundaryCoeffs();

    forAll(sourceInternal, patchI)
    {
        label curBPatch = boundaryProcAddressing[patchI];

        if (curBPatch == -1)
        {
            // Unknown patch. Do not change any values.
        }
        else
        {
            // Get addressing slice for this patch
            const fvMesh& sourceMesh = source.psi().mesh();

            const labelList::subList cp =
                sourceMesh.boundary()[patchI].patchSlice
                (
                    faceProcAddressing
                );

            const label curPatchStart =
                dest.psi().mesh().boundaryMesh()[curBPatch].start();

            labelList reverseAddressing(cp.size());

            forAll(cp, faceI)
            {
                // Subtract one to take into account offsets for
                // face direction.
                reverseAddressing[faceI] = cp[faceI] - 1 - curPatchStart;
            }
            dest.internalCoeffs()[curBPatch].rmap
            (
                sourceInternal[patchI],
                reverseAddressing
            );
            dest.boundaryCoeffs()[curBPatch].rmap
            (
                sourceBoundary[patchI],
                reverseAddressing
            );
        }
    }
}


// ************************************************************************* //
