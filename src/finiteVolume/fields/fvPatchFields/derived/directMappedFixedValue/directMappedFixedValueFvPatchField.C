/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "directMappedFixedValueFvPatchField.H"
#include "directMappedFvPatch.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    setAverage_(false),
    average_(pTraits<Type>::zero)
{}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const directMappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{
    if (!isType<directMappedFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "directMappedFixedValueFvPatchField<Type>::"
            "directMappedFixedValueFvPatchField\n"
            "(\n"
            "    const directMappedFixedValueFvPatchField<Type>&,\n"
            "    const fvPatch&,\n"
            "    const Field<Type>&,\n"
            "    const fvPatchFieldMapper&\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    setAverage_(readBool(dict.lookup("setAverage"))),
    average_(pTraits<Type>(dict.lookup("average")))
{
    if (!isType<directMappedFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "directMappedFixedValueFvPatchField<Type>::"
            "directMappedFixedValueFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const directMappedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{}


template<class Type>
directMappedFixedValueFvPatchField<Type>::directMappedFixedValueFvPatchField
(
    const directMappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void directMappedFixedValueFvPatchField<Type>::getNewValues
(
    const directMappedPolyPatch& mpp,
    const Field<Type>& sendValues,
    Field<Type>& newValues
) const
{
    // Get the scheduling information
    const List<labelPair>& schedule = mpp.schedule();
    const labelListList& sendLabels = mpp.sendLabels();
    const labelListList& receiveFaceLabels = mpp.receiveFaceLabels();

    forAll(schedule, i)
    {
        const labelPair& twoProcs = schedule[i];
        label sendProc = twoProcs[0];
        label recvProc = twoProcs[1];

        if (Pstream::myProcNo() == sendProc)
        {
            OPstream toProc(Pstream::scheduled, recvProc);
            toProc<< IndirectList<Type>(sendValues, sendLabels[recvProc])();
        }
        else
        {
            // I am receiver. Receive from sendProc.
            IPstream fromProc(Pstream::scheduled, sendProc);

            Field<Type> fromFld(fromProc);

            // Destination faces
            const labelList& faceLabels = receiveFaceLabels[sendProc];

            forAll(fromFld, i)
            {
                label patchFaceI = faceLabels[i];

                newValues[patchFaceI] = fromFld[i];
            }
        }
    }

    // Do data from myself
    {
        IndirectList<Type> fromFld(sendValues, sendLabels[Pstream::myProcNo()]);

        // Destination faces
        const labelList& faceLabels = receiveFaceLabels[Pstream::myProcNo()];

        forAll(fromFld, i)
        {
            label patchFaceI = faceLabels[i];

            newValues[patchFaceI] = fromFld[i];
        }
    }
}


template<class Type>
void directMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Get the directMappedPolyPatch
    const directMappedPolyPatch& mpp = refCast<const directMappedPolyPatch>
    (
        directMappedFixedValueFvPatchField<Type>::patch().patch()
    );

    Field<Type> newValues(this->size());

    switch (mpp.mode())
    {
        case directMappedPolyPatch::NEARESTCELL:
        {
            getNewValues(mpp, this->internalField(), newValues);

            break;
        }
        case directMappedPolyPatch::NEARESTPATCHFACE:
        {
            const label patchID =
                 this->patch().patch().boundaryMesh().findPatchID
                 (
                    mpp.samplePatch()
                 );
            if (patchID < 0)
            {
                FatalErrorIn
                (
                    "void directMappedFixedValueFvPatchField<Type>::"
                    "updateCoeffs()"
                )<< "Unable to find sample patch " << mpp.samplePatch()
                 << " for patch " << this->patch().name() << nl
                 << abort(FatalError);
            }
            typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
            const word& fieldName = this->dimensionedInternalField().name();
            const fieldType& sendField =
                this->db().objectRegistry::lookupObject<fieldType>(fieldName);

            getNewValues(mpp, sendField.boundaryField()[patchID], newValues);

            break;
        }
        case directMappedPolyPatch::NEARESTFACE:
        {
            typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
            const word& fieldName = this->dimensionedInternalField().name();
            const fieldType& sendField =
                this->db().objectRegistry::lookupObject<fieldType>(fieldName);

            Field<Type> allValues
            (
                this->patch().patch().boundaryMesh().mesh().nFaces(),
                pTraits<Type>::zero
            );

            forAll(sendField.boundaryField(), patchI)
            {
                const fvPatchField<Type>& pf =
                    sendField.boundaryField()[patchI];
                label faceStart = pf.patch().patch().start();

                forAll(pf, faceI)
                {
                    allValues[faceStart++] = pf[faceI];
                }
            }

            getNewValues(mpp, allValues, newValues);

            newValues = this->patch().patchSlice(newValues);

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "directMappedFixedValueFvPatchField<Type>::updateCoeffs()"
            )<< "Unknown sampling mode: " << mpp.mode()
             << nl << abort(FatalError);
        }
    }

    if (setAverage_)
    {
        Type averagePsi =
            gSum(this->patch().magSf()*newValues)
           /gSum(this->patch().magSf());

        if (mag(averagePsi)/mag(average_) > 0.5)
        {
            newValues *= mag(average_)/mag(averagePsi);
        }
        else
        {
            newValues += (average_ - averagePsi);
        }
    }

    this->operator==(newValues);

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void directMappedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;
    os.writeKeyword("average") << average_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
