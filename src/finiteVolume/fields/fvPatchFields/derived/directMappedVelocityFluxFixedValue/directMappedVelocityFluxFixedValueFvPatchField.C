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

#include "directMappedVelocityFluxFixedValueFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "directMappedFvPatch.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("undefinedPhi")
{}


directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const directMappedVelocityFluxFixedValueFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{
    if (!isType<directMappedFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "directMappedVelocityFluxFixedValueFvPatchField::"
            "directMappedFixedValueFvPatchField\n"
            "(\n"
            "    const directMappedVelocityFluxFixedValueFvPatchField&,\n"
            "    const fvPatch&,\n"
            "    const DimensionedField<vector, volMesh>&,\n"
            "    const fvPatchFieldMapper&\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.lookup("phi"))
{
    if (!isType<directMappedFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "directMappedVelocityFluxFixedValueFvPatchField::"
            "directMappedFixedValueFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const directMappedVelocityFluxFixedValueFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    phiName_(ptf.phiName_)
{}


directMappedVelocityFluxFixedValueFvPatchField::
directMappedVelocityFluxFixedValueFvPatchField
(
    const directMappedVelocityFluxFixedValueFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_)
{}


void directMappedVelocityFluxFixedValueFvPatchField::getNewValues
(
    const directMappedPolyPatch& mpp,
    const vectorField& sendUValues,
    const scalarField& sendPhiValues,
    vectorField& newUValues,
    scalarField& newPhiValues
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
            toProc<< IndirectList<vector>(sendUValues, sendLabels[recvProc])();
            toProc<< IndirectList<scalar>
                (
                    sendPhiValues,
                    sendLabels[recvProc]
                )();
        }
        else
        {
            // I am receiver. Receive from sendProc.
            IPstream fromProc(Pstream::scheduled, sendProc);

            vectorField fromUFld(fromProc);
            scalarField fromPhiFld(fromProc);

            // Destination faces
            const labelList& faceLabels = receiveFaceLabels[sendProc];

            forAll(fromUFld, i)
            {
                label patchFaceI = faceLabels[i];

                newUValues[patchFaceI] = fromUFld[i];
                newPhiValues[patchFaceI] = fromPhiFld[i];
            }
        }
    }

    // Do data from myself
    {
        IndirectList<vector> fromUFld
            (
                sendUValues,
                sendLabels[Pstream::myProcNo()]
            );

        IndirectList<scalar> fromPhiFld
            (
                sendPhiValues,
                sendLabels[Pstream::myProcNo()]
            );

        // Destination faces
        const labelList& faceLabels = receiveFaceLabels[Pstream::myProcNo()];

        forAll(fromUFld, i)
        {
            label patchFaceI = faceLabels[i];

            newUValues[patchFaceI] = fromUFld[i];
            newPhiValues[patchFaceI] = fromPhiFld[i];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void directMappedVelocityFluxFixedValueFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the directMappedPolyPatch
    const directMappedPolyPatch& mpp = refCast<const directMappedPolyPatch>
    (
        directMappedVelocityFluxFixedValueFvPatchField::patch().patch()
    );

    vectorField newUValues(size());
    scalarField newPhiValues(size());

    const word& fieldName = dimensionedInternalField().name();
    const volVectorField& UField = db().lookupObject<volVectorField>(fieldName);
    surfaceScalarField& phiField =
    const_cast<surfaceScalarField&>
        (
            db().lookupObject<surfaceScalarField>(phiName_)
        );

    switch (mpp.mode())
    {
        case directMappedPolyPatch::NEARESTFACE:
        {
            vectorField allUValues
            (
                patch().patch().boundaryMesh().mesh().nFaces(),
                vector::zero
            );
            scalarField allPhiValues
            (
                patch().patch().boundaryMesh().mesh().nFaces(),
                0.0
            );

            forAll(UField.boundaryField(), patchI)
            {
                const fvPatchVectorField& Upf = UField.boundaryField()[patchI];
                const scalarField& phipf = phiField.boundaryField()[patchI];

                label faceStart = Upf.patch().patch().start();

                forAll(Upf, faceI)
                {
                    allUValues[faceStart++] = Upf[faceI];
                    allPhiValues[faceStart] = phipf[faceI];
                }
            }

            getNewValues
            (
                mpp,
                allUValues,
                allPhiValues,
                newUValues,
                newPhiValues
            );

            newUValues = patch().patchSlice(newUValues);
            newPhiValues = patch().patchSlice(newPhiValues);

            break;
        }
        case directMappedPolyPatch::NEARESTPATCHFACE:
        {
            const label patchID =
                patch().patch().boundaryMesh().findPatchID
                (
                    mpp.samplePatch()
                );

            getNewValues
            (
                mpp,
                UField.boundaryField()[patchID],
                phiField.boundaryField()[patchID],
                newUValues,
                newPhiValues
            );

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "directMappedVelocityFluxFixedValueFvPatchField::updateCoeffs()"
            )<< "patch can only be used in NEARESTPATCHFACE or NEARESTFACE "
             << "mode" << nl << abort(FatalError);
        }
    }

    operator==(newUValues);
    phiField.boundaryField()[patch().index()] == newPhiValues;

    fixedValueFvPatchVectorField::updateCoeffs();
}


void directMappedVelocityFluxFixedValueFvPatchField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    directMappedVelocityFluxFixedValueFvPatchField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
