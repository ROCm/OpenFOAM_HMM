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

#include "directMappedFieldFvPatchField.H"

#include "volFields.H"
#include "interpolationCell.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
directMappedFieldFvPatchField<Type>::directMappedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    directMappedPatchBase(p.patch()),
    fixedValueFvPatchField<Type>(p, iF),
    fieldName_(iF.name()),
    setAverage_(false),
    average_(pTraits<Type>::zero)
{}


template<class Type>
directMappedFieldFvPatchField<Type>::directMappedFieldFvPatchField
(
    const directMappedFieldFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directMappedPatchBase(p.patch(), ptf),
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{}


template<class Type>
directMappedFieldFvPatchField<Type>::directMappedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    directMappedPatchBase(p.patch(), dict),
    fixedValueFvPatchField<Type>(p, iF, dict),
    fieldName_(dict.template lookupOrDefault<word>("fieldName", iF.name())),
    setAverage_(readBool(dict.lookup("setAverage"))),
    average_(pTraits<Type>(dict.lookup("average")))
{}


template<class Type>
directMappedFieldFvPatchField<Type>::directMappedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,

    // directMappedPatchBase
    const word& sampleRegion,
    const sampleMode sampleMode,
    const word& samplePatch,
    const scalar distance,

    // My settings
    const word& fieldName,
    const bool setAverage,
    const Type average
)
:
    directMappedPatchBase
    (
        p.patch(),
        sampleRegion,
        sampleMode,
        samplePatch,
        distance
    ),
    fixedValueFvPatchField<Type>(p, iF),
    fieldName_(fieldName),
    setAverage_(setAverage),
    average_(average)
{}


template<class Type>
directMappedFieldFvPatchField<Type>::directMappedFieldFvPatchField
(
    const directMappedFieldFvPatchField<Type>& ptf
)
:
    directMappedPatchBase(ptf.patch().patch(), ptf),
    fixedValueFvPatchField<Type>(ptf),
    fieldName_(ptf.fieldName_),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{}


template<class Type>
directMappedFieldFvPatchField<Type>::directMappedFieldFvPatchField
(
    const directMappedFieldFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    directMappedPatchBase(ptf.patch().patch(), ptf),
    fixedValueFvPatchField<Type>(ptf, iF),
    fieldName_(ptf.fieldName_),
    setAverage_(ptf.setAverage_),
    average_(ptf.average_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const GeometricField<Type, fvPatchField, volMesh>&
directMappedFieldFvPatchField<Type>::sampleField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& nbrMesh = refCast<const fvMesh>(sampleMesh());

    if (sameRegion())
    {
        if (fieldName_ == this->dimensionedInternalField().name())
        {
            // Optimisation: bypass field lookup
            return
                dynamic_cast<const fieldType&>
                (
                    this->dimensionedInternalField()
                );
        }
        else
        {
            const fvMesh& thisMesh = this->patch().boundaryMesh().mesh();
            return thisMesh.template lookupObject<fieldType>(fieldName_);
        }
    }
    else
    {
        return nbrMesh.template lookupObject<fieldType>(fieldName_);
    }
}


template<class Type>
void directMappedFieldFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const fvMesh& nbrMesh = refCast<const fvMesh>(sampleMesh());
    const mapDistribute& distMap = directMappedPatchBase::map();

    // Result of obtaining remote values
    Field<Type> newValues;

    switch (mode())
    {
        case NEARESTCELL:
        {
            newValues = sampleField();

            distMap.distribute(newValues);

            break;
        }
        case NEARESTPATCHFACE:
        {
            const label nbrPatchID = nbrMesh.boundaryMesh().findPatchID
            (
                samplePatch()
            );
            if (nbrPatchID < 0)
            {
                FatalErrorIn
                (
                    "void directMappedFieldFvPatchField<Type>::updateCoeffs()"
                )<< "Unable to find sample patch " << samplePatch()
                 << " in region " << sampleRegion()
                 << " for patch " << this->patch().name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues = nbrField.boundaryField()[nbrPatchID];
            distMap.distribute(newValues);

            break;
        }
        case NEARESTFACE:
        {
            Field<Type> allValues(nbrMesh.nFaces(), pTraits<Type>::zero);

            const fieldType& nbrField = sampleField();

            forAll(nbrField.boundaryField(), patchI)
            {
                const fvPatchField<Type>& pf =
                    nbrField.boundaryField()[patchI];
                label faceStart = pf.patch().patch().start();

                forAll(pf, faceI)
                {
                    allValues[faceStart++] = pf[faceI];
                }
            }

            distMap.distribute(allValues);

            newValues.transfer(allValues);

            break;
        }
        default:
        {
            FatalErrorIn("directMappedFieldFvPatchField<Type>::updateCoeffs()")
                << "Unknown sampling mode: " << mode()
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

    if (debug)
    {
        Info<< "operating on field:" << this->dimensionedInternalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void directMappedFieldFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    directMappedPatchBase::write(os);
    os.writeKeyword("fieldName") << fieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("setAverage") << setAverage_ << token::END_STATEMENT << nl;
    os.writeKeyword("average") << average_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
