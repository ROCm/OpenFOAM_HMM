/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd
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

#include "faceSource.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "sampledSurface.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::fieldValues::faceSource::validField(const word& fieldName) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (source_ != stSampledSurface && obr_.foundObject<sf>(fieldName))
    {
        return true;
    }
    else if (obr_.foundObject<vf>(fieldName))
    {
        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::faceSource::setFieldValues
(
    const word& fieldName,
    const bool mustGet,
    const bool applyOrientation
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (source_ != stSampledSurface && obr_.foundObject<sf>(fieldName))
    {
        return filterField(obr_.lookupObject<sf>(fieldName), applyOrientation);
    }
    else if (obr_.foundObject<vf>(fieldName))
    {
        const vf& fld = obr_.lookupObject<vf>(fieldName);

        if (surfacePtr_.valid())
        {
            if (surfacePtr_().interpolate())
            {
                const interpolationCellPoint<Type> interp(fld);
                tmp<Field<Type> > tintFld(surfacePtr_().interpolate(interp));
                const Field<Type>& intFld = tintFld();

                // Average
                const faceList& faces = surfacePtr_().faces();
                tmp<Field<Type> > tavg
                (
                    new Field<Type>(faces.size(), pTraits<Type>::zero)
                );
                Field<Type>& avg = tavg();

                forAll(faces, faceI)
                {
                    const face& f = faces[faceI];
                    forAll(f, fp)
                    {
                        avg[faceI] += intFld[f[fp]];
                    }
                    avg[faceI] /= f.size();
                }

                return tavg;
            }
            else
            {
                return surfacePtr_().sample(fld);
            }
        }
        else
        {
            return filterField(fld, applyOrientation);
        }
    }

    if (mustGet)
    {
        FatalErrorIn
        (
            "Foam::tmp<Foam::Field<Type> > "
            "Foam::fieldValues::faceSource::setFieldValues"
            "("
                "const word&, "
                "const bool, "
                "const bool"
            ") const"
        )   << "Field " << fieldName << " not found in database"
            << abort(FatalError);
    }

    return tmp<Field<Type> >(new Field<Type>(0));
}


template<class Type>
Type Foam::fieldValues::faceSource::processSameTypeValues
(
    const Field<Type>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    Type result = pTraits<Type>::zero;
    switch (operation_)
    {
        case opSum:
        {
            result = gSum(values);
            break;
        }
        case opSumMag:
        {
            result = gSum(cmptMag(values));
            break;
        }
        case opSumDirection:
        {
            FatalErrorIn
            (
                "template<class Type>"
                "Type Foam::fieldValues::faceSource::processSameTypeValues"
                "("
                    "const Field<Type>&, "
                    "const vectorField&, "
                    "const scalarField&"
                ") const"
            )
                << "Operation " << operationTypeNames_[operation_]
                << " not available for values of type "
                << pTraits<Type>::typeName
                << exit(FatalError);

            result = pTraits<Type>::zero;
            break;
        }
        case opSumDirectionBalance:
        {
            FatalErrorIn
            (
                "template<class Type>"
                "Type Foam::fieldValues::faceSource::processSameTypeValues"
                "("
                    "const Field<Type>&, "
                    "const vectorField&, "
                    "const scalarField&"
                ") const"
            )
                << "Operation " << operationTypeNames_[operation_]
                << " not available for values of type "
                << pTraits<Type>::typeName
                << exit(FatalError);

            result = pTraits<Type>::zero;
            break;
        }
        case opAverage:
        {
            label n = returnReduce(values.size(), sumOp<label>());
            result = gSum(values)/(scalar(n) + ROOTVSMALL);
            break;
        }
        case opWeightedAverage:
        {
            label wSize = returnReduce(weightField.size(), sumOp<label>());

            if (wSize > 0)
            {
                result =
                    gSum(weightField*values)/(gSum(weightField) + ROOTVSMALL);
            }
            else
            {
                label n = returnReduce(values.size(), sumOp<label>());
                result = gSum(values)/(scalar(n) + ROOTVSMALL);
            }
            break;
        }
        case opAreaAverage:
        {
            const scalarField magSf(mag(Sf));

            result = gSum(magSf*values)/gSum(magSf);
            break;
        }
        case opWeightedAreaAverage:
        {
            const scalarField magSf(mag(Sf));
            label wSize = returnReduce(weightField.size(), sumOp<label>());

            if (wSize > 0)
            {
                result = gSum(weightField*magSf*values)/gSum(magSf*weightField);
            }
            else
            {
                result = gSum(magSf*values)/gSum(magSf);
            }
            break;
        }
        case opAreaIntegrate:
        {
            const scalarField magSf(mag(Sf));

            result = gSum(magSf*values);
            break;
        }
        case opMin:
        {
            result = gMin(values);
            break;
        }
        case opMax:
        {
            result = gMax(values);
            break;
        }
        case opCoV:
        {
            const scalarField magSf(mag(Sf));

            const scalar gSumMagSf = gSum(magSf);

            Type meanValue = gSum(values*magSf)/gSumMagSf;

            const label nComp = pTraits<Type>::nComponents;

            for (direction d=0; d<nComp; ++d)
            {
                scalarField vals(values.component(d));
                scalar mean = component(meanValue, d);
                scalar& res = setComponent(result, d);

                res =
                    sqrt(gSum(magSf*sqr(vals - mean))/gSumMagSf)
                   /(mean + ROOTVSMALL);
            }

            break;
        }
        default:
        {
            // Do nothing
        }
    }

    return result;
}


template<class Type>
Type Foam::fieldValues::faceSource::processValues
(
    const Field<Type>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    return processSameTypeValues(values, Sf, weightField);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fieldValues::faceSource::writeValues
(
    const word& fieldName,
    const scalarField& weightField,
    const bool orient
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values(setFieldValues<Type>(fieldName, true, orient));

        vectorField Sf;
        if (surfacePtr_.valid())
        {
            // Get oriented Sf
            Sf = surfacePtr_().Sf();
        }
        else
        {
            // Get oriented Sf
            Sf = filterField(mesh().Sf(), true);
        }

        // Write raw values on surface if specified
        if (surfaceWriterPtr_.valid())
        {
            Field<Type> allValues(values);
            combineFields(allValues);

            faceList faces;
            pointField points;

            if (surfacePtr_.valid())
            {
                combineSurfaceGeometry(faces, points);
            }
            else
            {
                combineMeshGeometry(faces, points);
            }

            if (Pstream::master())
            {
                fileName outputDir =
                    baseFileDir()/name_/"surface"/obr_.time().timeName();

                surfaceWriterPtr_->write
                (
                    outputDir,
                    word(sourceTypeNames_[source_]) + "_" + sourceName_,
                    points,
                    faces,
                    fieldName,
                    allValues,
                    false
                );
            }
        }

        // Apply scale factor
        values *= scaleFactor_;

        Type result = processValues(values, Sf, weightField);

        file()<< tab << result;

        if (log_)
        {
            Info<< "    " << operationTypeNames_[operation_]
                << "(" << sourceName_ << ") for " << fieldName
                <<  " = " << result << endl;
        }

        // Write state/results information
        const word& opName = operationTypeNames_[operation_];
        word resultName = opName + '(' + sourceName_ + ',' + fieldName + ')';
        this->setResult(resultName, result);
    }

    return ok;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::faceSource::filterField
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const bool applyOrientation
) const
{
    tmp<Field<Type> > tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues();

    forAll(values, i)
    {
        label faceI = faceId_[i];
        label patchI = facePatchId_[i];
        if (patchI >= 0)
        {
            values[i] = field.boundaryField()[patchI][faceI];
        }
        else
        {
            FatalErrorIn
            (
                "fieldValues::faceSource::filterField"
                "("
                    "const GeometricField<Type, fvPatchField, volMesh>&, "
                    "const bool"
                ") const"
            )   << type() << " " << name_ << ": "
                << sourceTypeNames_[source_] << "(" << sourceName_ << "):"
                << nl
                << "    Unable to process internal faces for volume field "
                << field.name() << nl << abort(FatalError);
        }
    }

    if (applyOrientation)
    {
        forAll(values, i)
        {
            values[i] *= faceSign_[i];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValues::faceSource::filterField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
    const bool applyOrientation
) const
{
    tmp<Field<Type> > tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues();

    forAll(values, i)
    {
        label faceI = faceId_[i];
        label patchI = facePatchId_[i];
        if (patchI >= 0)
        {
            values[i] = field.boundaryField()[patchI][faceI];
        }
        else
        {
            values[i] = field[faceI];
        }
    }

    if (applyOrientation)
    {
        forAll(values, i)
        {
            values[i] *= faceSign_[i];
        }
    }

    return tvalues;
}


// ************************************************************************* //
