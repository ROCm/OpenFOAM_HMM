/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "surfaceFieldValue.H"
#include "surfaceFields.H"
#include "surfFields.H"
#include "volFields.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class WeightType>
inline bool Foam::functionObjects::fieldValues::surfaceFieldValue::canWeight
(
    const Field<WeightType>& weightField
) const
{
    return
    (
        usesWeight()
     && returnReduce(!weightField.empty(), orOp<bool>()) // On some processor
    );
}


template<class Type>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::validField
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;
    typedef DimensionedField<Type, surfGeoMesh> smt;

    return
    (
        foundObject<smt>(fieldName)
     || foundObject<vf>(fieldName)
     || (regionType_ != stSampledSurface && foundObject<sf>(fieldName))
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::getFieldValues
(
    const word& fieldName,
    const bool mustGet
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;
    typedef DimensionedField<Type, surfGeoMesh> smt;

    if (foundObject<smt>(fieldName))
    {
        return lookupObject<smt>(fieldName);
    }
    else if (regionType_ != stSampledSurface && foundObject<sf>(fieldName))
    {
        return filterField(lookupObject<sf>(fieldName));
    }
    else if (foundObject<vf>(fieldName))
    {
        const vf& fld = lookupObject<vf>(fieldName);

        if (surfacePtr_.valid())
        {
            if (surfacePtr_().interpolate())
            {
                const interpolationCellPoint<Type> interp(fld);
                tmp<Field<Type>> tintFld(surfacePtr_().interpolate(interp));
                const Field<Type>& intFld = tintFld();

                // Average
                const faceList& faces = surfacePtr_().faces();
                tmp<Field<Type>> tavg
                (
                    new Field<Type>(faces.size(), Zero)
                );
                Field<Type>& avg = tavg.ref();

                forAll(faces, facei)
                {
                    const face& f = faces[facei];
                    forAll(f, fp)
                    {
                        avg[facei] += intFld[f[fp]];
                    }
                    avg[facei] /= f.size();
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
            return filterField(fld);
        }
    }

    if (mustGet)
    {
        FatalErrorInFunction
            << "Field " << fieldName << " not found in database"
            << abort(FatalError);
    }

    return tmp<Field<Type>>(new Field<Type>(0));
}


template<class Type, class WeightType>
Type Foam::functionObjects::fieldValues::surfaceFieldValue::
processSameTypeValues
(
    const Field<Type>& values,
    const vectorField& Sf,
    const Field<WeightType>& weightField
) const
{
    Type result = Zero;
    switch (operation_)
    {
        case opNone:
        {
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
        case opSumMag:
        {
            result = gSum(cmptMag(values));
            break;
        }
        case opSum:
        case opWeightedSum:
        case opAbsWeightedSum:
        {
            if (canWeight(weightField))
            {
                tmp<scalarField> weight(weightingFactor(weightField));

                result = gSum(weight*values);
            }
            else
            {
                // Unweighted form
                result = gSum(values);
            }
            break;
        }
        case opSumDirection:
        case opSumDirectionBalance:
        {
            FatalErrorInFunction
                << "Operation " << operationTypeNames_[operation_]
                << " not available for values of type "
                << pTraits<Type>::typeName
                << exit(FatalError);

            break;
        }
        case opAverage:
        case opWeightedAverage:
        case opAbsWeightedAverage:
        {
            if (canWeight(weightField))
            {
                const scalarField factor(weightingFactor(weightField));

                result = gSum(factor*values)/(gSum(factor) + ROOTVSMALL);
            }
            else
            {
                // Unweighted form
                const label n = returnReduce(values.size(), sumOp<label>());
                result = gSum(values)/(scalar(n) + ROOTVSMALL);
            }
            break;
        }
        case opAreaAverage:
        case opWeightedAreaAverage:
        case opAbsWeightedAreaAverage:
        {
            if (canWeight(weightField))
            {
                const scalarField factor(weightingFactor(weightField, Sf));

                result = gSum(factor*values)/gSum(factor + ROOTVSMALL);
            }
            else
            {
                // Unweighted form
                const scalarField factor(mag(Sf));
                result = gSum(factor*values)/gSum(factor);
            }
            break;
        }
        case opAreaIntegrate:
        case opWeightedAreaIntegrate:
        case opAbsWeightedAreaIntegrate:
        {
            if (canWeight(weightField))
            {
                tmp<scalarField> factor(weightingFactor(weightField, Sf));
                result = gSum(factor*values);
            }
            else
            {
                // Unweighted form
                tmp<scalarField> factor(mag(Sf));
                result = gSum(factor*values);
            }
            break;
        }
        case opCoV:
        {
            const scalarField magSf(mag(Sf));
            const scalar gSumMagSf = gSum(magSf);

            Type meanValue = gSum(values*magSf)/gSumMagSf;

            for (direction d=0; d < pTraits<Type>::nComponents; ++d)
            {
                tmp<scalarField> vals(values.component(d));
                const scalar mean = component(meanValue, d);
                scalar& res = setComponent(result, d);

                res =
                    sqrt(gSum(magSf*sqr(vals - mean))/gSumMagSf)
                   /(mean + ROOTVSMALL);
            }

            break;
        }

        case opAreaNormalAverage:
        case opAreaNormalIntegrate:
        {
            // Handled in specializations only
            break;
        }
    }

    return result;
}


template<class Type, class WeightType>
Type Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<Type>& values,
    const vectorField& Sf,
    const Field<WeightType>& weightField
) const
{
    return processSameTypeValues(values, Sf, weightField);
}


template<class WeightType>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<WeightType>& weightField
) const
{
    // The scalar form is specialized.
    // For other types always need mag() to generate a scalar field.
    return mag(weightField);
}


template<class WeightType>
Foam::label Foam::functionObjects::fieldValues::surfaceFieldValue::writeAll
(
    const vectorField& Sf,
    const Field<WeightType>& weightField,
    const meshedSurf& surfToWrite
)
{
    label nProcessed = 0;

    for (const word& fieldName : fields_)
    {
        if
        (
            writeValues<scalar>(fieldName, Sf, weightField, surfToWrite)
         || writeValues<vector>(fieldName, Sf, weightField, surfToWrite)
         || writeValues<sphericalTensor>
            (
                fieldName, Sf, weightField, surfToWrite
            )
         || writeValues<symmTensor>(fieldName, Sf, weightField, surfToWrite)
         || writeValues<tensor>(fieldName, Sf, weightField, surfToWrite)
        )
        {
            ++nProcessed;
        }
        else
        {
            WarningInFunction
                << "Requested field " << fieldName
                << " not found in database and not processed"
                << endl;
        }
    }

    return nProcessed;
}


template<class Type, class WeightType>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::writeValues
(
    const word& fieldName,
    const vectorField& Sf,
    const Field<WeightType>& weightField,
    const meshedSurf& surfToWrite
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values(getFieldValues<Type>(fieldName, true));

        // Write raw values on surface if specified
        if (surfaceWriterPtr_.valid())
        {
            Field<Type> allValues(values);
            combineFields(allValues);

            if (Pstream::master())
            {
                surfaceWriterPtr_->write
                (
                    outputDir(),
                    regionTypeNames_[regionType_] + ("_" + regionName_),
                    surfToWrite,
                    fieldName,
                    allValues,
                    false
                );
            }
        }

        if (operation_ != opNone)
        {
            // Apply scale factor
            values *= scaleFactor_;

            Type result = processValues(values, Sf, weightField);

            switch (postOperation_)
            {
                case postOpNone:
                {
                    break;
                }
                case postOpSqrt:
                {
                    // sqrt: component-wise - doesn't change the type
                    for (direction d=0; d < pTraits<Type>::nComponents; ++d)
                    {
                        setComponent(result,  d)
                            = sqrt(mag(component(result,  d)));
                    }
                    break;
                }
            }

            file()<< tab << result;

            // Write state/results information
            word prefix, suffix;
            {
                if (postOperation_ != postOpNone)
                {
                    // Adjust result name to include post-operation
                    prefix += postOperationTypeNames_[postOperation_];
                    prefix += '(';
                    suffix += ')';
                }

                prefix += operationTypeNames_[operation_];
                prefix += '(';
                suffix += ')';
            }

            Log << "    " << prefix << regionName_ << suffix
                << " of " << fieldName
                <<  " = " << result << endl;

            // Write state/results information
            word resultName = prefix + regionName_ + ',' + fieldName + suffix;
            this->setResult(resultName, result);
        }
    }

    return ok;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::filterField
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    tmp<Field<Type>> tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(values, i)
    {
        const label facei = faceId_[i];
        const label patchi = facePatchId_[i];
        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << regionTypeNames_[regionType_] << "(" << regionName_ << "):"
                << nl
                << "    Unable to process internal faces for volume field "
                << field.name() << nl << abort(FatalError);
        }
    }

    // No need to flip values - all boundary faces point outwards

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::filterField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
) const
{
    tmp<Field<Type>> tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(values, i)
    {
        const label facei = faceId_[i];
        const label patchi = facePatchId_[i];
        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            values[i] = field[facei];
        }
    }

    if (debug)
    {
        Pout<< "field " << field.name() << " oriented: "
            << field.oriented()() << endl;
    }

    if (field.oriented()())
    {
        forAll(values, i)
        {
            if (faceFlip_[i])
            {
                values[i] *= -1;
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
