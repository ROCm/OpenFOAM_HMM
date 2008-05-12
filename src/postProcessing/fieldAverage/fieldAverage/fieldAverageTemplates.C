/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

#include "fieldAverageItem.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fieldAverage::addMeanFields
(
    const label fieldi,
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& fieldList
)
{
    if (faItems_[fieldi].mean())
    {
        typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const word& fieldName = faItems_[fieldi].fieldName();
        const fieldType& baseField = mesh.lookupObject<fieldType>(fieldName);
        const word meanFieldName = fieldName + EXT_MEAN;

        Info<< "Reading/calculating field " << meanFieldName << nl << endl;
        fieldList.set
        (
            fieldi,
            new fieldType
            (
                IOobject
                (
                    meanFieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                baseField
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::addPrime2MeanFields
(
    const label fieldi,
    PtrList<GeometricField<Type2, fvPatchField, volMesh> >& fieldList
)
{
    if (faItems_[fieldi].mean())
    {
        typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
        typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const word& fieldName = faItems_[fieldi].fieldName();
        const fieldType1& baseField = mesh.lookupObject<fieldType1>(fieldName);
        const fieldType1& meanField =
            mesh.lookupObject<fieldType1>(fieldName + EXT_MEAN);
        const word meanFieldName = fieldName + EXT_PRIME2MEAN;

        Info<< "Reading/calculating field " << meanFieldName << nl << endl;
        fieldList.set
        (
            fieldi,
            new fieldType2
            (
                IOobject
                (
                    meanFieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                sqr(baseField) - sqr(meanField)
            )
        );
    }
}


template<class Type>
void Foam::fieldAverage::calculateMeanFields
(
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& fieldList
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const scalar dt = obr_.time().deltaT().value();

    forAll(faItems_, i)
    {
        if (fieldList.set(i))
        {
            if (faItems_[i].mean())
            {
                const fvMesh& mesh = refCast<const fvMesh>(obr_);
                const word& fieldName = faItems_[i].fieldName();
                const fieldType& baseField =
                    mesh.lookupObject<fieldType>(fieldName);
                fieldType& meanField = fieldList[i];

                scalar alpha = 0.0;
                scalar beta = 0.0;
                if (faItems_[i].timeBase())
                {
                    alpha = (totalTime_[i] - dt)/totalTime_[i];
                    beta = dt/totalTime_[i];
                }
                else
                {
                    alpha = scalar(totalIter_[i] - 1)/scalar(totalIter_[i]);
                    beta = 1.0/scalar(totalIter_[i]);
                }

                meanField = alpha*meanField + beta*baseField;
            }
        }
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::calculatePrime2MeanFields
(
    PtrList<GeometricField<Type2, fvPatchField, volMesh> >& fieldList
)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

    const scalar dt = obr_.time().deltaT().value();

    forAll(faItems_, i)
    {
        if (fieldList.set(i))
        {
            if (faItems_[i].prime2Mean())
            {
                const fvMesh& mesh = refCast<const fvMesh>(obr_);
                const word& fieldName = faItems_[i].fieldName();
                const fieldType1& baseField =
                    mesh.lookupObject<fieldType1>(fieldName);
                const fieldType1& meanField =
                    mesh.lookupObject<fieldType1>(fieldName + EXT_MEAN);
                fieldType2& prime2MeanField = fieldList[i];

                scalar alpha = 0.0;
                scalar beta = 0.0;
                if (faItems_[i].timeBase())
                {
                    alpha = (totalTime_[i] - dt)/totalTime_[i];
                    beta = dt/totalTime_[i];
                }
                else
                {
                    alpha = scalar(totalIter_[i] - 1)/scalar(totalIter_[i]);
                    beta = 1.0/scalar(totalIter_[i]);
                }

                prime2MeanField =
                    alpha*prime2MeanField
                  + beta*sqr(baseField)
                  - sqr(meanField);
            }
        }
    }
}


template<class Type1, class Type2>
void Foam::fieldAverage::addMeanSqrToPrime2Mean
(
    PtrList<GeometricField<Type2, fvPatchField, volMesh> >& fieldList
)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> fieldType1;
    typedef GeometricField<Type2, fvPatchField, volMesh> fieldType2;

    forAll(faItems_, i)
    {
        if (fieldList.set(i))
        {
            if (faItems_[i].prime2Mean())
            {
                const fvMesh& mesh = refCast<const fvMesh>(obr_);
                const word& fieldName = faItems_[i].fieldName();
                const fieldType1& meanField =
                    mesh.lookupObject<fieldType1>(fieldName + EXT_MEAN);
                fieldType2& prime2MeanField = fieldList[i];

                prime2MeanField += sqr(meanField);
            }
        }
    }
}


template<class Type>
void Foam::fieldAverage::writeFieldList
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& fieldList
) const
{
    forAll(fieldList, i)
    {
        if (fieldList.set(i))
        {
            fieldList[i].write();
        }
    }
}


// ************************************************************************* //
