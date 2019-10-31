/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "fieldAverageItem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "polySurfaceFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldAverage::addMeanFieldType
(
    fieldAverageItem& item
)
{
    const word& fieldName = item.fieldName();

    if (!foundObject<Type>(fieldName))
    {
        return;
    }

    // Field has been found, so set active flag to true
    item.active() = true;

    const word& meanFieldName = item.meanFieldName();

    Log << "    Reading/initialising field " << meanFieldName << endl;

    if (foundObject<Type>(meanFieldName))
    {}
    else if (obr().found(meanFieldName))
    {
        Log << "    Cannot allocate average field " << meanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        item.mean() = false;
    }
    else
    {
        const Type& baseField = lookupObject<Type>(fieldName);

        // Store on registry
        obr().store
        (
            new Type
            (
                IOobject
                (
                    meanFieldName,
                    obr().time().timeName(obr().time().startTime().value()),
                    obr(),
                    (
                        restartOnOutput_
                      ? IOobject::NO_READ
                      : IOobject::READ_IF_PRESENT
                    ),
                    IOobject::NO_WRITE
                ),
                1*baseField
            )
        );
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::addMeanField
(
    fieldAverageItem& item
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    if (item.mean())
    {
        addMeanFieldType<VolFieldType>(item);
        addMeanFieldType<SurfaceFieldType>(item);
        addMeanFieldType<SurfFieldType>(item);
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::restoreWindowFieldsType
(
    const fieldAverageItem& item
)
{
    if (restartOnOutput_)
    {
        return;
    }

    const word& fieldName = item.fieldName();

    const Type* fieldPtr = findObject<Type>(fieldName);

    if (!fieldPtr)
    {
        return;
    }

    const FIFOStack<word>& fieldNames = item.windowFieldNames();

    forAllConstIters(fieldNames, fieldIter)
    {
        const word& name = fieldIter();

        IOobject io
        (
            name,
            obr().time().timeName(obr().time().startTime().value()),
            obr(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (io.typeHeaderOk<Type>(true))
        {
            DebugInfo << "Read and store: " << name << endl;
            obr().store(new Type(io, fieldPtr->mesh()));
        }
        else
        {
            WarningInFunction
                << "Unable to read window " << Type::typeName << " " << name
                << ".  Averaging restart behaviour may be compromised"
                << endl;
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::restoreWindowFields
(
    const fieldAverageItem& item
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    if (item.window() > 0)
    {
        restoreWindowFieldsType<VolFieldType>(item);
        restoreWindowFieldsType<SurfaceFieldType>(item);
        restoreWindowFieldsType<SurfFieldType>(item);
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addPrime2MeanFieldType
(
    fieldAverageItem& item
)
{
    const word& fieldName = item.fieldName();

    if (!foundObject<Type1>(fieldName))
    {
        return;
    }

    const word& meanFieldName = item.meanFieldName();
    const word& prime2MeanFieldName = item.prime2MeanFieldName();

    Log << "    Reading/initialising field " << prime2MeanFieldName << nl;

    if (foundObject<Type2>(prime2MeanFieldName))
    {}
    else if (obr().found(prime2MeanFieldName))
    {
        Log << "    Cannot allocate average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        item.prime2Mean() = false;
    }
    else
    {
        const Type1& baseField = lookupObject<Type1>(fieldName);
        const Type1& meanField = lookupObject<Type1>(meanFieldName);

        // Store on registry
        obr().store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr().time().timeName(obr().time().startTime().value()),
                    obr(),
                    restartOnOutput_?
                        IOobject::NO_READ
                      : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                sqr(baseField) - sqr(meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addPrime2MeanField
(
    fieldAverageItem& item
)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, polySurfaceGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, polySurfaceGeoMesh> SurfFieldType2;

    if (item.prime2Mean())
    {
        if (!item.mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << item.fieldName() << nl << exit(FatalError);
        }

        addPrime2MeanFieldType<VolFieldType1, VolFieldType2>(item);
        addPrime2MeanFieldType<SurfaceFieldType1, SurfaceFieldType2>(item);
        addPrime2MeanFieldType<SurfFieldType1, SurfFieldType2>(item);
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::storeWindowFieldType
(
    fieldAverageItem& item
)
{
    const word& fieldName = item.fieldName();
    if (!foundObject<Type>(fieldName))
    {
        return;
    }

    const Type& baseField = lookupObject<Type>(fieldName);

    const word windowFieldName = item.windowFieldName(this->name());

    // Store on registry
    obr().store
    (
        new Type
        (
            IOobject
            (
                windowFieldName,
                obr().time().timeName(obr().time().startTime().value()),
                obr(),
                restartOnOutput_ ?
                    IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            1*baseField
        )
    );

    DebugInfo << "Create and store: " << windowFieldName << endl;

    item.addToWindow(windowFieldName, obr().time().deltaTValue());
}


template<class Type>
void Foam::functionObjects::fieldAverage::storeWindowFields()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    for (fieldAverageItem& item : faItems_)
    {
        if (item.storeWindowFields())
        {
            storeWindowFieldType<VolFieldType>(item);
            storeWindowFieldType<SurfaceFieldType>(item);
            storeWindowFieldType<SurfFieldType>(item);
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    for (const fieldAverageItem& item : faItems_)
    {
        item.calculateMeanField<VolFieldType>(obr());
        item.calculateMeanField<SurfaceFieldType>(obr());
        item.calculateMeanField<SurfFieldType>(obr());
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculatePrime2MeanFields() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, polySurfaceGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, polySurfaceGeoMesh> SurfFieldType2;

    for (const fieldAverageItem& item : faItems_)
    {
        item.calculatePrime2MeanField<VolFieldType1, VolFieldType2>(obr());
        item.calculatePrime2MeanField<SurfaceFieldType1, SurfaceFieldType2>
        (
            obr()
        );
        item.calculatePrime2MeanField<SurfFieldType1, SurfFieldType2>(obr());
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2MeanType
(
    const fieldAverageItem& item
) const
{
    const word& fieldName = item.fieldName();

    if (!foundObject<Type1>(fieldName))
    {
        return;
    }

    const Type1& meanField = lookupObject<Type1>(item.meanFieldName());

    Type2& prime2MeanField = lookupObjectRef<Type2>(item.prime2MeanFieldName());

    prime2MeanField += sqr(meanField);
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2Mean() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, polySurfaceGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, polySurfaceGeoMesh> SurfFieldType2;

    for (const fieldAverageItem& item : faItems_)
    {
        if (item.prime2Mean())
        {
            addMeanSqrToPrime2MeanType<VolFieldType1, VolFieldType2>(item);
            addMeanSqrToPrime2MeanType<SurfaceFieldType1, SurfaceFieldType2>
            (
                item
            );
            addMeanSqrToPrime2MeanType<SurfFieldType1, SurfFieldType2>(item);
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFieldType
(
    const word& fieldName
) const
{
    if (foundObject<Type>(fieldName))
    {
        const Type& f = lookupObject<Type>(fieldName);
        f.write();
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    for (const fieldAverageItem& item : faItems_)
    {
        if (item.mean())
        {
            const word& fieldName = item.meanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
            writeFieldType<SurfFieldType>(fieldName);
        }

        if (item.prime2Mean())
        {
            const word& fieldName = item.prime2MeanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
            writeFieldType<SurfFieldType>(fieldName);
        }

        if (item.writeWindowFields())
        {
            FIFOStack<word> fieldNames = item.windowFieldNames();
            forAllConstIters(fieldNames, fieldNameIter)
            {
                const word& fieldName = fieldNameIter();
                writeFieldType<VolFieldType>(fieldName);
                writeFieldType<SurfaceFieldType>(fieldName);
                writeFieldType<SurfFieldType>(fieldName);
            }
        }
    }
}


// ************************************************************************* //
