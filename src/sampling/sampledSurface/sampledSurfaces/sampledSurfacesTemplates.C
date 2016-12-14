/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "sampledSurfaces.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ListListOps.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::sampledSurfaces::writeSurface
(
    const Field<Type>& values,
    const label surfI,
    const word& fieldName,
    const fileName& outputDir
)
{
    const sampledSurface& s = operator[](surfI);

    if (Pstream::parRun())
    {
        // Collect values from all processors
        List<Field<Type>> gatheredValues(Pstream::nProcs());
        gatheredValues[Pstream::myProcNo()] = values;
        Pstream::gatherList(gatheredValues);

        fileName sampleFile;
        if (Pstream::master())
        {
            // Combine values into single field
            Field<Type> allValues
            (
                ListListOps::combine<Field<Type>>
                (
                    gatheredValues,
                    accessOp<Field<Type>>()
                )
            );

            // Renumber (point data) to correspond to merged points
            if (mergedList_[surfI].pointsMap().size() == allValues.size())
            {
                inplaceReorder(mergedList_[surfI].pointsMap(), allValues);
                allValues.setSize(mergedList_[surfI].points().size());
            }

            // Write to time directory under outputPath_
            // skip surface without faces (eg, a failed cut-plane)
            if (mergedList_[surfI].size())
            {
                sampleFile = formatter_->write
                (
                    outputDir,
                    s.name(),
                    mergedList_[surfI],
                    fieldName,
                    allValues,
                    s.interpolate()
                );
            }
        }

        Pstream::scatter(sampleFile);
        if (sampleFile.size())
        {
            dictionary propsDict;
            propsDict.add("file", sampleFile);
            setProperty(fieldName, propsDict);
        }
    }
    else
    {
        // Write to time directory under outputPath_
        // skip surface without faces (eg, a failed cut-plane)
        if (s.faces().size())
        {
            fileName fName = formatter_->write
            (
                outputDir,
                s.name(),
                s,
                fieldName,
                values,
                s.interpolate()
            );

            dictionary propsDict;
            propsDict.add("file", fName);
            setProperty(fieldName, propsDict);
        }
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    // interpolator for this field
    autoPtr<interpolation<Type>> interpolatorPtr;

    const word& fieldName = vField.name();
    const fileName outputDir = outputPath_/vField.time().timeName();

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);

        Field<Type> values;

        if (s.interpolate())
        {
            if (interpolatorPtr.empty())
            {
                interpolatorPtr = interpolation<Type>::New
                (
                    interpolationScheme_,
                    vField
                );
            }

            values = s.interpolate(interpolatorPtr());
        }
        else
        {
            values = s.sample(vField);
        }

        writeSurface<Type>(values, surfI, fieldName, outputDir);
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
)
{
    const word& fieldName = sField.name();
    const fileName outputDir = outputPath_/sField.time().timeName();

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);
        Field<Type> values(s.sample(sField));
        writeSurface<Type>(values, surfI, fieldName, outputDir);
    }
}


template<class GeoField>
void Foam::sampledSurfaces::sampleAndWrite(const IOobjectList& objects)
{
    wordList fieldNames;
    if (loadFromFiles_)
    {
        fieldNames = objects.sortedNames(GeoField::typeName, fieldSelection_);
    }
    else
    {
        fieldNames = mesh_.thisDb().sortedNames<GeoField>(fieldSelection_);

        writeOriginalIds();
    }

    forAll(fieldNames, fieldi)
    {
        const word& fieldName = fieldNames[fieldi];

        if ((Pstream::master()) && verbose_)
        {
            Pout<< "sampleAndWrite: " << fieldName << endl;
        }

        if (loadFromFiles_)
        {
            const GeoField fld
            (
                IOobject
                (
                    fieldName,
                    time_.timeName(),
                    mesh_,
                    IOobject::MUST_READ
                ),
                mesh_
            );

            sampleAndWrite(fld);
        }
        else
        {
            sampleAndWrite
            (
                mesh_.thisDb().lookupObject<GeoField>(fieldName)
            );
        }
    }
}


// ************************************************************************* //
