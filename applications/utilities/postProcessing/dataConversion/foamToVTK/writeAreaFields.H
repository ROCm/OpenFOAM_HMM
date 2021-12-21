/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Read finite-area fields from disk
    and write with vtk::uindirectPatchGeoFieldsWriter

SourceFiles
    writeAreaFields.H

\*---------------------------------------------------------------------------*/

#ifndef writeAreaFields_H
#define writeAreaFields_H

#include "readFields.H"
#include "foamVtkUIndPatchGeoFieldsWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Writer type for finite-area mesh + fields
typedef vtk::uindirectPatchGeoFieldsWriter vtkWriterType_areaMesh;

template<class GeoField>
bool writeAreaField
(
    vtkWriterType_areaMesh& writer,
    const tmp<GeoField>& tfield
)
{
    if (tfield.valid())
    {
        writer.write(tfield());
        tfield.clear();

        return true;
    }

    return false;
}


template<class GeoField>
label writeAreaFields
(
    vtkWriterType_areaMesh& writer,
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const bool syncPar
)
{
    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        if
        (
            writeAreaField<GeoField>
            (
                writer,
                getField<GeoField>(mesh, objects, fieldName, syncPar)
            )
        )
        {
            ++count;
        }
    }

    return count;
}


label writeAllAreaFields
(
    vtkWriterType_areaMesh& writer,
    const faMesh& mesh,
    const IOobjectList& objects,
    const bool syncPar
)
{
    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writeAreaFields<FieldType>              \
        (                                       \
            writer,                             \
            mesh,                               \
            objects,                            \
            syncPar                             \
        )

    label count = 0;
    count += foamToVtk_WRITE_FIELD(areaScalarField);
    count += foamToVtk_WRITE_FIELD(areaVectorField);
    count += foamToVtk_WRITE_FIELD(areaSphericalTensorField);
    count += foamToVtk_WRITE_FIELD(areaSymmTensorField);
    count += foamToVtk_WRITE_FIELD(areaTensorField);

    #undef foamToVtk_WRITE_FIELD
    return count;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
