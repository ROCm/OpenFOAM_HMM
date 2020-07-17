/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Read volume fields from disk
    and write with vtk::internalWriter and vtk::patchWriter

SourceFiles
    writeVolFields.H

\*---------------------------------------------------------------------------*/

#ifndef writeVolFields_H
#define writeVolFields_H

#include "readFields.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkPatchWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class GeoField>
bool writeVolField
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,
    const tmp<GeoField>& tfield
)
{
    if (!tfield.valid())
    {
        return false;
    }

    const auto& field = tfield();

    // Internal
    if (internalWriter)
    {
        internalWriter->write(field);
    }

    // Boundary
    for (vtk::patchWriter& writer : patchWriters)
    {
        writer.write(field);
    }

    tfield.clear();
    return true;
}


template<class GeoField>
bool writeVolField
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,

    UPtrList<vtk::patchWriter>& patchWriters,
    const UPtrList<PrimitivePatchInterpolation<primitivePatch>>& patchInterps,

    const tmp<GeoField>& tfield
)
{
    if (!tfield.valid())
    {
        return false;
    }

    const auto& field = tfield();

    // Internal
    if (internalWriter && pInterp)
    {
        internalWriter->write(field, *pInterp);
    }

    // Boundary
    label writeri = 0;
    for (vtk::patchWriter& writer : patchWriters)
    {
        if (writeri < patchInterps.size() && patchInterps.set(writeri))
        {
            writer.write(field, patchInterps[writeri]);
        }
        ++writeri;
    }

    tfield.clear();
    return true;
}


template<class GeoField>
label writeVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,
    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar
)
{
    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        if
        (
            writeVolField<GeoField>
            (
                internalWriter,
                patchWriters,
                getField<GeoField>(proxy, objects, fieldName, syncPar)
            )
        )
        {
            ++count;
        }
    }

    return count;
}


template<class GeoField>
label writeVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,

    UPtrList<vtk::patchWriter>& patchWriters,
    const UPtrList<PrimitivePatchInterpolation<primitivePatch>>& patchInterps,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar
)
{
    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        if
        (
            writeVolField<GeoField>
            (
                internalWriter, pInterp,
                patchWriters,   patchInterps,
                getField<GeoField>(proxy, objects, fieldName, syncPar)
            )
        )
        {
            ++count;
        }
    }

    return count;
}


label writeAllVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar
)
{
    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writeVolFields<FieldType>               \
        (                                       \
            internalWriter,                     \
            patchWriters,                       \
            proxy,                              \
            objects,                            \
            syncPar                             \
        )

    label count = 0;
    count += foamToVtk_WRITE_FIELD(volScalarField);
    count += foamToVtk_WRITE_FIELD(volVectorField);
    count += foamToVtk_WRITE_FIELD(volSphericalTensorField);
    count += foamToVtk_WRITE_FIELD(volSymmTensorField);
    count += foamToVtk_WRITE_FIELD(volTensorField);

    #undef foamToVTK_WRITE_FIELD
    return count;
}


label writeAllVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,

    UPtrList<vtk::patchWriter>& patchWriters,
    const UPtrList<PrimitivePatchInterpolation<primitivePatch>>& patchInterps,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar
)
{
    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writeVolFields<FieldType>               \
        (                                       \
            internalWriter, pInterp,            \
            patchWriters,   patchInterps,       \
            proxy,                              \
            objects,                            \
            syncPar                             \
        )


    label count = 0;
    count += foamToVtk_WRITE_FIELD(volScalarField);
    count += foamToVtk_WRITE_FIELD(volVectorField);
    count += foamToVtk_WRITE_FIELD(volSphericalTensorField);
    count += foamToVtk_WRITE_FIELD(volSymmTensorField);
    count += foamToVtk_WRITE_FIELD(volTensorField);

    #undef foamToVtk_WRITE_FIELD
    return count;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
