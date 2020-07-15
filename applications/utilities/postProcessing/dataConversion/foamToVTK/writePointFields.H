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
    Read point fields from disk
    and write with vtk::internalWriter and vtk::patchWriter

SourceFiles
    writePointFields.H

\*---------------------------------------------------------------------------*/

#ifndef writePointFields_H
#define writePointFields_H

#include "readFields.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkPatchWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class GeoField>
bool writePointField
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,

    const tmp<GeoField>& tfield,
    const fvMeshSubsetProxy& proxy
)
{
    if (!tfield.valid())
    {
        return false;
    }

    tmp<GeoField> tproxied;
    if (proxy.useSubMesh())
    {
        tproxied = proxy.interpolate(tfield());
        tfield.clear();
    }
    else
    {
        tproxied = tfield;
    }

    if (!tproxied.valid())
    {
        // Or Error?
        return false;
    }


    const auto& field = tproxied();

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


    tproxied.clear();

    return true;
}


template<class GeoField>
label writePointFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,

    const fvMeshSubsetProxy& proxy,
    const typename GeoField::Mesh& ptMesh,
    const IOobjectList& objects,
    const bool syncPar
)
{
    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        if
        (
            writePointField<GeoField>
            (
                internalWriter,
                patchWriters,
                getField<GeoField>(ptMesh, objects, fieldName, syncPar),
                proxy
            )
        )
        {
            ++count;
        }
    }

    return count;
}


label writeAllPointFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar
)
{
    const pointMesh& ptMesh = pointMesh::New(proxy.baseMesh());

    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writePointFields<FieldType>             \
        (                                       \
            internalWriter,                     \
            patchWriters,                       \
            proxy, ptMesh,                      \
            objects,                            \
            syncPar                             \
        )

    label count = 0;
    count += foamToVtk_WRITE_FIELD(pointScalarField);
    count += foamToVtk_WRITE_FIELD(pointVectorField);
    count += foamToVtk_WRITE_FIELD(pointSphericalTensorField);
    count += foamToVtk_WRITE_FIELD(pointSymmTensorField);
    count += foamToVtk_WRITE_FIELD(pointTensorField);

    #undef foamToVTK_WRITE_FIELD
    return count;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
