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
    Read dimensioned fields from disk
    and write with vtk::internalWriter

SourceFiles
    writeDimFields.H

\*---------------------------------------------------------------------------*/

#ifndef writeDimFields_H
#define writeDimFields_H

#include "readFields.H"
#include "foamVtkInternalWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class GeoField>
bool writeDimField
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const tmp<GeoField>& tfield
)
{
    if (!tfield.valid())
    {
        return false;
    }

    const auto& field = tfield();

    if (internalWriter)
    {
        internalWriter->write(field);
    }

    tfield.clear();
    return true;
}


template<class GeoField>
bool writeDimField
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,
    const tmp<GeoField>& tfield
)
{
    if (!tfield.valid())
    {
        return false;
    }

    const auto& field = tfield();

    if (internalWriter && pInterp)
    {
        internalWriter->write(field, *pInterp);
    }

    tfield.clear();
    return true;
}


template<class GeoField>
label writeDimFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
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
            writeDimField<GeoField>
            (
                internalWriter,
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
label writeDimFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,

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
            writeDimField<GeoField>
            (
                internalWriter, pInterp,
                getField<GeoField>(proxy, objects, fieldName, syncPar)
            )
        )
        {
            ++count;
        }
    }

    return count;
}


label writeAllDimFields
(
    autoPtr<vtk::internalWriter>& internalWriter,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar
)
{
    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writeDimFields<FieldType>               \
        (                                       \
            internalWriter,                     \
            proxy,                              \
            objects,                            \
            syncPar                             \
        )

    label count = 0;
    count += foamToVtk_WRITE_FIELD(volScalarField::Internal);
    count += foamToVtk_WRITE_FIELD(volVectorField::Internal);
    count += foamToVtk_WRITE_FIELD(volSphericalTensorField::Internal);
    count += foamToVtk_WRITE_FIELD(volSymmTensorField::Internal);
    count += foamToVtk_WRITE_FIELD(volTensorField::Internal);

    #undef foamToVTK_WRITE_FIELD
    return count;
}


label writeAllDimFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar
)
{
    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writeDimFields<FieldType>               \
        (                                       \
            internalWriter, pInterp,            \
            proxy,                              \
            objects,                            \
            syncPar                             \
        )

    label count = 0;
    count += foamToVtk_WRITE_FIELD(volScalarField::Internal);
    count += foamToVtk_WRITE_FIELD(volVectorField::Internal);
    count += foamToVtk_WRITE_FIELD(volSphericalTensorField::Internal);
    count += foamToVtk_WRITE_FIELD(volSymmTensorField::Internal);
    count += foamToVtk_WRITE_FIELD(volTensorField::Internal);

    #undef foamToVTK_WRITE_FIELD
    return count;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
