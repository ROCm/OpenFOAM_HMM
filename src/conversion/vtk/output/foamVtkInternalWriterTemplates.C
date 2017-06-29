/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "foamVtkInternalWriter.H"
#include "foamVtkOutput.H"
#include "volPointInterpolation.H"
#include "interpolatePointToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::internalWriter::write
(
    const DimensionedField<Type, volMesh>& field
)
{
    const labelList& cellMap = vtuCells_.cellMap();

    const int nCmpt(pTraits<Type>::nComponents);
    // const uint64_t payLoad(cellMap.size() * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os(), field.name(), nCmpt, cellMap.size());
    }
    else
    {
        format().openDataArray<float, nCmpt>(field.name())
            .closeTag();
    }

    // writeField includes payload size, and flush
    vtk::writeField(format(), field, cellMap);

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type, template<class> class PatchField>
void Foam::vtk::internalWriter::write
(
    const GeometricField<Type, PatchField, volMesh>& field
)
{
    write(field.internalField());
}


template<class Type, template<class> class PatchField>
void Foam::vtk::internalWriter::write
(
    const GeometricField<Type, PatchField, pointMesh>& field
)
{
    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();

    const int nCmpt(pTraits<Type>::nComponents);
    const label nVals(vtuCells_.nFieldPoints());

    // Only needed for non-legacy
    const uint64_t payLoad(nVals * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os(), field.name(), nCmpt, nVals);
    }
    else
    {
        format().openDataArray<float, nCmpt>(field.name())
            .closeTag();
    }

    format().writeSize(payLoad);
    vtk::writeList(format(), field);

    for (const label cellId : addPointCellLabels)
    {
        const Type val = interpolatePointToCell(field, cellId);
        vtk::write(format(), val);
    }

    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const volPointInterpolation& pInterp,
    const DimensionedField<Type, volMesh>& vfield
)
{
    typedef DimensionedField<Type, pointMesh> PointFieldType;

    // Use tmp intermediate. Compiler sometimes weird otherwise.
    tmp<PointFieldType> tfield = pInterp.interpolate(vfield);
    const PointFieldType& pfield = tfield();

    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();

    const int nCmpt(pTraits<Type>::nComponents);
    const label nVals(vtuCells_.nFieldPoints());

    // Only needed for non-legacy
    const uint64_t payLoad(nVals * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os(), vfield.name(), nCmpt, nVals);
    }
    else
    {
        format().openDataArray<float, nCmpt>(vfield.name())
            .closeTag();
    }

    format().writeSize(payLoad);
    vtk::writeList(format(), pfield);
    vtk::writeList(format(), vfield, addPointCellLabels);
    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const volPointInterpolation& pInterp,
    const GeometricField<Type, fvPatchField, volMesh>& vfield
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> PointFieldType;

    // Use tmp intermediate. Compiler sometimes weird otherwise.
    tmp<PointFieldType> tfield = pInterp.interpolate(vfield);
    const PointFieldType& pfield = tfield();

    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();

    const int nCmpt(pTraits<Type>::nComponents);
    const label nVals(vtuCells_.nFieldPoints());

    // Only needed for non-legacy
    const uint64_t payLoad(nVals * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os(), vfield.name(), nCmpt, nVals);
    }
    else
    {
        format().openDataArray<float, nCmpt>(vfield.name())
            .closeTag();
    }

    format().writeSize(payLoad);
    vtk::writeList(format(), pfield);
    vtk::writeList(format(), vfield, addPointCellLabels);
    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const UPtrList<const DimensionedField<Type, volMesh>>& flds
)
{
    for (const auto& field : flds)
    {
        write(field);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::vtk::internalWriter::write
(
    const UPtrList<const GeometricField<Type, PatchField, GeoMesh>>& flds
)
{
    for (const auto& field : flds)
    {
        write(field);
    }
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const volPointInterpolation& pInterp,
    const UPtrList<const DimensionedField<Type, volMesh>>& flds
)
{
    for (const auto& field : flds)
    {
        write(pInterp, field);
    }
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const volPointInterpolation& pInterp,
    const UPtrList<const GeometricField<Type, fvPatchField, volMesh>>& flds
)
{
    for (const auto& field : flds)
    {
        write(pInterp, field);
    }
}


// ************************************************************************* //
