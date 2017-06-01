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
void Foam::foamVtkOutput::internalWriter::write
(
    const UPtrList<const DimensionedField<Type, volMesh>>& flds
)
{
    const labelList& cellMap = vtkCells_.cellMap();

    const int nCmpt(pTraits<Type>::nComponents);
    // const uint64_t payLoad(cellMap.size() * nCmpt * sizeof(float));

    forAll(flds, i)
    {
        const auto& fld = flds[i];

        if (legacy_)
        {
            legacy::floatField(os(), fld.name(), nCmpt, cellMap.size());
        }
        else
        {
            format().openDataArray<float, nCmpt>(fld.name())
                .closeTag();
        }

        // writeField includes payload size
        foamVtkOutput::writeField(format(), fld, cellMap);

        if (!legacy_)
        {
            format().endDataArray();
        }
    }
}


template<class Type, template<class> class PatchField>
void Foam::foamVtkOutput::internalWriter::write
(
    const UPtrList<const GeometricField<Type, PatchField, volMesh>>& flds
)
{
    const labelList& cellMap = vtkCells_.cellMap();

    const int nCmpt(pTraits<Type>::nComponents);
    // const uint64_t payLoad(cellMap.size() * nCmpt * sizeof(float));

    forAll(flds, i)
    {
        const auto& fld = flds[i];

        if (legacy_)
        {
            legacy::floatField(os(), fld.name(), nCmpt, cellMap.size());
        }
        else
        {
            format().openDataArray<float, nCmpt>(fld.name())
                .closeTag();
        }

        // writeField includes payload size
        foamVtkOutput::writeField(format(), fld, cellMap);

        if (!legacy_)
        {
            format().endDataArray();
        }
    }
}


template<class Type, template<class> class PatchField>
void Foam::foamVtkOutput::internalWriter::write
(
    const UPtrList<const GeometricField<Type, PatchField, pointMesh>>& flds
)
{
    const labelList& addPointCellLabels = vtkCells_.addPointCellLabels();

    const int nCmpt(pTraits<Type>::nComponents);
    const int nVals(vtkCells_.nFieldPoints());

    // Only needed for non-legacy
    const uint64_t payLoad(nVals * nCmpt * sizeof(float));

    forAll(flds, i)
    {
        const auto& fld = flds[i];

        if (legacy_)
        {
            legacy::floatField(os(), fld.name(), nCmpt, nVals);
        }
        else
        {
            format().openDataArray<float, nCmpt>(fld.name())
                .closeTag();
        }

        format().writeSize(payLoad);
        foamVtkOutput::writeList(format(), fld);

        forAll(addPointCellLabels, i)
        {
            const Type val = interpolatePointToCell(fld, addPointCellLabels[i]);
            foamVtkOutput::write(format(), val);
        }

        format().flush();

        if (!legacy_)
        {
            format().endDataArray();
        }
    }
}


template<class Type>
void Foam::foamVtkOutput::internalWriter::write
(
    const volPointInterpolation& pInterp,
    const UPtrList<const DimensionedField<Type, volMesh>>& flds
)
{
    const labelList& addPointCellLabels = vtkCells_.addPointCellLabels();

    const int nCmpt(pTraits<Type>::nComponents);
    const int nVals(vtkCells_.nFieldPoints());

    // Only needed for non-legacy
    const uint64_t payLoad(nVals * nCmpt * sizeof(float));

    forAll(flds, i)
    {
        const auto& vfield = flds[i];
        const auto& pfield = pInterp.interpolate(vfield)();

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
        foamVtkOutput::writeList(format(), pfield);
        foamVtkOutput::writeList(format(), vfield, addPointCellLabels);
        format().flush();

        if (!legacy_)
        {
            format().endDataArray();
        }
    }
}


template<class Type, template<class> class PatchField>
void Foam::foamVtkOutput::internalWriter::write
(
    const volPointInterpolation& pInterp,
    const UPtrList<const GeometricField<Type, PatchField, volMesh>>& flds
)
{
    const labelList& addPointCellLabels = vtkCells_.addPointCellLabels();

    const int nCmpt(pTraits<Type>::nComponents);
    const int nVals(vtkCells_.nFieldPoints());

    // Only needed for non-legacy
    const uint64_t payLoad(nVals * nCmpt * sizeof(float));

    forAll(flds, i)
    {
        const auto& vfield = flds[i];
        const auto& pfield = pInterp.interpolate(vfield)();

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
        foamVtkOutput::writeList(format(), pfield);
        foamVtkOutput::writeList(format(), vfield, addPointCellLabels);
        format().flush();

        if (!legacy_)
        {
            format().endDataArray();
        }
    }
}


// ************************************************************************* //
