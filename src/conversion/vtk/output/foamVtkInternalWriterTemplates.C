/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
void Foam::vtk::internalWriter::writeUniform
(
    const word& fieldName,
    const Type& val
)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
        vtk::fileWriter::writeUniform<Type>(fieldName, val, numberOfCells_);
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
        vtk::fileWriter::writeUniform<Type>(fieldName, val, numberOfPoints_);
    }
    else
    {
        WarningInFunction
            << "Ignore bad writer state (" << stateNames[state_]
            << ") for field " << fieldName << nl << endl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField>
void Foam::vtk::internalWriter::write
(
    const GeometricField<Type, PatchField, pointMesh>& field
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::POINT_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    const direction nCmpt(pTraits<Type>::nComponents);

    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();

    if (format_)
    {
        if (legacy())
        {
            legacy::floatField<nCmpt>(format(), field.name(), numberOfPoints_);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float, nCmpt>(numberOfPoints_);

            format().beginDataArray<float, nCmpt>(field.name());
            format().writeSize(payLoad);
        }
    }

    if (parallel_)
    {
        List<Type> addedValues(addPointCellLabels.size());
        label outi = 0;

        for (const label cellId : addPointCellLabels)
        {
            addedValues[outi++] = interpolatePointToCell(field, cellId);
        }

        vtk::writeListsParallel(format_.ref(), field, addedValues);
    }
    else
    {
        vtk::writeList(format(), field);

        for (const label cellId : addPointCellLabels)
        {
            const Type val = interpolatePointToCell(field, cellId);
            vtk::write(format(), val);
        }
    }

    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const DimensionedField<Type, volMesh>& field
)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    const direction nCmpt(pTraits<Type>::nComponents);

    const labelList& cellMap = vtuCells_.cellMap();

    if (format_)
    {
        if (legacy())
        {
            legacy::floatField<nCmpt>(format(), field.name(), numberOfCells_);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float, nCmpt>(numberOfCells_);

            format().beginDataArray<float, nCmpt>(field.name());
            format().writeSize(payLoad);
        }
    }

    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), field, cellMap);
    }
    else
    {
        vtk::writeList(format(), field, cellMap);
    }

    if (format_)
    {
        format().flush();
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


template<class Type>
void Foam::vtk::internalWriter::write
(
    const DimensionedField<Type, volMesh>& vfield,
    const volPointInterpolation& pInterp
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::POINT_DATA]
            << ") for field " << vfield.name() << nl << endl
            << exit(FatalError);
    }

    const direction nCmpt(pTraits<Type>::nComponents);

    typedef DimensionedField<Type, pointMesh> PointFieldType;

    // Use tmp intermediate. Compiler sometimes weird otherwise.
    tmp<PointFieldType> tfield = pInterp.interpolate(vfield);
    const PointFieldType& pfield = tfield();

    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();

    if (format_)
    {
        if (legacy())
        {
            legacy::floatField<nCmpt>(format(), vfield.name(), numberOfPoints_);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float, nCmpt>(numberOfPoints_);

            format().beginDataArray<float, nCmpt>(vfield.name());
            format().writeSize(payLoad);
        }
    }

    if (parallel_)
    {
        vtk::writeListsParallel
        (
            format_.ref(),
            pfield,
            vfield,
            addPointCellLabels
        );
    }
    else
    {
        vtk::writeLists(format(), pfield, vfield, addPointCellLabels);
    }

    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::internalWriter::write
(
    const GeometricField<Type, fvPatchField, volMesh>& vfield,
    const volPointInterpolation& pInterp
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::POINT_DATA]
            << ") for field " << vfield.name() << nl << endl
            << exit(FatalError);
    }

    const direction nCmpt(pTraits<Type>::nComponents);

    typedef GeometricField<Type, pointPatchField, pointMesh> PointFieldType;

    // Use tmp intermediate. Compiler sometimes weird otherwise.
    tmp<PointFieldType> tfield = pInterp.interpolate(vfield);
    const PointFieldType& pfield = tfield();

    const labelList& addPointCellLabels = vtuCells_.addPointCellLabels();

    if (format_)
    {
        if (legacy())
        {
            legacy::floatField<nCmpt>(format(), vfield.name(), numberOfPoints_);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float, nCmpt>(numberOfPoints_);

            format().beginDataArray<float, nCmpt>(vfield.name());
            format().writeSize(payLoad);
        }
    }

    if (parallel_)
    {
        vtk::writeListsParallel
        (
            format_.ref(),
            pfield,
            vfield,
            addPointCellLabels
        );
    }
    else
    {
        vtk::writeLists(format(), pfield, vfield, addPointCellLabels);
    }

    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


// ************************************************************************* //
