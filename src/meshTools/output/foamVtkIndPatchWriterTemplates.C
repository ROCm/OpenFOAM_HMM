/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::vtk::indirectPatchWriter::writeUniform
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
            << "Bad writer state (" << stateNames[state_]
            << ") for field " << fieldName << nl << endl
            << exit(FatalError);
    }
}


template<class Type>
void Foam::vtk::indirectPatchWriter::write
(
    const word& fieldName,
    const UList<Type>& field
)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") or (" << stateNames[outputState::POINT_DATA]
            << ") for field " << fieldName << nl << endl
            << exit(FatalError);
    }

    static_assert
    (
        (
            std::is_same<label, typename pTraits<Type>::cmptType>::value
         || std::is_floating_point<typename pTraits<Type>::cmptType>::value
        ),
        "Label and Floating-point vector space only"
    );

    const direction nCmpt(pTraits<Type>::nComponents);

    label nValues = field.size();

    // Could check sizes:
    //     nValues == nLocalFaces (CELL_DATA)
    //     nValues == nLocalPoints (POINT_DATA)

    if (parallel_)
    {
        reduce(nValues, sumOp<label>());
    }

    if (format_)
    {
        if (std::is_same<label, typename pTraits<Type>::cmptType>::value)
        {
            if (legacy())
            {
                legacy::intField<nCmpt>(format(), fieldName, nValues);
            }
            else
            {
                const uint64_t payLoad = vtk::sizeofData<label, nCmpt>(nValues);

                format().beginDataArray<label, nCmpt>(fieldName);
                format().writeSize(payLoad);
            }
        }
        else
        {
            if (legacy())
            {
                legacy::floatField<nCmpt>(format(), fieldName, nValues);
            }
            else
            {
                const uint64_t payLoad = vtk::sizeofData<float, nCmpt>(nValues);

                format().beginDataArray<float, nCmpt>(fieldName);
                format().writeSize(payLoad);
            }
        }
    }


    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), field);
    }
    else
    {
        vtk::writeList(format(), field);
    }


    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


// ************************************************************************* //
