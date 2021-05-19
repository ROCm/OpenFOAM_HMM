/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "foamVtkInternalMeshWriter.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::internalMeshWriter::writeUniform
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
        reportBadState
        (
            FatalErrorInFunction,
            outputState::CELL_DATA,
            outputState::POINT_DATA
        )   << " for field " << fieldName << nl << endl
            << exit(FatalError);
    }
}


template<class Type>
void Foam::vtk::internalMeshWriter::writeCellData
(
    const word& fieldName,
    const UList<Type>& field
)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else
    {
        reportBadState(FatalErrorInFunction, outputState::CELL_DATA)
            << " for field " << fieldName << nl << endl
            << exit(FatalError);
    }

    const labelList& cellMap = vtuCells_.cellMap();

    this->beginDataArray<Type>(fieldName, numberOfCells_);

    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), field, cellMap);
    }
    else
    {
        vtk::writeList(format(), field, cellMap);
    }

    this->endDataArray();
}


// ************************************************************************* //
