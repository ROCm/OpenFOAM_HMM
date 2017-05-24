/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#ifndef foamVtkAdaptors_H
#define foamVtkAdaptors_H

// OpenFOAM includes
#include "labelList.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkIdTypeArray.h"
#include "vtkSmartPointer.h"
#include "vtkUnsignedCharArray.h"
#include "vtkAOSDataArrayTemplate.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Attach a smart pointer, or generate a non-null one.
    template<class T>
    inline vtkSmartPointer<T> nonNullSmartPointer(T* ptr)
    {
        return vtkSmartPointer<T>(ptr ? ptr : T::New());
    }


    //- Helper to wrap vtkUnsignedCharArray as a UList
    inline UList<uint8_t> vtkUList
    (
        vtkUnsignedCharArray* array,
        const label size
    )
    {
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(size);

        UList<uint8_t> list
        (
            array->WritePointer(0, size),
            size
        );

        return list;
    }


    //- Helper to wrap vtkIdTypeArray as a UList
    inline UList<vtkIdType> vtkUList
    (
        vtkIdTypeArray* array,
        const label size
    )
    {
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(size);

        UList<vtkIdType> list
        (
            array->WritePointer(0, size),
            size
        );

        return list;
    }


    //- Special helper to wrap vtkCellArray as a UList
    inline UList<vtkIdType> vtkUList
    (
        vtkCellArray* cells,
        const label nCells,
        const label size
    )
    {
        cells->GetData()->SetNumberOfTuples(size);

        UList<vtkIdType> list
        (
            cells->WritePointer(nCells, size),
            size
        );

        return list;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
