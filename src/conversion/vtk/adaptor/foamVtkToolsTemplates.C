/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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

InClass
    Foam::vtk::Tools

\*---------------------------------------------------------------------------*/

#ifndef foamVtkToolsTemplates_C
#define foamVtkToolsTemplates_C

// OpenFOAM includes
#include "error.H"

// VTK includes
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class PatchType>
vtkSmartPointer<vtkPoints>
Foam::vtk::Tools::Patch::points(const PatchType& p)
{
    // Local patch points to vtkPoints
    const pointField& pts = p.localPoints();

    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    vtkpoints->SetNumberOfPoints(pts.size());
    vtkIdType pointId = 0;

    for (const point& p : pts)
    {
        vtkpoints->SetPoint(pointId, p.v_);
        ++pointId;
    }

    return vtkpoints;
}


template<class PatchType>
vtkSmartPointer<vtkCellArray>
Foam::vtk::Tools::Patch::faces(const PatchType& p)
{
    // Faces as polygons
    const faceList& fcs = p.localFaces();

    label nAlloc = fcs.size();
    for (const face& f : fcs)
    {
        nAlloc += f.size();
    }

    auto vtkcells = vtkSmartPointer<vtkCellArray>::New();

    UList<vtkIdType> list = asUList(vtkcells, fcs.size(), nAlloc);

    // Cell connectivity for polygons
    // [size, verts..., size, verts... ]
    auto iter = list.begin();
    for (const face& f : fcs)
    {
        *(iter++) = f.size();

        for (const label verti : f)
        {
            *(iter++) = verti;
        }
    }

    return vtkcells;
}


template<class PatchType>
vtkSmartPointer<vtkPolyData>
Foam::vtk::Tools::Patch::mesh(const PatchType& p)
{
    auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

    vtkmesh->SetPoints(points(p));
    vtkmesh->SetPolys(faces(p));

    return vtkmesh;
}


template<class PatchType>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::Tools::Patch::faceNormals(const PatchType& p)
{
    auto array = vtkSmartPointer<vtkFloatArray>::New();

    array->SetNumberOfComponents(3);
    array->SetNumberOfTuples(p.size());

    // Unit normals for patch faces.
    // If not already cached, could be more memory efficient to loop over
    // the individual faces instead.

    const vectorField& norms = p.faceNormals();

    vtkIdType faceId = 0;
    for (const vector& n : norms)
    {
        array->SetTuple(faceId, n.v_);
        ++faceId;
    }

    return array;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Low-Level conversions
//

template<class Type>
Foam::label Foam::vtk::Tools::transcribeFloatData
(
    vtkFloatArray* array,
    const UList<Type>& input,
    const label start
)
{
    const int nComp(pTraits<Type>::nComponents);

    if (array->GetNumberOfComponents() != nComp)
    {
        FatalErrorInFunction
            << "vtk array '" << array->GetName()
            << "' has mismatch in number of components for type '"
            << pTraits<Type>::typeName
            << "' : target array has " << array->GetNumberOfComponents()
            << " components instead of " << nComp
            << nl;
    }

    const vtkIdType maxSize = array->GetNumberOfTuples();
    const vtkIdType endPos = vtkIdType(start) + vtkIdType(input.size());

    if (!maxSize)
    {
        // no-op
        return 0;
    }
    else if (start < 0 || vtkIdType(start) >= maxSize)
    {
        WarningInFunction
            << "vtk array '" << array->GetName()
            << "' copy with out-of-range [0," << long(maxSize) << ")"
            << " starting at " << start
            << nl;

        return 0;
    }
    else if (endPos > maxSize)
    {
        WarningInFunction
            << "vtk array '" << array->GetName()
            << "' copy ends out-of-range (" << long(maxSize) << ")"
            << " using sizing (start,size) = ("
            << start << "," << input.size() << ")"
            << nl;

        return 0;
    }

    float scratch[nComp];
    forAll(input, idx)
    {
        const Type& t = input[idx];
        for (direction d=0; d<nComp; ++d)
        {
            scratch[d] = component(t, d);
        }
        remapTuple<Type>(scratch);

        array->SetTuple(start+idx, scratch);
    }

    return input.size();
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::Tools::zeroField
(
    const word& name,
    const label size
)
{
    auto data = vtkSmartPointer<vtkFloatArray>::New();

    data->SetName(name.c_str());
    data->SetNumberOfComponents(int(pTraits<Type>::nComponents));
    data->SetNumberOfTuples(size);

    data->Fill(0);

    return data;
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::Tools::convertFieldToVTK
(
    const word& name,
    const UList<Type>& fld
)
{
    auto data = vtkSmartPointer<vtkFloatArray>::New();

    data->SetName(name.c_str());
    data->SetNumberOfComponents(int(pTraits<Type>::nComponents));
    data->SetNumberOfTuples(fld.size());

    transcribeFloatData(data, fld);

    return data;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
