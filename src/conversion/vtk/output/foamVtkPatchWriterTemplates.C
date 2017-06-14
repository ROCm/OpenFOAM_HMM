/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "foamVtkPatchWriter.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField>
void Foam::vtk::patchWriter::write
(
    const GeometricField<Type, PatchField, volMesh>& field
)
{
    const int nCmpt(pTraits<Type>::nComponents);
    const uint64_t payLoad(nFaces_ * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os_, field.name(), nCmpt, nFaces_);
    }
    else
    {
        format().openDataArray<float, nCmpt>(field.name())
            .closeTag();
    }

    format().writeSize(payLoad);

    for (const label patchId : patchIDs_)
    {
        const auto& pfld = field.boundaryField()[patchId];

        if (nearCellValue_)
        {
            vtk::writeList(format(), pfld.patchInternalField()());
        }
        else
        {
            vtk::writeList(format(), pfld);
        }
    }

    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type, template<class> class PatchField>
void Foam::vtk::patchWriter::write
(
    const GeometricField<Type, PatchField, pointMesh>& field
)
{
    const int nCmpt(pTraits<Type>::nComponents);
    const uint64_t payLoad(nPoints_ * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os_, field.name(), nCmpt, nPoints_);
    }
    else
    {
        format().openDataArray<float, nCmpt>(field.name())
            .closeTag();
    }

    format().writeSize(payLoad);

    for (const label patchId : patchIDs_)
    {
        const auto& pfld = field.boundaryField()[patchId];

        vtk::writeList(format(), pfld.patchInternalField()());
    }

    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::patchWriter::write
(
    const PrimitivePatchInterpolation<primitivePatch>& pInter,
    const GeometricField<Type, fvPatchField, volMesh>& field
)
{
    const int nCmpt(pTraits<Type>::nComponents);
    const uint64_t payLoad(nPoints_ * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os_, field.name(), nCmpt, nPoints_);
    }
    else
    {
        format().openDataArray<float, nCmpt>(field.name())
            .closeTag();
    }

    format().writeSize(payLoad);

    for (const label patchId : patchIDs_)
    {
        const auto& pfld = field.boundaryField()[patchId];

        if (nearCellValue_)
        {
            auto tfield =
                pInter.faceToPointInterpolate(pfld.patchInternalField()());

            vtk::writeList(format(), tfield());
        }
        else
        {
            auto tfield = pInter.faceToPointInterpolate(pfld);

            vtk::writeList(format(), tfield());
        }
    }

    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::vtk::patchWriter::write
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
void Foam::vtk::patchWriter::write
(
    const PrimitivePatchInterpolation<primitivePatch>& pInter,
    const UPtrList<const GeometricField<Type, fvPatchField, volMesh>>& flds
)
{
    for (const auto& field : flds)
    {
        write(pInter, field);
    }
}


// ************************************************************************* //
