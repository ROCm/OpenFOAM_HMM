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

#include "foamVtkSurfaceMeshWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::vtk::surfaceMeshWriter::getFaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld
) const
{
    const polyBoundaryMesh& patches = sfld.mesh().boundaryMesh();

    tmp<Field<Type>> tfld(new Field<Type>(pp_.size()));
    Field<Type>& fld = tfld.ref();

    forAll(pp_.addressing(), i)
    {
        const label facei = pp_.addressing()[i];
        const label patchi = patches.whichPatch(facei);

        if (patchi == -1)
        {
            fld[i] = sfld[facei];
        }
        else
        {
            const label localFacei = facei - patches[patchi].start();
            fld[i] = sfld.boundaryField()[patchi][localFacei];
        }
    }

    return tfld;
}


template<class Type>
void Foam::vtk::surfaceMeshWriter::write
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
)
{
    const int nCmpt(pTraits<Type>::nComponents);
    const uint64_t payLoad(pp_.size() * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os(), field.name(), nCmpt, pp_.size());
    }
    else
    {
        format().openDataArray<float, nCmpt>(field.name()).closeTag();
    }

    format().writeSize(payLoad);
    vtk::writeList(format(), getFaceField(field)());

    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::surfaceMeshWriter::write
(
    const GeometricField<Type, faPatchField, areaMesh>& field
)
{
    const int nCmpt(pTraits<Type>::nComponents);
    const uint64_t payLoad(pp_.size() * nCmpt * sizeof(float));

    if (legacy_)
    {
        legacy::floatField(os(), field.name(), nCmpt, pp_.size());
    }
    else
    {
        format().openDataArray<float, nCmpt>(field.name()).closeTag();
    }

    format().writeSize(payLoad);
    vtk::writeList(format(), field.primitiveField());

    format().flush();

    if (!legacy_)
    {
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::surfaceMeshWriter::write
(
    const UPtrList
    <
        const GeometricField<Type, fvsPatchField, surfaceMesh>
    >& sflds
)
{
    for (const auto& field : sflds)
    {
        write(field);
    }
}


template<class Type>
void Foam::vtk::surfaceMeshWriter::write
(
    const UPtrList
    <
        const GeometricField<Type, faPatchField, areaMesh>
    >& sflds
)
{
    for (const auto& field : sflds)
    {
        write(field);
    }
}


// ************************************************************************* //
