/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoFieldType>
bool Foam::functionObjects::STDMD::getSnapshot()
{
    if (!initialised_)
    {
        init();
    }

    // Move previous-time snapshot into previous-time slot in z_
    // Effectively moves the lower half of z_ to its upper half
    std::rotate(z_.begin(), z_.begin() + nSnap_, z_.end());

    // Copy new current-time snapshot into current-time slot in z_
    // Effectively copies the new field elements into the lower half of z_
    const GeoFieldType& Field = lookupObject<GeoFieldType>(fieldName_);
    const label nField = Field.size();

    for (direction dir = 0; dir < nComps_; ++dir)
    {
        z_.subColumn(0, nSnap_ + dir*nField, nField) = Field.component(dir);
    }

    return true;
}


template<class Type>
bool Foam::functionObjects::STDMD::getSnapshotType()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (foundObject<VolFieldType>(fieldName_))
    {
        return getSnapshot<VolFieldType>();
    }
    else if (foundObject<SurfaceFieldType>(fieldName_))
    {
        return getSnapshot<SurfaceFieldType>();
    }

    return false;
}


template<class Type>
bool Foam::functionObjects::STDMD::getComps()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (foundObject<VolFieldType>(fieldName_))
    {
        nComps_ = pTraits<typename VolFieldType::value_type>::nComponents;
        return true;
    }
    else if (foundObject<SurfaceFieldType>(fieldName_))
    {
        nComps_ = pTraits<typename SurfaceFieldType::value_type>::nComponents;
        return true;
    }

    return false;
}


template<class Type>
void Foam::functionObjects::STDMD::filterIndexed
(
    List<Type>& lst,
    const UList<label>& indices
)
{
    // Elems within [a, b]
    List<Type> lstWithin(indices.size());

    // Copy if frequency of elem is within [a, b]
    label j = 0;
    for (const auto& i : indices)
    {
        lstWithin[j] = lst[i];
        ++j;
    }
    lst.transfer(lstWithin);
}


template<class MatrixType>
void Foam::functionObjects::STDMD::filterIndexed
(
    MatrixType& mat,
    const UList<label>& indices
)
{
    // Elems within [a, b]
    MatrixType matWithin(labelPair(mat.m(), indices.size()));

    // Copy if frequency of elem is within [a, b]
    label j = 0;
    for (const auto& i : indices)
    {
        matWithin.subColumn(j) = mat.subColumn(i);
        ++j;
    }
    mat.transfer(matWithin);
}


// ************************************************************************* //
