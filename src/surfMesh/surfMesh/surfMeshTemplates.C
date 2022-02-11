/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "surfMesh.H"
#include "surfFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMeshType>
void Foam::surfMesh::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    const Field<Type>& values
)
{
    const objectRegistry& fieldDb = *this;

    auto* dimfield =
        fieldDb.getObjectPtr<DimensionedField<Type, GeoMeshType>>(fieldName);

    if (dimfield)
    {
        dimfield->dimensions().reset(dims);  // Dimensions may have changed
        dimfield->field() = values;
    }
    else
    {
        dimfield = new DimensionedField<Type, GeoMeshType>
        (
            IOobject
            (
                fieldName,
                fieldDb,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dims,
            values
        );

        dimfield->store();
    }
}


template<class Type, class GeoMeshType>
void Foam::surfMesh::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    Field<Type>&& values
)
{
    const objectRegistry& fieldDb = *this;

    auto* dimfield =
        fieldDb.getObjectPtr<DimensionedField<Type, GeoMeshType>>(fieldName);

    if (dimfield)
    {
        dimfield->dimensions().reset(dims);  // Dimensions may have changed
        dimfield->field() = std::move(values);
    }
    else
    {
        dimfield = new DimensionedField<Type, GeoMeshType>
        (
            IOobject
            (
                fieldName,
                fieldDb,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dims,
            std::move(values)
        );

        dimfield->store();
    }
}


// ************************************************************************* //
