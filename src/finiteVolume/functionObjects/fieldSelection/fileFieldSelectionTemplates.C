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

#include "IOobjectList.H"
#include "GeometricField.H"
#include "fvMesh.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fileFieldSelection::addFromFile
(
    const IOobjectList& allFileObjects,
    wordHashSet& set
) const
{
    DynamicList<word> names(this->size());

    for (const wordRe& fieldName : *this)
    {
        names.append(allFileObjects.names(Type::typeName, fieldName));
    }

    set.insert(names);
}


template<template<class> class PatchType, class MeshType>
void Foam::functionObjects::fileFieldSelection::addGeoFieldTypes
(
    wordHashSet& set
) const
{
    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);

    const IOobjectList allObjects(mesh, mesh.time().timeName());

    addFromFile<GeometricField<scalar, PatchType, MeshType>>(allObjects, set);
    addFromFile<GeometricField<vector, PatchType, MeshType>>(allObjects, set);
    addFromFile<GeometricField<sphericalTensor, PatchType, MeshType>>
    (
        allObjects,
        set
    );
    addFromFile<GeometricField<symmTensor, PatchType, MeshType>>
    (
        allObjects,
        set
    );
    addFromFile<GeometricField<tensor, PatchType, MeshType>>(allObjects, set);
}


// ************************************************************************* //
