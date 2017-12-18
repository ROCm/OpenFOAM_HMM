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

#include "fileFieldSelection.H"
#include "objectRegistry.H"
#include "volMesh.H"
#include "fvPatchField.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "UniformDimensionedField.H"

void Foam::functionObjects::fileFieldSelection::addInternalFieldTypes
(
    wordHashSet& set
) const
{
    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);

    const IOobjectList allObjects(mesh, mesh.time().timeName());

    addFromFile<DimensionedField<scalar, volMesh>>(allObjects, set);
    addFromFile<DimensionedField<vector, volMesh>>(allObjects, set);
    addFromFile<DimensionedField<sphericalTensor, volMesh>>(allObjects, set);
    addFromFile<DimensionedField<symmTensor, volMesh>>(allObjects, set);
    addFromFile<DimensionedField<tensor, volMesh>>(allObjects, set);
}


void Foam::functionObjects::fileFieldSelection::addUniformFieldTypes
(
    wordHashSet& set
) const
{
    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);

    const IOobjectList allObjects(mesh, mesh.time().timeName());

    addFromFile<UniformDimensionedField<scalar>>(allObjects, set);
    addFromFile<UniformDimensionedField<vector>>(allObjects, set);
    addFromFile<UniformDimensionedField<sphericalTensor>>(allObjects, set);
    addFromFile<UniformDimensionedField<symmTensor>>(allObjects, set);
    addFromFile<UniformDimensionedField<tensor>>(allObjects, set);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fileFieldSelection::fileFieldSelection
(
    const objectRegistry& obr
)
:
    fieldSelection(obr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fileFieldSelection::~fileFieldSelection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fileFieldSelection::updateSelection()
{
    wordHashSet oldSet;
    oldSet.swap(selection_);

    // Geometric fields
    addGeoFieldTypes<fvPatchField, volMesh>(selection_);
    addGeoFieldTypes<fvsPatchField, surfaceMesh>(selection_);
    addGeoFieldTypes<pointPatchField, pointMesh>(selection_);

    // Internal fields
    addInternalFieldTypes(selection_);

    // Uniform fields
    addUniformFieldTypes(selection_);

    return selection_ != oldSet;
}


// ************************************************************************* //
