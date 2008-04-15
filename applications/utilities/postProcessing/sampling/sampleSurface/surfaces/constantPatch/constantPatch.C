/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "constantPatch.H"
#include "meshSearch.H"
#include "polyMesh.H"
#include "interpolation.H"
#include "dictionary.H"
#include "polyPatch.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(constantPatch, 0);

addToRunTimeSelectionTable(surface, constantPatch, word);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::constantPatch::makeTriangles()
{
    if (patchIndex_ != -1)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchIndex_];

        // Count triangles
        label nTris = 0;

        const faceList& localFaces = patch.localFaces();

        forAll(localFaces, patchFaceI)
        {
            const face& f = localFaces[patchFaceI];

            nTris += f.nTriangles(patch.localPoints());
        }

        // Triangulation is done using all localPoints
        points_ = patch.localPoints();

        faces_.setSize(nTris);
        patchFaceLabels_.setSize(nTris);

        label triI = 0;
        label oldTriI = 0;

        forAll(localFaces, patchFaceI)
        {
            const face& f = localFaces[patchFaceI];

            f.triangles(patch.localPoints(), triI, faces_);

            for(label i = oldTriI; i < triI; i++)
            {
                patchFaceLabels_[i] = patchFaceI;
            }

            oldTriI = triI;
        }
    }
}


void Foam::constantPatch::copyFaces()
{
    if (patchIndex_ != -1)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchIndex_];

        points_ = patch.localPoints();
        faces_ = patch.localFaces();
        patchFaceLabels_.setSize(faces_.size());
        forAll(patchFaceLabels_, i)
        {
            patchFaceLabels_[i] = i;
        }
    }
}


void Foam::constantPatch::createGeometry()
{
    if (triangulate_)
    {
        makeTriangles();
    }
    else
    {
        copyFaces();
    }

    Pout<< "Created " << name() << " :"
        << "  patch:" << patchName_
        << "  faces:" << faces_.size()
        << "  points:" << points_.size() << endl;
}


template <class Type>
Foam::tmp<Foam::Field<Type> > Foam::constantPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<Type>& cache
) const
{
    // One value per face
    tmp<Field<Type> > tvalues(new Field<Type>(patchFaceLabels_.size()));
    Field<Type>& values = tvalues();

    if (patchIndex_ != -1)
    {
        const GeometricField<Type, fvPatchField, volMesh>& vField =
            *cache[fieldName];

        const Field<Type>& bField = vField.boundaryField()[patchIndex_];

        forAll(patchFaceLabels_, elemI)
        {
            values[elemI] = bField[patchFaceLabels_[elemI]];
        }
    }

    return tvalues;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantPatch::constantPatch
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& name,
    const word& patchName,
    const bool triangulate
)
:
    surface(mesh, searchEngine, name),
    patchName_(patchName),
    patchIndex_(mesh.boundaryMesh().findPatchID(patchName_)),
    triangulate_(triangulate),
    points_(0),
    faces_(0),
    patchFaceLabels_(0)
{
    createGeometry();
}


Foam::constantPatch::constantPatch
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict          
)
:
    surface(mesh, searchEngine, dict),
    patchName_(dict.lookup("patchName")),
    patchIndex_(mesh.boundaryMesh().findPatchID(patchName_)),
    triangulate_(getBool(dict, "triangulate", true)),
    points_(0),
    faces_(0),
    patchFaceLabels_(0)
{
    createGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantPatch::~constantPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantPatch::correct
(
    const bool meshChanged,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes,
    const fieldsCache<scalar>& scalarCache
)
{
    if (meshChanged)
    {
        createGeometry();
    }
}


Foam::tmp<Foam::scalarField> Foam::constantPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<scalar>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<scalar>(fieldName, cache);
}


Foam::tmp<Foam::vectorField> Foam::constantPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<vector>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<vector>(fieldName, cache);
}


Foam::tmp<Foam::sphericalTensorField> Foam::constantPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<sphericalTensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<sphericalTensor>(fieldName, cache);
}


Foam::tmp<Foam::symmTensorField> Foam::constantPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<symmTensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<symmTensor>(fieldName, cache);
}


Foam::tmp<Foam::tensorField> Foam::constantPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<tensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<tensor>(fieldName, cache);
}


// ************************************************************************* //
