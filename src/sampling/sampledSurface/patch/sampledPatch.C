/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "sampledPatch.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "volFields.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledPatch, 0);
    addNamedToRunTimeSelectionTable(sampledSurface, sampledPatch, word, patch);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledPatch::createGeometry()
{
    sampledSurface::clearGeom();
    MeshStorage::clear();
    patchFaceLabels_.clear();

    if (patchIndex_ != -1)
    {
        const polyPatch& p = mesh().boundaryMesh()[patchIndex_];
        this->storedPoints() = p.localPoints();
        this->storedFaces()  = p.localFaces();

        // an identity map
        patchFaceLabels_.setSize(faces().size());
        forAll(patchFaceLabels_, i)
        {
            patchFaceLabels_[i] = i;
        }

        // triangulate uses remapFaces()
        // - this is somewhat less efficient since it recopies the faces
        // that we just created, but we probably don't want to do this
        // too often anyhow.
        if (triangulate_)
        {
            MeshStorage::triangulate();
        }
    }

    if (debug)
    {
        print(Pout);
        Pout << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledPatch::sampledPatch
(
    const word& name,
    const polyMesh& mesh,
    const word& patchName,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    patchName_(patchName),
    patchIndex_(mesh.boundaryMesh().findPatchID(patchName_)),
    triangulate_(triangulate),
    patchFaceLabels_(0)
{
    createGeometry();
}


Foam::sampledPatch::sampledPatch
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    patchName_(dict.lookup("patchName")),
    patchIndex_(mesh.boundaryMesh().findPatchID(patchName_)),
    triangulate_(dict.lookupOrDefault("triangulate", false)),
    patchFaceLabels_(0)
{
    createGeometry();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledPatch::~sampledPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledPatch::correct(const bool meshChanged)
{
    if (meshChanged)
    {
        createGeometry();
    }
}


// remap action on triangulation
void Foam::sampledPatch::remapFaces
(
    const UList<label>& faceMap
)
{
    // recalculate the cells cut
    if (&faceMap && faceMap.size())
    {
        MeshStorage::remapFaces(faceMap);
//
//        List<label> newCutCells(faceMap.size());
//        forAll(faceMap, faceI)
//        {
//            newCutCells[faceI] = cutCells_[faceMap[faceI]];
//        }
//        cutCells_.transfer(newCutCells);
    }
}



Foam::tmp<Foam::scalarField>
Foam::sampledPatch::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledPatch::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledPatch::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledPatch::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledPatch::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledPatch::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledPatch::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledPatch::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledPatch::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledPatch::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledPatch::print(Ostream& os) const
{
    os  << "sampledPatch: " << name() << " :"
        << "  patch:" << patchName()
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
