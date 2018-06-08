/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "surfMeshSampleDiscrete.H"
#include "MeshedSurfaces.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshSampleDiscrete, 0);

    // Add under name "sampledTriSurfaceMesh"
    // for symmetry with regular sampledSurface
    addNamedToRunTimeSelectionTable
    (
        surfMeshSample,
        surfMeshSampleDiscrete,
        word,
        sampledTriSurfaceMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMeshSampleDiscrete::transferContent()
{
    getOrCreateSurfMesh().transfer
    (
        static_cast<SurfaceSource&>(*this)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSampleDiscrete::surfMeshSampleDiscrete
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName,
    const discreteSurface::samplingSource sampleSource
)
:
    surfMeshSample(name, mesh),
    SurfaceSource(mesh, surfaceName, sampleSource, false) // no interpolate
{}


Foam::surfMeshSampleDiscrete::surfMeshSampleDiscrete
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    surfMeshSample(name, mesh),
    SurfaceSource(mesh, dict, false) // no interpolate
{}


Foam::surfMeshSampleDiscrete::surfMeshSampleDiscrete
(
    const word& name,
    const polyMesh& mesh,
    const triSurface& surface,
    const word& sampleSourceName
)
:
    surfMeshSample(name, mesh),
    SurfaceSource(name, mesh, surface, sampleSourceName, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfMeshSampleDiscrete::needsUpdate() const
{
    return SurfaceSource::needsUpdate();
}


bool Foam::surfMeshSampleDiscrete::expire()
{
    return SurfaceSource::expire();
}


bool Foam::surfMeshSampleDiscrete::update()
{
    if (SurfaceSource::update())
    {
        transferContent();
        return true;
    }

    return false;
}


bool Foam::surfMeshSampleDiscrete::update(const treeBoundBox& bb)
{
    if (SurfaceSource::update(bb))
    {
        transferContent();
        return true;
    }

    return false;
}


bool Foam::surfMeshSampleDiscrete::sample
(
    const word& fieldName,
    const word& sampleScheme
) const
{
    return
    (
        sampleType<scalar>(fieldName, sampleScheme)
     || sampleType<vector>(fieldName, sampleScheme)
     || sampleType<sphericalTensor>(fieldName, sampleScheme)
     || sampleType<symmTensor>(fieldName, sampleScheme)
     || sampleType<tensor>(fieldName, sampleScheme)
    );
}


// ************************************************************************* //
