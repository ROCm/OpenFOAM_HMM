/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "surfMeshDiscreteSampler.H"
#include "MeshedSurfaces.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshDiscreteSampler, 0);

    // Add under name "sampledTriSurfaceMesh"
    // for symmetry with normal sampledSurface
    addNamedToRunTimeSelectionTable
    (
        surfMeshSampler,
        surfMeshDiscreteSampler,
        word,
        sampledTriSurfaceMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfMeshDiscreteSampler::transferContent()
{
    SurfaceSource& src = static_cast<SurfaceSource&>(*this);
    surfMesh& dst = getOrCreateSurfMesh();

    dst.transfer(src);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshDiscreteSampler::surfMeshDiscreteSampler
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName,
    const discreteSurface::samplingSource sampleSource
)
:
    surfMeshSampler(name, mesh),
    SurfaceSource(mesh, surfaceName, sampleSource, false) // no interpolate
{}


Foam::surfMeshDiscreteSampler::surfMeshDiscreteSampler
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    surfMeshSampler(name, mesh),
    SurfaceSource(mesh, dict, false) // no interpolate
{}


Foam::surfMeshDiscreteSampler::surfMeshDiscreteSampler
(
    const word& name,
    const polyMesh& mesh,
    const triSurface& surface,
    const word& sampleSourceName
)
:
    surfMeshSampler(name, mesh),
    SurfaceSource(name, mesh, surface, sampleSourceName, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMeshDiscreteSampler::~surfMeshDiscreteSampler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfMeshDiscreteSampler::needsUpdate() const
{
    return SurfaceSource::needsUpdate();
}


bool Foam::surfMeshDiscreteSampler::expire()
{
    return SurfaceSource::expire();
}


bool Foam::surfMeshDiscreteSampler::update()
{
    if (SurfaceSource::update())
    {
        transferContent();
        return true;
    }

    return false;
}


bool Foam::surfMeshDiscreteSampler::update(const treeBoundBox& bb)
{
    if (SurfaceSource::update(bb))
    {
        transferContent();
        return true;
    }

    return false;
}


bool Foam::surfMeshDiscreteSampler::sample
(
    const word& fieldName
) const
{
    return
    (
        sampleType<scalar>(fieldName)
     || sampleType<vector>(fieldName)
     || sampleType<sphericalTensor>(fieldName)
     || sampleType<symmTensor>(fieldName)
     || sampleType<tensor>(fieldName)
    );
}


// ************************************************************************* //
