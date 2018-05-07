/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "surfMeshSampleDistanceSurface.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"
#include "coordinateSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshSampleDistanceSurface, 0);
    addNamedToRunTimeSelectionTable
    (
        surfMeshSample,
        surfMeshSampleDistanceSurface,
        word,
        distanceSurface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSampleDistanceSurface::surfMeshSampleDistanceSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    surfMeshSample(name, mesh, dict),
    SurfaceSource(name, mesh, dict),
    needsUpdate_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfMeshSampleDistanceSurface::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::surfMeshSampleDistanceSurface::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::surfMeshSampleDistanceSurface::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    SurfaceSource::createGeometry();

    // Transfer content
    getOrCreateSurfMesh().transfer
    (
        SurfaceSource::surface()
    );

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


bool Foam::surfMeshSampleDistanceSurface::sample
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


void Foam::surfMeshSampleDistanceSurface::print(Ostream& os) const
{
    os  << "distanceSurface: " << name() << " :"
        << "  surface:" << surfaceName()
        << "  distance:" << distance()
        << "  faces:" << surface().faces().size()
        << "  points:" << surface().points().size();
}


// ************************************************************************* //
