/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "sampledDiscreteSurface.H"
#include "meshSearch.H"
#include "Tuple2.H"
#include "globalIndex.H"
#include "treeDataCell.H"
#include "treeDataFace.H"
#include "meshTools.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledDiscreteSurface, 0);
    addToRunTimeSelectionTable
    (
        sampledSurface,
        sampledDiscreteSurface,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledDiscreteSurface::sampledDiscreteSurface
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName,
    const discreteSurface::samplingSource sampleSource
)
:
    sampledSurface(name, mesh),
    SurfaceSource(mesh, surfaceName, sampleSource)
{}


Foam::sampledDiscreteSurface::sampledDiscreteSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    SurfaceSource(mesh, dict)
{}


Foam::sampledDiscreteSurface::sampledDiscreteSurface
(
    const word& name,
    const polyMesh& mesh,
    const triSurface& surface,
    const word& sampleSourceName
)
:
    sampledSurface(name, mesh),
    SurfaceSource(name, mesh, surface, sampleSourceName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledDiscreteSurface::~sampledDiscreteSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledDiscreteSurface::needsUpdate() const
{
    return SurfaceSource::needsUpdate();
}


bool Foam::sampledDiscreteSurface::expire()
{
    if (SurfaceSource::expire())
    {
        // merged information etc
        sampledSurface::clearGeom();

        return true;
    }

    return false;
}


bool Foam::sampledDiscreteSurface::update()
{
    return SurfaceSource::update();
}


bool Foam::sampledDiscreteSurface::update(const treeBoundBox& bb)
{
    return SurfaceSource::update(bb);
}


bool Foam::sampledDiscreteSurface::sampleAndStore
(
    const objectRegistry& store,
    const word& fieldName
) const
{
    return SurfaceSource::sampleAndStore(store, fieldName);
}


Foam::tmp<Foam::scalarField> Foam::sampledDiscreteSurface::sample
(
    const volScalarField& vField
) const
{
    return SurfaceSource::sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledDiscreteSurface::sample
(
    const volVectorField& vField
) const
{
    return SurfaceSource::sampleField(vField);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledDiscreteSurface::sample
(
    const volSphericalTensorField& vField
) const
{
    return SurfaceSource::sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledDiscreteSurface::sample
(
    const volSymmTensorField& vField
) const
{
    return SurfaceSource::sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledDiscreteSurface::sample
(
    const volTensorField& vField
) const
{
    return SurfaceSource::sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledDiscreteSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return SurfaceSource::interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledDiscreteSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return SurfaceSource::interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledDiscreteSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return SurfaceSource::interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledDiscreteSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return SurfaceSource::interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledDiscreteSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return SurfaceSource::interpolateField(interpolator);
}


void Foam::sampledDiscreteSurface::print(Ostream& os) const
{
    os  << "sampledDiscreteSurface: " << name() << " :";
    SurfaceSource::print(os);
}


// ************************************************************************* //
