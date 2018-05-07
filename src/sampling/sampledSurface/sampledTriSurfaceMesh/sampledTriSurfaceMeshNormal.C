/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "sampledTriSurfaceMeshNormal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledTriSurfaceMeshNormal, 0);
    addToRunTimeSelectionTable
    (
        sampledSurface,
        sampledTriSurfaceMeshNormal,
        word
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledTriSurfaceMeshNormal::sampledTriSurfaceMeshNormal
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName,
    const samplingSource sampleSource
)
:
    sampledTriSurfaceMesh(name, mesh, surfaceName, sampleSource)
{}


Foam::sampledTriSurfaceMeshNormal::sampledTriSurfaceMeshNormal
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledTriSurfaceMesh(name, mesh, dict)
{}


Foam::sampledTriSurfaceMeshNormal::sampledTriSurfaceMeshNormal
(
    const word& name,
    const polyMesh& mesh,
    const triSurface& surface,
    const word& sampleSourceName
)
:
    sampledTriSurfaceMesh(name, mesh, surface, sampleSourceName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledTriSurfaceMeshNormal::sample
(
    const interpolation<vector>& sampler
) const
{
    auto tvalues = tmp<Field<vector>>::New(size(), Zero);

    tvalues.ref().replace
    (
        0,
        meshedSurface::faceNormals()
       &sampledTriSurfaceMesh::sample(sampler)
    );

    return tvalues;
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledTriSurfaceMeshNormal::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    auto tvalues = tmp<Field<vector>>::New(points().size(), Zero);

    pointField allNormals(points().size(), Zero);
    UIndirectList<vector>(allNormals, meshPoints()) = pointNormals();

    tvalues.ref().replace
    (
        0,
        allNormals
       &sampledTriSurfaceMesh::interpolate(interpolator)
    );

    return tvalues;
}


// ************************************************************************* //
