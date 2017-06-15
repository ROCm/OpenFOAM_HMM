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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledTriSurfaceMeshNormal::~sampledTriSurfaceMeshNormal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledTriSurfaceMeshNormal::sample
(
    const GeometricField<vector, fvPatchField, volMesh>& vField
) const
{
    tmp<Field<vector>> tfld(new Field<vector>(size(), vector::zero));
    tfld.ref().replace
    (
        0,
        meshedSurface::faceNormals()
       &sampledTriSurfaceMesh::sample(vField)
    );
    return tfld;
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::sampledTriSurfaceMeshNormal::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    // One value per vertex
    tmp<vectorField> tn
    (
        new vectorField
        (
            points().size(),
            vector::zero
        )
    );

    pointField allNormals(tn().size(), vector::zero);
    UIndirectList<vector>(allNormals, meshPoints()) = pointNormals();

    tn.ref().replace
    (
        0,
        allNormals
       &sampledTriSurfaceMesh::interpolate(interpolator)
    );
    return tn;
}


// ************************************************************************* //
