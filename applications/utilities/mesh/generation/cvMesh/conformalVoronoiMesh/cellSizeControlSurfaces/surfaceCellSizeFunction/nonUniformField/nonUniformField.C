/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "nonUniformField.H"
#include "triSurfaceMesh.H"
#include "searchableSurface.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonUniformField, 0);
    addToRunTimeSelectionTable
    (
        surfaceCellSizeFunction,
        nonUniformField,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonUniformField::nonUniformField
(
    const dictionary& cellSizeFunctionDict,
    const searchableSurface& surface,
    const scalar& defaultCellSize
)
:
    surfaceCellSizeFunction
    (
        typeName,
        cellSizeFunctionDict,
        surface,
        defaultCellSize
    ),
    surfaceTriMesh_(refCast<const triSurfaceMesh>(surface)),
    surfaceCellSize_
    (
        IOobject
        (
            "surfaceCellSize.cellSize",
            surface.time().constant(),
            "triSurface",
            surface.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surfaceTriMesh_,
        dimLength,
        false
    ),
    cellSizeCalculationType_
    (
        cellSizeCalculationType::New
        (
            coeffsDict(),
            surfaceTriMesh_,
            defaultCellSize
        )
    ),
    pointCellSize_(),
    patch_()
{
    Info<< incrIndent;

    surfaceCellSize_ = cellSizeCalculationType_().load();

    Info<< indent << "Cell size field statistics:" << nl
        << indent << "    Minimum: " << min(surfaceCellSize_).value() << nl
        << indent << "    Average: " << average(surfaceCellSize_).value() << nl
        << indent << "    Maximum: " << max(surfaceCellSize_).value() << endl;

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalar& Foam::nonUniformField::surfaceSize(const label index) const
{
    return surfaceCellSize_[index];
}


const Foam::scalar& Foam::nonUniformField::refineSurfaceSize(const label index)
{
    surfaceCellSize_[index] *= refinementFactor_;

    return surfaceCellSize_[index];
}


Foam::scalar Foam::nonUniformField::interpolate
(
    const point& pt,
    const label index
) const
{
    if (patch_.empty())
    {
        recalculateInterpolation();
    }

    const face& faceHitByPt = surfaceTriMesh_.triSurface::operator[](index);

    const pointField& pts = surfaceTriMesh_.points();
    const Map<label>& pMap = surfaceTriMesh_.meshPointMap();

    triPointRef tri
    (
        pts[faceHitByPt[0]],
        pts[faceHitByPt[1]],
        pts[faceHitByPt[2]]
    );

    scalarList bary(3, 0.0);

    tri.barycentric(pt, bary);

    return pointCellSize_()[pMap[faceHitByPt[0]]]*bary[0]
         + pointCellSize_()[pMap[faceHitByPt[1]]]*bary[1]
         + pointCellSize_()[pMap[faceHitByPt[2]]]*bary[2];
}


void Foam::nonUniformField::recalculateInterpolation() const
{
    patch_.reset(new primitivePatchInterpolation(surfaceTriMesh_));

    pointCellSize_.reset
    (
        new triSurfacePointScalarField
        (
            IOobject
            (
                "pointCellSize.cellSize",
                surfaceTriMesh_.searchableSurface::time().constant(),
                "triSurface",
                surfaceTriMesh_.searchableSurface::time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            surfaceTriMesh_,
            dimLength,
            patch_().faceToPointInterpolate(surfaceCellSize_)
        )
    );
}


// ************************************************************************* //
