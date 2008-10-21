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

#include "isoSurface.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"
#include "mergePoints.H"
#include "tetMatcher.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(isoSurface, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isoSurface::isoSurface
(
    const polyMesh& mesh,
    const scalarField& cellValues,
    const scalarField& pointValues,
    const scalar iso
)
:
    mesh_(mesh),
    cellValues_(cellValues),
    pointValues_(pointValues),
    iso_(iso)
{
    const pointField& cellCentres = mesh.cellCentres();

    tetMatcher tet;

    DynamicList<point> triPoints;
    DynamicList<label> triMeshCells;

    forAll(mesh.cells(), cellI)
    {
        label oldNPoints = triPoints.size();

        const cell& cFaces = mesh.cells()[cellI];

        if (tet.isA(mesh, cellI))
        {
            // For tets don't do cell-centre decomposition, just use the
            // tet points and values

            const face& f0 = mesh.faces()[cFaces[0]];

            // Get the other point
            const face& f1 = mesh.faces()[cFaces[1]];
            label oppositeI = -1;
            forAll(f1, fp)
            {
                oppositeI = f1[fp];

                if (findIndex(f0, oppositeI) == -1)
                {
                    break;
                }
            }

            vertexInterp
            (
                iso,
                pointValues[f0[0]],
                pointValues[f0[1]],
                pointValues[f0[2]],
                pointValues[oppositeI],

                mesh.points()[f0[0]],
                mesh.points()[f0[1]],
                mesh.points()[f0[2]],
                mesh.points()[oppositeI],

                triPoints
            );
        }
        else
        {
            const cell& cFaces = mesh.cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                label faceI = cFaces[cFaceI];
                const face& f = mesh.faces()[faceI];

                // Do a tetrahedrisation. Each face to cc becomes pyr.
                // Each pyr gets split into two tets by diagionalisation
                // of face. So
                // - f[0], f[1], f[2], cc
                // - f[0], f[2], f[3], cc

                for(label fp = 1; fp < f.size() - 1; fp++)
                {
                    vertexInterp
                    (
                        iso,
                        pointValues[f[0]],
                        pointValues[f[fp]],
                        pointValues[f[f.fcIndex(fp)]],
                        cellValues[cellI],

                        mesh.points()[f[0]],
                        mesh.points()[f[fp]],
                        mesh.points()[f[f.fcIndex(fp)]],
                        cellCentres[cellI],

                        triPoints
                    );
                }
            }
        }


        // Every three triPoints is a cell
        label nCells = (triPoints.size()-oldNPoints)/3;
        for (label i = 0; i < nCells; i++)
        {
            triMeshCells.append(cellI);
        }
    }

    triPoints.shrink();
    triMeshCells.shrink();
    meshCells_.transfer(triMeshCells);

    pointField newPoints;
    mergePoints(triPoints, SMALL, false, triPointMergeMap_, newPoints);

    DynamicList<labelledTri> tris(meshCells_.size());
    forAll(meshCells_, triI)
    {
        tris.append
        (
            labelledTri
            (
                triPointMergeMap_[3*triI],
                triPointMergeMap_[3*triI+1],
                triPointMergeMap_[3*triI+2],
                0
            )
        );
    }

    //- 1.5.x:
    //tris.shrink();
    //triSurface::operator=
    //(
    //    triSurface(tris, geometricSurfacePatchList(0), newPoints)
    //);

    triSurface::operator=
    (
        triSurface(tris, geometricSurfacePatchList(0), newPoints, true)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//Foam::tmp<Foam::scalarField>
//Foam::isoSurface::sample
//(
//    const volScalarField& vField
//) const
//{
//    return sampleField(vField);
//}
//
//
//Foam::tmp<Foam::vectorField>
//Foam::isoSurface::sample
//(
//    const volVectorField& vField
//) const
//{
//    return sampleField(vField);
//}
//
//Foam::tmp<Foam::sphericalTensorField>
//Foam::isoSurface::sample
//(
//    const volSphericalTensorField& vField
//) const
//{
//    return sampleField(vField);
//}
//
//
//Foam::tmp<Foam::symmTensorField>
//Foam::isoSurface::sample
//(
//    const volSymmTensorField& vField
//) const
//{
//    return sampleField(vField);
//}
//
//
//Foam::tmp<Foam::tensorField>
//Foam::isoSurface::sample
//(
//    const volTensorField& vField
//) const
//{
//    return sampleField(vField);
//}
//
//
//Foam::tmp<Foam::scalarField>
//Foam::isoSurface::interpolate
//(
//    const interpolation<scalar>& interpolator
//) const
//{
//    return interpolateField(interpolator);
//}
//
//
//Foam::tmp<Foam::vectorField>
//Foam::isoSurface::interpolate
//(
//    const interpolation<vector>& interpolator
//) const
//{
//    return interpolateField(interpolator);
//}
//
//Foam::tmp<Foam::sphericalTensorField>
//Foam::isoSurface::interpolate
//(
//    const interpolation<sphericalTensor>& interpolator
//) const
//{
//    return interpolateField(interpolator);
//}
//
//
//Foam::tmp<Foam::symmTensorField>
//Foam::isoSurface::interpolate
//(
//    const interpolation<symmTensor>& interpolator
//) const
//{
//    return interpolateField(interpolator);
//}
//
//
//Foam::tmp<Foam::tensorField>
//Foam::isoSurface::interpolate
//(
//    const interpolation<tensor>& interpolator
//) const
//{
//    return interpolateField(interpolator);
//}


// ************************************************************************* //
