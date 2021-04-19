/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "cellDistFuncs.H"
#include "polyMesh.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellDistFuncs, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellDistFuncs::cellDistFuncs(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelHashSet Foam::cellDistFuncs::getPatchIDs
(
    const UList<wordRe>& patchNames
) const
{
    return mesh().boundaryMesh().patchSet(patchNames, false);
}


// Return smallest true distance from p to any of wallFaces.
// Note that even if normal hits face we still check other faces.
// Note that wallFaces is untruncated and we explicitly pass in size.
Foam::scalar Foam::cellDistFuncs::smallestDist
(
    const point& p,
    const polyPatch& patch,
    const labelUList& wallFaces,
    label& minFacei
) const
{
    const pointField& points = patch.points();

    scalar minDist = GREAT;
    minFacei = -1;

    for (const label patchFacei : wallFaces)
    {
        const pointHit curHit = patch[patchFacei].nearestPoint(p, points);

        if (curHit.distance() < minDist)
        {
            minDist = curHit.distance();
            minFacei = patch.start() + patchFacei;
        }
    }

    return minDist;
}


// Get point neighbours of facei (including facei). Returns number of faces.
// Note: does not allocate storage but does use linear search to determine
// uniqueness. For polygonal faces this might be quite inefficient.
void Foam::cellDistFuncs::getPointNeighbours
(
    const primitivePatch& patch,
    const label patchFacei,
    DynamicList<label>& neighbours
) const
{
    neighbours.clear();

    // Add myself
    neighbours.append(patchFacei);

    // Add all face neighbours
    const labelList& faceNeighbours = patch.faceFaces()[patchFacei];

    for (const label nbr : faceNeighbours)
    {
        neighbours.appendUniq(nbr);
    }

    // Add all point-only neighbours by linear searching in edge neighbours.
    // Assumes that point-only neighbours are not using multiple points on
    // face.

    const face& f = patch.localFaces()[patchFacei];

    forAll(f, fp)
    {
        label pointi = f[fp];

        const labelList& pointNbs = patch.pointFaces()[pointi];

        for (const label facei : pointNbs)
        {
            // Check for facei in edge-neighbours part of neighbours
            neighbours.appendUniq(facei);
        }
    }


    if (debug)
    {
        // Check for duplicates

        // Use hashSet to determine nbs.
        labelHashSet nbs(4*f.size());

        forAll(f, fp)
        {
            const labelList& pointNbs = patch.pointFaces()[f[fp]];
            nbs.insert(pointNbs);
        }

        // Subtract ours.
        for (const label nb : neighbours)
        {
            if (!nbs.found(nb))
            {
                SeriousErrorInFunction
                    << "getPointNeighbours : patchFacei:" << patchFacei
                    << " verts:" << f << endl;

                forAll(f, fp)
                {
                    SeriousErrorInFunction
                        << "point:" << f[fp] << " pointFaces:"
                        << patch.pointFaces()[f[fp]] << endl;
                }

                for (const label facei : neighbours)
                {
                    SeriousErrorInFunction
                        << "fast nbr:" << facei
                        << endl;
                }

                FatalErrorInFunction
                    << "Problem: fast pointNeighbours routine included " << nb
                    << " which is not in proper neighbour list " << nbs.toc()
                    << abort(FatalError);
            }
            nbs.erase(nb);
        }

        if (nbs.size())
        {
            FatalErrorInFunction
                << "Problem: fast pointNeighbours routine did not find "
                << nbs.toc() << abort(FatalError);
        }
    }
}


// size of largest patch (out of supplied subset of patches)
Foam::label Foam::cellDistFuncs::maxPatchSize
(
    const labelHashSet& patchIDs
) const
{
    label maxSize = 0;

    forAll(mesh().boundaryMesh(), patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchi];

            maxSize = Foam::max(maxSize, patch.size());
        }
    }
    return maxSize;
}


// sum of patch sizes (out of supplied subset of patches)
Foam::label Foam::cellDistFuncs::sumPatchSize
(
    const labelHashSet& patchIDs
)
const
{
    label sum = 0;

    forAll(mesh().boundaryMesh(), patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchi];

            sum += patch.size();
        }
    }
    return sum;
}


// Gets nearest wall for cells next to wall
void Foam::cellDistFuncs::correctBoundaryFaceCells
(
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    Map<label>& nearestFace
) const
{
    // Size neighbours array for maximum possible (= size of largest patch)
    DynamicList<label> neighbours(maxPatchSize(patchIDs));

    // Correct all cells with face on wall
    const vectorField& cellCentres = mesh().cellCentres();
    const labelList& faceOwner = mesh().faceOwner();

    forAll(mesh().boundaryMesh(), patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchi];

            // Check cells with face on wall
            forAll(patch, patchFacei)
            {
                getPointNeighbours(patch, patchFacei, neighbours);

                label celli = faceOwner[patch.start() + patchFacei];

                label minFacei = -1;

                wallDistCorrected[celli] = smallestDist
                (
                    cellCentres[celli],
                    patch,
                    neighbours,
                    minFacei
                );

                // Store wallCell and its nearest neighbour
                nearestFace.insert(celli, minFacei);
            }
        }
    }
}


// Correct all cells connected to wall (via point) and not in nearestFace
void Foam::cellDistFuncs::correctBoundaryPointCells
(
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    Map<label>& nearestFace
) const
{
    // Correct all (non-visited) cells with point on wall

    const vectorField& cellCentres = mesh().cellCentres();

    forAll(mesh().boundaryMesh(), patchi)
    {
        if (patchIDs.found(patchi))
        {
            const polyPatch& patch = mesh().boundaryMesh()[patchi];

            const labelList& meshPoints = patch.meshPoints();
            const labelListList& pointFaces = patch.pointFaces();

            forAll(meshPoints, meshPointi)
            {
                const label vertI = meshPoints[meshPointi];

                const labelList& neighbours = mesh().pointCells(vertI);

                for (const label celli : neighbours)
                {
                    if (!nearestFace.found(celli))
                    {
                        const labelList& wallFaces = pointFaces[meshPointi];

                        label minFacei = -1;

                        wallDistCorrected[celli] = smallestDist
                        (
                            cellCentres[celli],
                            patch,
                            wallFaces,
                            minFacei
                        );

                        // Store wallCell and its nearest neighbour
                        nearestFace.insert(celli, minFacei);
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
