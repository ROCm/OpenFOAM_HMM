/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "duplicatePoints.H"
#include "localPointRegion.H"
#include "polyTopoChange.H"
#include "polyAddPoint.H"
#include "polyModifyFace.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(duplicatePoints, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::duplicatePoints::duplicatePoints(const polyMesh& mesh)
:
    mesh_(mesh),
    duplicates_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::duplicatePoints::setRefinement
(
    const labelList& nonManifPoints,
    const localPointRegion& regionSide,
    polyTopoChange& meshMod
)
{
    const labelList& faceRegion = regionSide.faceRegion();
    const labelList& meshFaces = regionSide.meshFaces();
    const Map<label>& localFaces = regionSide.localFaces();

    // New faces
    faceList newFaces(faceRegion.size());

    // Initialise to current faces
    forAll(meshFaces, localFaceI)
    {
        newFaces[localFaceI] = mesh_.faces()[meshFaces[localFaceI]];
    }


    duplicates_.setSize(regionSide.pointRegions().size());

    // Create new point for point with more than one region
    forAll(regionSide.pointRegions(), localPointI)
    {
        label pointI = nonManifPoints[localPointI];
        const labelList& regions = regionSide.pointRegions()[localPointI];

        if (regions.size() > 1)
        {
            // First region for point gets the original point label
            duplicates_[localPointI].setSize(regions.size());
            duplicates_[localPointI][0] = pointI;

            for (label i = 1; i < regions.size(); i++)
            {
                // Add a point for the point for all faces with the same region
                label addedPointI = meshMod.setAction
                (
                    polyAddPoint
                    (
                        mesh_.points()[pointI], // point
                        pointI,                 // master point
                        -1,                     // zone for point
                        true                    // supports a cell
                    )
                );

                // Store added point
                duplicates_[localPointI][i] = addedPointI;

                const labelList& pFaces = mesh_.pointFaces()[pointI];

                // Replace all the vertices with the same region with the new
                // point label.
                forAll(pFaces, pFaceI)
                {
                    label faceI = pFaces[pFaceI];
                    label localFaceI = localFaces[faceI];

                    if (faceRegion[localFaceI] == regions[i])
                    {
                        const face& f = mesh_.faces()[faceI];

                        forAll(f, fp)
                        {
                            if (f[fp] == pointI)
                            {
                                newFaces[localFaceI][fp] = addedPointI;
                            }
                        }
                    }
                }
            }
        }
    }

    // Modify the faces

    forAll(meshFaces, localFaceI)
    {
        label faceI = meshFaces[localFaceI];

        label own = mesh_.faceOwner()[faceI];
        label nei = -1;
        label patchID = -1;
        if (mesh_.isInternalFace(faceI))
        {
            nei = mesh_.faceNeighbour()[faceI];
        }
        else
        {
            patchID = mesh_.boundaryMesh().whichPatch(faceI);
        }

        // Get current zone info
        label zoneID = mesh_.faceZones().whichZone(faceI);
        bool zoneFlip = false;
        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];
            zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
        }

        meshMod.setAction
        (
            polyModifyFace
            (
                newFaces[localFaceI],   // modified face
                faceI,                  // label of face being modified
                own,                    // owner
                nei,                    // neighbour
                false,                  // face flip
                patchID,                // patch for face
                false,                  // remove from zone
                zoneID,                 // zone for face
                zoneFlip                // face flip in zone
            )
        );
    }
}


void Foam::duplicatePoints::updateMesh(const mapPolyMesh& map)
{
    forAll(duplicates_, masterI)
    {
        inplaceRenumber(map.reversePointMap(), duplicates_[masterI]);
    }
}


// ************************************************************************* //
