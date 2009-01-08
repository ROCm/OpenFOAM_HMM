/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "PrimitivePatchExtra.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Finds area, starting at faceI, delimited by borderEdge. Marks all visited
// faces (from face-edge-face walk) with currentZone.
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::markZone
(
    const UList<bool>& borderEdge,
    const label faceI,
    const label currentZone,
    labelList& faceZone
) const
{
    // List of faces whose faceZone has been set.
    labelList changedFaces(1, faceI);

    const labelListList& faceEs = this->faceEdges();
    const labelListList& eFaces = this->edgeFaces();

    while (true)
    {
        // Pick up neighbours of changedFaces
        DynamicList<label> newChangedFaces(2*changedFaces.size());

        forAll(changedFaces, i)
        {
            label faceI = changedFaces[i];

            const labelList& fEdges = faceEs[faceI];

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (!borderEdge[edgeI])
                {
                    const labelList& eFaceLst = eFaces[edgeI];

                    forAll(eFaceLst, j)
                    {
                        label nbrFaceI = eFaceLst[j];

                        if (faceZone[nbrFaceI] == -1)
                        {
                            faceZone[nbrFaceI] = currentZone;
                            newChangedFaces.append(nbrFaceI);
                        }
                        else if (faceZone[nbrFaceI] != currentZone)
                        {
                            FatalErrorIn
                            (
                                "PrimitivePatchExtra<Face, FaceList, PointField>::markZone"
                                "(const boolList&, const label, const label, labelList&) const"
                            )
                                << "Zones " << faceZone[nbrFaceI]
                                << " at face " << nbrFaceI
                                << " connects to zone " << currentZone
                                << " at face " << faceI
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        if (newChangedFaces.size() == 0)
        {
            break;
        }

        // New dynamicList: can leave dynamicList unshrunk
        changedFaces.transfer(newChangedFaces);
    }
}


// Finds areas delimited by borderEdge (or 'real' edges).
// Fills faceZone accordingly
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::label Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
markZones
(
    const UList<bool>& borderEdge,
    labelList& faceZone
) const
{
    const label numEdges = this->nEdges();
    const label numFaces = this->size();

    if (borderEdge.size() != numEdges)
    {
        FatalErrorIn
        (
            "PrimitivePatchExtra<Face, FaceList, PointField>::markZones"
            "(const boolList&, labelList&)"
        )
            << "borderEdge boolList not same size as number of edges" << endl
            << "borderEdge:" << borderEdge.size() << endl
            << "nEdges    :" << numEdges
            << exit(FatalError);
    }

    faceZone.setSize(numFaces);
    faceZone = -1;

    label zoneI = 0;
    label startFaceI = 0;

    while (true)
    {
        // Find first non-visited face
        for (; startFaceI < numFaces; startFaceI++)
        {
            if (faceZone[startFaceI] == -1)
            {
                faceZone[startFaceI] = zoneI;
                markZone(borderEdge, startFaceI, zoneI, faceZone);
                break;
            }
        }

        if (startFaceI >= numFaces)
        {
            // Finished
            break;
        }

        zoneI++;
    }

    return zoneI;
}



// Finds areas delimited by borderEdge (or 'real' edges).
// Fills faceZone accordingly
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
subsetMap
(
    const UList<bool>& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const List<Face>& locFaces = this->localFaces();
    const label numPoints = this->nPoints();

    label faceI = 0;
    label pointI = 0;

    faceMap.setSize(locFaces.size());
    pointMap.setSize(numPoints);

    boolList pointHad(numPoints, false);

    forAll(include, oldFaceI)
    {
        if (include[oldFaceI])
        {
            // Store new faces compact
            faceMap[faceI++] = oldFaceI;

            // Renumber labels for face
            const Face& f = locFaces[oldFaceI];

            forAll(f, fp)
            {
                const label ptLabel = f[fp];
                if (!pointHad[ptLabel])
                {
                    pointHad[ptLabel] = true;
                    pointMap[pointI++] = ptLabel;
                }
            }
        }
    }

    // Trim
    faceMap.setSize(faceI);
    pointMap.setSize(pointI);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
