/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "tetOverlapVolume.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::tetOverlapVolume, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetOverlapVolume::tetOverlapVolume()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetOverlapVolume::~tetOverlapVolume()
{}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::tetOverlapVolume::planeIntersection
(
    const FixedList<scalar, 4>& d,
    const tetPoints& t,
    const label negI,
    const label posI
)
{
    return (d[posI]*t[negI] - d[negI]*t[posI])/(-d[negI]+d[posI]);
}


template <class TetOp>
inline void Foam::tetOverlapVolume::decomposePrism
(
    const FixedList<point, 6>& points,
    TetOp& op
)
{
    op(tetPoints(points[1], points[3], points[2], points[0]));
    op(tetPoints(points[1], points[2], points[3], points[4]));
    op(tetPoints(points[4], points[2], points[3], points[5]));
}


template <class AboveTetOp, class BelowTetOp>
inline void Foam::tetOverlapVolume::tetSliceWithPlane
(
    const tetPoints& tet,
    const plane& pl,

    AboveTetOp& aboveOp,
    BelowTetOp& belowOp
)
{
    // Distance to plane
    FixedList<scalar, 4> d;
    label nPos = 0;
    forAll(tet, i)
    {
        d[i] = ((tet[i]-pl.refPoint()) & pl.normal());
        if (d[i] > 0)
        {
            nPos++;
        }
    }

    if (nPos == 4)
    {
        aboveOp(tet);
    }
    else if (nPos == 3)
    {
        // Sliced into below tet and above prism. Prism gets split into
        // two tets.

        // Find the below tet
        label i0 = -1;
        forAll(d, i)
        {
            if (d[i] <= 0)
            {
                i0 = i;
                break;
            }
        }

        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);
        label i3 = d.fcIndex(i2);

        point p01 = planeIntersection(d, tet, i0, i1);
        point p02 = planeIntersection(d, tet, i0, i2);
        point p03 = planeIntersection(d, tet, i0, i3);

        // i0 = tetCell vertex 0: p01,p02,p03 outwards pointing triad
        //          ,,         1 :     ,,     inwards pointing triad
        //          ,,         2 :     ,,     outwards pointing triad
        //          ,,         3 :     ,,     inwards pointing triad

        //Pout<< "Split 3pos tet " << tet << " d:" << d << " into" << nl;

        if (i0 == 0 || i0 == 2)
        {
            tetPoints t(tet[i0], p01, p02, p03);
            //Pout<< "    belowtet:" << t << " around i0:" << i0 << endl;
            //checkTet(t, "nPos 3, belowTet i0==0 or 2");
            belowOp(t);

            // Prism
            FixedList<point, 6> p;
            p[0] = tet[i1];
            p[1] = tet[i3];
            p[2] = tet[i2];
            p[3] = p01;
            p[4] = p03;
            p[5] = p02;
            //Pout<< "    aboveprism:" << p << endl;
            decomposePrism(p, aboveOp);
        }
        else
        {
            tetPoints t(p01, p02, p03, tet[i0]);
            //Pout<< "    belowtet:" << t << " around i0:" << i0 << endl;
            //checkTet(t, "nPos 3, belowTet i0==1 or 3");
            belowOp(t);

            // Prism
            FixedList<point, 6> p;
            p[0] = tet[i3];
            p[1] = tet[i1];
            p[2] = tet[i2];
            p[3] = p03;
            p[4] = p01;
            p[5] = p02;
            //Pout<< "    aboveprism:" << p << endl;
            decomposePrism(p, aboveOp);
        }
    }
    else if (nPos == 2)
    {
        // Tet cut into two prisms. Determine the positive one.
        label pos0 = -1;
        label pos1 = -1;
        label neg0 = -1;
        label neg1 = -1;
        forAll(d, i)
        {
            if (d[i] > 0)
            {
                if (pos0 == -1)
                {
                    pos0 = i;
                }
                else
                {
                    pos1 = i;
                }
            }
            else
            {
                if (neg0 == -1)
                {
                    neg0 = i;
                }
                else
                {
                    neg1 = i;
                }
            }
        }

        //Pout<< "Split 2pos tet " << tet << " d:" << d
        //    << " around pos0:" << pos0 << " pos1:" << pos1
        //    << " neg0:" << neg0 << " neg1:" << neg1 << " into" << nl;

        const edge posEdge(pos0, pos1);

        if (posEdge == edge(0, 1))
        {
            point p02 = planeIntersection(d, tet, 0, 2);
            point p03 = planeIntersection(d, tet, 0, 3);
            point p12 = planeIntersection(d, tet, 1, 2);
            point p13 = planeIntersection(d, tet, 1, 3);
            // Split the resulting prism
            {
                FixedList<point, 6> p;
                p[0] = tet[0];
                p[1] = p02;
                p[2] = p03;
                p[3] = tet[1];
                p[4] = p12;
                p[5] = p13;
                //Pout<< "    01 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p;
                p[0] = tet[2];
                p[1] = p02;
                p[2] = p12;
                p[3] = tet[3];
                p[4] = p03;
                p[5] = p13;
                //Pout<< "    01 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(1, 2))
        {
            point p01 = planeIntersection(d, tet, 0, 1);
            point p13 = planeIntersection(d, tet, 1, 3);
            point p02 = planeIntersection(d, tet, 0, 2);
            point p23 = planeIntersection(d, tet, 2, 3);
            // Split the resulting prism
            {
                FixedList<point, 6> p;
                p[0] = tet[1];
                p[1] = p01;
                p[2] = p13;
                p[3] = tet[2];
                p[4] = p02;
                p[5] = p23;
                //Pout<< "    12 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p;
                p[0] = tet[3];
                p[1] = p23;
                p[2] = p13;
                p[3] = tet[0];
                p[4] = p02;
                p[5] = p01;
                //Pout<< "    12 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(2, 0))
        {
            point p01 = planeIntersection(d, tet, 0, 1);
            point p03 = planeIntersection(d, tet, 0, 3);
            point p12 = planeIntersection(d, tet, 1, 2);
            point p23 = planeIntersection(d, tet, 2, 3);
            // Split the resulting prism
            {
                FixedList<point, 6> p;
                p[0] = tet[2];
                p[1] = p12;
                p[2] = p23;
                p[3] = tet[0];
                p[4] = p01;
                p[5] = p03;
                //Pout<< "    20 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p;
                p[0] = tet[1];
                p[1] = p12;
                p[2] = p01;
                p[3] = tet[3];
                p[4] = p23;
                p[5] = p03;
                //Pout<< "    20 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(0, 3))
        {
            point p01 = planeIntersection(d, tet, 0, 1);
            point p02 = planeIntersection(d, tet, 0, 2);
            point p13 = planeIntersection(d, tet, 1, 3);
            point p23 = planeIntersection(d, tet, 2, 3);
            // Split the resulting prism
            {
                FixedList<point, 6> p;
                p[0] = tet[3];
                p[1] = p23;
                p[2] = p13;
                p[3] = tet[0];
                p[4] = p02;
                p[5] = p01;
                //Pout<< "    03 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p;
                p[0] = tet[2];
                p[1] = p23;
                p[2] = p02;
                p[3] = tet[1];
                p[4] = p13;
                p[5] = p01;
                //Pout<< "    03 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(1, 3))
        {
            point p01 = planeIntersection(d, tet, 0, 1);
            point p12 = planeIntersection(d, tet, 1, 2);
            point p03 = planeIntersection(d, tet, 0, 3);
            point p23 = planeIntersection(d, tet, 2, 3);
            // Split the resulting prism
            {
                FixedList<point, 6> p;
                p[0] = tet[1];
                p[1] = p12;
                p[2] = p01;
                p[3] = tet[3];
                p[4] = p23;
                p[5] = p03;
                //Pout<< "    13 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p;
                p[0] = tet[2];
                p[1] = p12;
                p[2] = p23;
                p[3] = tet[0];
                p[4] = p01;
                p[5] = p03;
                //Pout<< "    13 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else if (posEdge == edge(2, 3))
        {
            point p02 = planeIntersection(d, tet, 0, 2);
            point p12 = planeIntersection(d, tet, 1, 2);
            point p03 = planeIntersection(d, tet, 0, 3);
            point p13 = planeIntersection(d, tet, 1, 3);
            // Split the resulting prism
            {
                FixedList<point, 6> p;
                p[0] = tet[2];
                p[1] = p02;
                p[2] = p12;
                p[3] = tet[3];
                p[4] = p03;
                p[5] = p13;
                //Pout<< "    23 aboveprism:" << p << endl;
                decomposePrism(p, aboveOp);
            }
            {
                FixedList<point, 6> p;
                p[0] = tet[0];
                p[1] = p02;
                p[2] = p03;
                p[3] = tet[1];
                p[4] = p12;
                p[5] = p13;
                //Pout<< "    23 belowprism:" << p << endl;
                decomposePrism(p, belowOp);
            }
        }
        else
        {
            FatalErrorIn("tetSliceWithPlane(..)") << "Missed edge:" << posEdge
                << abort(FatalError);
        }
    }
    else if (nPos == 1)
    {
        // Find the positive tet
        label i0 = -1;
        forAll(d, i)
        {
            if (d[i] > 0)
            {
                i0 = i;
                break;
            }
        }

        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);
        label i3 = d.fcIndex(i2);

        point p01 = planeIntersection(d, tet, i0, i1);
        point p02 = planeIntersection(d, tet, i0, i2);
        point p03 = planeIntersection(d, tet, i0, i3);

        //Pout<< "Split 1pos tet " << tet << " d:" << d << " into" << nl;

        if (i0 == 0 || i0 == 2)
        {
            tetPoints t(tet[i0], p01, p02, p03);
            //Pout<< "    abovetet:" << t << " around i0:" << i0 << endl;
            //checkTet(t, "nPos 1, aboveTets i0==0 or 2");
            aboveOp(t);

            // Prism
            FixedList<point, 6> p;
            p[0] = tet[i1];
            p[1] = tet[i3];
            p[2] = tet[i2];
            p[3] = p01;
            p[4] = p03;
            p[5] = p02;
            //Pout<< "    belowprism:" << p << endl;
            decomposePrism(p, belowOp);
        }
        else
        {
            tetPoints t(p01, p02, p03, tet[i0]);
            //Pout<< "    abovetet:" << t << " around i0:" << i0 << endl;
            //checkTet(t, "nPos 1, aboveTets i0==1 or 3");
            aboveOp(t);

            // Prism
            FixedList<point, 6> p;
            p[0] = tet[i3];
            p[1] = tet[i1];
            p[2] = tet[i2];
            p[3] = p03;
            p[4] = p01;
            p[5] = p02;
            //Pout<< "    belowprism:" << p << endl;
            decomposePrism(p, belowOp);
        }
    }
    else    // nPos == 0
    {
        belowOp(tet);
    }
}


void Foam::tetOverlapVolume::tetTetOverlap
(
    const tetPoints& tetA,
    const tetPoints& tetB,
    FixedList<tetPoints, 200>& insideTets,
    label& nInside,
    FixedList<tetPoints, 200>& outsideTets,
    label& nOutside
)
{
    // Work storage
    FixedList<tetPoints, 200> cutInsideTets;
    label nCutInside = 0;

    storeTetOp inside(insideTets, nInside);
    storeTetOp cutInside(cutInsideTets, nCutInside);
    dummyTetOp outside;



    // Cut tetA with all inwards pointing faces of tetB. Any tets remaining
    // in aboveTets are inside tetB.

    {
        // face0
        plane pl0(tetB[1], tetB[3], tetB[2]);

        // Cut and insert subtets into cutInsideTets (either by getting
        // an index from freeSlots or by appending to insideTets) or
        // insert into outsideTets
        tetSliceWithPlane(tetA, pl0, cutInside, outside);
    }

    if (nCutInside == 0)
    {
        nInside = nCutInside;
        return;
    }

    {
        // face1
        plane pl1(tetB[0], tetB[2], tetB[3]);

        nInside = 0;

        for (label i = 0; i < nCutInside; i++)
        {
            tetSliceWithPlane(cutInsideTets[i], pl1, inside, outside);
        }

        if (nInside == 0)
        {
            return;
        }
    }

    {
        // face2
        plane pl2(tetB[0], tetB[3], tetB[1]);

        nCutInside = 0;

        for (label i = 0; i < nInside; i++)
        {
            tetSliceWithPlane(insideTets[i], pl2, cutInside, outside);
        }

        if (nCutInside == 0)
        {
            nInside = nCutInside;
            return;
        }
    }

    {
        // face3
        plane pl3(tetB[0], tetB[1], tetB[2]);

        nInside = 0;

        for (label i = 0; i < nCutInside; i++)
        {
            tetSliceWithPlane(cutInsideTets[i], pl3, inside, outside);
        }
    }
}


inline Foam::scalar
Foam::tetOverlapVolume::tetTetOverlapVol
(
    const tetPoints& tetA,
    const tetPoints& tetB
)
{
    FixedList<tetPoints, 200> insideTets;
    label nInside = 0;
    FixedList<tetPoints, 200> cutInsideTets;
    label nCutInside = 0;

    storeTetOp inside(insideTets, nInside);
    storeTetOp cutInside(cutInsideTets, nCutInside);
    sumTetVolOp volInside;
    dummyTetOp outside;

    // face0
    plane pl0(tetB[1], tetB[3], tetB[2]);
    tetA.tet().sliceWithPlane(pl0, cutInside, outside);
    if (nCutInside == 0)
    {
        return 0.0;
    }

    // face1
    plane pl1(tetB[0], tetB[2], tetB[3]);
    nInside = 0;
    for (label i = 0; i < nCutInside; i++)
    {
        const tetPointRef t = cutInsideTets[i].tet();
        t.sliceWithPlane(pl1, inside, outside);
    }
    if (nInside == 0)
    {
        return 0.0;
    }

    // face2
    plane pl2(tetB[0], tetB[3], tetB[1]);
    nCutInside = 0;
    for (label i = 0; i < nInside; i++)
    {
        const tetPointRef t = insideTets[i].tet();
        t.sliceWithPlane(pl2, cutInside, outside);
    }
    if (nCutInside == 0)
    {
        return 0.0;
    }

    // face3
    plane pl3(tetB[0], tetB[1], tetB[2]);
    for (label i = 0; i < nCutInside; i++)
    {
        const tetPointRef t = cutInsideTets[i].tet();
        t.sliceWithPlane(pl3, volInside, outside);
    }

    return volInside.vol_;
}


inline const Foam::treeBoundBox Foam::tetOverlapVolume::pyrBb
(
    const pointField& points,
    const face& f,
    const point& fc
)
{
    treeBoundBox bb(fc, fc);
    forAll(f, fp)
    {
        const point& pt = points[f[fp]];
        bb.min() = min(bb.min(), pt);
        bb.max() = max(bb.max(), pt);
    }
    return bb;
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::tetOverlapVolume::cellCellOverlapVolumeMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,

    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB
)
{
    const cell& cFacesA = meshA.cells()[cellAI];
    const point& ccA = meshA.cellCentres()[cellAI];

    const cell& cFacesB = meshB.cells()[cellBI];
    const point& ccB = meshB.cellCentres()[cellBI];

    scalar vol = 0.0;

    forAll(cFacesA, cFA)
    {
        label faceAI = cFacesA[cFA];

        const face& fA = meshA.faces()[faceAI];
        const treeBoundBox pyrA = pyrBb(meshA.points(), fA, ccA);
        if (!pyrA.overlaps(cellBbB))
        {
            continue;
        }

        bool ownA = (meshA.faceOwner()[faceAI] == cellAI);

        label tetBasePtAI = 0;

        const point& tetBasePtA = meshA.points()[fA[tetBasePtAI]];

        for (label tetPtI = 1; tetPtI < fA.size() - 1; tetPtI++)
        {
            label facePtAI = (tetPtI + tetBasePtAI) % fA.size();
            label otherFacePtAI = fA.fcIndex(facePtAI);

            label pt0I = -1;
            label pt1I = -1;

            if (ownA)
            {
                pt0I = fA[facePtAI];
                pt1I = fA[otherFacePtAI];
            }
            else
            {
                pt0I = fA[otherFacePtAI];
                pt1I = fA[facePtAI];
            }

            const tetPoints tetA
            (
                ccA,
                tetBasePtA,
                meshA.points()[pt0I],
                meshA.points()[pt1I]
            );
            const treeBoundBox tetABb(tetA.bounds());


            // Loop over tets of cellB
            forAll(cFacesB, cFB)
            {
                label faceBI = cFacesB[cFB];

                const face& fB = meshB.faces()[faceBI];
                const treeBoundBox pyrB = pyrBb(meshB.points(), fB, ccB);
                if (!pyrB.overlaps(pyrA))
                {
                    continue;
                }

                bool ownB = (meshB.faceOwner()[faceBI] == cellBI);

                label tetBasePtBI = 0;

                const point& tetBasePtB = meshB.points()[fB[tetBasePtBI]];

                for (label tetPtI = 1; tetPtI < fB.size() - 1; tetPtI++)
                {
                    label facePtBI = (tetPtI + tetBasePtBI) % fB.size();
                    label otherFacePtBI = fB.fcIndex(facePtBI);

                    label pt0I = -1;
                    label pt1I = -1;

                    if (ownB)
                    {
                        pt0I = fB[facePtBI];
                        pt1I = fB[otherFacePtBI];
                    }
                    else
                    {
                        pt0I = fB[otherFacePtBI];
                        pt1I = fB[facePtBI];
                    }

                    const tetPoints tetB
                    (
                        ccB,
                        tetBasePtB,
                        meshB.points()[pt0I],
                        meshB.points()[pt1I]
                    );
                    if (!tetB.bounds().overlaps(tetABb))
                    {
                        continue;
                    }

                    vol += tetTetOverlapVol(tetA, tetB);
                }
            }
        }
    }

    return vol;
}


Foam::labelList Foam::tetOverlapVolume::overlappingCells
(
    const fvMesh& fromMesh,
    const fvMesh& toMesh,
    const label iTo
) const
{
    const indexedOctree<treeDataCell>& treeA = fromMesh.cellTree();

    treeBoundBox bbB
    (
        pointField(toMesh.points(), toMesh.cellPoints()[iTo])
    );

    return treeA.findBox(bbB);
}


/*

        forAll(cellsA, i)
        {
            label cellAI = cellsA[i];
            treeBoundBox bbA
            (
                pointField(meshA.points(), meshA.cellPoints()[cellAI])
            );

            scalar v = cellCellOverlapVolumeMinDecomp
            (
                meshA,
                cellAI,
                bbA,
                meshB,
                cellBI,
                bbB
            );

            overlapVol += v;
            nOverlapTests++;
        }

*/

// ************************************************************************* //
