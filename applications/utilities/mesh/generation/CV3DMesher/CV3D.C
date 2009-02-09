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

#include "CV3D.H"
#include "Random.H"
#include "uint.H"
#include "ulong.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV3D::insertBoundingBox()
{
    Info<< "insertBoundingBox: creating bounding mesh" << nl << endl;
    scalar bigSpan = 10*tols_.span;
    insertPoint(point(-bigSpan, -bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan, -bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan,  bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan,  bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan, -bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan, -bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan,  bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan,  bigSpan , bigSpan), Vb::FAR_POINT);
}


void Foam::CV3D::reinsertPoints(const pointField& points)
{
    Info<< nl << "Reinserting points after motion. ";

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    // Add the points and index them
    forAll(points, i)
    {
        const point& p = points[i];

        insert(toPoint(p))->index() = nVert++;
    }

    Info<< nVert << " vertices reinserted" << endl;
}


void Foam::CV3D::setVertexAlignmentDirections()
{
    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            List<vector>& alignmentDirections(vit->alignmentDirections());

            point vert(topoint(vit->point()));

            pointIndexHit pHit = qSurf_.tree().findNearest
            (
                vert,
                controls_.nearWallAlignedDist2
            );

            if (pHit.hit())
            {
                // Primary alignment
                vector np = qSurf_.faceNormals()[pHit.index()];

                // Generate equally spaced 'spokes' in a circle normal to the
                // direction from the vertex to the closest point on the surface
                // and look for a secondary intersection.

                vector d = pHit.hitPoint() - vert;

                tensor R = rotationTensor(vector(0,0,1), np);

                label s = 36;

                scalar closestSpokeHitDistance = GREAT;

                point closestSpokeHitPoint = point(GREAT,GREAT,GREAT);

                label closestSpokeHitDistanceIndex = -1;

                for(label i = 0; i < s; i++)
                {
                    vector spoke
                    (
                        Foam::cos(i*mathematicalConstant::twoPi/s),
                        Foam::sin(i*mathematicalConstant::twoPi/s),
                        0
                    );

                    spoke *= controls_.nearWallAlignedDist;

                    spoke = R & spoke;

                    pointIndexHit spokeHit;

                    // internal spoke

                    spokeHit = qSurf_.tree().findLine
                    (
                        vert,
                        vert + spoke
                    );

                    if (spokeHit.hit())
                    {
                        scalar spokeHitDistance = mag
                        (
                            spokeHit.hitPoint() - vert
                        );

                        if (spokeHitDistance < closestSpokeHitDistance)
                        {
                            closestSpokeHitPoint =  spokeHit.hitPoint();

                            closestSpokeHitDistanceIndex = spokeHit.index();
                        }
                    }

                    //external spoke

                    point mirrorVert = vert + 2*d;

                    spokeHit = qSurf_.tree().findLine
                    (
                        mirrorVert,
                        mirrorVert + spoke
                    );

                    if (spokeHit.hit())
                    {
                        scalar spokeHitDistance = mag
                        (
                            spokeHit.hitPoint() - mirrorVert
                        );

                        if (spokeHitDistance < closestSpokeHitDistance)
                        {
                            closestSpokeHitPoint =  spokeHit.hitPoint();

                            closestSpokeHitDistanceIndex = spokeHit.index();
                        }
                    }
                }

                if (closestSpokeHitDistanceIndex > -1)
                {
                    // Auxiliary alignment generated by spoke intersection
                    // normal.

                    vector na =
                    qSurf_.faceNormals()[closestSpokeHitDistanceIndex];

                    // Secondary alignment
                    vector ns = np ^ na;

                    if (mag(ns) <  SMALL)
                    {
                        FatalErrorIn("Foam::CV3D::setVertexAlignmentDirections")
                            << "Parallel normals detected in spoke search."
                            << nl << exit(FatalError);
                    }

                    ns /= mag(ns);

                    // Tertiary alignment
                    vector nt = ns ^ np;

                    // this normalisation is not necessary if np and ns are
                    // perpendicular unit vectors.

                    nt /= mag(nt);

                    alignmentDirections.setSize(3);

                    alignmentDirections[0] = np;

                    alignmentDirections[1] = ns;

                    alignmentDirections[2] = nt;

                    // Info<< "internal " << vit->internalPoint()
                    //     << nl << alignmentDirections
                    //     << nl << "v " << vert + alignmentDirections[0]
                    //     << nl << "v " << vert + alignmentDirections[1]
                    //     << nl << "v " << vert + alignmentDirections[2]
                    //     << nl << "v " << vert
                    //     << nl << "v " << pHit.hitPoint()
                    //     << nl << "v " << closestSpokeHitPoint
                    //     << nl << "f 4 1"
                    //     << nl << "f 4 2"
                    //     << nl << "f 4 3"
                    //     << nl << endl;
                }
                else
                {
                    // Using only the primary alignment

                    alignmentDirections.setSize(1);

                    alignmentDirections[0] = np;
                }
            }
            else
            {
                alignmentDirections.setSize(0);
            }
        }
    }
}


Foam::scalar Foam::CV3D::alignmentDistanceWeight(scalar dist) const
{
    scalar w;

    scalar x = dist/controls_.nearWallAlignedDist;

    if (x < 0.5)
    {
        w = 1;
    }

    else if (x < 1)
    {
        w = 2*(1 - x);
    }
    else
    {
        w = 0;
    }

    return w;
}


Foam::scalar Foam::CV3D::faceAreaWeight(scalar faceArea) const
{
    scalar fl2 = 0.5;

    scalar fu2 = 1.0;

    scalar m2 = controls_.minCellSize2;

    if (faceArea < fl2*m2)
    {
        return 0;
    }
    else if (faceArea < fu2*m2)
    {
        return faceArea/((fu2 - fl2)*m2) - 1/((fu2/fl2) - 1);
    }
    else
    {
        return 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CV3D::CV3D
(
    const Time& runTime,
    const querySurface& qSurf
)
:
    HTriangulation(),
    qSurf_(qSurf),
    runTime_(runTime),
    controls_
    (
        IOdictionary
        (
            IOobject
            (
                "CV3DMesherDict",
                runTime_.system(),
                runTime_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    ),
    tols_
    (
        IOdictionary
        (
            IOobject
            (
                "CV3DMesherDict",
                runTime_.system(),
                runTime_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ),
        controls_.minCellSize,
        qSurf.bb()
    ),
    startOfInternalPoints_(0),
    startOfSurfacePointPairs_(0),
    featureConstrainingVertices_(0)
{
    // insertBoundingBox();
    insertFeaturePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CV3D::~CV3D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CV3D::insertPoints
(
    const pointField& points,
    const scalar nearness
)
{
    Info<< "insertInitialPoints(const pointField& points): ";

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    // Add the points and index them
    forAll(points, i)
    {
        const point& p = points[i];

        if (qSurf_.wellInside(p, nearness))
        {
            insert(toPoint(p))->index() = nVert++;
        }
        else
        {
            Warning
                << "Rejecting point " << p << " outside surface" << endl;
        }
    }

    Info<< nVert << " vertices inserted" << endl;

    if (controls_.writeInitialTriangulation)
    {
        // Checking validity of triangulation
        assert(is_valid());

        writeTriangles("initial_triangles.obj", true);
    }
}


void Foam::CV3D::insertPoints(const fileName& pointFileName)
{
    pointIOField points
    (
        IOobject
        (
            pointFileName.name(),
            pointFileName.path(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    insertPoints(points, 0.5*controls_.minCellSize2);
}


void Foam::CV3D::insertGrid()
{
    Info<< nl << "Inserting initial grid." << endl;

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    Info<< nl << nVert << " existing vertices." << endl;

    scalar x0 = qSurf_.bb().min().x();
    scalar xR = qSurf_.bb().max().x() - x0;
    int ni = int(xR/controls_.minCellSize) + 1;

    scalar y0 = qSurf_.bb().min().y();
    scalar yR = qSurf_.bb().max().y() - y0;
    int nj = int(yR/controls_.minCellSize) + 1;

    scalar z0 = qSurf_.bb().min().z();
    scalar zR = qSurf_.bb().max().z() - z0;
    int nk = int(zR/controls_.minCellSize) + 1;

    vector delta(xR/ni, yR/nj, zR/nk);

    delta *= pow((1.0),-(1.0/3.0));

    Random rndGen(1321);
    scalar pert = controls_.randomPerturbation*cmptMin(delta);

    std::vector<Point> initialPoints;

    for (int i=0; i<ni; i++)
    {
        for (int j=0; j<nj; j++)
        {
            for (int k=0; k<nk; k++)
            {
                point p
                (
                    x0 + i*delta.x(),
                    y0 + j*delta.y(),
                    z0 + k*delta.z()
                );

                if (controls_.randomiseInitialGrid)
                {
                    p.x() += pert*(rndGen.scalar01() - 0.5);
                    p.y() += pert*(rndGen.scalar01() - 0.5);
                    p.z() += pert*(rndGen.scalar01() - 0.5);
                }

                if (qSurf_.wellInside(p, 0.5*controls_.minCellSize2))
                {
                    initialPoints.push_back(Point(p.x(), p.y(), p.z()));
                    // insert(Point(p.x(), p.y(), p.z()))->index() = nVert++;
                }
            }
        }
    }

    Info<< nl << initialPoints.size() << " vertices to insert." << endl;

    // using the range insert (it is faster than inserting points one by one)
    insert(initialPoints.begin(), initialPoints.end());

    Info<< nl << number_of_vertices() - startOfInternalPoints_
        << " vertices inserted." << endl;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->uninitialised())
        {
            vit->index() = nVert++;
        }
    }

    if (controls_.writeInitialTriangulation)
    {
        assert(is_valid());

        writePoints("initial_points.obj", true);
        writeTriangles("initial_triangles.obj", true);
    }
}


void Foam::CV3D::relaxPoints(const scalar relaxation)
{
    Info<< "Calculating new points: " << endl;

    pointField dualVertices(number_of_cells());

    pointField displacementAccumulator(startOfSurfacePointPairs_, vector::zero);

    scalarField weightAccumulator(startOfSurfacePointPairs_, 0);

    label dualVerti = 0;

    // Find the dual point of each tetrahedron and assign it an index.

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        cit->cellIndex() = -1;

        if
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
        {
            cit->cellIndex() = dualVerti;

            // To output Delaunay tet which causes CGAL assertion failure.
            // Info<< nl << topoint(cit->vertex(0)->point())
            //     << nl << topoint(cit->vertex(1)->point())
            //     << nl << topoint(cit->vertex(2)->point())
            //     << nl << topoint(cit->vertex(3)->point())
            //     << endl;

            dualVertices[dualVerti] = topoint(dual(cit));

            dualVerti++;
        }
    }

    setVertexAlignmentDirections();

    dualVertices.setSize(dualVerti);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop around the Delaunay edges to construct the dual faces.
    // Find the face-centre and use it to calculate the displacement vector
    // contribution to the Delaunay vertices (Dv) attached to the edge.  Add the
    // contribution to the running displacement vector of each Dv.

    // for
    // (
    //     Triangulation::Finite_edges_iterator eit = finite_edges_begin();
    //     eit != finite_edges_end();
    //     ++eit
    // )
    // {
    //     if
    //     (
    //         eit->first->vertex(eit->second)->internalOrBoundaryPoint()
    //      && eit->first->vertex(eit->third)->internalOrBoundaryPoint()
    //     )
    //     {
    //         Cell_circulator ccStart = incident_cells(*eit);
    //         Cell_circulator cc = ccStart;

    //         DynamicList<label> verticesOnFace;

    //         do
    //         {
    //             if (!is_infinite(cc))
    //             {
    //                 if (cc->cellIndex() < 0)
    //                 {
    //                     FatalErrorIn("Foam::CV3D::relaxPoints")
    //                         << "Dual face uses circumcenter defined by a "
    //                         << " Delaunay tetrahedron with no internal "
    //                         << "or boundary points."
    //                         << exit(FatalError);
    //                 }

    //                 verticesOnFace.append(cc->cellIndex());
    //             }
    //         } while (++cc != ccStart);

    //         verticesOnFace.shrink();

    //         face dualFace(verticesOnFace);

    //         Cell_handle c = eit->first;
    //         Vertex_handle vA = c->vertex(eit->second);
    //         Vertex_handle vB = c->vertex(eit->third);

    //         point dVA = topoint(vA->point());
    //         point dVB = topoint(vB->point());

    //         point dualFaceCentre(dualFace.centre(dualVertices));

    //         vector rAB = dVA - dVB;

    //         scalar rABMag = mag(rAB);

    //         scalar faceArea = dualFace.mag(dualVertices);

    //         scalar directStiffness = 2.0;

    //         scalar transverseStiffness = 0.0001;

    //         scalar r0 = 0.9*controls_.minCellSize;

    //         vector dA = -directStiffness*(1 - r0/rABMag)
    //         *faceAreaWeight(faceArea)*rAB;

    //         vector dT = transverseStiffness*faceAreaWeight(faceArea)
    //         *(dualFaceCentre - 0.5*(dVA - dVB));

    //         if (vA->internalPoint())
    //         {
    //             displacementAccumulator[vA->index()] += (dA + dT);
    //         }
    //         if (vB->internalPoint())
    //         {
    //             displacementAccumulator[vB->index()] += (-dA + dT);
    //         }
    //     }
    // }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Rotate faces that are sufficiently large and well enough aligned with the
    // cell alignment direction(s)

    // vector n2 = vector(1,1,1);

    // n2 /= mag(n2);

    // tensor R = rotationTensor(vector(1,0,0), n2);

    // List<vector> globalAlignmentDirs(3);

    // globalAlignmentDirs[0] = R & vector(1,0,0);
    // globalAlignmentDirs[1] = R & vector(0,1,0);
    // globalAlignmentDirs[2] = R & vector(0,0,1);

    // Info<< "globalAlignmentDirs " << globalAlignmentDirs << endl;

    // for
    // (
    //     Triangulation::Finite_edges_iterator eit = finite_edges_begin();
    //     eit != finite_edges_end();
    //     ++eit
    // )
    // {
    //     if
    //     (
    //         eit->first->vertex(eit->second)->internalOrBoundaryPoint()
    //      && eit->first->vertex(eit->third)->internalOrBoundaryPoint()
    //     )
    //     {
    //         Cell_circulator ccStart = incident_cells(*eit);
    //         Cell_circulator cc = ccStart;

    //         DynamicList<label> verticesOnFace;

    //         do
    //         {
    //             if (!is_infinite(cc))
    //             {
    //                 if (cc->cellIndex() < 0)
    //                 {
    //                     FatalErrorIn("Foam::CV3D::relaxPoints")
    //                     << "Dual face uses circumcenter defined by a "
    //                         << " Delaunay tetrahedron with no internal "
    //                         << "or boundary points."
    //                         << exit(FatalError);
    //                 }

    //                 verticesOnFace.append(cc->cellIndex());
    //             }
    //         } while (++cc != ccStart);

    //         verticesOnFace.shrink();

    //         face dualFace(verticesOnFace);

    //         Cell_handle c = eit->first;
    //         Vertex_handle vA = c->vertex(eit->second);
    //         Vertex_handle vB = c->vertex(eit->third);

    //         point dVA = topoint(vA->point());
    //         point dVB = topoint(vB->point());

    //         Field<vector> alignmentDirs;

    //         if
    //         (
    //             vA->alignmentDirections().size() > 0
    //          || vB->alignmentDirections().size() > 0
    //         )
    //         {
    //             if
    //             (
    //                 vA->alignmentDirections().size() == 3
    //              || vB->alignmentDirections().size() == 3
    //             )
    //             {
    //                 alignmentDirs.setSize(3);

    //                 if (vB->alignmentDirections().size() == 0)
    //                 {
    //                     alignmentDirs = vA->alignmentDirections();
    //                 }
    //                 else if (vA->alignmentDirections().size() == 0)
    //                 {
    //                     alignmentDirs = vB->alignmentDirections();
    //                 }
    //                 else if
    //                 (
    //                     vA->alignmentDirections().size() == 3
    //                  && vB->alignmentDirections().size() == 3
    //                 )
    //                 {
    //                     forAll(vA->alignmentDirections(), aA)
    //                     {
    //                         const vector& a(vA->alignmentDirections()[aA]);

    //                         scalar maxDotProduct = 0.0;

    //                         forAll(vB->alignmentDirections(), aB)
    //                         {
    //                             const vector& b(vB->alignmentDirections()[aB]);

    //                             if (mag(a & b) > maxDotProduct)
    //                             {
    //                                 maxDotProduct = mag(a & b);

    //                                 alignmentDirs[aA] = a + sign(a & b)*b;

    //                                 alignmentDirs[aA] /= mag(alignmentDirs[aA]);
    //                             }
    //                         }
    //                     }
    //                 }
    //                 else
    //                 {
    //                     // One of the vertices has 3 alignments and the other
    //                     // has 1

    //                     vector otherAlignment;

    //                     if (vA->alignmentDirections().size() == 3)
    //                     {
    //                         alignmentDirs = vA->alignmentDirections();

    //                         otherAlignment = vB->alignmentDirections()[0];
    //                     }
    //                     else
    //                     {
    //                         alignmentDirs = vB->alignmentDirections();

    //                         otherAlignment = vA->alignmentDirections()[0];
    //                     }

    //                     label matchingDirection = -1;

    //                     scalar maxDotProduct = 0.0;

    //                     forAll(alignmentDirs, aD)
    //                     {
    //                         const vector& a(alignmentDirs[aD]);

    //                         if (mag(a & otherAlignment) > maxDotProduct)
    //                         {
    //                             maxDotProduct = mag(a & otherAlignment);

    //                             matchingDirection = aD;
    //                         }
    //                     }

    //                     vector& matchingAlignment =
    //                     alignmentDirs[matchingDirection];

    //                     matchingAlignment = matchingAlignment
    //                     + sign(matchingAlignment & otherAlignment)
    //                     *otherAlignment;

    //                     matchingAlignment /= mag(matchingAlignment);
    //                 }

    //                 // vector midpoint = 0.5*(dVA + dVB);

    //                 // Info<< "midpoint " << midpoint
    //                 //     << nl << alignmentDirs
    //                 //     << nl << "v " << midpoint + alignmentDirs[0]
    //                 //     << nl << "v " << midpoint + alignmentDirs[1]
    //                 //     << nl << "v " << midpoint + alignmentDirs[2]
    //                 //     << nl << "v " << midpoint
    //                 //     << nl << "f 4 1"
    //                 //     << nl << "f 4 2"
    //                 //     << nl << "f 4 3"
    //                 //     << nl << endl;
    //             }
    //             else
    //             {
    //                 alignmentDirs.setSize(1);

    //                 vector& alignmentDir = alignmentDirs[0];

    //                 if
    //                 (
    //                     vA->alignmentDirections().size() > 0
    //                  && vB->alignmentDirections().size() == 0
    //                 )
    //                 {
    //                     alignmentDir = vA->alignmentDirections()[0];
    //                 }
    //                 else if
    //                 (
    //                     vA->alignmentDirections().size() == 0
    //                  && vB->alignmentDirections().size() > 0
    //                 )
    //                 {
    //                     alignmentDir = vB->alignmentDirections()[0];
    //                 }
    //                 else
    //                 {
    //                     // Both vertices have an alignment

    //                     const vector& a(vA->alignmentDirections()[0]);

    //                     const vector& b(vB->alignmentDirections()[0]);

    //                     if (mag(a & b) < 1 - SMALL)
    //                     {
    //                         alignmentDirs.setSize(3);

    //                         alignmentDirs[0] = a + b;

    //                         alignmentDirs[1] = a - b;

    //                         alignmentDirs[2] =
    //                         alignmentDirs[0] ^ alignmentDirs[1];

    //                         alignmentDirs /= mag(alignmentDirs);
    //                     }
    //                     else
    //                     {
    //                         alignmentDir = a + sign(a & b)*b;

    //                         alignmentDir /= mag(alignmentDir);
    //                     }
    //                 }
    //                 if (alignmentDirs.size() ==1)
    //                 {
    //                     // Use the least similar of globalAlignmentDirs as the
    //                     // 2nd alignment and then generate the third.

    //                     scalar minDotProduct = 1+SMALL;

    //                     alignmentDirs.setSize(3);

    //                     forAll(alignmentDirs, aD)
    //                     {
    //                         if
    //                         (
    //                             mag(alignmentDir & globalAlignmentDirs[aD])
    //                           < minDotProduct
    //                         )
    //                         {
    //                             minDotProduct = mag
    //                             (
    //                                 alignmentDir & globalAlignmentDirs[aD]
    //                             );

    //                             alignmentDirs[1] = globalAlignmentDirs[aD];
    //                         }
    //                     }

    //                     alignmentDirs[2] = alignmentDirs[0] ^ alignmentDirs[1];

    //                     alignmentDirs[2] /= mag(alignmentDirs[2]);
    //                 }
    //             }
    //         }
    //         else
    //         {
    //             alignmentDirs = globalAlignmentDirs;
    //         }

    //         // alignmentDirs found, now apply them

    //         vector rAB = dVA - dVB;

    //         scalar rABMag = mag(rAB);

    //         forAll(alignmentDirs, aD)
    //         {
    //             vector& alignmentDir = alignmentDirs[aD];

    //             if ((rAB & alignmentDir) < 0)
    //             {
    //                 // swap the direction of the alignment so that has the
    //                 // same sense as rAB
    //                 alignmentDir *= -1;
    //             }

    //             //scalar cosAcceptanceAngle = 0.743;
    //             scalar cosAcceptanceAngle = 0.67;

    //             if (((rAB/rABMag) & alignmentDir) > cosAcceptanceAngle)
    //             {
    //                 alignmentDir *= 0.5*controls_.minCellSize;

    //                 vector delta = alignmentDir - 0.5*rAB;

    //                 scalar faceArea = dualFace.mag(dualVertices);

    //                 delta *= faceAreaWeight(faceArea);

    //                 if (vA->internalPoint())
    //                 {
    //                     displacementAccumulator[vA->index()] += delta;
    //                 }
    //                 if (vB->internalPoint())
    //                 {
    //                     displacementAccumulator[vB->index()] += -delta;
    //                 }
    //             }
    //         }
    //     }
    // }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Simple isotropic forcing using tension between the Delaunay vertex and
    // and the dual face centre

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        if
        (
            eit->first->vertex(eit->second)->internalOrBoundaryPoint()
         && eit->first->vertex(eit->third)->internalOrBoundaryPoint()
        )
        {
            Cell_circulator ccStart = incident_cells(*eit);
            Cell_circulator cc = ccStart;

            DynamicList<label> verticesOnFace;

            do
            {
                if (!is_infinite(cc))
                {
                    if (cc->cellIndex() < 0)
                    {
                        FatalErrorIn("Foam::CV3D::relaxPoints")
                            << "Dual face uses circumcenter defined by a "
                            << " Delaunay tetrahedron with no internal "
                            << "or boundary points."
                            << exit(FatalError);
                    }

                    verticesOnFace.append(cc->cellIndex());
                }
            } while (++cc != ccStart);

            verticesOnFace.shrink();

            face dualFace(verticesOnFace);

            Cell_handle c = eit->first;
            Vertex_handle vA = c->vertex(eit->second);
            Vertex_handle vB = c->vertex(eit->third);

            point dVA = topoint(vA->point());
            point dVB = topoint(vB->point());

            scalar faceArea = dualFace.mag(dualVertices);

            // point dualFaceCentre(dualFace.centre(dualVertices));
            point dEMid = 0.5*(dVA + dVB);

            if (vA->internalPoint())
            {
                displacementAccumulator[vA->index()] +=
                faceArea*(dEMid - dVA);
                //faceArea*(dualFaceCentre - dVA);

                weightAccumulator[vA->index()] += faceArea;
            }
            if (vB->internalPoint())
            {
                displacementAccumulator[vB->index()] +=
                faceArea*(dEMid - dVB);
                //faceArea*(dualFaceCentre - dVB);

                weightAccumulator[vB->index()] += faceArea;
            }
        }
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Cell based looping.

    // Loop over all Delaunay vertices (Dual cells)

    // for
    // (
    //     Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
    //     vit != finite_vertices_end();
    //     vit++
    // )
    // {
    //     if (vit->internalOrBoundaryPoint())
    //     {
    //         std::list<Vertex_handle> incidentVertices;

    //         incident_vertices(vit, std::back_inserter(incidentVertices));

    //         for
    //         (
    //             std::list<Vertex_handle>::iterator ivit =
    //                 incidentVertices.begin();
    //             ivit != incidentVertices.end();
    //             ++ivit
    //         )
    //         {
    //             // Pick up all incident vertices to the Dv - check the number to
    //             // see if this is a reasonable result

    //             // Form the edge from the current Dv and iterate around the
    //             // incident cells, finding their circumcentres, then form (and
    //             // store?) the dual face.

    //             // Wait!

    //             // There is no sensible way in CGAL-3.3.1 to turn the Delaunay
    //             // vertex in question and the incident vertex into an edge.
    //             // There is, however, an incident_edges function in CGAL-3.4
    //         }
    //     }
    // }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    displacementAccumulator /= weightAccumulator;

    vector totalDisp = sum(displacementAccumulator);
    scalar totalDist = sum(mag(displacementAccumulator));

    Info<< "Total displacement = " << totalDisp
        << nl << "Total distance = " << totalDist << endl;

    displacementAccumulator *= relaxation;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            displacementAccumulator[vit->index()] += topoint(vit->point());
        }
    }

    pointField internalDelaunayVertices = SubField<point>
    (
        displacementAccumulator,
        displacementAccumulator.size() - startOfInternalPoints_,
        startOfInternalPoints_
    );

    // Write the mesh before clearing it
    if (runTime_.outputTime())
    {
        writeMesh(true);
    }

    // Remove the entire triangulation
    this->clear();

    reinsertFeaturePoints();

    reinsertPoints(internalDelaunayVertices);
}


void Foam::CV3D::insertSurfacePointPairs()
{
    startOfSurfacePointPairs_ = number_of_vertices();

    if (controls_.insertSurfaceNearestPointPairs)
    {
        insertSurfaceNearestPointPairs();
    }

    if (controls_.writeNearestTriangulation)
    {
        // writeFaces("near_allFaces.obj", false);
        // writeFaces("near_faces.obj", true);
        writeTriangles("near_triangles.obj", true);
    }

    if (controls_.insertSurfaceNearPointPairs)
    {
        insertSurfaceNearPointPairs();
    }

    startOfBoundaryConformPointPairs_ = number_of_vertices();
}

void Foam::CV3D::boundaryConform()
{
}


void Foam::CV3D::removeSurfacePointPairs()
{
    Info<< "Removing surface point pairs." << nl << endl;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->index() >= startOfSurfacePointPairs_)
        {
            remove(vit);
        }
    }
}


void Foam::CV3D::write()
{
    if (controls_.writeFinalTriangulation)
    {
        writePoints("allPoints.obj", false);
        writePoints("points.obj", true);
        writeTriangles("allTriangles.obj", false);
        writeTriangles("triangles.obj", true);
    }

    writeMesh();
}

// ************************************************************************* //
