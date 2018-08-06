/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "cuttingPlane.H"
#include "fvMesh.H"
#include "volFields.H"
#include "linePointRef.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::cuttingPlane::debug(Foam::debug::debugSwitch("cuttingPlane", 0));


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Check edge/plane intersection based on crossings ... trivial check.
    inline bool intersectsEdge(const PackedList<2>& sides, const edge& e)
    {
        return (sides[e.first()] != sides[e.last()]);
    }


    // Check for face/plane intersection based on crossings
    // Took (-1,0,+1) from plane::sign and packed as (0,1,2).
    // Now use for left shift to obtain (1,2,4).
    //
    // Test accumulated value for an intersection with the plane.
    inline bool intersectsFace
    (
        const PackedList<2>& sides,
        const labelUList& indices
    )
    {
        unsigned accum = 0u;

        for (const label pointi : indices)
        {
            accum |= (1u << sides[pointi]);
        }

        // Accumulated value 3,5,6,7 indicates an intersection
        return (accum == 3 || accum >= 5);
    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::PackedList<2> Foam::cuttingPlane::classifySides
(
    const plane& pln,
    const pointField& pts
)
{
    const label len = pts.size();

    PackedList<2> output(len);

    // From signed (-1,0,+1) to (0,1,2) for PackedList
    for (label i=0; i < len; ++i)
    {
        output.set(i, unsigned(1 + pln.sign(pts[i], SMALL)));
    }

    return output;
}


Foam::label Foam::cuttingPlane::calcCellCuts
(
    const primitiveMesh& mesh,
    const PackedList<2>& sides,
    bitSet& cellCuts
)
{
    const faceList& faces = mesh.faces();

    const label nCells = mesh.nCells();
    const label nFaces = mesh.nFaces();
    const label nInternFaces = mesh.nInternalFaces();

    // Detect cells cuts from the face cuts

    label nFaceCuts = 0;

    // 1st face cut of cell
    bitSet hasCut1(nCells);

    // 2nd face cut of cell
    bitSet hasCut2(nCells);

    for (label facei = 0; facei < nInternFaces; ++facei)
    {
        if (intersectsFace(sides, faces[facei]))
        {
            const label own = mesh.faceOwner()[facei];
            const label nei = mesh.faceNeighbour()[facei];

            ++nFaceCuts;

            if (!hasCut1.set(own))
            {
                hasCut2.set(own);
            }
            if (!hasCut1.set(nei))
            {
                hasCut2.set(nei);
            }
        }
    }

    for (label facei = nInternFaces; facei < nFaces; ++facei)
    {
        if (intersectsFace(sides, faces[facei]))
        {
            const label own = mesh.faceOwner()[facei];

            ++nFaceCuts;

            if (!hasCut1.set(own))
            {
                hasCut2.set(own);
            }
        }
    }

    hasCut1.clearStorage();   // Not needed now

    if (cellCuts.size())
    {
        cellCuts.resize(nCells);    // safety
        cellCuts &= hasCut2;        // restrict to cell subset

        if (debug)
        {
            Pout<<"detected " << cellCuts.count() << "/" << nCells
                << " cells cut, subsetted from "
                << hasCut2.count() << "/" << nCells << " cells." << endl;
        }
    }
    else
    {
        cellCuts = std::move(hasCut2);

        if (debug)
        {
            Pout<<"detected " << cellCuts.count() << "/" << nCells
                << " cells cut." << endl;
        }
    }


    if (debug && isA<fvMesh>(mesh))
    {
        const fvMesh& fvm = dynamicCast<const fvMesh&>(mesh);

        volScalarField debugField
        (
            IOobject
            (
                "cuttingPlane.cellCuts",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvm,
            dimensionedScalar(dimless, Zero)
        );

        auto& debugFld = debugField.primitiveFieldRef();

        for (const label celli : cellCuts)
        {
            debugFld[celli] = 1;
        }

        Pout<< "Writing cut types:"
            << debugField.objectPath() << endl;

        debugField.write();
    }


    return nFaceCuts;
}


void Foam::cuttingPlane::intersectEdges
(
    const primitiveMesh& mesh,
    const PackedList<2>& sides,
    const bitSet& cellCuts,
    List<label>& edgePoint
)
{
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    // Per edge -1 or the label of the intersection point
    edgePoint.resize(edges.size());

    DynamicList<point> dynCutPoints(4*cellCuts.count());

    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];

        if (intersectsEdge(sides, points, e))
        {
            // Edge is cut
            edgePoint[edgei] = dynCutPoints.size();

            const point& p0 = points[e[0]];
            const point& p1 = points[e[1]];

            const scalar alpha = lineIntersect(linePointRef(p0, p1));

            if (alpha < SMALL)
            {
                dynCutPoints.append(p0);
            }
            else if (alpha >= 1.0)
            {
                dynCutPoints.append(p1);
            }
            else
            {
                dynCutPoints.append((1-alpha)*p0 + alpha*p1);
            }
        }
        else
        {
            edgePoint[edgei] = -1;
        }
    }

    this->storedPoints().transfer(dynCutPoints);
}


bool Foam::cuttingPlane::walkCell
(
    const primitiveMesh& mesh,
    const labelUList& edgePoint,
    const label celli,
    const label startEdgei,
    DynamicList<label>& faceVerts
)
{
    label facei = -1;
    label edgei = startEdgei;

    label nIter = 0;

    faceVerts.clear();
    do
    {
        faceVerts.append(edgePoint[edgei]);

        // Cross edge to other face
        facei = meshTools::otherFace(mesh, celli, facei, edgei);

        // Find next cut edge on face.
        const labelList& fEdges = mesh.faceEdges()[facei];

        label nextEdgei = -1;

        //Note: here is where we should check for whether there are more
        // than 2 intersections with the face (warped/non-convex face).
        // If so should e.g. decompose the cells on both faces and redo
        // the calculation.

        for (const label edge2i : fEdges)
        {
            if (edge2i != edgei && edgePoint[edge2i] != -1)
            {
                nextEdgei = edge2i;
                break;
            }
        }

        if (nextEdgei == -1)
        {
            // Did not find another cut edge on facei. Do what?
            WarningInFunction
                << "Did not find closed walk along surface of cell " << celli
                << " starting from edge " << startEdgei
                << " in " << nIter << " iterations." << nl
                << "Collected cutPoints so far:" << faceVerts
                << endl;

            return false;
        }

        edgei = nextEdgei;

        ++nIter;

        if (nIter > 1000)
        {
            WarningInFunction
                << "Did not find closed walk along surface of cell " << celli
                << " at " << mesh.cellCentres()[celli]
                << " starting from edge " << startEdgei
                << " in " << nIter << " iterations."
                << endl;
            return false;
        }

    } while (edgei != startEdgei);


    if (faceVerts.size() >= 3)
    {
        return true;
    }

    WarningInFunction
        << "Did not find closed walk along surface of cell " << celli
        << " starting from edge " << startEdgei << nl
        << "Collected cutPoints so far:" << faceVerts
        << endl;

    return false;
}


void Foam::cuttingPlane::walkCellCuts
(
    const primitiveMesh& mesh,
    const bitSet& cellCuts,
    const labelUList& edgePoint,
    const bool triangulate
)
{
    const pointField& cutPoints = this->points();

    // Dynamic lists to handle triangulation and/or missed cuts

    DynamicList<face>  dynCutFaces(cellCuts.count());
    DynamicList<label> dynCutCells(cellCuts.count());

    // Scratch space for calculating the face vertices
    DynamicList<label> faceVerts(16);

    for (const label celli : cellCuts)
    {
        // Find the starting edge to walk from.
        const labelList& cEdges = mesh.cellEdges()[celli];

        label startEdgei = -1;

        for (const label edgei : cEdges)
        {
            if (edgePoint[edgei] != -1)
            {
                startEdgei = edgei;
                break;
            }
        }

        // Check for the unexpected ...
        if (startEdgei == -1)
        {
            FatalErrorInFunction
                << "Cannot find cut edge for cut cell " << celli
                << abort(FatalError);
        }

        // Walk from starting edge around the circumference of the cell.
        bool okCut = walkCell
        (
            mesh,
            edgePoint,
            celli,
            startEdgei,
            faceVerts
        );

        if (okCut)
        {
            face f(faceVerts);

            // Orient face to point in the same direction as the plane normal
            if ((f.normal(cutPoints) & normal()) < 0)
            {
                f.flip();
            }

            // The cut faces can be quite ugly, so optionally triangulate
            if (triangulate)
            {
                label nTri = f.triangles(cutPoints, dynCutFaces);
                while (nTri--)
                {
                    dynCutCells.append(celli);
                }
            }
            else
            {
                dynCutFaces.append(f);
                dynCutCells.append(celli);
            }
        }
    }

    // No cuts? Then no need for any of this information
    if (dynCutCells.empty())
    {
        this->storedPoints().clear();
        this->storedFaces().clear();
        meshCells_.clear();
    }
    else
    {
        this->storedFaces().transfer(dynCutFaces);
        meshCells_.transfer(dynCutCells);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cuttingPlane::cuttingPlane(const plane& pln)
:
    plane(pln)
{}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    const bitSet& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    bitSet&& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    const labelUList& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cuttingPlane::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    bitSet&& cellIdLabels
)
{
    MeshStorage::clear();
    meshCells_.clear();

    // Pre-populate with restriction
    bitSet cellCuts(std::move(cellIdLabels));

    if (cellCuts.size())
    {
        cellCuts.resize(mesh.nCells());
    }

    // Classification of points with respect to the plane
    PackedList<2> sides = classifySides(*this, mesh.points());

    // Determine cells that are (likely) cut
    // - some ambiguity when plane is exactly between cells
    calcCellCuts(mesh, sides, cellCuts);

    // Determine cutPoints and return list of edge cuts.
    // per edge -1 or the label of the intersection point
    labelList edgePoint;
    intersectEdges(mesh, sides, cellCuts, edgePoint);

    // Do topological walk around cell to find closed loop.
    walkCellCuts(mesh, cellCuts, edgePoint, triangulate);
}


void Foam::cuttingPlane::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    const bitSet& cellIdLabels
)
{
    bitSet subsetCells(cellIdLabels);

    performCut(mesh, triangulate, std::move(subsetCells));
}


void Foam::cuttingPlane::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    const labelUList& cellIdLabels
)
{
    bitSet subsetCells;

    if (notNull(cellIdLabels))
    {
        // Pre-populate with restriction
        subsetCells.resize(mesh.nCells());
        subsetCells.set(cellIdLabels);
    }

    performCut(mesh, triangulate, std::move(subsetCells));
}


void Foam::cuttingPlane::remapFaces
(
    const labelUList& faceMap
)
{
    if (notNull(faceMap) && !faceMap.empty())
    {
        MeshStorage::remapFaces(faceMap);

        List<label> newCutCells(faceMap.size());
        forAll(faceMap, facei)
        {
            newCutCells[facei] = meshCells_[faceMap[facei]];
        }
        meshCells_.transfer(newCutCells);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cuttingPlane::operator=(const cuttingPlane& rhs)
{
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    static_cast<MeshStorage&>(*this) = rhs;
    static_cast<plane&>(*this) = rhs;
    meshCells_ = rhs.meshCells();
}


// ************************************************************************* //
