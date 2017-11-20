/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "triSurface.H"
#include "demandDrivenData.H"
#include "Time.H"
#include "surfZoneList.H"
#include "MeshedSurface.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(triSurface, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::fileName Foam::triSurface::triSurfInstance(const Time& d)
{
    fileName foamName(d.caseName() + ".ftr");

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = d.times();
    label i;

    for (i=ts.size()-1; i>=0; i--)
    {
        if (ts[i].value() <= d.timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (i>=0)
    {
        for (label j=i; j>=0; j--)
        {
            if (isFile(d.path()/ts[j].name()/typeName/foamName))
            {
                if (debug)
                {
                    Pout<< " triSurface::triSurfInstance(const Time& d)"
                        << "reading " << foamName
                        << " from " << ts[j].name()/typeName
                        << endl;
                }

                return ts[j].name();
            }
        }
    }

    if (debug)
    {
        Pout<< " triSurface::triSurfInstance(const Time& d)"
            << "reading " << foamName
            << " from constant/" << endl;
    }
    return d.constant();
}


Foam::List<Foam::labelledTri> Foam::triSurface::convertToTri
(
    const faceList& faces,
    const label defaultRegion
)
{
    List<labelledTri> triFaces(faces.size());

    forAll(triFaces, facei)
    {
        const face& f = faces[facei];

        if (f.size() != 3)
        {
            FatalErrorInFunction
                << "Face at position " << facei
                << " does not have three vertices:" << f
                << abort(FatalError);
        }

        labelledTri& tri = triFaces[facei];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = defaultRegion;
    }

    return triFaces;
}


Foam::List<Foam::labelledTri> Foam::triSurface::convertToTri
(
    const triFaceList& faces,
    const label defaultRegion
)
{
    List<labelledTri> triFaces(faces.size());

    forAll(triFaces, facei)
    {
        const triFace& f = faces[facei];

        labelledTri& tri = triFaces[facei];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = defaultRegion;
    }

    return triFaces;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Remove non-triangles, double triangles.
void Foam::triSurface::checkTriangles(const bool verbose)
{
    // Simple check on indices ok.
    const label maxPointi = points().size() - 1;

    for (const triSurface::FaceType& f : *this)
    {
        for (const label verti : f)
        {
            if (verti < 0 || verti > maxPointi)
            {
                FatalErrorInFunction
                    << "triangle " << f
                    << " uses point indices outside point range 0.."
                    << maxPointi
                    << exit(FatalError);
            }
        }
    }

    // Two phase process
    //   1. mark invalid faces
    //   2. pack
    // Done to keep numbering constant in phase 1

    // List of valid triangles
    boolList valid(size(), true);
    bool hasInvalid = false;

    forAll(*this, facei)
    {
        const labelledTri& f = (*this)[facei];

        if ((f[0] == f[1]) || (f[0] == f[2]) || (f[1] == f[2]))
        {
            // 'degenerate' triangle check
            valid[facei] = false;
            hasInvalid = true;

            if (verbose)
            {
                WarningInFunction
                    << "triangle " << facei
                    << " does not have three unique vertices:\n";
                printTriangle(Warning, "    ", f, points());
            }
        }
        else
        {
            // duplicate triangle check
            const labelList& fEdges = faceEdges()[facei];

            // Check if faceNeighbours use same points as this face.
            // Note: discards normal information - sides of baffle are merged.

            for (const label edgei : fEdges)
            {
                const labelList& eFaces = edgeFaces()[edgei];

                for (const label neighbour : eFaces)
                {
                    if (neighbour > facei)
                    {
                        // lower numbered faces already checked
                        const labelledTri& n = (*this)[neighbour];

                        if
                        (
                            ((f[0] == n[0]) || (f[0] == n[1]) || (f[0] == n[2]))
                         && ((f[1] == n[0]) || (f[1] == n[1]) || (f[1] == n[2]))
                         && ((f[2] == n[0]) || (f[2] == n[1]) || (f[2] == n[2]))
                        )
                        {
                            valid[facei] = false;
                            hasInvalid = true;

                            if (verbose)
                            {
                                WarningInFunction
                                    << "triangles share the same vertices:\n"
                                    << "    face 1 :" << facei << endl;
                                printTriangle(Warning, "    ", f, points());

                                Warning
                                    << endl
                                    << "    face 2 :"
                                    << neighbour << endl;
                                printTriangle(Warning, "    ", n, points());
                            }

                            break;
                        }
                    }
                }
            }
        }
    }

    if (hasInvalid)
    {
        // Pack
        label newFacei = 0;
        forAll(*this, facei)
        {
            if (valid[facei])
            {
                const labelledTri& f = (*this)[facei];
                (*this)[newFacei++] = f;
            }
        }

        if (verbose)
        {
            WarningInFunction
                << "Removing " << size() - newFacei
                << " illegal faces." << endl;
        }
        (*this).setSize(newFacei);

        // Topology can change because of renumbering
        clearOut();
    }
}


// Check/fix edges with more than two triangles
void Foam::triSurface::checkEdges(const bool verbose)
{
    const labelListList& eFaces = edgeFaces();

    forAll(eFaces, edgei)
    {
        const labelList& myFaces = eFaces[edgei];

        if (myFaces.empty())
        {
            FatalErrorInFunction
                << "Edge " << edgei << " with vertices " << edges()[edgei]
                << " has no edgeFaces"
                << exit(FatalError);
        }
        else if (myFaces.size() > 2 && verbose)
        {
            WarningInFunction
                << "Edge " << edgei << " with vertices " << edges()[edgei]
                << " has more than 2 faces connected to it : " << myFaces
                << endl;
        }
    }
}


// Returns patch info. Sets faceMap to the indexing according to patch
// numbers. Patch numbers start at 0.
Foam::surfacePatchList
Foam::triSurface::calcPatches(labelList& faceMap) const
{
    // Determine the sorted order:
    // use sortedOrder directly (the intermediate list is discared anyhow)

    List<label> regions(size());
    forAll(regions, facei)
    {
        regions[facei] = operator[](facei).region();
    }
    sortedOrder(regions, faceMap);
    regions.clear();

    // Extend regions
    label maxRegion = patches_.size()-1;    // for non-compacted regions

    if (faceMap.size())
    {
        maxRegion = max
        (
            maxRegion,
            operator[](faceMap.last()).region()
        );
    }

    // Get new region list
    surfacePatchList newPatches(maxRegion + 1);

    // Fill patch sizes
    forAll(*this, facei)
    {
        label region = operator[](facei).region();

        newPatches[region].size()++;
    }


    // Fill rest of patch info

    label startFacei = 0;
    forAll(newPatches, newPatchi)
    {
        surfacePatch& newPatch = newPatches[newPatchi];

        newPatch.index() = newPatchi;
        newPatch.start() = startFacei;

        // Take over any information from existing patches
        if
        (
            newPatchi < patches_.size()
         && !patches_[newPatchi].name().empty()
        )
        {
            newPatch.name() = patches_[newPatchi].name();
        }
        else
        {
            newPatch.name() = word("patch") + Foam::name(newPatchi);
        }

        if
        (
            newPatchi < patches_.size()
         && !patches_[newPatchi].geometricType().empty()
        )
        {
            newPatch.geometricType() = patches_[newPatchi].geometricType();
        }
        else
        {
            newPatch.geometricType() = geometricSurfacePatch::emptyType;
        }

        startFacei += newPatch.size();
    }

    return newPatches;
}


void Foam::triSurface::setDefaultPatches()
{
    labelList faceMap;

    // Get names, types and sizes
    surfacePatchList newPatches(calcPatches(faceMap));

    // Take over names and types (but not size)
    patches_.setSize(newPatches.size());

    forAll(newPatches, patchi)
    {
        patches_[patchi].index() = patchi;
        patches_[patchi].name() = newPatches[patchi].name();
        patches_[patchi].geometricType() = newPatches[patchi].geometricType();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurface::triSurface()
:
    ParentType(List<Face>(), pointField()),
    patches_(0),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface(const triSurface& surf)
:
    ParentType(surf, surf.points()),
    patches_(surf.patches()),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface(triSurface&& surf)
:
    triSurface()
{
    FaceListType::operator=(std::move(static_cast<FaceListType&>(surf)));
    storedPoints() = std::move(surf.storedPoints());
    patches_ = std::move(surf.patches());
}


Foam::triSurface::triSurface
(
    const List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    const pointField& points
)
:
    ParentType(triangles, points),
    patches_(patches),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface
(
    List<labelledTri>& triangles,
    const geometricSurfacePatchList& patches,
    pointField& points,
    const bool reuse
)
:
    ParentType(triangles, points, reuse),
    patches_(patches),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface
(
    const Xfer<List<labelledTri>>& triangles,
    const geometricSurfacePatchList& patches,
    const Xfer<List<point>>& points
)
:
    ParentType(triangles, points),
    patches_(patches),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{}


Foam::triSurface::triSurface
(
    const List<labelledTri>& triangles,
    const pointField& points
)
:
    ParentType(triangles, points),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    setDefaultPatches();
}


Foam::triSurface::triSurface
(
    const triFaceList& triangles,
    const pointField& points
)
:
    ParentType(convertToTri(triangles, 0), points),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    setDefaultPatches();
}


Foam::triSurface::triSurface(const fileName& name, const scalar scaleFactor)
:
    ParentType(List<Face>(), pointField()),
    patches_(),
    sortedEdgeFacesPtr_(nullptr),
    edgeOwnerPtr_(nullptr)
{
    const word ext = name.ext();

    read(name, ext);

    scalePoints(scaleFactor);

    setDefaultPatches();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurface::~triSurface()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::clearTopology()
{
    ParentType::clearTopology();
    deleteDemandDrivenData(sortedEdgeFacesPtr_);
    deleteDemandDrivenData(edgeOwnerPtr_);
}


void Foam::triSurface::clearPatchMeshAddr()
{
    ParentType::clearPatchMeshAddr();
}


void Foam::triSurface::clearOut()
{
    ParentType::clearOut();
    clearTopology();
    clearPatchMeshAddr();
}


void Foam::triSurface::swap(triSurface& surf)
{
    clearOut();
    surf.clearOut();

    FaceListType::swap(static_cast<FaceListType&>(surf));
    storedPoints().swap(surf.storedPoints());
    patches_.swap(surf.patches());
}


const Foam::labelListList& Foam::triSurface::sortedEdgeFaces() const
{
    if (!sortedEdgeFacesPtr_)
    {
        calcSortedEdgeFaces();
    }

    return *sortedEdgeFacesPtr_;
}


const Foam::labelList& Foam::triSurface::edgeOwner() const
{
    if (!edgeOwnerPtr_)
    {
        calcEdgeOwner();
    }

    return *edgeOwnerPtr_;
}


void Foam::triSurface::movePoints(const pointField& newPoints)
{
    // Remove all geometry dependent data
    deleteDemandDrivenData(sortedEdgeFacesPtr_);

    // Adapt for new point position
    ParentType::movePoints(newPoints);

    // Copy new points
    storedPoints() = newPoints;
}


void Foam::triSurface::scalePoints(const scalar scaleFactor)
{
    // Avoid bad scaling
    if (scaleFactor > VSMALL && scaleFactor != 1.0)
    {
        // Remove all geometry dependent data
        clearTopology();

        // Adapt for new point position
        ParentType::movePoints(pointField());

        storedPoints() *= scaleFactor;
    }
}


// Remove non-triangles, double triangles.
void Foam::triSurface::cleanup(const bool verbose)
{
    // Merge points (already done for STL, TRI)
    stitchTriangles(SMALL, verbose);

    // Merging points might have changed geometric factors
    clearOut();

    checkTriangles(verbose);

    checkEdges(verbose);
}


Foam::List<Foam::surfZone>
Foam::triSurface::sortedZones(labelList& faceMap) const
{
    surfacePatchList patches(calcPatches(faceMap));

    surfZoneList zones(patches.size());
    forAll(patches, patchi)
    {
        zones[patchi] = surfZone(patches[patchi]);
    }

    return zones;
}


void Foam::triSurface::triFaceFaces(List<face>& plainFaces) const
{
    plainFaces.setSize(size());

    forAll(*this, facei)
    {
        plainFaces[facei] = this->operator[](facei);
    }
}


Foam::Xfer<Foam::List<Foam::labelledTri>>
Foam::triSurface::xferFaces()
{
    // Topology changed because of transfer
    clearOut();

    return this->storedFaces().xfer();
}


Foam::Xfer<Foam::List<Foam::point>>
Foam::triSurface::xferPoints()
{
    // Topology changed because of transfer
    clearOut();

    return this->storedPoints().xfer();
}


// Finds area, starting at facei, delimited by borderEdge. Marks all visited
// faces (from face-edge-face walk) with currentZone.
void Foam::triSurface::markZone
(
    const boolList& borderEdge,
    const label facei,
    const label currentZone,
    labelList& faceZone
) const
{
    // List of faces whose faceZone has been set.
    labelList changedFaces(1, facei);

    while (true)
    {
        // Pick up neighbours of changedFaces
        DynamicList<label> newChangedFaces(2*changedFaces.size());

        for (const label facei : changedFaces)
        {
            const labelList& fEdges = faceEdges()[facei];

            for (const label edgei : fEdges)
            {
                if (!borderEdge[edgei])
                {
                    const labelList& eFaces = edgeFaces()[edgei];

                    for (const label nbrFacei : eFaces)
                    {
                        if (faceZone[nbrFacei] == -1)
                        {
                            faceZone[nbrFacei] = currentZone;
                            newChangedFaces.append(nbrFacei);
                        }
                        else if (faceZone[nbrFacei] != currentZone)
                        {
                            FatalErrorInFunction
                                << "Zones " << faceZone[nbrFacei]
                                << " at face " << nbrFacei
                                << " connects to zone " << currentZone
                                << " at face " << facei
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        if (newChangedFaces.empty())
        {
            break;
        }

        changedFaces.transfer(newChangedFaces);
    }
}


// Finds areas delimited by borderEdge (or 'real' edges).
// Fills faceZone accordingly
Foam::label Foam::triSurface::markZones
(
    const boolList& borderEdge,
    labelList& faceZone
) const
{
    faceZone.setSize(size());
    faceZone = -1;

    if (borderEdge.size() != nEdges())
    {
        FatalErrorInFunction
            << "borderEdge boolList not same size as number of edges" << endl
            << "borderEdge:" << borderEdge.size() << endl
            << "nEdges    :" << nEdges()
            << exit(FatalError);
    }

    label zoneI = 0;

    for (label startFacei = 0;; zoneI++)
    {
        // Find first non-coloured face
        for (; startFacei < size(); startFacei++)
        {
            if (faceZone[startFacei] == -1)
            {
                break;
            }
        }

        if (startFacei >= size())
        {
            break;
        }

        faceZone[startFacei] = zoneI;

        markZone(borderEdge, startFacei, zoneI, faceZone);
    }

    return zoneI;
}


void Foam::triSurface::subsetMeshMap
(
    const boolList& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const List<labelledTri>& locFaces = localFaces();

    label facei = 0;
    label pointi = 0;

    faceMap.setSize(locFaces.size());

    pointMap.setSize(nPoints());

    boolList pointHad(nPoints(), false);

    forAll(include, oldFacei)
    {
        if (include[oldFacei])
        {
            // Store new faces compact
            faceMap[facei++] = oldFacei;

            // Renumber labels for face
            const triSurface::FaceType& f = locFaces[oldFacei];

            for (const label verti : f)
            {
                if (!pointHad[verti])
                {
                    pointHad[verti] = true;
                    pointMap[pointi++] = verti;
                }
            }
        }
    }

    // Trim
    faceMap.setSize(facei);
    pointMap.setSize(pointi);
}


Foam::triSurface Foam::triSurface::subsetMesh
(
    const boolList& include,
    labelList& pointMap,
    labelList& faceMap
) const
{
    const pointField& locPoints = localPoints();
    const List<labelledTri>& locFaces = localFaces();

    // Fill pointMap, faceMap
    subsetMeshMap(include, pointMap, faceMap);


    // Create compact coordinate list and forward mapping array
    pointField newPoints(pointMap.size());
    labelList oldToNew(locPoints.size());
    forAll(pointMap, pointi)
    {
        newPoints[pointi] = locPoints[pointMap[pointi]];
        oldToNew[pointMap[pointi]] = pointi;
    }

    // Renumber triangle node labels and compact
    List<labelledTri> newTriangles(faceMap.size());

    forAll(faceMap, facei)
    {
        // Get old vertex labels
        const labelledTri& tri = locFaces[faceMap[facei]];

        newTriangles[facei][0] = oldToNew[tri[0]];
        newTriangles[facei][1] = oldToNew[tri[1]];
        newTriangles[facei][2] = oldToNew[tri[2]];
        newTriangles[facei].region() = tri.region();
    }

    // Construct subsurface
    return triSurface(newTriangles, patches(), newPoints, true);
}


void Foam::triSurface::reset
(
    const Xfer<List<labelledTri>>& triangles,
    const geometricSurfacePatchList& patches,
    const Xfer<List<point>>& points
)
{
    clearOut();

    this->storedFaces().transfer(triangles());
    this->storedPoints().transfer(points());

    patches_ = patches;
}


void Foam::triSurface::reset(MeshedSurface<labelledTri>& input)
{
    const auto& zones = input.surfZones();

    geometricSurfacePatchList patches(zones.size());
    forAll(zones, zonei)
    {
        patches[zonei] = geometricSurfacePatch(zones[zonei]);
    }

    this->reset(input.xferFaces(), patches, input.xferPoints());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::triSurface::operator=(const triSurface& surf)
{
    clearOut();

    FaceListType::operator=(static_cast<const FaceListType&>(surf));
    storedPoints() = surf.points();
    patches_ = surf.patches();
}


void Foam::triSurface::operator=(triSurface&& surf)
{
    clearOut();

    FaceListType::operator=(std::move(static_cast<FaceListType&>(surf)));
    storedPoints() = std::move(surf.storedPoints());
    patches_ = std::move(surf.patches());
}


// ************************************************************************* //
