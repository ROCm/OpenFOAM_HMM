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

#include "triSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "EdgeMap.H"
#include "triSurfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(triSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, triSurfaceMesh, dict);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Check file existence
const Foam::fileName& Foam::triSurfaceMesh::checkFile
(
    const fileName& fName,
    const fileName& objectName
)
{
    if (fName.empty())
    {
        FatalErrorIn
        (
            "triSurfaceMesh::checkFile(const fileName&, const fileName&)"
        )   << "Cannot find triSurfaceMesh starting from "
            << objectName << exit(FatalError);
    }
    return fName;
}


bool Foam::triSurfaceMesh::isSurfaceClosed() const
{
    // Construct pointFaces. Let's hope surface has compact point
    // numbering ...
    labelListList pointFaces;
    invertManyToMany(points().size(), *this, pointFaces);

    // Loop over all faces surrounding point. Count edges emanating from point.
    // Every edge should be used by two faces exactly.
    // To prevent doing work twice per edge only look at edges to higher
    // point
    EdgeMap<label> facesPerEdge(100);
    forAll(pointFaces, pointI)
    {
        const labelList& pFaces = pointFaces[pointI];

        facesPerEdge.clear();
        forAll(pFaces, i)
        {
            const labelledTri& f = triSurface::operator[](pFaces[i]);
            label fp = findIndex(f, pointI);

            // Forward edge
            {
                label p1 = f[f.fcIndex(fp)];

                if (p1 > pointI)
                {
                    const edge e(pointI, p1);
                    EdgeMap<label>::iterator eFnd = facesPerEdge.find(e);
                    if (eFnd != facesPerEdge.end())
                    {
                        if (eFnd() == 2)
                        {
                            return false;
                        }
                        eFnd()++;
                    }
                    else
                    {
                        facesPerEdge.insert(e, 1);
                    }
                }
            }
            // Reverse edge
            {
                label p1 = f[f.rcIndex(fp)];

                if (p1 > pointI)
                {
                    const edge e(pointI, p1);
                    EdgeMap<label>::iterator eFnd = facesPerEdge.find(e);
                    if (eFnd != facesPerEdge.end())
                    {
                        if (eFnd() == 2)
                        {
                            return false;
                        }
                        eFnd()++;
                    }
                    else
                    {
                        facesPerEdge.insert(e, 1);
                    }
                }
            }
        }

        // Check for any edges used only once.
        forAllConstIter(EdgeMap<label>, facesPerEdge, iter)
        {
            if (iter() != 2)
            {
                return false;
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const triSurface& s)
:
    searchableSurface(io),
    objectRegistry(io),
    triSurface(s),
    surfaceClosed_(-1)
{}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io)
:
    searchableSurface(io),
    objectRegistry(io),
    triSurface
    (
        checkFile
        (
            searchableSurface::filePath(),
            searchableSurface::objectPath()
        )
    ),
    surfaceClosed_(-1)
{}


Foam::triSurfaceMesh::triSurfaceMesh
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    objectRegistry(io),
    triSurface
    (
        checkFile
        (
            searchableSurface::filePath(),
            searchableSurface::objectPath()
        )
    ),
    surfaceClosed_(-1)
{
    scalar scaleFactor = 0;

    // allow rescaling of the surface points
    // eg, CAD geometries are often done in millimeters
    if (dict.readIfPresent("scale", scaleFactor) && scaleFactor > 0)
    {
        triSurface::scalePoints(scaleFactor);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::~triSurfaceMesh()
{
    clearOut();
}


void Foam::triSurfaceMesh::clearOut()
{
    tree_.clear();
    edgeTree_.clear();
    triSurface::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurfaceMesh::movePoints(const pointField& newPoints)
{
    tree_.clear();
    edgeTree_.clear();
    triSurface::movePoints(newPoints);
}


const Foam::indexedOctree<Foam::treeDataTriSurface>&
    Foam::triSurfaceMesh::tree() const
{
    if (tree_.empty())
    {
        treeBoundBox bb(points(), meshPoints());

        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        tree_.reset
        (
            new indexedOctree<treeDataTriSurface>
            (
                treeDataTriSurface(*this),
                bb.extend(rndGen, 1E-4),    // slightly randomize bb
                10,     // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return tree_();
}


const Foam::indexedOctree<Foam::treeDataEdge>&
    Foam::triSurfaceMesh::edgeTree() const
{
    if (edgeTree_.empty())
    {
        treeBoundBox bb(localPoints());

        // Boundary edges
        labelList bEdges
        (
            identity
            (
                nEdges()
               -nInternalEdges()
            )
          + nInternalEdges()
        );

        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    localPoints(),  // points
                    bEdges          // selected edges
                ),
                bb.extend(rndGen, 1E-4),    // slightly randomize bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }
    return edgeTree_();
}


const Foam::wordList& Foam::triSurfaceMesh::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(patches().size());
        forAll(regions_, regionI)
        {
            regions_[regionI] = patches()[regionI].name();
        }
    }
    return regions_;
}


// Find out if surface is closed.
bool Foam::triSurfaceMesh::hasVolumeType() const
{
    if (surfaceClosed_ == -1)
    {
        if (isSurfaceClosed())
        {
            surfaceClosed_ = 1;
        }
        else
        {
            surfaceClosed_ = 0;
        }
    }

    return surfaceClosed_ == 1;
}


void Foam::triSurfaceMesh::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(samples.size());

    forAll(samples, i)
    {
        static_cast<pointIndexHit&>(info[i]) =
            octree.findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::triSurfaceMesh::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    forAll(start, i)
    {
        static_cast<pointIndexHit&>(info[i]) = octree.findLine
        (
            start[i],
            end[i]
        );
    }
}


void Foam::triSurfaceMesh::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    forAll(start, i)
    {
        static_cast<pointIndexHit&>(info[i]) =
            octree.findLineAny(start[i], end[i]);
    }
}


void Foam::triSurfaceMesh::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit> >& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    info.setSize(start.size());

    // Work array
    DynamicList<pointIndexHit, 1, 1> hits;

    // Tolerances:
    // To find all intersections we add a small vector to the last intersection
    // This is chosen such that
    // - it is significant (SMALL is smallest representative relative tolerance;
    //   we need something bigger since we're doing calculations)
    // - if the start-end vector is zero we still progress
    const vectorField dirVec(end-start);
    const scalarField magSqrDirVec(magSqr(dirVec));
    const vectorField smallVec
    (
        Foam::sqrt(SMALL)*dirVec
      + vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL)
    );

    forAll(start, pointI)
    {
        // See if any intersection between pt and end
        pointIndexHit inter = octree.findLine(start[pointI], end[pointI]);

        if (inter.hit())
        {
            hits.clear();
            hits.append(inter);

            point pt = inter.hitPoint() + smallVec[pointI];

            while (((pt-start[pointI])&dirVec[pointI]) <= magSqrDirVec[pointI])
            {
                // See if any intersection between pt and end
                pointIndexHit inter = octree.findLine(pt, end[pointI]);

                // Check for not hit or hit same triangle as before (can happen
                // if vector along surface of triangle)
                if
                (
                    !inter.hit()
                 || (inter.index() == hits[hits.size()-1].index())
                )
                {
                    break;
                }
                hits.append(inter);

                pt = inter.hitPoint() + smallVec[pointI];
            }

            info[pointI].transfer(hits);
        }
        else
        {
            info[pointI].clear();
        }
    }
}


void Foam::triSurfaceMesh::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    forAll(info, i)
    {
        if (info[i].hit())
        {
            region[i] = triSurface::operator[](info[i].index()).region();
        }
        else
        {
            region[i] = -1;
        }
    }
}


void Foam::triSurfaceMesh::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());

    forAll(info, i)
    {
        if (info[i].hit())
        {
            normal[i] = faceNormals()[info[i].index()];
        }
        else
        {
            // Set to what?
            normal[i] = vector::zero;
        }
    }
}


void Foam::triSurfaceMesh::getField
(
    const word& fieldName,
    const List<pointIndexHit>& info,
    labelList& values
) const
{
    const triSurfaceLabelField& fld = lookupObject<triSurfaceLabelField>
    (
        fieldName
    );

    values.setSize(info.size());
    forAll(info, i)
    {
        if (info[i].hit())
        {
            values[i] = fld[info[i].index()];
        }
    }
}


void Foam::triSurfaceMesh::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());

    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        // - use cached volume type per each tree node
        // - cheat conversion since same values
        volType[pointI] = static_cast<volumeType>(tree().getVolumeType(pt));
    }
}


//- Write using given format, version and compression
bool Foam::triSurfaceMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    fileName fullPath(searchableSurface::objectPath());

    if (!mkDir(fullPath.path()))
    {
        return false;
    }

    triSurface::write(fullPath);

    if (!isFile(fullPath))
    {
        return false;
    }

    //return objectRegistry::writeObject(fmt, ver, cmp);
    return true;
}


// ************************************************************************* //
