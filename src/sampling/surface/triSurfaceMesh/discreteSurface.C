/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "discreteSurface.H"
#include "meshSearch.H"
#include "Tuple2.H"
#include "globalIndex.H"
#include "treeDataCell.H"
#include "treeDataFace.H"
#include "meshTools.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::discreteSurface::samplingSource
>
Foam::discreteSurface::samplingSourceNames_
{
    { samplingSource::cells, "cells" },
    { samplingSource::insideCells, "insideCells" },
    { samplingSource::boundaryFaces, "boundaryFaces" },
};


namespace Foam
{
    defineTypeNameAndDebug(discreteSurface, 0);

    //- Private class for finding nearest
    //  Comprising:
    //  - global index
    //  - sqr(distance)
    typedef Tuple2<scalar, label> nearInfo;

    class nearestEqOp
    {
    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first() < x.first())
            {
                x = y;
            }
        }
    };
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::discreteSurface::setZoneMap
(
    const surfZoneList& zoneLst,
    labelList& zoneIds
)
{
    label sz = 0;
    forAll(zoneLst, zonei)
    {
        const surfZone& zn = zoneLst[zonei];
        sz += zn.size();
    }

    zoneIds.setSize(sz);
    forAll(zoneLst, zonei)
    {
        const surfZone& zn = zoneLst[zonei];

        // Assign sub-zone Ids
        SubList<label>(zoneIds, zn.size(), zn.start()) = zonei;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataFace>&
Foam::discreteSurface::nonCoupledboundaryTree() const
{
    // Variant of meshSearch::boundaryTree() that only does non-coupled
    // boundary faces.

    if (!boundaryTreePtr_.valid())
    {
        // all non-coupled boundary faces (not just walls)
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        labelList bndFaces(mesh().nFaces()-mesh().nInternalFaces());
        label bndI = 0;
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];
            if (!pp.coupled())
            {
                forAll(pp, i)
                {
                    bndFaces[bndI++] = pp.start()+i;
                }
            }
        }
        bndFaces.setSize(bndI);


        treeBoundBox overallBb(mesh().points());
        Random rndGen(123456);
        // Extend slightly and make 3D
        overallBb = overallBb.extend(rndGen, 1e-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        boundaryTreePtr_.reset
        (
            new indexedOctree<treeDataFace>
            (
                treeDataFace    // all information needed to search faces
                (
                    false,                      // do not cache bb
                    mesh(),
                    bndFaces                    // boundary faces only
                ),
                overallBb,                      // overall search domain
                8,                              // maxLevel
                10,                             // leafsize
                3.0                             // duplicity
            )
        );
    }

    return boundaryTreePtr_();
}


bool Foam::discreteSurface::update(const meshSearch& meshSearcher)
{
    // Find the cells the triangles of the surface are in.
    // Does approximation by looking at the face centres only
    const pointField& fc = surface_.faceCentres();

    List<nearInfo> nearest(fc.size());

    // Global numbering for cells/faces - only used to uniquely identify local
    // elements
    globalIndex globalCells(onBoundary() ? mesh().nFaces() : mesh().nCells());

    forAll(nearest, i)
    {
        nearest[i].first()  = GREAT;
        nearest[i].second() = labelMax;
    }

    if (sampleSource_ == cells)
    {
        // Search for nearest cell

        const indexedOctree<treeDataCell>& cellTree = meshSearcher.cellTree();

        forAll(fc, triI)
        {
            pointIndexHit nearInfo = cellTree.findNearest
            (
                fc[triI],
                sqr(GREAT)
            );
            if (nearInfo.hit())
            {
                nearest[triI].first()  = magSqr(nearInfo.hitPoint()-fc[triI]);
                nearest[triI].second() = globalCells.toGlobal(nearInfo.index());
            }
        }
    }
    else if (sampleSource_ == insideCells)
    {
        // Search for cell containing point

        const indexedOctree<treeDataCell>& cellTree = meshSearcher.cellTree();

        forAll(fc, triI)
        {
            if (cellTree.bb().contains(fc[triI]))
            {
                const label index = cellTree.findInside(fc[triI]);
                if (index != -1)
                {
                    nearest[triI].first()  = 0.0;
                    nearest[triI].second() = globalCells.toGlobal(index);
                }
            }
        }
    }
    else
    {
        // Search for nearest boundaryFace

        ////- Search on all (including coupled) boundary faces
        //const indexedOctree<treeDataFace>& bTree = meshSearcher.boundaryTree()
        //- Search on all non-coupled boundary faces
        const indexedOctree<treeDataFace>& bTree = nonCoupledboundaryTree();

        forAll(fc, triI)
        {
            pointIndexHit nearInfo = bTree.findNearest
            (
                fc[triI],
                sqr(GREAT)
            );
            if (nearInfo.hit())
            {
                nearest[triI].first()  = magSqr(nearInfo.hitPoint()-fc[triI]);
                nearest[triI].second() = globalCells.toGlobal
                (
                    bTree.shapes().faceLabels()[nearInfo.index()]
                );
            }
        }
    }


    // See which processor has the nearest. Mark and subset
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Pstream::listCombineGather(nearest, nearestEqOp());
    Pstream::listCombineScatter(nearest);

    labelList cellOrFaceLabels(fc.size(), -1);

    label nFound = 0;
    forAll(nearest, triI)
    {
        if (nearest[triI].second() == labelMax)
        {
            // Not found on any processor. How to map?
        }
        else if (globalCells.isLocal(nearest[triI].second()))
        {
            cellOrFaceLabels[triI] = globalCells.toLocal
            (
                nearest[triI].second()
            );
            nFound++;
        }
    }


    if (debug)
    {
        Pout<< "Local out of faces:" << cellOrFaceLabels.size()
            << " keeping:" << nFound << endl;
    }

    // Now subset the surface. Do not use triSurface::subsetMesh since requires
    // original surface to be in compact numbering.

    const triSurface& s = surface_;

    // Compact to original triangle
    labelList faceMap(s.size());
    // Compact to original points
    labelList pointMap(s.points().size());
    // From original point to compact points
    labelList reversePointMap(s.points().size(), -1);

    // Handle region-wise sorting (makes things slightly more complicated)
    zoneIds_.setSize(s.size(), -1);

    // Better not to use triSurface::sortedZones here,
    // since we'll sort ourselves

    // Get zone/region sizes used, store under the original region Id
    Map<label> zoneSizes;

    // Recover region names from the input surface
    Map<word> zoneNames;
    {
        const geometricSurfacePatchList& patches = s.patches();

        forAll(patches, patchi)
        {
            zoneNames.set
            (
                patchi,
                (
                    patches[patchi].name().empty()
                  ? Foam::name("patch%d", patchi)
                  : patches[patchi].name()
                )
            );

            zoneSizes.set(patchi, 0);
        }
    }


    {
        label newPointi = 0;
        label newFacei = 0;

        forAll(s, facei)
        {
            if (cellOrFaceLabels[facei] != -1)
            {
                const triSurface::FaceType& f = s[facei];
                const label regionid = f.region();

                Map<label>::iterator fnd = zoneSizes.find(regionid);
                if (fnd != zoneSizes.end())
                {
                    fnd()++;
                }
                else
                {
                    // This shouldn't happen
                    zoneSizes.insert(regionid, 1);
                    zoneNames.set
                    (
                        regionid,
                        Foam::name("patch%d", regionid)
                    );
                }

                // Store new faces compact
                faceMap[newFacei] = facei;
                zoneIds_[newFacei] = regionid;
                ++newFacei;

                // Renumber face points
                forAll(f, fp)
                {
                    const label labI = f[fp];

                    if (reversePointMap[labI] == -1)
                    {
                        pointMap[newPointi] = labI;
                        reversePointMap[labI] = newPointi++;
                    }
                }
            }
        }

        // Trim
        faceMap.setSize(newFacei);
        zoneIds_.setSize(newFacei);
        pointMap.setSize(newPointi);
    }


    // Assign start/size (and name) to the newZones
    // re-use the lookup to map (zoneId => zoneI)
    surfZoneList zoneLst(zoneSizes.size());
    label start = 0;
    label zoneI = 0;
    forAllIter(Map<label>, zoneSizes, iter)
    {
        // No negative regionids, so Map<label> sorts properly
        const label regionid = iter.key();

        word name;
        Map<word>::const_iterator fnd = zoneNames.find(regionid);
        if (fnd != zoneNames.end())
        {
            name = fnd();
        }
        if (name.empty())
        {
            name = ::Foam::name("patch%d", regionid);
        }

        zoneLst[zoneI] = surfZone
        (
            name,
            0,           // initialize with zero size
            start,
            zoneI
        );

        // Adjust start for the next zone and save (zoneId => zoneI) mapping
        start += iter();
        iter() = zoneI++;
    }


    // At this stage:
    // - faceMap to map the (unsorted) compact to original triangle
    // - zoneIds for the next sorting
    // - zoneSizes contains region -> count information

    // Rebuild the faceMap for the sorted order
    labelList sortedFaceMap(faceMap.size());

    forAll(zoneIds_, facei)
    {
        label zonei = zoneIds_[facei];
        label sortedFacei = zoneLst[zonei].start() + zoneLst[zonei].size()++;
        sortedFaceMap[sortedFacei] = faceMap[facei];
    }

    // zoneIds are now simply flat values
    setZoneMap(zoneLst, zoneIds_);

    // Replace the faceMap with the properly sorted face map
    faceMap.transfer(sortedFaceMap);

    if (keepIds_)
    {
        originalIds_ = faceMap;
    }

    // Subset cellOrFaceLabels (for compact faces)
    cellOrFaceLabels = labelUIndList(cellOrFaceLabels, faceMap)();

    // Store any face per point (without using pointFaces())
    labelList pointToFace(pointMap.size());

    // Create faces and points for subsetted surface
    faceList& surfFaces = this->storedFaces();
    surfFaces.setSize(faceMap.size());

    this->storedZones().transfer(zoneLst);

    forAll(faceMap, facei)
    {
        const labelledTri& origF = s[faceMap[facei]];
        face& f = surfFaces[facei];

        f = triFace
        (
            reversePointMap[origF[0]],
            reversePointMap[origF[1]],
            reversePointMap[origF[2]]
        );

        forAll(f, fp)
        {
            pointToFace[f[fp]] = facei;
        }
    }

    this->storedPoints() = pointField(s.points(), pointMap);

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }


    // Collect the samplePoints and sampleElements
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (interpolate())
    {
        samplePoints_.setSize(pointMap.size());
        sampleElements_.setSize(pointMap.size(), -1);

        if (sampleSource_ == cells)
        {
            // samplePoints_   : per surface point a location inside the cell
            // sampleElements_ : per surface point the cell

            forAll(points(), pointi)
            {
                const point& pt = points()[pointi];
                label celli = cellOrFaceLabels[pointToFace[pointi]];
                sampleElements_[pointi] = celli;

                // Check if point inside cell
                if
                (
                    mesh().pointInCell
                    (
                        pt,
                        sampleElements_[pointi],
                        meshSearcher.decompMode()
                    )
                )
                {
                    samplePoints_[pointi] = pt;
                }
                else
                {
                    // Find nearest point on faces of cell
                    const cell& cFaces = mesh().cells()[celli];

                    scalar minDistSqr = VGREAT;

                    forAll(cFaces, i)
                    {
                        const face& f = mesh().faces()[cFaces[i]];
                        pointHit info = f.nearestPoint(pt, mesh().points());
                        if (info.distance() < minDistSqr)
                        {
                            minDistSqr = info.distance();
                            samplePoints_[pointi] = info.rawPoint();
                        }
                    }
                }
            }
        }
        else if (sampleSource_ == insideCells)
        {
            // samplePoints_   : per surface point a location inside the cell
            // sampleElements_ : per surface point the cell

            forAll(points(), pointi)
            {
                const point& pt = points()[pointi];
                label celli = cellOrFaceLabels[pointToFace[pointi]];
                sampleElements_[pointi] = celli;
                samplePoints_[pointi] = pt;
            }
        }
        else
        {
            // samplePoints_   : per surface point a location on the boundary
            // sampleElements_ : per surface point the boundary face containing
            //                   the location

            forAll(points(), pointi)
            {
                const point& pt = points()[pointi];
                label facei = cellOrFaceLabels[pointToFace[pointi]];
                sampleElements_[pointi] = facei;
                samplePoints_[pointi] = mesh().faces()[facei].nearestPoint
                (
                    pt,
                    mesh().points()
                ).rawPoint();
            }
        }
    }
    else
    {
        // if sampleSource_ == cells:
        //      sampleElements_ : per surface triangle the cell
        //      samplePoints_   : n/a
        // if sampleSource_ == insideCells:
        //      sampleElements_ : -1 or per surface triangle the cell
        //      samplePoints_   : n/a
        // else:
        //      sampleElements_ : per surface triangle the boundary face
        //      samplePoints_   : n/a
        sampleElements_.transfer(cellOrFaceLabels);
        samplePoints_.clear();
    }


    if (debug)
    {
        this->clearOut();
        OFstream str(mesh().time().path()/"surfaceToMesh.obj");
        Info<< "Dumping correspondence from local surface (points or faces)"
            << " to mesh (cells or faces) to " << str.name() << endl;

        const vectorField& centres =
        (
            onBoundary()
          ? mesh().faceCentres()
          : mesh().cellCentres()
        );

        if (interpolate())
        {
            label vertI = 0;
            forAll(samplePoints_, pointi)
            {
                meshTools::writeOBJ(str, points()[pointi]);
                vertI++;

                meshTools::writeOBJ(str, samplePoints_[pointi]);
                vertI++;

                label elemi = sampleElements_[pointi];
                meshTools::writeOBJ(str, centres[elemi]);
                vertI++;
                str << "l " << vertI-2 << ' ' << vertI-1 << ' ' << vertI << nl;
            }
        }
        else
        {
            label vertI = 0;
            forAll(sampleElements_, triI)
            {
                meshTools::writeOBJ(str, faceCentres()[triI]);
                vertI++;

                label elemi = sampleElements_[triI];
                meshTools::writeOBJ(str, centres[elemi]);
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }

    needsUpdate_ = false;
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::discreteSurface::discreteSurface
(
    const polyMesh& mesh,
    const word& surfaceName,
    const samplingSource sampleSource,
    const bool allowInterpolate
)
:
    MeshStorage(),
    mesh_(mesh),
    allowInterpolate_(allowInterpolate),
    interpolate_(false),
    surface_
    (
        IOobject
        (
            surfaceName,
            mesh.time().constant(), // instance
            "triSurface",           // local
            mesh,                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    sampleSource_(sampleSource),
    needsUpdate_(true),
    keepIds_(false),
    originalIds_(),
    zoneIds_(),
    sampleElements_(0),
    samplePoints_(0)
{}


Foam::discreteSurface::discreteSurface
(
    const polyMesh& mesh,
    const dictionary& dict,
    const bool allowInterpolate
)
:
    MeshStorage(),
    mesh_(mesh),
    allowInterpolate_(allowInterpolate),
    interpolate_
    (
        allowInterpolate
     && dict.lookupOrDefault<Switch>("interpolate", false)
    ),
    surface_
    (
        IOobject
        (
            dict.lookup("surface"),
            mesh.time().constant(), // instance
            "triSurface",           // local
            mesh,                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    sampleSource_(samplingSourceNames_.lookup("source", dict)),
    needsUpdate_(true),
    keepIds_(dict.lookupOrDefault<Switch>("keepIds", false)),
    originalIds_(),
    zoneIds_(),
    sampleElements_(0),
    samplePoints_(0)
{}


Foam::discreteSurface::discreteSurface
(
    const word& name,
    const polyMesh& mesh,
    const triSurface& surface,
    const word& sampleSourceName,
    const bool allowInterpolate
)
:
    MeshStorage(),
    mesh_(mesh),
    allowInterpolate_(allowInterpolate),
    interpolate_(false),
    surface_
    (
        IOobject
        (
            name,
            mesh.time().constant(), // instance
            "triSurface",           // local
            mesh,                   // registry
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        surface
    ),
    sampleSource_(samplingSourceNames_[sampleSourceName]),
    needsUpdate_(true),
    keepIds_(false),
    originalIds_(),
    zoneIds_(),
    sampleElements_(0),
    samplePoints_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::discreteSurface::~discreteSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::discreteSurface::interpolate() const
{
    return interpolate_;
}

bool Foam::discreteSurface::interpolate(bool b)
{
    if (interpolate_ == b)
    {
        return false;
    }

    if (b && !allowInterpolate_)
    {
        return false;
    }

    // Value changed, old sampling information is invalid
    interpolate_ = b;
    expire();

    return true;
}


bool Foam::discreteSurface::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::discreteSurface::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    MeshStorage::clear();
    zoneIds_.clear();

    originalIds_.clear();
    boundaryTreePtr_.clear();
    sampleElements_.clear();
    samplePoints_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::discreteSurface::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Calculate surface and mesh overlap bounding box
    treeBoundBox bb
    (
        surface_.triSurface::points(),
        surface_.triSurface::meshPoints()
    );
    bb.min() = max(bb.min(), mesh().bounds().min());
    bb.max() = min(bb.max(), mesh().bounds().max());

    // Extend a bit
    const vector span(bb.span());

    bb.min() -= 0.5*span;
    bb.max() += 0.5*span;

    bb.inflate(1e-6);

    // Mesh search engine, no triangulation of faces.
    meshSearch meshSearcher(mesh(), bb, polyMesh::FACE_PLANES);

    return update(meshSearcher);
}


bool Foam::discreteSurface::update(const treeBoundBox& bb)
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Mesh search engine on subset, no triangulation of faces.
    meshSearch meshSearcher(mesh(), bb, polyMesh::FACE_PLANES);

    return update(meshSearcher);
}


bool Foam::discreteSurface::sampleAndStore
(
    const objectRegistry& store,
    const word& fieldName
) const
{
    return
    (
        sampleType<scalar>(store, fieldName)
     || sampleType<vector>(store, fieldName)
     || sampleType<sphericalTensor>(store, fieldName)
     || sampleType<symmTensor>(store, fieldName)
     || sampleType<tensor>(store, fieldName)
    );
}


void Foam::discreteSurface::print(Ostream& os) const
{
    os  << "discreteSurface: "
        << " surface:" << surface_.objectRegistry::name()
        << " faces:"   << this->surfFaces().size()
        << " points:"  << this->points().size()
        << " zoneids:" << this->zoneIds().size();
}


// ************************************************************************* //
