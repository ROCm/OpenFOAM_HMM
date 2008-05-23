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

Class
    autoHexMeshDriver

\*----------------------------------------------------------------------------*/

#include "autoHexMeshDriver.H"
#include "fvMesh.H"
#include "Time.H"
#include "boundBox.H"
#include "globalIndex.H"
#include "wallPolyPatch.H"
#include "mapDistributePolyMesh.H"
#include "cellSet.H"
#include "syncTools.H"
#include "motionSmoother.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(autoHexMeshDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Check writing tolerance before doing any serious work
Foam::scalar Foam::autoHexMeshDriver::getMergeDistance() const
{
    const scalar mergeTol = readScalar(dict_.lookup("mergeTolerance"));
    const boundBox& meshBb = mesh_.bounds();
    scalar mergeDist = mergeTol*mag(meshBb.max() - meshBb.min());
    scalar writeTol = std::pow
    (
        scalar(10.0),
       -scalar(IOstream::defaultPrecision())
    );

    Info<< "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    if (mesh_.time().writeFormat() == IOstream::ASCII && mergeTol < writeTol)
    {
        FatalErrorIn("autoHexMeshDriver::getMergeDistance() const")
            << "Your current settings specify ASCII writing with "
            << IOstream::defaultPrecision() << " digits precision." << endl
            << "Your merging tolerance (" << mergeTol << ") is finer than this."
            << endl
            << "Please change your writeFormat to binary"
            << " or increase the writePrecision" << endl
            << "or adjust the merge tolerance (-mergeTol)."
            << exit(FatalError);
    }

    return mergeDist;
}


// Return per keeppoint -1 or the local cell label the point is in. Guaranteed
// to be only on one processor.
Foam::labelList Foam::autoHexMeshDriver::findCells(const pointField& keepPoints)
 const
{
    // Global calculation engine
    globalIndex globalCells(mesh_.nCells());

    // Cell label per point
    labelList cellLabels(keepPoints.size());

    forAll(keepPoints, i)
    {
        const point& keepPoint = keepPoints[i];

        label localCellI = mesh_.findCell(keepPoint);

        label globalCellI = -1;

        if (localCellI != -1)
        {
            Pout<< "Found point " << keepPoint << " in cell " << localCellI
                << " on processor " << Pstream::myProcNo() << endl;
            globalCellI = globalCells.toGlobal(localCellI);
        }

        reduce(globalCellI, maxOp<label>());

        if (globalCellI == -1)
        {
            FatalErrorIn
            (
                "autoHexMeshDriver::findCells(const pointField&) const"
            )   << "Point " << keepPoint << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }

        if (globalCells.isLocal(globalCellI))
        {
            cellLabels[i] = localCellI;
        }
        else
        {
            cellLabels[i] = -1;
        }
    }
    return cellLabels;
}


// Specifically orient using a calculated point outside
void Foam::autoHexMeshDriver::orientOutside(PtrList<searchableSurface>& shells)
{
    // Determine outside point.
    boundBox overallBb
    (
        point(GREAT, GREAT, GREAT),
        point(-GREAT, -GREAT, -GREAT)
    );

    bool hasSurface = false;

    forAll(shells, shellI)
    {
        if (isA<triSurfaceMesh>(shells[shellI]))
        {
            const triSurfaceMesh& shell =
                refCast<const triSurfaceMesh>(shells[shellI]);

            hasSurface = true;

            boundBox shellBb(shell.localPoints(), false);

            overallBb.min() = min(overallBb.min(), shellBb.min());
            overallBb.max() = max(overallBb.max(), shellBb.max());
        }
    }

    if (hasSurface)
    {
        const point outsidePt(2*overallBb.max() - overallBb.min());

        //Info<< "Using point " << outsidePt << " to orient shells" << endl;

        forAll(shells, shellI)
        {
            if (isA<triSurfaceMesh>(shells[shellI]))
            {
                triSurfaceMesh& shell = refCast<triSurfaceMesh>(shells[shellI]);

                if (!refinementSurfaces::isSurfaceClosed(shell))
                {
                    FatalErrorIn("orientOutside(PtrList<searchableSurface>&)")
                        << "Refinement shell "
                        << shell.searchableSurface::name()
                        << " is not closed." << exit(FatalError);
                }

                refinementSurfaces::orientSurface(outsidePt, shell);
            }
        }
    }
}


// Check that face zones are synced
void Foam::autoHexMeshDriver::checkCoupledFaceZones() const
{
    const faceZoneMesh& fZones = mesh_.faceZones();

    // Check any zones are present anywhere and in same order

    {
        List<wordList> zoneNames(Pstream::nProcs());
        zoneNames[Pstream::myProcNo()] = fZones.names();
        Pstream::gatherList(zoneNames);
        Pstream::scatterList(zoneNames);
        // All have same data now. Check.
        forAll(zoneNames, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                if (zoneNames[procI] != zoneNames[Pstream::myProcNo()])
                {
                    FatalErrorIn
                    (
                        "autoHexMeshDriver::checkCoupledFaceZones() const"
                    )   << "faceZones are not synchronised on processors." << nl
                        << "Processor " << procI << " has faceZones "
                        << zoneNames[procI] << nl
                        << "Processor " << Pstream::myProcNo()
                        << " has faceZones "
                        << zoneNames[Pstream::myProcNo()] << nl
                        << exit(FatalError);
                }
            }
        }
    }

    // Check that coupled faces are present on both sides.

    labelList faceToZone(mesh_.nFaces()-mesh_.nInternalFaces(), -1);

    forAll(fZones, zoneI)
    {
        const faceZone& fZone = fZones[zoneI];

        forAll(fZone, i)
        {
            label bFaceI = fZone[i]-mesh_.nInternalFaces();

            if (bFaceI >= 0)
            {
                if (faceToZone[bFaceI] == -1)
                {
                    faceToZone[bFaceI] = zoneI;
                }
                else if (faceToZone[bFaceI] == zoneI)
                {
                    FatalErrorIn
                    (
                        "autoHexMeshDriver::checkCoupledFaceZones()"
                    )   << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is twice in zone!"
                        << abort(FatalError);
                }
                else
                {
                    FatalErrorIn
                    (
                        "autoHexMeshDriver::checkCoupledFaceZones()"
                    )   << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is also in zone "
                        << fZones[faceToZone[bFaceI]].name()
                        << abort(FatalError);
                }
            }
        }
    }

    labelList neiFaceToZone(faceToZone);
    syncTools::swapBoundaryFaceList(mesh_, neiFaceToZone, false);

    forAll(faceToZone, i)
    {
        if (faceToZone[i] != neiFaceToZone[i])
        {
            FatalErrorIn
            (
                "autoHexMeshDriver::checkCoupledFaceZones()"
            )   << "Face " << mesh_.nInternalFaces()+i
                << " is in zone " << faceToZone[i]
                << ", its coupled face is in zone " << neiFaceToZone[i]
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::autoHexMeshDriver::autoHexMeshDriver
(
    fvMesh& mesh,
    const dictionary& dict,
    const dictionary& decomposeDict
)
:
    mesh_(mesh),
    dict_(dict),
    debug_(readLabel(dict_.lookup("debug"))),
    maxGlobalCells_(readLabel(dict_.lookup("cellLimit"))),
    maxLocalCells_(readLabel(dict_.lookup("procCellLimit"))),
    minRefineCells_(readLabel(dict_.lookup("minimumRefine"))),
    curvature_(readScalar(dict_.lookup("curvature"))),
    nBufferLayers_(readLabel(dict_.lookup("nBufferLayers"))),
    keepPoints_(dict_.lookup("keepPoints")),
    mergeDist_(getMergeDistance())
{

    if (debug_ > 0)
    {
        meshRefinement::debug = debug_;
        autoHexMeshDriver::debug = debug_;
    }

    Info<< "Overall cell limit                         : " << maxGlobalCells_
        << endl;
    Info<< "Per processor cell limit                   : " << maxLocalCells_
        << endl;
    Info<< "Minimum number of cells to refine          : " << minRefineCells_
        << endl;
    Info<< "Curvature                                  : " << curvature_
        << nl << endl;
    Info<< "Layers between different refinement levels : " << nBufferLayers_
        << endl;

    // Check keepPoints are sensible
    findCells(keepPoints_);


    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    {
        Info<< "Reading refinement shells." << endl;

        PtrList<dictionary> shellDicts(dict_.lookup("refinementShells"));

        shells_.setSize(shellDicts.size());
        shellLevels_.setSize(shellDicts.size());
        shellRefineInside_.setSize(shellDicts.size());

        forAll(shellDicts, i)
        {
            const dictionary& dict = shellDicts[i];

            shells_.set
            (
                i,
                searchableSurface::New
                (
                    dict.lookup("type"),
                    dict.lookup("name"),
                    mesh_.time(),
                    dict
                )
            );
            shellLevels_[i] = readLabel(dict.lookup("level"));
            shellRefineInside_[i] = Switch(dict.lookup("refineInside"));

            if (shellRefineInside_[i])
            {
                Info<< "Refinement level " << shellLevels_[i]
                    << " for all cells inside " << shells_[i].name() << endl;
            }
            else
            {
                Info<< "Refinement level " << shellLevels_[i]
                    << " for all cells outside " << shells_[i].name() << endl;
            }
        }

        Info<< "Read refinement shells in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // Orient shell surfaces before any searching is done.
        Info<< "Orienting triSurface shells so point far away is outside."
            << endl;
        orientOutside(shells_);
        Info<< "Oriented shells in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;
    }


    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    {
        Info<< "Reading surfaces and constructing search trees." << endl;

        surfacesPtr_.reset
        (
            new refinementSurfaces
            (
                IOobject
                (
                    "",                                 // dummy name
                    mesh_.time().constant(),            // directory
                    "triSurface",                       // instance
                    mesh_.time(),                       // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                dict_.lookup("surfaces")
            )
        );
        Info<< "Read surfaces in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // Orient surfaces (if they're closed) before any searching is done.
        Info<< "Orienting (closed) surfaces so keepPoint is outside." << endl;
        forAll(surfaces(), i)
        {
            if (refinementSurfaces::isSurfaceClosed(surfaces()[i]))
            {
                refinementSurfaces::orientSurface
                (
                    keepPoints_[0],
                    surfacesPtr_()[i]
                );
            }
        }
        Info<< "Oriented closed surfaces in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        Info<< "Setting refinement level of surface to be consistent"
            << " with shells." << endl;
        surfacesPtr_().setMinLevelFields
        (
            shells_,
            shellLevels_,
            shellRefineInside_
        );
        Info<< "Checked shell refinement in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;
    }

    // Check faceZones are synchronised
    checkCoupledFaceZones();


    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        Info<< nl
            << "Adding patches for surface regions" << nl
            << "----------------------------------" << nl
            << endl;

        // From global region number to mesh patch.
        globalToPatch_.setSize(surfaces().nRegions(), -1);

        Info<< "Patch\tRegion" << nl
            << "-----\t------"
            << endl;

        forAll(surfaces(), surfI)
        {
            const triSurfaceMesh& s = surfaces()[surfI];

            Info<< surfaces().names()[surfI] << ':' << nl << nl;

            const geometricSurfacePatchList& regions = s.patches();

            labelList nTrisPerRegion(surfaces().countRegions(s));

            forAll(regions, i)
            {
                if (nTrisPerRegion[i] > 0)
                {
                    label patchI = meshRefinement::addPatch
                    (
                        mesh,
                        //s.searchableSurface::name() + '_' + regions[i].name(),
                        surfaces().names()[surfI] + '_' + regions[i].name(),
                        wallPolyPatch::typeName
                    );

                    Info<< patchI << '\t' << regions[i].name() << nl;

                    globalToPatch_[surfaces().globalRegion(surfI, i)] = patchI;
                }
            }

            Info<< nl;
        }
        Info<< "Added patches in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;
    }


    //// Add cyclics for any named faceZones
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //// (these cyclics are used later on to temporarily put the faceZones
    ////  in when snapping)
    //
    //labelList namedSurfaces(surfaces().getNamedSurfaces());
    //if (namedSurfaces.size() > 0)
    //{
    //    Info<< nl
    //        << "Introducing cyclics for faceZones" << nl
    //        << "---------------------------------" << nl
    //        << endl;
    //
    //    // From surface to cyclic patch
    //    surfaceToCyclicPatch_.setSize(surfaces().size(), -1);
    //
    //    Info<< "Patch\tZone" << nl
    //        << "----\t-----"
    //        << endl;
    //
    //    forAll(namedSurfaces, i)
    //    {
    //        label surfI = namedSurfaces[i];
    //
    //        surfaceToCyclicPatch_[surfI] = meshRefinement::addPatch
    //        (
    //            mesh,
    //            surfaces().faceZoneNames()[surfI],
    //            cyclicPolyPatch::typeName
    //        );
    //
    //        Info<< surfaceToCyclicPatch_[surfI] << '\t'
    //            << surfaces().faceZoneNames()[surfI] << nl << endl;
    //    }
    //    Info<< "Added cyclic patches in = "
    //        << mesh_.time().cpuTimeIncrement() << " s" << endl;
    //}


    // Parallel
    // ~~~~~~~~

    {
        // Decomposition
        decomposerPtr_ = decompositionMethod::New
        (
            decomposeDict,
            mesh_
        );
        decompositionMethod& decomposer = decomposerPtr_();


        if (Pstream::parRun() && !decomposer.parallelAware())
        {
            FatalErrorIn("autoHexMeshDriver::autoHexMeshDriver(const IOobject&, fvMesh&)")
                << "You have selected decomposition method "
                << decomposer.typeName
                << " which is not parallel aware." << endl
                << "Please select one that is (parMetis, hierarchical)"
                << exit(FatalError);
        }

        // Mesh distribution engine (uses tolerance to reconstruct meshes)
        distributorPtr_.reset(new fvMeshDistribute(mesh_, mergeDist_));
    }


    // Refinement engine
    // ~~~~~~~~~~~~~~~~~

    {
        Info<< nl
            << "Determining initial surface intersections" << nl
            << "-----------------------------------------" << nl
            << endl;

        // Main refinement engine
        meshRefinerPtr_.reset
        (
            new meshRefinement
            (
                mesh,
                mergeDist_,         // tolerance used in sorting coordinates
                surfaces()
            )
        );
        Info<< "Calculated surface intersections in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // Some stats
        meshRefinerPtr_().printMeshInfo(debug_, "Initial mesh");

        meshRefinerPtr_().write
        (
            debug_&meshRefinement::OBJINTERSECTIONS,
            mesh_.time().path()/mesh_.time().timeName()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Read explicit feature edges
Foam::label Foam::autoHexMeshDriver::readFeatureEdges()
{
    Info<< "Reading external feature lines." << endl;

    PtrList<dictionary> featDicts(dict_.lookup("features"));

    featureMeshes_.setSize(featDicts.size());
    featureLevels_.setSize(featDicts.size());

    forAll(featDicts, i)
    {
        const dictionary& dict = featDicts[i];

        fileName featFileName(dict.lookup("file"));

        featureMeshes_.set
        (
            i,
            new featureEdgeMesh
            (
                IOobject
                (
                    featFileName,           // name
                    mesh_.time().constant(),// directory
                    "triSurface",           // instance
                    mesh_.db(),             // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );

        featureMeshes_[i].mergePoints(mergeDist_);
        featureLevels_[i] = readLabel(dict.lookup("level"));

        Info<< "Refinement level " << featureLevels_[i]
            << " for all cells crossed by feature " << featFileName
            << " (" << featureMeshes_[i].points().size() << " points, "
            << featureMeshes_[i].edges().size() << ")." << endl;
    }

    Info<< "Read feature lines in = "
        << mesh_.time().cpuTimeIncrement() << " s" << endl;

    return featureMeshes_.size();
}


Foam::label Foam::autoHexMeshDriver::featureEdgeRefine
(
    const label maxIter,
    const label minRefine
)
{
    // Read explicit feature edges
    readFeatureEdges();


    meshRefinement& meshRefiner = meshRefinerPtr_();

    label iter = 0;

    if (featureMeshes_.size() > 0 && maxIter > 0)
    {
        for (; iter < maxIter; iter++)
        {
            Info<< nl
                << "Feature refinement iteration " << iter << nl
                << "------------------------------" << nl
                << endl;

            labelList candidateCells
            (
                meshRefiner.refineCandidates
                (
                    keepPoints_[0],     // For now only use one.
                    globalToPatch_,
                    curvature_,

                    featureMeshes_,
                    featureLevels_,

                    shells_,
                    shellLevels_,
                    shellRefineInside_,

                    true,               // featureRefinement
                    false,              // internalRefinement
                    false,              // surfaceRefinement
                    false,              // curvatureRefinement
                    maxGlobalCells_,
                    maxLocalCells_
                )
            );
            labelList cellsToRefine
            (
                meshRefiner.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );
            Info<< "Determined cells to refine in = "
                << mesh_.time().cpuTimeIncrement() << " s" << endl;



            label nCellsToRefine = cellsToRefine.size();
            reduce(nCellsToRefine, sumOp<label>());

            Info<< "Selected for feature refinement : " << nCellsToRefine
                << " cells (out of " << mesh_.globalData().nTotalCells()
                << ')' << endl;

            if (nCellsToRefine <= minRefine)
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }


            if (debug_ > 0)
            {
                const_cast<Time&>(mesh_.time())++;
            }

            meshRefiner.refineAndBalance
            (
                "feature refinement iteration " + name(iter),
                decomposerPtr_(),
                distributorPtr_(),
                cellsToRefine
            );
        }
    }
    return iter;
}


Foam::label Foam::autoHexMeshDriver::surfaceOnlyRefine(const label maxIter)
{
    meshRefinement& meshRefiner = meshRefinerPtr_();

    // Determine the maximum refinement level over all surfaces. This
    // determines the minumum number of surface refinement iterations.
    label overallMaxLevel = max(surfaces().maxLevel());

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Surface refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        labelList candidateCells
        (
            meshRefiner.refineCandidates
            (
                keepPoints_[0],
                globalToPatch_,
                curvature_,

                featureMeshes_,
                featureLevels_,

                shells_,
                shellLevels_,
                shellRefineInside_,

                false,              // featureRefinement
                false,              // internalRefinement
                true,               // surfaceRefinement
                true,               // curvatureRefinement
                maxGlobalCells_,
                maxLocalCells_
            )
        );
        labelList cellsToRefine
        (
            meshRefiner.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh_.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum nessecary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxLevel
             && nCellsToRefine <= minRefineCells_
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug_)
        {
            const_cast<Time&>(mesh_.time())++;
        }

        meshRefiner.refineAndBalance
        (
            "surface refinement iteration " + name(iter),
            decomposerPtr_(),
            distributorPtr_(),
            cellsToRefine
        );
    }
    return iter;
}


void Foam::autoHexMeshDriver::removeInsideCells(const label nBufferLayers)
{

    meshRefinement& meshRefiner = meshRefinerPtr_();

    Info<< nl
        << "Removing mesh beyond surface intersections" << nl
        << "------------------------------------------" << nl
        << endl;

    if (debug_)
    {
       const_cast<Time&>(mesh_.time())++;
    }

    meshRefiner.splitMesh
    (
        nBufferLayers,                  // nBufferLayers
        globalToPatch_,
        keepPoints_[0]
    );

    if (debug_)
    {
        Pout<< "Writing subsetted mesh to time "
            << mesh_.time().timeName() << '.' << endl;
        meshRefiner.write(debug_, mesh_.time().path()/mesh_.time().timeName());
        Pout<< "Dumped mesh in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


Foam::label Foam::autoHexMeshDriver::shellRefine(const label maxIter)
{
    meshRefinement& meshRefiner = meshRefinerPtr_();


    // Mark current boundary faces with 0. Have meshRefiner maintain them.
    meshRefiner.userFaceData().setSize(1);

    // mark list to remove any refined faces
    meshRefiner.userFaceData()[0].first() = meshRefinement::REMOVE;
    meshRefiner.userFaceData()[0].second() = createWithValues<labelList>
    (
        mesh_.nFaces(),
        -1,
        meshRefiner.intersectedFaces(),
        0
    );

    // Determine the maximum refinement level over all volume refinement
    // regions. This determines the minumum number of shell refinement
    // iterations.
    label overallMaxShellLevel;
    if (shellLevels_.size() == 0)
    {
        overallMaxShellLevel = 0;
    }
    else
    {
        overallMaxShellLevel = max(shellLevels_);
    }

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Shell refinement iteration " << iter << nl
            << "----------------------------" << nl
            << endl;

        labelList candidateCells
        (
            meshRefiner.refineCandidates
            (
                keepPoints_[0],
                globalToPatch_,
                curvature_,

                featureMeshes_,
                featureLevels_,

                shells_,
                shellLevels_,
                shellRefineInside_,

                false,              // featureRefinement
                true,               // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                maxGlobalCells_,
                maxLocalCells_
            )
        );

        if (debug_)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromShells." << endl;

            cellSet(mesh_, "candidateCellsFromShells", candidateCells).write();
        }

        // Problem choosing starting faces for bufferlayers (bFaces)
        //  - we can't use the current intersected boundary faces
        //    (intersectedFaces) since this grows indefinitely
        //  - if we use 0 faces we don't satisfy bufferLayers from the
        //    surface.
        //  - possibly we want to have bFaces only the initial set of faces
        //    and maintain the list while doing the refinement.
        labelList bFaces
        (
            findIndices(meshRefiner.userFaceData()[0].second(), 0)
        );

        Info<< "Collected boundary faces : "
            << returnReduce(bFaces.size(), sumOp<label>()) << endl;

        labelList cellsToRefine;

        if (nBufferLayers_ <= 2)
        {
            cellsToRefine = meshRefiner.meshCutter().consistentSlowRefinement
            (
                nBufferLayers_,
                candidateCells,                 // cells to refine
                bFaces,                         // faces for nBufferLayers
                1,                              // point difference
                meshRefiner.intersectedPoints() // points to check
            );
        }
        else
        {
            cellsToRefine = meshRefiner.meshCutter().consistentSlowRefinement2
            (
                nBufferLayers_,
                candidateCells,                 // cells to refine
                bFaces                          // faces for nBufferLayers
            );
        }

        Info<< "Determined cells to refine in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for internal refinement : " << nCellsToRefine
            << " cells (out of " << mesh_.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum nessecary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxShellLevel
             && nCellsToRefine <= minRefineCells_
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug_)
        {
            const_cast<Time&>(mesh_.time())++;
        }

        meshRefiner.refineAndBalance
        (
            "shell refinement iteration " + name(iter),
            decomposerPtr_(),
            distributorPtr_(),
            cellsToRefine
        );
    }
    meshRefiner.userFaceData().clear();

    return iter;
}


void Foam::autoHexMeshDriver::baffleAndSplitMesh(const bool handleSnapProblems)
{
    meshRefinement& meshRefiner = meshRefinerPtr_();

    Info<< nl
        << "Splitting mesh at surface intersections" << nl
        << "---------------------------------------" << nl
        << endl;

    // Introduce baffles at surface intersections. Note:
    // meshRefiment::surfaceIndex() will
    // be like boundary face from now on so not coupled anymore.
    meshRefiner.baffleAndSplitMesh
    (
        handleSnapProblems,
        !handleSnapProblems,            // merge free standing baffles?
        const_cast<Time&>(mesh_.time()),
        globalToPatch_,
        keepPoints_[0]
    );
}


void Foam::autoHexMeshDriver::zonify()
{
    // Mesh is at its finest. Do zoning
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This puts all faces with intersection across a zoneable surface
    // into that surface's faceZone. All cells inside faceZone get given the
    // same cellZone.


    meshRefinement& meshRefiner = meshRefinerPtr_();

    if (surfaces().getNamedSurfaces().size() > 0)
    {
        Info<< nl
            << "Introducing zones for interfaces" << nl
            << "--------------------------------" << nl
            << endl;

        if (debug_)
        {
            const_cast<Time&>(mesh_.time())++;
        }

        meshRefiner.zonify(keepPoints_[0]);

        if (debug_)
        {
            Pout<< "Writing zoned mesh to time "
                << mesh_.time().timeName() << '.' << endl;
            meshRefiner.write
            (
                debug_,
                mesh_.time().path()/mesh_.time().timeName()
            );
        }

        // Check that all faces are synced
        checkCoupledFaceZones();
    }
}


void Foam::autoHexMeshDriver::splitAndMergeBaffles
(
    const bool handleSnapProblems
)
{
    Info<< nl
        << "Handling cells with snap problems" << nl
        << "---------------------------------" << nl
        << endl;

    // Introduce baffles and split mesh
    if (debug_)
    {
        const_cast<Time&>(mesh_.time())++;
    }

    meshRefinement& meshRefiner = meshRefinerPtr_();

    meshRefiner.baffleAndSplitMesh
    (
        handleSnapProblems,
        false,                  // merge free standing baffles?
        const_cast<Time&>(mesh_.time()),
        globalToPatch_,
        keepPoints_[0]
    );

    if (debug_)
    {
        const_cast<Time&>(mesh_.time())++;
    }

    // Duplicate points on baffles that are on more than one cell
    // region. This will help snapping pull them to separate surfaces.
    meshRefiner.dupNonManifoldPoints();


    // Merge all baffles that are still remaining after duplicating points.
    List<labelPair> couples
    (
        meshRefiner.getDuplicateFaces   // get all baffles
        (
            identity(mesh_.nFaces()-mesh_.nInternalFaces())
           +mesh_.nInternalFaces()
        )
    );

    label nCouples = returnReduce(couples.size(), sumOp<label>());

    Info<< "Detected unsplittable baffles : "
        << nCouples << endl;

    if (nCouples > 0)
    {
        // Actually merge baffles. Note: not exactly parallellized. Should
        // convert baffle faces into processor faces if they resulted
        // from them.
        meshRefiner.mergeBaffles(couples);

        if (debug_)
        {
            // Debug:test all is still synced across proc patches
            meshRefiner.checkData();
        }
        Info<< "Merged free-standing baffles in = "
            << mesh_.time().cpuTimeIncrement() << " s." << endl;
    }

    if (debug_)
    {
        Pout<< "Writing handleProblemCells mesh to time "
            << mesh_.time().timeName() << '.' << endl;
        meshRefiner.write(debug_, mesh_.time().path()/mesh_.time().timeName());
    }
}


void Foam::autoHexMeshDriver::mergePatchFaces()
{
    Info<< nl
        << "Merge refined boundary faces" << nl
        << "----------------------------" << nl
        << endl;

    if (debug_)
    {
        const_cast<Time&>(mesh_.time())++;
    }

    meshRefinement& meshRefiner = meshRefinerPtr_();

    meshRefiner.mergePatchFaces
    (
        Foam::cos(45*mathematicalConstant::pi/180.0),
        Foam::cos(45*mathematicalConstant::pi/180.0),
        meshRefinement::addedPatches(globalToPatch_)
    );

    if (debug_)
    {
        meshRefiner.checkData();
    }

    meshRefiner.mergeEdges(Foam::cos(45*mathematicalConstant::pi/180.0));

    if (debug_)
    {
        meshRefiner.checkData();
    }
}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::autoHexMeshDriver::balance
(
    const bool keepZoneFaces,
    const bool keepBaffles
)
{
    autoPtr<mapDistributePolyMesh> map;

    if (Pstream::parRun())
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;

        if (debug_)
        {
            const_cast<Time&>(mesh_.time())++;
        }

        meshRefinement& meshRefiner = meshRefinerPtr_();
        decompositionMethod& decomposer = decomposerPtr_();
        fvMeshDistribute& distributor = distributorPtr_();

        // Wanted distribution
        labelList distribution;

        if (keepZoneFaces || keepBaffles)
        {
            // Faces where owner and neighbour are not 'connected' so can
            // go to different processors.
            boolList blockedFace(mesh_.nFaces(), true);
            // Pairs of baffles
            List<labelPair> couples;

            if (keepZoneFaces)
            {
                label nNamed = surfaces().getNamedSurfaces().size();

                Info<< "Found " << nNamed << " surfaces with faceZones."
                    << " Applying special decomposition to keep those together."
                    << endl;

                // Determine decomposition to keep/move surface zones
                // on one processor. The reason is that snapping will make these
                // into baffles, move and convert them back so if they were
                // proc boundaries after baffling&moving the points might be no
                // longer snychronised so recoupling will fail. To prevent this
                // keep owner&neighbour of such a surface zone on the same
                // processor.

                const wordList& fzNames = surfaces().faceZoneNames();
                const faceZoneMesh& fZones = mesh_.faceZones();

                // Get faces whose owner and neighbour should stay together, i.e.
                // they are not 'blocked'.

                label nZoned = 0;

                forAll(fzNames, surfI)
                {
                    if (fzNames[surfI].size() > 0)
                    {
                        // Get zone
                        label zoneI = fZones.findZoneID(fzNames[surfI]);

                        const faceZone& fZone = fZones[zoneI];

                        forAll(fZone, i)
                        {
                            if (blockedFace[fZone[i]])
                            {
                                blockedFace[fZone[i]] = false;
                                nZoned++;
                            }
                        }
                    }
                }
                Info<< "Found " << returnReduce(nZoned, sumOp<label>())
                    << " zoned faces to keep together."
                    << endl;
            }

            if (keepBaffles)
            {
                // Get boundary baffles that need to stay together.
                couples =  meshRefiner.getDuplicateFaces   // all baffles
                (
                    identity(mesh_.nFaces()-mesh_.nInternalFaces())
                   +mesh_.nInternalFaces()
                );

                Info<< "Found " << returnReduce(couples.size(), sumOp<label>())
                    << " baffles to keep together."
                    << endl;
            }

            distribution = meshRefiner.decomposeCombineRegions
            (
                blockedFace,
                couples,
                decomposer
            );

            labelList nProcCells(distributor.countCells(distribution));
            Pstream::listCombineGather(nProcCells, plusEqOp<label>());
            Pstream::listCombineScatter(nProcCells);

            Info<< "Calculated decomposition:" << endl;
            forAll(nProcCells, procI)
            {
                Info<< "    " << procI << '\t' << nProcCells[procI] << endl;
            }
            Info<< endl;
        }
        else
        {
            // Normal decomposition
            distribution = decomposer.decompose(mesh_.cellCentres());
        }

        if (debug_)
        {
            Pout<< "Wanted distribution:"
                << distributor.countCells(distribution)
                << endl;
        }
        // Do actual sending/receiving of mesh
        map = distributor.distribute(distribution);

        // Update numbering of meshRefiner
        meshRefiner.distribute(map);
    }
    return map;
}


void Foam::autoHexMeshDriver::writeMesh(const string& msg) const
{
    const meshRefinement& meshRefiner = meshRefinerPtr_();

    meshRefiner.printMeshInfo(debug_, msg);
    Info<< "Writing mesh to time " << mesh_.time().timeName() << endl;

    meshRefiner.write(meshRefinement::MESH|meshRefinement::SCALARLEVELS, "");
    if (debug_ & meshRefinement::OBJINTERSECTIONS)
    {
        meshRefiner.write
        (
            meshRefinement::OBJINTERSECTIONS,
            mesh_.time().path()/mesh_.time().timeName()
        );
    }
    Info<< "Written mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s." << endl;
}


void Foam::autoHexMeshDriver::doMesh()
{
    Switch doRefine(dict_.lookup("doRefine"));
    Switch doSnap(dict_.lookup("doSnap"));
    Switch doLayers(dict_.lookup("doLayers"));

    Info<< "Do refinement : " << doRefine << nl
        << "Do snapping   : " << doSnap << nl
        << "Do layers     : " << doLayers << nl
        << endl;

    if (doRefine)
    {
        Info<< nl
            << "Refinement phase" << nl
            << "----------------" << nl
            << endl;

        const_cast<Time&>(mesh_.time())++;

        // Refine around feature edges
        featureEdgeRefine
        (
            100,    // maxIter
            0       // min cells to refine
        );

        // Refine based on surface
        surfaceOnlyRefine
        (
            100     // maxIter
        );

        // Remove cells (a certain distance) beyond surface intersections
        removeInsideCells
        (
            1       // nBufferLayers
        );

        // Internal mesh refinement
        shellRefine
        (
            100    // maxIter
        );

        // Introduce baffles at surface intersections
        baffleAndSplitMesh(doSnap);

        // Mesh is at its finest. Do optional zoning.
        zonify();

        // Pull baffles apart
        splitAndMergeBaffles(doSnap);

        // Do something about cells with refined faces on the boundary
        if (doSnap)
        {
            mergePatchFaces();
        }

        // Do final balancing. Keep zoned faces on one processor.
        balance(true, false);

        // Write mesh
        writeMesh("Refined mesh");
    }


    if (doSnap)
    {
        Info<< nl
            << "Morphing phase" << nl
            << "--------------" << nl
            << endl;

        const_cast<Time&>(mesh_.time())++;

        // Get the labels of added patches.
        labelList adaptPatchIDs(getSurfacePatches());

        // Create baffles (pairs of faces that share the same points)
        // Baffles stored as owner and neighbour face that have been created.
        List<labelPair> baffles;
        createZoneBaffles(baffles);

        {
            autoPtr<indirectPrimitivePatch> ppPtr
            (
                meshRefinement::makePatch
                (
                    mesh_,
                    adaptPatchIDs
                )
            );
            indirectPrimitivePatch& pp = ppPtr();

            // Distance to attact to nearest feature on surface
            const scalarField snapDist(calcSnapDistance(pp));


            // Construct iterative mesh mover.
            Info<< "Constructing mesh displacer ..." << endl;
            const dictionary& motionDict = dict_.subDict("motionDict");
            Info<< "Using mesh parameters " << motionDict << nl << endl;

            pointMesh pMesh(mesh_);

            motionSmoother meshMover
            (
                mesh_,
                pp,
                adaptPatchIDs,
                meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs),
                motionDict
            );


            // Check initial mesh
            Info<< "Checking initial mesh ..." << endl;
            labelHashSet wrongFaces(mesh_.nFaces()/100);
            motionSmoother::checkMesh(false, mesh_, motionDict, wrongFaces);
            const label nInitErrors = returnReduce
            (
                wrongFaces.size(),
                sumOp<label>()
            );

            Info<< "Detected " << nInitErrors << " illegal faces"
                << " (concave, zero area or negative cell pyramid volume)"
                << endl;


            Info<< "Checked initial mesh in = "
                << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;

            // Pre-smooth patch vertices (so before determining nearest)
            preSmoothPatch(nInitErrors, baffles, meshMover);

            // Calculate displacement at every patch point. Insert into
            // meshMover.
            calcNearestSurface(snapDist, meshMover);

            // Get smoothly varying internal displacement field.
            smoothDisplacement(meshMover);

            // Apply internal displacement to mesh.
            scaleMesh(nInitErrors, baffles, meshMover);
        }

        // Merge any introduced baffles.
        mergeZoneBaffles(baffles);

        // Write mesh.
        writeMesh("Snapped mesh");
    }


    if (doLayers)
    {
        Info<< nl
            << "Shrinking and layer addition phase" << nl
            << "----------------------------------" << nl
            << endl;

        const_cast<Time&>(mesh_.time())++;

        // Merge coplanar boundary faces
        mergePatchFacesUndo();

        // Per global region the number of layers (0 if no layer)
        labelList nLayers(readNumLayers());

        // Patches that need to get a layer
        DynamicList<label> patchIDs(nLayers.size());
        label nFacesWithLayers = 0;
        forAll(nLayers, region)
        {
            if (nLayers[region] > 0 && globalToPatch()[region] != -1)
            {
                label patchI = globalToPatch()[region];
                patchIDs.append(patchI);
                nFacesWithLayers += mesh_.boundaryMesh()[patchI].size();
            }
        }
        patchIDs.shrink();

        if (returnReduce(nFacesWithLayers, sumOp<label>()) == 0)
        {
            Info<< nl << "No layers to generate ..." << endl;
        }
        else
        {
            autoPtr<indirectPrimitivePatch> ppPtr
            (
                meshRefinement::makePatch
                (
                    mesh_,
                    patchIDs
                )
            );
            indirectPrimitivePatch& pp = ppPtr();

            // Construct iterative mesh mover.
            Info<< "Constructing mesh displacer ..." << endl;
            const dictionary& motionDict = dict_.subDict("motionDict");
            Info<< "Using mesh parameters " << motionDict << nl << endl;

            {
                pointMesh pMesh(mesh_);

                motionSmoother meshMover
                (
                    mesh_,
                    pp,
                    patchIDs,
                    meshRefinement::makeDisplacementField(pMesh, patchIDs),
                    motionDict
                );

                // Check that outside of mesh is not multiply connected.
                checkMeshManifold();

                // Check initial mesh
                Info<< "Checking initial mesh ..." << endl;
                labelHashSet wrongFaces(mesh_.nFaces()/100);
                motionSmoother::checkMesh(false, mesh_, motionDict, wrongFaces);
                const label nInitErrors = returnReduce
                (
                    wrongFaces.size(),
                    sumOp<label>()
                );

                Info<< "Detected " << nInitErrors << " illegal faces"
                    << " (concave, zero area or negative cell pyramid volume)"
                    << endl;


                // Do all topo changes
                addLayers(nInitErrors, meshMover);
            }

            // Write mesh.
            writeMesh("Layer mesh");
        }
    }
}


// ************************************************************************* //
