/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

Application
    snappyHexMesh

Group
    grpMeshGenerationUtilities

Description
    Automatic split hex mesher. Refines and snaps to surface.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "snappyRefineDriver.H"
#include "snappySnapDriver.H"
#include "snappyLayerDriver.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "decompositionMethod.H"
#include "noDecomp.H"
#include "fvMeshDistribute.H"
#include "wallPolyPatch.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "layerParameters.H"
#include "vtkSetWriter.H"
#include "faceSet.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "uindirectPrimitivePatch.H"
#include "surfZoneIdentifierList.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurface.H"
#include "globalIndex.H"
#include "IOmanip.H"
#include "decompositionModel.H"
#include "fvMeshTools.H"
#include "profiling.H"
#include "processorMeshes.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convert size (as fraction of defaultCellSize) to refinement level
label sizeCoeffToRefinement
(
    const scalar level0Coeff,   // ratio of hex cell size v.s. defaultCellSize
    const scalar sizeCoeff
)
{
     return round(::log(level0Coeff/sizeCoeff)/::log(2));
}


autoPtr<refinementSurfaces> createRefinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& surfacesDict,
    const dictionary& shapeControlDict,
    const label gapLevelIncrement,
    const scalar level0Coeff
)
{
    autoPtr<refinementSurfaces> surfacePtr;

    // Count number of surfaces.
    label surfi = 0;
    forAll(allGeometry.names(), geomi)
    {
        const word& geomName = allGeometry.names()[geomi];

        if (surfacesDict.found(geomName))
        {
            surfi++;
        }
    }

    labelList surfaces(surfi);
    wordList names(surfi);
    PtrList<surfaceZonesInfo> surfZones(surfi);

    labelList regionOffset(surfi);

    labelList globalMinLevel(surfi, 0);
    labelList globalMaxLevel(surfi, 0);
    labelList globalLevelIncr(surfi, 0);
    PtrList<dictionary> globalPatchInfo(surfi);
    List<Map<label>> regionMinLevel(surfi);
    List<Map<label>> regionMaxLevel(surfi);
    List<Map<label>> regionLevelIncr(surfi);
    List<Map<scalar>> regionAngle(surfi);
    List<Map<autoPtr<dictionary>>> regionPatchInfo(surfi);

    HashSet<word> unmatchedKeys(surfacesDict.toc());

    surfi = 0;
    forAll(allGeometry.names(), geomi)
    {
        const word& geomName = allGeometry.names()[geomi];

        const entry* ePtr = surfacesDict.lookupEntryPtr(geomName, false, true);

        if (ePtr)
        {
            const dictionary& shapeDict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            names[surfi] = geomName;
            surfaces[surfi] = geomi;

            const searchableSurface& surface = allGeometry[geomi];

            // Find the index in shapeControlDict
            // Invert surfaceCellSize to get the refinementLevel

            const word scsFuncName =
                shapeDict.lookup("surfaceCellSizeFunction");
            const dictionary& scsDict =
                shapeDict.optionalSubDict(scsFuncName + "Coeffs");

            const scalar surfaceCellSize =
                readScalar(scsDict.lookup("surfaceCellSizeCoeff"));

            const label refLevel = sizeCoeffToRefinement
            (
                level0Coeff,
                surfaceCellSize
            );

            globalMinLevel[surfi] = refLevel;
            globalMaxLevel[surfi] = refLevel;
            globalLevelIncr[surfi] = gapLevelIncrement;

            // Surface zones
            surfZones.set(surfi, new surfaceZonesInfo(surface, shapeDict));


            // Global perpendicular angle
            if (shapeDict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfi,
                    shapeDict.subDict("patchInfo").clone()
                );
            }


            // Per region override of patchInfo

            if (shapeDict.found("regions"))
            {
                const dictionary& regionsDict = shapeDict.subDict("regions");
                const wordList& regionNames =
                    allGeometry[surfaces[surfi]].regions();

                forAll(regionNames, regioni)
                {
                    if (regionsDict.found(regionNames[regioni]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regioni]
                        );

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfi].insert
                            (
                                regioni,
                                regionDict.subDict("patchInfo").clone()
                            );
                        }
                    }
                }
            }

            // Per region override of cellSize
            if (shapeDict.found("regions"))
            {
                const dictionary& shapeControlRegionsDict =
                    shapeDict.subDict("regions");
                const wordList& regionNames =
                    allGeometry[surfaces[surfi]].regions();

                forAll(regionNames, regioni)
                {
                    if (shapeControlRegionsDict.found(regionNames[regioni]))
                    {
                        const dictionary& shapeControlRegionDict =
                            shapeControlRegionsDict.subDict
                            (
                                regionNames[regioni]
                            );

                        const word scsFuncName =
                            shapeControlRegionDict.lookup
                            (
                                "surfaceCellSizeFunction"
                            );
                        const dictionary& scsDict =
                            shapeControlRegionDict.subDict
                            (
                                scsFuncName + "Coeffs"
                            );

                        const scalar surfaceCellSize =
                            readScalar
                            (
                                scsDict.lookup("surfaceCellSizeCoeff")
                            );

                        const label refLevel = sizeCoeffToRefinement
                        (
                            level0Coeff,
                            surfaceCellSize
                        );

                        regionMinLevel[surfi].insert(regioni, refLevel);
                        regionMaxLevel[surfi].insert(regioni, refLevel);
                        regionLevelIncr[surfi].insert(regioni, 0);
                    }
                }
            }

            surfi++;
        }
    }

    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces, surfi)
    {
        regionOffset[surfi] = nRegions;
        nRegions += allGeometry[surfaces[surfi]].regions().size();
    }

    // Rework surface specific information into information per global region
    labelList minLevel(nRegions, 0);
    labelList maxLevel(nRegions, 0);
    labelList gapLevel(nRegions, -1);
    PtrList<dictionary> patchInfo(nRegions);

    forAll(globalMinLevel, surfi)
    {
        label nRegions = allGeometry[surfaces[surfi]].regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegioni = regionOffset[surfi] + i;
            minLevel[globalRegioni] = globalMinLevel[surfi];
            maxLevel[globalRegioni] = globalMaxLevel[surfi];
            gapLevel[globalRegioni] =
                maxLevel[globalRegioni]
              + globalLevelIncr[surfi];

            if (globalPatchInfo.set(surfi))
            {
                patchInfo.set
                (
                    globalRegioni,
                    globalPatchInfo[surfi].clone()
                );
            }
        }

        // Overwrite with region specific information
        forAllConstIter(Map<label>, regionMinLevel[surfi], iter)
        {
            label globalRegioni = regionOffset[surfi] + iter.key();

            minLevel[globalRegioni] = iter();
            maxLevel[globalRegioni] = regionMaxLevel[surfi][iter.key()];
            gapLevel[globalRegioni] =
                maxLevel[globalRegioni]
              + regionLevelIncr[surfi][iter.key()];
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfi];
        forAllConstIter(Map<autoPtr<dictionary>>, localInfo, iter)
        {
            label globalRegioni = regionOffset[surfi] + iter.key();
            patchInfo.set(globalRegioni, iter()().clone());
        }
    }

    surfacePtr.reset
    (
        new refinementSurfaces
        (
            allGeometry,
            surfaces,
            names,
            surfZones,
            regionOffset,
            minLevel,
            maxLevel,
            gapLevel,
            scalarField(nRegions, -GREAT),  //perpendicularAngle,
            patchInfo
        )
    );


    const refinementSurfaces& rf = surfacePtr();

    // Determine maximum region name length
    label maxLen = 0;
    forAll(rf.surfaces(), surfi)
    {
        label geomi = rf.surfaces()[surfi];
        const wordList& regionNames = allGeometry.regionNames()[geomi];
        forAll(regionNames, regioni)
        {
            maxLen = Foam::max(maxLen, label(regionNames[regioni].size()));
        }
    }


    Info<< setw(maxLen) << "Region"
        << setw(10) << "Min Level"
        << setw(10) << "Max Level"
        << setw(10) << "Gap Level" << nl
        << setw(maxLen) << "------"
        << setw(10) << "---------"
        << setw(10) << "---------"
        << setw(10) << "---------" << endl;

    forAll(rf.surfaces(), surfi)
    {
        label geomi = rf.surfaces()[surfi];

        Info<< rf.names()[surfi] << ':' << nl;

        const wordList& regionNames = allGeometry.regionNames()[geomi];

        forAll(regionNames, regioni)
        {
            label globali = rf.globalRegion(surfi, regioni);

            Info<< setw(maxLen) << regionNames[regioni]
                << setw(10) << rf.minLevel()[globali]
                << setw(10) << rf.maxLevel()[globali]
                << setw(10) << rf.gapLevel()[globali] << endl;
        }
    }


    return surfacePtr;
}


void extractSurface
(
    const polyMesh& mesh,
    const Time& runTime,
    const labelHashSet& includePatches,
    const fileName& outFileName
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // Collect sizes. Hash on names to handle local-only patches (e.g.
    //  processor patches)
    HashTable<label> patchSize(1024);
    label nFaces = 0;
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        patchSize.insert(pp.name(), pp.size());
        nFaces += pp.size();
    }
    Pstream::mapCombineGather(patchSize, plusEqOp<label>());


    // Allocate zone/patch for all patches
    HashTable<label> compactZoneID(1024);
    forAllConstIter(HashTable<label>, patchSize, iter)
    {
        label sz = compactZoneID.size();
        compactZoneID.insert(iter.key(), sz);
    }
    Pstream::mapCombineScatter(compactZoneID);


    // Rework HashTable into labelList just for speed of conversion
    labelList patchToCompactZone(bMesh.size(), -1);
    forAllConstIter(HashTable<label>, compactZoneID, iter)
    {
        label patchi = bMesh.findPatchID(iter.key());
        if (patchi != -1)
        {
            patchToCompactZone[patchi] = iter();
        }
    }

    // Collect faces on zones
    DynamicList<label> faceLabels(nFaces);
    DynamicList<label> compactZones(nFaces);
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        forAll(pp, i)
        {
            faceLabels.append(pp.start()+i);
            compactZones.append(patchToCompactZone[pp.index()]);
        }
    }

    // Addressing engine for all faces
    uindirectPrimitivePatch allBoundary
    (
        UIndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );


    // Find correspondence to master points
    labelList pointToGlobal;
    labelList uniqueMeshPoints;
    autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
    (
        allBoundary.meshPoints(),
        allBoundary.meshPointMap(),
        pointToGlobal,
        uniqueMeshPoints
    );

    // Gather all unique points on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()] = pointField
    (
        mesh.points(),
        uniqueMeshPoints
    );
    Pstream::gatherList(gatheredPoints);

    // Gather all faces
    List<faceList> gatheredFaces(Pstream::nProcs());
    gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
    forAll(gatheredFaces[Pstream::myProcNo()], i)
    {
        inplaceRenumber(pointToGlobal, gatheredFaces[Pstream::myProcNo()][i]);
    }
    Pstream::gatherList(gatheredFaces);

    // Gather all ZoneIDs
    List<labelList> gatheredZones(Pstream::nProcs());
    gatheredZones[Pstream::myProcNo()] = compactZones.xfer();
    Pstream::gatherList(gatheredZones);

    // On master combine all points, faces, zones
    if (Pstream::master())
    {
        pointField allPoints = ListListOps::combine<pointField>
        (
            gatheredPoints,
            accessOp<pointField>()
        );
        gatheredPoints.clear();

        faceList allFaces = ListListOps::combine<faceList>
        (
            gatheredFaces,
            accessOp<faceList>()
        );
        gatheredFaces.clear();

        labelList allZones = ListListOps::combine<labelList>
        (
            gatheredZones,
            accessOp<labelList>()
        );
        gatheredZones.clear();


        // Zones
        surfZoneIdentifierList surfZones(compactZoneID.size());
        forAllConstIter(HashTable<label>, compactZoneID, iter)
        {
            surfZones[iter()] = surfZoneIdentifier(iter.key(), iter());
            Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
                << endl;
        }

        UnsortedMeshedSurface<face> unsortedFace
        (
            xferMove(allPoints),
            xferMove(allFaces),
            xferMove(allZones),
            xferMove(surfZones)
        );


        MeshedSurface<face> sortedFace(unsortedFace);

        fileName globalCasePath
        (
            runTime.processorCase()
          ? runTime.path()/".."/outFileName
          : runTime.path()/outFileName
        );
        globalCasePath.clean();

        Info<< "Writing merged surface to " << globalCasePath << endl;

        sortedFace.write(globalCasePath);
    }
}


// Check writing tolerance before doing any serious work
scalar getMergeDistance(const polyMesh& mesh, const scalar mergeTol)
{
    const boundBox& meshBb = mesh.bounds();
    scalar mergeDist = mergeTol * meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    // check writing tolerance
    if (mesh.time().writeFormat() == IOstream::ASCII)
    {
        const scalar writeTol = std::pow
        (
            scalar(10.0),
            -scalar(IOstream::defaultPrecision())
        );

        if (mergeTol < writeTol)
        {
            FatalErrorInFunction
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << nl
                << "Your merging tolerance (" << mergeTol
                << ") is finer than this." << nl
                << "Change to binary writeFormat, "
                << "or increase the writePrecision" << endl
                << "or adjust the merge tolerance (mergeTol)."
                << exit(FatalError);
        }
    }

    return mergeDist;
}


void removeZeroSizedPatches(fvMesh& mesh)
{
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduce(pp.size(), sumOp<label>())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
        Info<< endl;
    }
}


// Write mesh and additional information
void writeMesh
(
    const string& msg,
    const meshRefinement& meshRefiner,
    const meshRefinement::debugType debugLevel,
    const meshRefinement::writeType writeLevel
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    meshRefiner.printMeshInfo(debugLevel, msg);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    processorMeshes::removeFiles(mesh);
    if (!debugLevel && !(writeLevel&meshRefinement::WRITELAYERSETS))
    {
        topoSet::removeFiles(mesh);
    }
    refinementHistory::removeFiles(mesh);

    meshRefiner.write
    (
        debugLevel,
        meshRefinement::writeType(writeLevel | meshRefinement::WRITEMESH),
        mesh.time().path()/meshRefiner.timeName()
    );
    Info<< "Wrote mesh in = "
        << mesh.time().cpuTimeIncrement() << " s." << endl;
}


int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "addOverwriteOption.H"
    Foam::argList::addBoolOption
    (
        "checkGeometry",
        "check all surface geometry for quality"
    );
    Foam::argList::addOption
    (
        "surfaceSimplify",
        "boundBox",
        "simplify the surface using snappyHexMesh starting from a boundBox"
    );
    Foam::argList::addOption
    (
        "patches",
        "(patch0 .. patchN)",
        "only triangulate selected patches (wildcards supported)"
    );
    Foam::argList::addOption
    (
        "outFile",
        "file",
        "name of the file to save the simplified surface to"
    );
    #include "addProfilingOption.H"
    #include "addDictOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    const bool overwrite = args.optionFound("overwrite");
    const bool checkGeometry = args.optionFound("checkGeometry");
    const bool surfaceSimplify = args.optionFound("surfaceSimplify");

    autoPtr<fvMesh> meshPtr;

    {
        word regionName;
        if (args.optionReadIfPresent("region", regionName))
        {
            Info<< "Create mesh " << regionName << " for time = "
                << runTime.timeName() << nl << endl;
        }
        else
        {
            regionName = fvMesh::defaultRegion;
            Info<< "Create mesh for time = "
                << runTime.timeName() << nl << endl;
        }

        meshPtr.reset
        (
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }

    fvMesh& mesh = meshPtr();

    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    // Check patches and faceZones are synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    meshRefinement::checkCoupledFaceZones(mesh);


    // Read meshing dictionary
    const word dictName("snappyHexMeshDict");
    #include "setSystemMeshDictionaryIO.H"
    const IOdictionary meshDict(dictIO);


    // all surface geometry
    const dictionary& geometryDict = meshDict.subDict("geometry");

    // refinement parameters
    const dictionary& refineDict = meshDict.subDict("castellatedMeshControls");

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    // snap-to-surface parameters
    const dictionary& snapDict = meshDict.subDict("snapControls");

    // layer addition parameters
    const dictionary& layerDict = meshDict.subDict("addLayersControls");

    // absolute merge distance
    const scalar mergeDist = getMergeDistance
    (
        mesh,
        readScalar(meshDict.lookup("mergeTolerance"))
    );

    const Switch keepPatches(meshDict.lookupOrDefault("keepPatches", false));


    // Read decomposePar dictionary
    dictionary decomposeDict;
    if (Pstream::parRun())
    {
        fileName decompDictFile;
        args.optionReadIfPresent("decomposeParDict", decompDictFile);

        // A demand-driven decompositionMethod can have issues finding
        // an alternative decomposeParDict location.

        IOdictionary* dictPtr = new IOdictionary
        (
            decompositionModel::selectIO
            (
                IOobject
                (
                    "decomposeParDict",
                    runTime.system(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                decompDictFile
            )
        );

        // Store it on the object registry, but to be found it must also
        // have the expected "decomposeParDict" name.

        dictPtr->rename("decomposeParDict");
        runTime.store(dictPtr);

        decomposeDict = *dictPtr;
    }
    else
    {
        decomposeDict.add("method", "none");
        decomposeDict.add("numberOfSubdomains", 1);
    }


    // Debug
    // ~~~~~

    // Set debug level
    meshRefinement::debugType debugLevel = meshRefinement::debugType
    (
        meshDict.lookupOrDefault<label>
        (
            "debug",
            0
        )
    );
    {
        wordList flags;
        if (meshDict.readIfPresent("debugFlags", flags))
        {
            debugLevel = meshRefinement::debugType
            (
                meshRefinement::readFlags
                (
                    meshRefinement::debugTypeNames,
                    flags
                )
            );
        }
    }
    if (debugLevel > 0)
    {
        meshRefinement::debug   = debugLevel;
        snappyRefineDriver::debug = debugLevel;
        snappySnapDriver::debug   = debugLevel;
        snappyLayerDriver::debug  = debugLevel;
    }

    // Set file writing level
    {
        wordList flags;
        if (meshDict.readIfPresent("writeFlags", flags))
        {
            meshRefinement::writeLevel
            (
                meshRefinement::writeType
                (
                    meshRefinement::readFlags
                    (
                        meshRefinement::writeTypeNames,
                        flags
                    )
                )
            );
        }
    }

    // Set output level
    {
        wordList flags;
        if (meshDict.readIfPresent("outputFlags", flags))
        {
            meshRefinement::outputLevel
            (
                meshRefinement::outputType
                (
                    meshRefinement::readFlags
                    (
                        meshRefinement::outputTypeNames,
                        flags
                    )
                )
            );
        }
    }

    // for the impatient who want to see some output files:
    profiling::writeNow();

    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            //mesh.time().findInstance("triSurface", word::null),// instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        meshDict.lookupOrDefault("singleRegionName", true)
    );


    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<refinementSurfaces> surfacesPtr;

    Info<< "Reading refinement surfaces." << endl;

    if (surfaceSimplify)
    {
        addProfiling(surfaceSimplify, "snappyHexMesh::surfaceSimplify");
        IOdictionary foamyHexMeshDict
        (
           IOobject
           (
                "foamyHexMeshDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
           )
        );

        const dictionary& conformationDict =
            foamyHexMeshDict.subDict("surfaceConformation").subDict
            (
                "geometryToConformTo"
            );

        const dictionary& motionDict =
            foamyHexMeshDict.subDict("motionControl");

        const dictionary& shapeControlDict =
            motionDict.subDict("shapeControlFunctions");

        // Calculate current ratio of hex cells v.s. wanted cell size
        const scalar defaultCellSize =
            readScalar(motionDict.lookup("defaultCellSize"));

        const scalar initialCellSize = ::pow(meshPtr().V()[0], 1.0/3.0);

        //Info<< "Wanted cell size  = " << defaultCellSize << endl;
        //Info<< "Current cell size = " << initialCellSize << endl;
        //Info<< "Fraction          = " << initialCellSize/defaultCellSize
        //    << endl;

        surfacesPtr =
            createRefinementSurfaces
            (
                allGeometry,
                conformationDict,
                shapeControlDict,
                refineDict.lookupOrDefault("gapLevelIncrement", 0),
                initialCellSize/defaultCellSize
            );

        profiling::writeNow();
    }
    else
    {
        surfacesPtr.reset
        (
            new refinementSurfaces
            (
                allGeometry,
                refineDict.subDict("refinementSurfaces"),
                refineDict.lookupOrDefault("gapLevelIncrement", 0)
            )
        );

        Info<< "Read refinement surfaces in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }

    refinementSurfaces& surfaces = surfacesPtr();


    // Checking only?

    if (checkGeometry)
    {
        // Extract patchInfo
        List<wordList> patchTypes(allGeometry.size());

        const PtrList<dictionary>& patchInfo = surfaces.patchInfo();
        const labelList& surfaceGeometry = surfaces.surfaces();
        forAll(surfaceGeometry, surfi)
        {
            label geomi = surfaceGeometry[surfi];
            const wordList& regNames = allGeometry.regionNames()[geomi];

            patchTypes[geomi].setSize(regNames.size());
            forAll(regNames, regioni)
            {
                label globalRegioni = surfaces.globalRegion(surfi, regioni);

                if (patchInfo.set(globalRegioni))
                {
                    patchTypes[geomi][regioni] =
                        word(patchInfo[globalRegioni].lookup("type"));
                }
                else
                {
                    patchTypes[geomi][regioni] = wallPolyPatch::typeName;
                }
            }
        }

        // Write some stats
        allGeometry.writeStats(patchTypes, Info);
        // Check topology
        allGeometry.checkTopology(true);
        // Check geometry
        allGeometry.checkGeometry
        (
            100.0,      // max size ratio
            1e-9,       // intersection tolerance
            autoPtr<writer<scalar>>(new vtkSetWriter<scalar>()),
            0.01,       // min triangle quality
            true
        );

        return 0;
    }



    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement shells." << endl;
    shellSurfaces shells
    (
        allGeometry,
        refineDict.subDict("refinementRegions")
    );
    Info<< "Read refinement shells in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    Info<< "Setting refinement level of surface to be consistent"
        << " with shells." << endl;
    surfaces.setMinLevelFields(shells);
    Info<< "Checked shell refinement in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    // Optionally read limit shells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const dictionary limitDict(refineDict.subOrEmptyDict("limitRegions"));

    if (!limitDict.empty())
    {
        Info<< "Reading limit shells." << endl;
    }

    shellSurfaces limitShells(allGeometry, limitDict);

    if (!limitDict.empty())
    {
        Info<< "Read refinement shells in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }



    // Read feature meshes
    // ~~~~~~~~~~~~~~~~~~~

    Info<< "Reading features." << endl;
    refinementFeatures features
    (
        mesh,
        refineDict.lookup("features")
    );
    Info<< "Read features in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;



    // Refinement engine
    // ~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Determining initial surface intersections" << nl
        << "-----------------------------------------" << nl
        << endl;

    // Main refinement engine
    meshRefinement meshRefiner
    (
        mesh,
        mergeDist,          // tolerance used in sorting coordinates
        overwrite,          // overwrite mesh files?
        surfaces,           // for surface intersection refinement
        features,           // for feature edges/point based refinement
        shells,             // for volume (inside/outside) refinement
        limitShells         // limit of volume refinement
    );
    Info<< "Calculated surface intersections in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;

    // Some stats
    meshRefiner.printMeshInfo(debugLevel, "Initial mesh");

    meshRefiner.write
    (
        meshRefinement::debugType(debugLevel&meshRefinement::OBJINTERSECTIONS),
        meshRefinement::writeType(0),
        mesh.time().path()/meshRefiner.timeName()
    );


    // Refinement parameters
    const refinementParameters refineParams(refineDict);

    // Snap parameters
    const snapParameters snapParams(snapDict);



    // Add all the cellZones and faceZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. cellZones relating to surface (faceZones added later)

    const labelList namedSurfaces
    (
        surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones())
    );

    labelList surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
    (
        surfaces.surfZones(),
        namedSurfaces,
        mesh
    );


    // 2. cellZones relating to locations

    refineParams.addCellZonesToMesh(mesh);



    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Global surface region to patch (non faceZone surface) or patches
    //  (faceZone surfaces)
    labelList globalToMasterPatch;
    labelList globalToSlavePatch;


    {
        Info<< nl
            << "Adding patches for surface regions" << nl
            << "----------------------------------" << nl
            << endl;

        // From global region number to mesh patch.
        globalToMasterPatch.setSize(surfaces.nRegions(), -1);
        globalToSlavePatch.setSize(surfaces.nRegions(), -1);

        Info<< setf(ios_base::left)
            << setw(6) << "Patch"
            << setw(20) << "Type"
            << setw(30) << "Region" << nl
            << setw(6) << "-----"
            << setw(20) << "----"
            << setw(30) << "------" << endl;

        const labelList& surfaceGeometry = surfaces.surfaces();
        const PtrList<dictionary>& surfacePatchInfo = surfaces.patchInfo();
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        forAll(surfaceGeometry, surfi)
        {
            label geomi = surfaceGeometry[surfi];

            const wordList& regNames = allGeometry.regionNames()[geomi];

            Info<< surfaces.names()[surfi] << ':' << nl << nl;

            const word& fzName = surfaces.surfZones()[surfi].faceZoneName();

            if (fzName.empty())
            {
                // 'Normal' surface
                forAll(regNames, i)
                {
                    label globalRegioni = surfaces.globalRegion(surfi, i);

                    label patchi;

                    if (surfacePatchInfo.set(globalRegioni))
                    {
                        patchi = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            surfacePatchInfo[globalRegioni]
                        );
                    }
                    else
                    {
                        dictionary patchInfo;
                        patchInfo.set("type", wallPolyPatch::typeName);

                        patchi = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            patchInfo
                        );
                    }

                    Info<< setf(ios_base::left)
                        << setw(6) << patchi
                        << setw(20) << pbm[patchi].type()
                        << setw(30) << regNames[i] << nl;

                    globalToMasterPatch[globalRegioni] = patchi;
                    globalToSlavePatch[globalRegioni] = patchi;
                }
            }
            else
            {
                // Zoned surface
                forAll(regNames, i)
                {
                    label globalRegioni = surfaces.globalRegion(surfi, i);

                    // Add master side patch
                    {
                        label patchi;

                        if (surfacePatchInfo.set(globalRegioni))
                        {
                            patchi = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                surfacePatchInfo[globalRegioni]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchi = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                patchInfo
                            );
                        }

                        Info<< setf(ios_base::left)
                            << setw(6) << patchi
                            << setw(20) << pbm[patchi].type()
                            << setw(30) << regNames[i] << nl;

                        globalToMasterPatch[globalRegioni] = patchi;
                    }
                    // Add slave side patch
                    {
                        const word slaveName = regNames[i] + "_slave";
                        label patchi;

                        if (surfacePatchInfo.set(globalRegioni))
                        {
                            patchi = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                surfacePatchInfo[globalRegioni]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchi = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                patchInfo
                            );
                        }

                        Info<< setf(ios_base::left)
                            << setw(6) << patchi
                            << setw(20) << pbm[patchi].type()
                            << setw(30) << slaveName << nl;

                        globalToSlavePatch[globalRegioni] = patchi;
                    }
                }

                // For now: have single faceZone per surface. Use first
                // region in surface for patch for zoneing
                if (regNames.size())
                {
                    label globalRegioni = surfaces.globalRegion(surfi, 0);

                    meshRefiner.addFaceZone
                    (
                        fzName,
                        pbm[globalToMasterPatch[globalRegioni]].name(),
                        pbm[globalToSlavePatch[globalRegioni]].name(),
                        surfaces.surfZones()[surfi].faceType()
                    );
                }
            }

            Info<< nl;
        }
        Info<< "Added patches in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }



    // Add all information for all the remaining faceZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HashTable<Pair<word>> faceZoneToPatches;
    forAll(mesh.faceZones(), zonei)
    {
        const word& fzName = mesh.faceZones()[zonei].name();

        label mpI, spI;
        surfaceZonesInfo::faceZoneType fzType;
        bool hasInfo = meshRefiner.getFaceZoneInfo(fzName, mpI, spI, fzType);

        if (!hasInfo)
        {
            // faceZone does not originate from a surface but presumably
            // from a cellZone pair instead
            string::size_type i = fzName.find("_to_");
            if (i != string::npos)
            {
                word cz0 = fzName.substr(0, i);
                word cz1 = fzName.substr(i+4, fzName.size()-i+4);
                word slaveName(cz1 + "_to_" + cz0);
                faceZoneToPatches.insert(fzName, Pair<word>(fzName, slaveName));
            }
            else
            {
                // Add as fzName + fzName_slave
                const word slaveName = fzName + "_slave";
                faceZoneToPatches.insert(fzName, Pair<word>(fzName, slaveName));
            }
        }
    }

    if (faceZoneToPatches.size())
    {
        snappyRefineDriver::addFaceZones
        (
            meshRefiner,
            refineParams,
            faceZoneToPatches
        );
    }



    // Re-do intersections on meshed boundaries since they use an extrapolated
    // other side
    {
        const labelList adaptPatchIDs(meshRefiner.meshedPatches());

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        label nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            nFaces += pbm[adaptPatchIDs[i]].size();
        }

        labelList faceLabels(nFaces);
        nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            const polyPatch& pp = pbm[adaptPatchIDs[i]];
            forAll(pp, i)
            {
                faceLabels[nFaces++] = pp.start()+i;
            }
        }
        meshRefiner.updateIntersections(faceLabels);
    }



    // Parallel
    // ~~~~~~~~

    // Construct decomposition engine. Note: cannot use decompositionModel
    // MeshObject since we're clearing out the mesh inside the mesh generation.
    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New
        (
            decomposeDict
        )
    );
    decompositionMethod& decomposer = decomposerPtr();

    if (Pstream::parRun() && !decomposer.parallelAware())
    {
        FatalErrorInFunction
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which is not parallel aware." << endl
            << "Please select one that is (hierarchical, ptscotch)"
            << exit(FatalError);
    }

    // Mesh distribution engine (uses tolerance to reconstruct meshes)
    fvMeshDistribute distributor(mesh, mergeDist);





    // Now do the real work -refinement -snapping -layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const Switch wantRefine(meshDict.lookup("castellatedMesh"));
    const Switch wantSnap(meshDict.lookup("snap"));
    const Switch wantLayers(meshDict.lookup("addLayers"));

    const Switch mergePatchFaces
    (
        meshDict.lookupOrDefault("mergePatchFaces", true)
    );

    if (!mergePatchFaces)
    {
        Info<< "Not merging patch-faces of cell to preserve"
            << " (split)hex cell shape."
            << nl << endl;
    }


    if (wantRefine)
    {
        cpuTime timer;

        snappyRefineDriver refineDriver
        (
            meshRefiner,
            decomposer,
            distributor,
            globalToMasterPatch,
            globalToSlavePatch
        );


        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }


        refineDriver.doRefine
        (
            refineDict,
            refineParams,
            snapParams,
            refineParams.handleSnapProblems(),
            mergePatchFaces,        // merge co-planar faces
            motionDict
        );

        // Remove zero sized patches originating from faceZones
        if (!keepPatches && !wantSnap && !wantLayers)
        {
            fvMeshTools::removeEmptyPatches(mesh, true);
        }

        writeMesh
        (
            "Refined mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Mesh refined in = "
            << timer.cpuTimeIncrement() << " s." << endl;

        profiling::writeNow();
    }

    if (wantSnap)
    {
        cpuTime timer;

        snappySnapDriver snapDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch
        );

        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Use the resolveFeatureAngle from the refinement parameters
        scalar curvature = refineParams.curvature();
        scalar planarAngle = refineParams.planarAngle();

        snapDriver.doSnap
        (
            snapDict,
            motionDict,
            mergePatchFaces,
            curvature,
            planarAngle,
            snapParams
        );

        // Remove zero sized patches originating from faceZones
        if (!keepPatches && !wantLayers)
        {
            fvMeshTools::removeEmptyPatches(mesh, true);
        }

        writeMesh
        (
            "Snapped mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Mesh snapped in = "
            << timer.cpuTimeIncrement() << " s." << endl;

        profiling::writeNow();
    }

    if (wantLayers)
    {
        cpuTime timer;

        // Layer addition parameters
        const layerParameters layerParams(layerDict, mesh.boundaryMesh());

        snappyLayerDriver layerDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Use the maxLocalCells from the refinement parameters
        bool preBalance = returnReduce
        (
            (mesh.nCells() >= refineParams.maxLocalCells()),
            orOp<bool>()
        );


        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }

        layerDriver.doLayers
        (
            layerDict,
            motionDict,
            layerParams,
            mergePatchFaces,
            preBalance,
            decomposer,
            distributor
        );

        // Remove zero sized patches originating from faceZones
        if (!keepPatches)
        {
            fvMeshTools::removeEmptyPatches(mesh, true);
        }

        writeMesh
        (
            "Layer mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Layers added in = "
            << timer.cpuTimeIncrement() << " s." << endl;

        profiling::writeNow();
    }


    {
        addProfiling(checkMesh, "snappyHexMesh::checkMesh");

        // Check final mesh
        Info<< "Checking final mesh ..." << endl;
        faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);
        const label nErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        if (nErrors > 0)
        {
            Info<< "Finished meshing with " << nErrors << " illegal faces"
                << " (concave, zero area or negative cell pyramid volume)"
                << endl;
            wrongFaces.write();
        }
        else
        {
            Info<< "Finished meshing without any errors" << endl;
        }

        profiling::writeNow();
    }


    if (surfaceSimplify)
    {
        addProfiling(surfaceSimplify, "snappyHexMesh::surfaceSimplify");

        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        labelHashSet includePatches(bMesh.size());

        if (args.optionFound("patches"))
        {
            includePatches = bMesh.patchSet
            (
                wordReList(args.optionLookup("patches")())
            );
        }
        else
        {
            forAll(bMesh, patchi)
            {
                const polyPatch& patch = bMesh[patchi];

                if (!isA<processorPolyPatch>(patch))
                {
                    includePatches.insert(patchi);
                }
            }
        }

        fileName outFileName
        (
            args.optionLookupOrDefault<fileName>
            (
                "outFile",
                "constant/triSurface/simplifiedSurface.stl"
            )
        );

        extractSurface
        (
            mesh,
            runTime,
            includePatches,
            outFileName
        );

        pointIOField cellCentres
        (
            IOobject
            (
                "internalCellCentres",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh.cellCentres()
        );

        cellCentres.write();
    }

    profiling::writeNow();

    Info<< "Finished meshing in = "
        << runTime.elapsedCpuTime() << " s." << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
