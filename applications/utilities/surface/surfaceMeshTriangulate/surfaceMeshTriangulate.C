/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    surfaceMeshTriangulate

Group
    grpSurfaceUtilities

Description
    Extract patch or faceZone surfaces from a polyMesh.
    Depending on output surface format triangulates faces.

    Region numbers on faces no guaranteed to be the same as the patch indices.

    Optionally only extracts named patches.

    If run in parallel, processor patches get filtered out by default and
    the mesh is merged (based on topology).

\*---------------------------------------------------------------------------*/

#include "MeshedSurface.H"
#include "UnsortedMeshedSurface.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "ListListOps.H"
#include "uindirectPrimitivePatch.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "timeSelector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Extract patch or faceZone surfaces from a polyMesh."
        " The name is historical, it only triangulates faces"
        " when the output format requires it."
    );
    timeSelector::addOptions();

    argList::addArgument("output", "The output surface file");

    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "excludeProcPatches",
        "Exclude processor patches"
    );
    argList::addOption
    (
        "patches",
        "wordRes"
        "Specify single patch or multiple patches to extract.\n"
        "Eg, 'top' or '( front \".*back\" )'"
    );
    argList::addOption
    (
        "faceZones",
        "wordRes",
        "Specify single or multiple faceZones to extract\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const fileName userOutFileName(args[1]);

    if (!userOutFileName.hasExt())
    {
        FatalErrorInFunction
            << "Missing extension on output name " << userOutFileName
            << exit(FatalError);
    }

    Info<< "Extracting surface from boundaryMesh ..." << nl << nl;

    const bool includeProcPatches =
       !(
            args.found("excludeProcPatches")
         || Pstream::parRun()
        );

    if (includeProcPatches)
    {
        Info<< "Including all processor patches." << nl << endl;
    }
    else if (Pstream::parRun())
    {
        Info<< "Excluding all processor patches." << nl << endl;
    }


    Info<< "Reading mesh from time " << runTime.value() << endl;

    #include "createNamedPolyMesh.H"

    // User specified times
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeIndex)
    {
        runTime.setTime(timeDirs[timeIndex], timeIndex);
        Info<< "Time [" << timeIndex << "] = " << runTime.timeName();

        fileName outFileName;
        if (timeDirs.size() == 1)
        {
            outFileName = userOutFileName;
        }
        else
        {
            polyMesh::readUpdateState meshState = mesh.readUpdate();
            if (timeIndex && meshState == polyMesh::UNCHANGED)
            {
                Info<<"  ... no mesh change." << nl;
                continue;
            }

            // The filename based on the original, but with additional
            // time information. The extension was previously checked that
            // it exists
            const auto dot = userOutFileName.rfind('.');

            outFileName =
                userOutFileName.substr(0, dot) + "_"
              + Foam::name(runTime.value()) + "."
              + userOutFileName.ext();
        }

        Info<< nl;

        // Create local surface from:
        // - explicitly named patches only (-patches (at your option)
        // - all patches (default in sequential mode)
        // - all non-processor patches (default in parallel mode)
        // - all non-processor patches (sequential mode, -excludeProcPatches
        //   (at your option)

        // Construct table of patches to include.
        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        labelList includePatches;

        if (args.found("patches"))
        {
            includePatches =
                bMesh.patchSet(args.getList<wordRe>("patches")).sortedToc();
        }
        else if (includeProcPatches)
        {
            includePatches = identity(bMesh.size());
        }
        else
        {
            includePatches = identity(bMesh.nNonProcessor());
        }


        labelList includeFaceZones;

        const faceZoneMesh& fzm = mesh.faceZones();

        if (args.found("faceZones"))
        {
            const wordList allZoneNames(fzm.names());

            const wordRes zoneNames(args.getList<wordRe>("faceZones"));

            labelHashSet hashed(2*fzm.size());

            for (const wordRe& zoneName : zoneNames)
            {
                labelList zoneIDs = findStrings(zoneName, allZoneNames);
                hashed.insert(zoneIDs);

                if (zoneIDs.empty())
                {
                    WarningInFunction
                        << "Cannot find any faceZone name matching "
                        << zoneName << endl;
                }
            }

            includeFaceZones = hashed.sortedToc();

            Info<< "Additionally extracting faceZones "
                << UIndirectList<word>(allZoneNames, includeFaceZones)
                << endl;
        }


        // From (name of) patch to compact 'zone' index
        HashTable<label> compactZoneID(1024);
        // Mesh face and compact zone indx
        DynamicList<label> faceLabels;
        DynamicList<label> compactZones;

        {
            // Collect sizes. Hash on names to handle local-only patches (e.g.
            //  processor patches)
            HashTable<label> patchSize(1024);
            label nFaces = 0;
            for (const label patchi : includePatches)
            {
                const polyPatch& pp = bMesh[patchi];
                patchSize.insert(pp.name(), pp.size());
                nFaces += pp.size();
            }

            HashTable<label> zoneSize(1024);
            for (const label zonei : includeFaceZones)
            {
                const faceZone& pp = fzm[zonei];
                zoneSize.insert(pp.name(), pp.size());
                nFaces += pp.size();
            }


            Pstream::mapCombineGather(patchSize, plusEqOp<label>());
            Pstream::mapCombineGather(zoneSize, plusEqOp<label>());


            // Allocate compact numbering for all patches/faceZones
            forAllConstIters(patchSize, iter)
            {
                compactZoneID.insert(iter.key(), compactZoneID.size());
            }

            forAllConstIters(zoneSize, iter)
            {
                compactZoneID.insert(iter.key(), compactZoneID.size());
            }


            Pstream::mapCombineScatter(compactZoneID);


            // Rework HashTable into labelList just for speed of conversion
            labelList patchToCompactZone(bMesh.size(), -1);
            labelList faceZoneToCompactZone(bMesh.size(), -1);
            forAllConstIters(compactZoneID, iter)
            {
                label patchi = bMesh.findPatchID(iter.key());
                if (patchi != -1)
                {
                    patchToCompactZone[patchi] = iter();
                }
                else
                {
                    label zoneI = fzm.findZoneID(iter.key());
                    faceZoneToCompactZone[zoneI] = iter();
                }
            }


            faceLabels.setCapacity(nFaces);
            compactZones.setCapacity(nFaces);

            // Collect faces on patches
            for (const label patchi : includePatches)
            {
                const polyPatch& pp = bMesh[patchi];
                forAll(pp, i)
                {
                    faceLabels.append(pp.start()+i);
                    compactZones.append(patchToCompactZone[pp.index()]);
                }
            }
            // Collect faces on faceZones
            for (const label zonei : includeFaceZones)
            {
                const faceZone& pp = fzm[zonei];
                forAll(pp, i)
                {
                    faceLabels.append(pp[i]);
                    compactZones.append(faceZoneToCompactZone[pp.index()]);
                }
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
            inplaceRenumber
            (
                pointToGlobal,
                gatheredFaces[Pstream::myProcNo()][i]
            );
        }
        Pstream::gatherList(gatheredFaces);

        // Gather all ZoneIDs
        List<labelList> gatheredZones(Pstream::nProcs());
        gatheredZones[Pstream::myProcNo()].transfer(compactZones);
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
                Info<< "surfZone " << iter()
                    <<  " : "      << surfZones[iter()].name()
                    << endl;
            }

            UnsortedMeshedSurface<face> unsortedFace
            (
                std::move(allPoints),
                std::move(allFaces),
                std::move(allZones),
                surfZones
            );


            MeshedSurface<face> sortedFace(unsortedFace);

            fileName globalCasePath
            (
                outFileName.isAbsolute()
              ? outFileName
              : (
                    runTime.processorCase()
                  ? runTime.globalPath()/outFileName
                  : runTime.path()/outFileName
                )
            );

            Info<< "Writing merged surface to " << globalCasePath << endl;

            sortedFace.write(globalCasePath);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
