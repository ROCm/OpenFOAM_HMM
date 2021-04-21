/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "checkTools.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "hexMatcher.H"
#include "wedgeMatcher.H"
#include "prismMatcher.H"
#include "pyrMatcher.H"
#include "tetWedgeMatcher.H"
#include "tetMatcher.H"
#include "IOmanip.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "Time.H"
#include "surfaceWriter.H"
#include "syncTools.H"
#include "globalIndex.H"
#include "PatchTools.H"
#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::printMeshStats(const polyMesh& mesh, const bool allTopology)
{
    if (mesh.name() == Foam::polyMesh::defaultRegion)
    {
        Info<< "Mesh stats" << nl;
    }
    else
    {
        Info<< "Mesh " << mesh.name() << " stats" << nl;
    }
    Info<< "    points:           "
        << returnReduce(mesh.points().size(), sumOp<label>()) << nl;


    // Count number of internal points (-1 if not sorted; 0 if no internal
    // points)
    const label minInt = returnReduce(mesh.nInternalPoints(), minOp<label>());
    const label maxInt = returnReduce(mesh.nInternalPoints(), maxOp<label>());

    if (minInt == -1 && maxInt > 0)
    {
        WarningInFunction
            << "Some processors have their points sorted into internal"
            << " and external and some do not." << endl
            << "    This can cause problems later on." << endl;
    }
    else if (minInt != -1)
    {
        // Assume all sorted
        label nInternalPoints = returnReduce
        (
            mesh.nInternalPoints(),
            sumOp<label>()
        );
        Info<< "    internal points:  " << nInternalPoints << nl;
    }

    if (allTopology && (minInt != -1))
    {
        label nEdges = returnReduce(mesh.nEdges(), sumOp<label>());
        label nInternalEdges = returnReduce
        (
            mesh.nInternalEdges(),
            sumOp<label>()
        );
        label nInternal1Edges = returnReduce
        (
            mesh.nInternal1Edges(),
            sumOp<label>()
        );
        label nInternal0Edges = returnReduce
        (
            mesh.nInternal0Edges(),
            sumOp<label>()
        );

        Info<< "    edges:            " << nEdges << nl
            << "    internal edges:   " << nInternalEdges << nl
            << "    internal edges using one boundary point:   "
            << nInternal1Edges-nInternal0Edges << nl
            << "    internal edges using two boundary points:  "
            << nInternalEdges-nInternal1Edges << nl;
    }

    label nFaces = returnReduce(mesh.faces().size(), sumOp<label>());
    label nIntFaces = returnReduce(mesh.faceNeighbour().size(), sumOp<label>());
    label nCells = returnReduce(mesh.cells().size(), sumOp<label>());
    label nPatches = mesh.boundaryMesh().size();

    Info<< "    faces:            " << nFaces << nl
        << "    internal faces:   " << nIntFaces << nl
        << "    cells:            " << nCells << nl
        << "    faces per cell:   "
        << (scalar(nFaces) + scalar(nIntFaces))/max(1, nCells) << nl
        << "    boundary patches: ";

    if (Pstream::parRun())
    {
        // Number of global patches and min-max range of total patches
        Info<< mesh.boundaryMesh().nNonProcessor() << ' '
            << returnReduce(labelMinMax(nPatches), minMaxOp<label>()) << nl;
    }
    else
    {
        Info<< nPatches << nl;
    }

    Info<< "    point zones:      " << mesh.pointZones().size() << nl
        << "    face zones:       " << mesh.faceZones().size() << nl
        << "    cell zones:       " << mesh.cellZones().size() << nl
        << endl;

    // Construct shape recognizers
    prismMatcher prism;
    wedgeMatcher wedge;
    tetWedgeMatcher tetWedge;

    // Counters for different cell types
    label nHex = 0;
    label nWedge = 0;
    label nPrism = 0;
    label nPyr = 0;
    label nTet = 0;
    label nTetWedge = 0;
    label nUnknown = 0;

    Map<label> polyhedralFaces;

    for (label celli = 0; celli < mesh.nCells(); celli++)
    {
        if (hexMatcher::test(mesh, celli))
        {
            nHex++;
        }
        else if (tetMatcher::test(mesh, celli))
        {
            nTet++;
        }
        else if (pyrMatcher::test(mesh, celli))
        {
            nPyr++;
        }
        else if (prism.isA(mesh, celli))
        {
            nPrism++;
        }
        else if (wedge.isA(mesh, celli))
        {
            nWedge++;
        }
        else if (tetWedge.isA(mesh, celli))
        {
            nTetWedge++;
        }
        else
        {
            nUnknown++;
            polyhedralFaces(mesh.cells()[celli].size())++;
        }
    }

    reduce(nHex,sumOp<label>());
    reduce(nPrism,sumOp<label>());
    reduce(nWedge,sumOp<label>());
    reduce(nPyr,sumOp<label>());
    reduce(nTetWedge,sumOp<label>());
    reduce(nTet,sumOp<label>());
    reduce(nUnknown,sumOp<label>());

    Info<< "Overall number of cells of each type:" << nl
        << "    hexahedra:     " << nHex << nl
        << "    prisms:        " << nPrism << nl
        << "    wedges:        " << nWedge << nl
        << "    pyramids:      " << nPyr << nl
        << "    tet wedges:    " << nTetWedge << nl
        << "    tetrahedra:    " << nTet << nl
        << "    polyhedra:     " << nUnknown
        << endl;

    if (nUnknown > 0)
    {
        Pstream::mapCombineGather(polyhedralFaces, plusEqOp<label>());

        Info<< "    Breakdown of polyhedra by number of faces:" << nl
            << "        faces" << "   number of cells" << endl;

        const labelList sortedKeys = polyhedralFaces.sortedToc();

        forAll(sortedKeys, keyi)
        {
            const label nFaces = sortedKeys[keyi];

            Info<< setf(std::ios::right) << setw(13)
                << nFaces << "   " << polyhedralFaces[nFaces] << nl;
        }
    }

    Info<< endl;
}


void Foam::mergeAndWrite
(
    const polyMesh& mesh,
    surfaceWriter& writer,
    const word& name,
    const indirectPrimitivePatch& setPatch,
    const fileName& outputDir
)
{
    if (Pstream::parRun())
    {
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPoints;
        autoPtr<globalIndex> globalFaces;
        faceList mergedFaces;
        pointField mergedPoints;
        Foam::PatchTools::gatherAndMerge
        (
            mesh,
            setPatch.localFaces(),
            setPatch.meshPoints(),
            setPatch.meshPointMap(),

            pointToGlobal,
            uniqueMeshPointLabels,
            globalPoints,
            globalFaces,

            mergedFaces,
            mergedPoints
        );

        // Write
        if (Pstream::master())
        {
            writer.open
            (
                mergedPoints,
                mergedFaces,
                (outputDir / name),
                false  // serial - already merged
            );

            writer.write();
            writer.clear();
        }
    }
    else
    {
        writer.open
        (
            setPatch.localPoints(),
            setPatch.localFaces(),
            (outputDir / name),
            false  // serial - already merged
        );

        writer.write();
        writer.clear();
    }
}


void Foam::mergeAndWrite
(
    surfaceWriter& writer,
    const faceSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());

    const indirectPrimitivePatch setPatch
    (
        IndirectList<face>(mesh.faces(), set.sortedToc()),
        mesh.points()
    );

    fileName outputDir
    (
        set.time().globalPath()
      / functionObject::outputPrefix
      / mesh.pointsInstance()
      / set.name()
    );
    outputDir.clean();  // Remove unneeded ".."

    mergeAndWrite(mesh, writer, set.name(), setPatch, outputDir);
}


void Foam::mergeAndWrite
(
    surfaceWriter& writer,
    const cellSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();


    // Determine faces on outside of cellSet
    bitSet isInSet(mesh.nCells());
    for (const label celli : set)
    {
        isInSet.set(celli);
    }


    boolList bndInSet(mesh.nBoundaryFaces());
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const labelList& fc = pp.faceCells();
        forAll(fc, i)
        {
            bndInSet[pp.start()+i-mesh.nInternalFaces()] = isInSet[fc[i]];
        }
    }
    syncTools::swapBoundaryFaceList(mesh, bndInSet);


    DynamicList<label> outsideFaces(3*set.size());
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        const bool ownVal = isInSet[mesh.faceOwner()[facei]];
        const bool neiVal = isInSet[mesh.faceNeighbour()[facei]];

        if (ownVal != neiVal)
        {
            outsideFaces.append(facei);
        }
    }


    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const labelList& fc = pp.faceCells();
        if (pp.coupled())
        {
            forAll(fc, i)
            {
                label facei = pp.start()+i;

                const bool neiVal = bndInSet[facei-mesh.nInternalFaces()];
                if (isInSet[fc[i]] && !neiVal)
                {
                    outsideFaces.append(facei);
                }
            }
        }
        else
        {
            forAll(fc, i)
            {
                if (isInSet[fc[i]])
                {
                    outsideFaces.append(pp.start()+i);
                }
            }
        }
    }


    const indirectPrimitivePatch setPatch
    (
        IndirectList<face>(mesh.faces(), outsideFaces),
        mesh.points()
    );

    fileName outputDir
    (
        set.time().globalPath()
      / functionObject::outputPrefix
      / mesh.pointsInstance()
      / set.name()
    );
    outputDir.clean();  // Remove unneeded ".."

    mergeAndWrite(mesh, writer, set.name(), setPatch, outputDir);
}


void Foam::mergeAndWrite
(
    const writer<scalar>& writer,
    const pointSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());

    pointField mergedPts;
    labelList mergedIDs;

    if (Pstream::parRun())
    {
        // Note: we explicitly do not merge the points
        // (mesh.globalData().mergePoints etc) since this might
        // hide any synchronisation problem

        globalIndex globalNumbering(mesh.nPoints());

        mergedPts.setSize(returnReduce(set.size(), sumOp<label>()));
        mergedIDs.setSize(mergedPts.size());

        labelList setPointIDs(set.sortedToc());

        // Get renumbered local data
        pointField myPoints(mesh.points(), setPointIDs);
        labelList myIDs(globalNumbering.toGlobal(setPointIDs));

        if (Pstream::master())
        {
            // Insert master data first
            label pOffset = 0;
            SubList<point>(mergedPts, myPoints.size(), pOffset) = myPoints;
            SubList<label>(mergedIDs, myIDs.size(), pOffset) = myIDs;
            pOffset += myPoints.size();

            // Receive slave ones
            for (const int slave : Pstream::subProcs())
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

                pointField slavePts(fromSlave);
                labelList slaveIDs(fromSlave);

                SubList<point>(mergedPts, slavePts.size(), pOffset) = slavePts;
                SubList<label>(mergedIDs, slaveIDs.size(), pOffset) = slaveIDs;
                pOffset += slaveIDs.size();
            }
        }
        else
        {
            // Construct processor stream with estimate of size. Could
            // be improved.
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                myPoints.byteSize() + myIDs.byteSize()
            );
            toMaster << myPoints << myIDs;
        }
    }
    else
    {
        mergedIDs = set.sortedToc();
        mergedPts = pointField(mesh.points(), mergedIDs);
    }


    // Write with scalar pointID
    if (Pstream::master())
    {
        scalarField scalarPointIDs(mergedIDs.size());
        forAll(mergedIDs, i)
        {
            scalarPointIDs[i] = 1.0*mergedIDs[i];
        }

        coordSet points(set.name(), "distance", mergedPts, mag(mergedPts));

        List<const scalarField*> flds(1, &scalarPointIDs);

        wordList fldNames(1, "pointID");

        // Output e.g. pointSet p0 to
        // postProcessing/<time>/p0.vtk
        fileName outputDir
        (
            set.time().globalPath()
          / functionObject::outputPrefix
          / mesh.pointsInstance()
          // set.name()
        );
        outputDir.clean();  // Remove unneeded ".."
        mkDir(outputDir);

        fileName outputFile(outputDir/writer.getFileName(points, wordList()));
        //fileName outputFile(outputDir/set.name());

        OFstream os(outputFile);

        writer.write(points, fldNames, flds, os);
    }
}


// ************************************************************************* //
