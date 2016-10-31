/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    surfaceCheck

Group
    grpSurfaceUtilities

Description
    Checks geometric and topological quality of a surface.

Usage
    - surfaceCheck surfaceFile [OPTION]

    \param -checkSelfIntersection \n
    Check for self-intersection.

    \param -splitNonManifold \n
    Split surface along non-manifold edges.

    \param -verbose \n
    Extra verbosity.

    \param -blockMesh \n
    Write vertices/blocks for tight-fitting 1 cell blockMeshDict.

    \param -outputThreshold \<num files\> \n
    Specifies upper limit for the number of files written. This is useful to
    prevent surfaces with lots of disconnected parts to write lots of files.
    Default is 10. A special case is 0 which prevents writing any files.

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "argList.H"
#include "OFstream.H"
#include "OBJstream.H"
#include "SortableList.H"
#include "PatchTools.H"
#include "vtkSurfaceWriter.H"

using namespace Foam;

// Does face use valid vertices?
bool validTri
(
    const bool verbose,
    const triSurface& surf,
    const label faceI
)
{
    // Simple check on indices ok.

    const labelledTri& f = surf[faceI];

    forAll(f, fp)
    {
        if (f[fp] < 0 || f[fp] >= surf.points().size())
        {
            WarningInFunction
                << "triangle " << faceI << " vertices " << f
                << " uses point indices outside point range 0.."
                << surf.points().size()-1 << endl;
            return false;
        }
    }

    if ((f[0] == f[1]) || (f[0] == f[2]) || (f[1] == f[2]))
    {
        WarningInFunction
            << "triangle " << faceI
            << " uses non-unique vertices " << f
            << " coords:" << f.points(surf.points())
            << endl;
        return false;
    }

    // duplicate triangle check

    const labelList& fFaces = surf.faceFaces()[faceI];

    // Check if faceNeighbours use same points as this face.
    // Note: discards normal information - sides of baffle are merged.
    forAll(fFaces, i)
    {
        label nbrFaceI = fFaces[i];

        if (nbrFaceI <= faceI)
        {
            // lower numbered faces already checked
            continue;
        }

        const labelledTri& nbrF = surf[nbrFaceI];

        if
        (
            ((f[0] == nbrF[0]) || (f[0] == nbrF[1]) || (f[0] == nbrF[2]))
         && ((f[1] == nbrF[0]) || (f[1] == nbrF[1]) || (f[1] == nbrF[2]))
         && ((f[2] == nbrF[0]) || (f[2] == nbrF[1]) || (f[2] == nbrF[2]))
        )
        {
            WarningInFunction
                << "triangle " << faceI << " vertices " << f
                << " has the same vertices as triangle " << nbrFaceI
                << " vertices " << nbrF
                << " coords:" << f.points(surf.points())
                << endl;

            return false;
        }
    }
    return true;
}


labelList countBins
(
    const scalar min,
    const scalar max,
    const label nBins,
    const scalarField& vals
)
{
    scalar dist = nBins/(max - min);

    labelList binCount(nBins, 0);

    forAll(vals, i)
    {
        scalar val = vals[i];

        label index = -1;

        if (Foam::mag(val - min) < SMALL)
        {
            index = 0;
        }
        else if (val >= max - SMALL)
        {
            index = nBins - 1;
        }
        else
        {
            index = label((val - min)*dist);

            if ((index < 0) || (index >= nBins))
            {
                WarningInFunction
                    << "value " << val << " at index " << i
                    << " outside range " << min << " .. " << max << endl;

                if (index < 0)
                {
                    index = 0;
                }
                else
                {
                    index = nBins - 1;
                }
            }
        }
        binCount[index]++;
    }

    return binCount;
}



void writeZoning
(
    const triSurface& surf,
    const labelList& faceZone,
    const word& fieldName,
    const fileName& surfFilePath,
    const fileName& surfFileNameBase
)
{
    Info<< "Writing zoning to "
        <<  fileName
            (
                surfFilePath
              / fieldName
              + '_'
              + surfFileNameBase
              + '.'
              + vtkSurfaceWriter::typeName
            )
        << "..." << endl << endl;

    // Convert data
    scalarField scalarFaceZone(faceZone.size());
    forAll(faceZone, i)
    {
        scalarFaceZone[i] = faceZone[i];
    }
    faceList faces(surf.size());
    forAll(surf, i)
    {
        faces[i] = surf[i].triFaceFace();
    }

    vtkSurfaceWriter().write
    (
        surfFilePath,
        surfFileNameBase,
        surf.points(),
        faces,
        fieldName,
        scalarFaceZone,
        false               // face based data
    );
}


void writeParts
(
    const triSurface& surf,
    const label nFaceZones,
    const labelList& faceZone,
    const fileName& surfFilePath,
    const fileName& surfFileNameBase
)
{
    for (label zone = 0; zone < nFaceZones; zone++)
    {
        boolList includeMap(surf.size(), false);

        forAll(faceZone, faceI)
        {
            if (faceZone[faceI] == zone)
            {
                includeMap[faceI] = true;
            }
        }

        labelList pointMap;
        labelList faceMap;

        triSurface subSurf
        (
            surf.subsetMesh
            (
                includeMap,
                pointMap,
                faceMap
            )
        );

        fileName subName
        (
            surfFilePath
           /surfFileNameBase + "_" + name(zone) + ".obj"
        );

        Info<< "writing part " << zone << " size " << subSurf.size()
            << " to " << subName << endl;

        subSurf.write(subName);
    }
}


void syncEdges(const triSurface& p, labelHashSet& markedEdges)
{
    // See comment below about having duplicate edges

    const edgeList& edges = p.edges();
    HashSet<edge, Hash<edge>> edgeSet(2*markedEdges.size());

    forAllConstIter(labelHashSet, markedEdges, iter)
    {
        edgeSet.insert(edges[iter.key()]);
    }

    forAll(edges, edgeI)
    {
        if (edgeSet.found(edges[edgeI]))
        {
            markedEdges.insert(edgeI);
        }
    }
}


void syncEdges(const triSurface& p, boolList& isMarkedEdge)
{
    // See comment below about having duplicate edges

    const edgeList& edges = p.edges();

    label n = 0;
    forAll(isMarkedEdge, edgeI)
    {
        if (isMarkedEdge[edgeI])
        {
            n++;
        }
    }

    HashSet<edge, Hash<edge>> edgeSet(2*n);

    forAll(isMarkedEdge, edgeI)
    {
        if (isMarkedEdge[edgeI])
        {
            edgeSet.insert(edges[edgeI]);
        }
    }

    forAll(edges, edgeI)
    {
        if (edgeSet.found(edges[edgeI]))
        {
            isMarkedEdge[edgeI] = true;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::addBoolOption
    (
        "checkSelfIntersection",
        "also check for self-intersection"
    );
    argList::addBoolOption
    (
        "splitNonManifold",
        "split surface along non-manifold edges"
        " (default split is fully disconnected)"
    );
    argList::addBoolOption
    (
        "verbose",
        "verbose operation"
    );
    argList::addBoolOption
    (
        "blockMesh",
        "write vertices/blocks for blockMeshDict"
    );

    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const bool checkSelfIntersect = args.optionFound("checkSelfIntersection");
    const bool verbose = args.optionFound("verbose");
    const bool splitNonManifold = args.optionFound("splitNonManifold");
    label outputThreshold = 10;
    args.optionReadIfPresent("outputThreshold", outputThreshold);

    Info<< "Reading surface from " << surfFileName << " ..." << nl << endl;


    // Read
    // ~~~~

    triSurface surf(surfFileName);


    Info<< "Statistics:" << endl;
    surf.writeStats(Info);
    Info<< endl;


    // Determine path and extension
    fileName surfFileNameBase(surfFileName.name());
    const word fileType = surfFileNameBase.ext();
    // Strip extension
    surfFileNameBase = surfFileNameBase.lessExt();
    // If extension was .gz strip original extension
    if (fileType == "gz")
    {
        surfFileNameBase = surfFileNameBase.lessExt();
    }
    const fileName surfFilePath(surfFileName.path());


    // write bounding box corners
    if (args.optionFound("blockMesh"))
    {
        pointField cornerPts(boundBox(surf.points(), false).points());

        Info<< "// blockMeshDict info" << nl << nl;

        Info<< "vertices\n(" << nl;
        forAll(cornerPts, ptI)
        {
            Info<< "    " << cornerPts[ptI] << nl;
        }

        // number of divisions needs adjustment later
        Info<< ");\n" << nl
            << "blocks\n"
            << "(\n"
            << "    hex (0 1 2 3 4 5 6 7) (10 10 10) simpleGrading (1 1 1)\n"
            << ");\n" << nl;

        Info<< "edges\n();" << nl
            << "patches\n();" << endl;

        Info<< nl << "// end blockMeshDict info" << nl << endl;
    }


    // Region sizes
    // ~~~~~~~~~~~~

    {
        labelList regionSize(surf.patches().size(), 0);

        forAll(surf, faceI)
        {
            label region = surf[faceI].region();

            if (region < 0 || region >= regionSize.size())
            {
                WarningInFunction
                    << "Triangle " << faceI << " vertices " << surf[faceI]
                    << " has region " << region << " which is outside the range"
                    << " of regions 0.." << surf.patches().size()-1
                    << endl;
            }
            else
            {
                regionSize[region]++;
            }
        }

        Info<< "Region\tSize" << nl
            << "------\t----" << nl;
        forAll(surf.patches(), patchI)
        {
            Info<< surf.patches()[patchI].name() << '\t'
                << regionSize[patchI] << nl;
        }
        Info<< nl << endl;
    }


    // Check triangles
    // ~~~~~~~~~~~~~~~

    {
        DynamicList<label> illegalFaces(surf.size()/100 + 1);

        forAll(surf, faceI)
        {
            if (!validTri(verbose, surf, faceI))
            {
                illegalFaces.append(faceI);
            }
        }

        if (illegalFaces.size())
        {
            Info<< "Surface has " << illegalFaces.size()
                << " illegal triangles." << endl;

            if (outputThreshold > 0)
            {
                OFstream str("illegalFaces");
                Info<< "Dumping conflicting face labels to " << str.name()
                    << endl
                    << "Paste this into the input for surfaceSubset" << endl;
                str << illegalFaces;
            }
        }
        else
        {
            Info<< "Surface has no illegal triangles." << endl;
        }
        Info<< endl;
    }



    // Triangle quality
    // ~~~~~~~~~~~~~~~~

    {
        scalarField triQ(surf.size(), 0);
        forAll(surf, faceI)
        {
            const labelledTri& f = surf[faceI];

            if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
            {
                //WarningInFunction
                //    << "Illegal triangle " << faceI << " vertices " << f
                //    << " coords " << f.points(surf.points()) << endl;
            }
            else
            {
                triQ[faceI] = triPointRef
                (
                    surf.points()[f[0]],
                    surf.points()[f[1]],
                    surf.points()[f[2]]
                ).quality();
            }
        }

        labelList binCount = countBins(0, 1, 20, triQ);

        Info<< "Triangle quality (equilateral=1, collapsed=0):"
            << endl;


        OSstream& os = Info;
        os.width(4);

        scalar dist = (1.0 - 0.0)/20.0;
        scalar min = 0;
        forAll(binCount, binI)
        {
            Info<< "    " << min << " .. " << min+dist << "  : "
                << 1.0/surf.size() * binCount[binI]
                << endl;
            min += dist;
        }
        Info<< endl;

        label minIndex = findMin(triQ);
        label maxIndex = findMax(triQ);

        Info<< "    min " << triQ[minIndex] << " for triangle " << minIndex
            << nl
            << "    max " << triQ[maxIndex] << " for triangle " << maxIndex
            << nl
            << endl;


        if (triQ[minIndex] < SMALL)
        {
            WarningInFunction
                << triQ[minIndex] << ". This might give problems in"
                << " self-intersection testing later on." << endl;
        }

        // Dump for subsetting
        if (outputThreshold > 0)
        {
            DynamicList<label> problemFaces(surf.size()/100+1);

            forAll(triQ, faceI)
            {
                if (triQ[faceI] < 1e-11)
                {
                    problemFaces.append(faceI);
                }
            }

            if (!problemFaces.empty())
            {
                OFstream str("badFaces");

                Info<< "Dumping bad quality faces to " << str.name() << endl
                    << "Paste this into the input for surfaceSubset" << nl
                    << nl << endl;

                str << problemFaces;
            }
        }
    }



    // Edges
    // ~~~~~
    {
        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();

        scalarField edgeMag(edges.size());

        forAll(edges, edgeI)
        {
            edgeMag[edgeI] = edges[edgeI].mag(localPoints);
        }

        label minEdgeI = findMin(edgeMag);
        label maxEdgeI = findMax(edgeMag);

        const edge& minE = edges[minEdgeI];
        const edge& maxE = edges[maxEdgeI];


        Info<< "Edges:" << nl
            << "    min " << edgeMag[minEdgeI] << " for edge " << minEdgeI
            << " points " << localPoints[minE[0]] << localPoints[minE[1]]
            << nl
            << "    max " << edgeMag[maxEdgeI] << " for edge " << maxEdgeI
            << " points " << localPoints[maxE[0]] << localPoints[maxE[1]]
            << nl
            << endl;
    }



    // Close points
    // ~~~~~~~~~~~~
    {
        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();

        const boundBox bb(localPoints);
        scalar smallDim = 1e-6 * bb.mag();

        Info<< "Checking for points less than 1e-6 of bounding box ("
            << bb.span() << " metre) apart."
            << endl;

        // Sort points
        SortableList<scalar> sortedMag(mag(localPoints));

        label nClose = 0;

        for (label i = 1; i < sortedMag.size(); i++)
        {
            label ptI = sortedMag.indices()[i];

            label prevPtI = sortedMag.indices()[i-1];

            if (mag(localPoints[ptI] - localPoints[prevPtI]) < smallDim)
            {
                // Check if neighbours.
                const labelList& pEdges = surf.pointEdges()[ptI];

                label edgeI = -1;

                forAll(pEdges, i)
                {
                    const edge& e = edges[pEdges[i]];

                    if (e[0] == prevPtI || e[1] == prevPtI)
                    {
                        // point1 and point0 are connected through edge.
                        edgeI = pEdges[i];

                        break;
                    }
                }

                nClose++;

                if (edgeI == -1)
                {
                    Info<< "    close unconnected points "
                        << ptI << ' ' << localPoints[ptI]
                        << " and " << prevPtI << ' '
                        << localPoints[prevPtI]
                        << " distance:"
                        << mag(localPoints[ptI] - localPoints[prevPtI])
                        << endl;
                }
                else
                {
                    Info<< "    small edge between points "
                        << ptI << ' ' << localPoints[ptI]
                        << " and " << prevPtI << ' '
                        << localPoints[prevPtI]
                        << " distance:"
                        << mag(localPoints[ptI] - localPoints[prevPtI])
                        << endl;
                }
            }
        }

        Info<< "Found " << nClose << " nearby points." << nl
            << endl;
    }



    // Check manifold
    // ~~~~~~~~~~~~~~

    DynamicList<label> problemFaces(surf.size()/100 + 1);

    const labelListList& edgeFaces = surf.edgeFaces();

    label nSingleEdges = 0;
    forAll(edgeFaces, edgeI)
    {
        const labelList& myFaces = edgeFaces[edgeI];

        if (myFaces.size() == 1)
        {
            problemFaces.append(myFaces[0]);

            nSingleEdges++;
        }
    }

    label nMultEdges = 0;
    forAll(edgeFaces, edgeI)
    {
        const labelList& myFaces = edgeFaces[edgeI];

        if (myFaces.size() > 2)
        {
            forAll(myFaces, myFaceI)
            {
                problemFaces.append(myFaces[myFaceI]);
            }

            nMultEdges++;
        }
    }
    problemFaces.shrink();

    if ((nSingleEdges != 0) || (nMultEdges != 0))
    {
        Info<< "Surface is not closed since not all edges ("
            << edgeFaces.size() << ") connected to "
            << "two faces:" << endl
            << "    connected to one face : " << nSingleEdges << endl
            << "    connected to >2 faces : " << nMultEdges << endl;

        Info<< "Conflicting face labels:" << problemFaces.size() << endl;

        if (outputThreshold > 0)
        {
            OFstream str("problemFaces");

            Info<< "Dumping conflicting face labels to " << str.name() << endl
                << "Paste this into the input for surfaceSubset" << endl;

            str << problemFaces;
        }
    }
    else
    {
        Info<< "Surface is closed. All edges connected to two faces." << endl;
    }
    Info<< endl;



    // Check singly connected domain
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        boolList borderEdge(surf.nEdges(), false);
        if (splitNonManifold)
        {
            forAll(edgeFaces, edgeI)
            {
                if (edgeFaces[edgeI].size() > 2)
                {
                    borderEdge[edgeI] = true;
                }
            }
            syncEdges(surf, borderEdge);
        }

        labelList faceZone;
        label numZones = surf.markZones(borderEdge, faceZone);

        Info<< "Number of unconnected parts : " << numZones << endl;

        if (numZones > 1 && outputThreshold > 0)
        {
            Info<< "Splitting surface into parts ..." << endl << endl;

            writeZoning(surf, faceZone, "zone", surfFilePath, surfFileNameBase);

            if (numZones > outputThreshold)
            {
                Info<< "Limiting number of files to " << outputThreshold
                    << endl;
            }
            writeParts
            (
                surf,
                min(outputThreshold, numZones),
                faceZone,
                surfFilePath,
                surfFileNameBase
            );
        }
    }



    // Check orientation
    // ~~~~~~~~~~~~~~~~~

    labelHashSet borderEdge(surf.size()/1000);
    PatchTools::checkOrientation(surf, false, &borderEdge);

    // Bit strange: if a triangle has two same vertices (illegal!) it will
    // still have three distinct edges (two of which have the same vertices).
    // In this case the faceEdges addressing is not symmetric, i.e. a
    // neighbouring, valid, triangle will have correct addressing so 3 distinct
    // edges so it will miss one of those two identical edges.
    // - we don't want to fix this in PrimitivePatch since it is too specific
    // - instead just make sure we mark all identical edges consistently
    //   when we use them for marking.

    syncEdges(surf, borderEdge);

    //
    // Colour all faces into zones using borderEdge
    //
    labelList normalZone;
    label numNormalZones = PatchTools::markZones(surf, borderEdge, normalZone);

    Info<< endl
        << "Number of zones (connected area with consistent normal) : "
        << numNormalZones << endl;

    if (numNormalZones > 1)
    {
        Info<< "More than one normal orientation." << endl;

        if (outputThreshold > 0)
        {
            writeZoning
            (
                surf,
                normalZone,
                "normal",
                surfFilePath,
                surfFileNameBase
            );

            if (numNormalZones > outputThreshold)
            {
                Info<< "Limiting number of files to " << outputThreshold
                    << endl;
            }
            writeParts
            (
                surf,
                min(outputThreshold, numNormalZones),
                normalZone,
                surfFilePath,
                surfFileNameBase + "_normal"
            );
        }
    }
    Info<< endl;



    // Check self-intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (checkSelfIntersect)
    {
        Info<< "Checking self-intersection." << endl;

        triSurfaceSearch querySurf(surf);

        const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

        autoPtr<OBJstream> intStreamPtr;
        if (outputThreshold > 0)
        {
            intStreamPtr.reset(new OBJstream("selfInterPoints.obj"));
        }

        label nInt = 0;

        forAll(surf.edges(), edgeI)
        {
            const edge& e = surf.edges()[edgeI];
            const point& start = surf.points()[surf.meshPoints()[e[0]]];
            const point& end = surf.points()[surf.meshPoints()[e[1]]];

            // Exclude hits of connected triangles
            treeDataTriSurface::findSelfIntersectOp exclOp(tree, edgeI);

            pointIndexHit hitInfo(tree.findLineAny(start, end, exclOp));

            if (hitInfo.hit())
            {
                nInt++;

                if (intStreamPtr.valid())
                {
                    intStreamPtr().write(hitInfo.hitPoint());
                }

                // Try and find from other side.
                pointIndexHit hitInfo2(tree.findLineAny(end, start, exclOp));

                if (hitInfo2.hit() && hitInfo.index() != hitInfo2.index())
                {
                    nInt++;

                    if (intStreamPtr.valid())
                    {
                        intStreamPtr().write(hitInfo2.hitPoint());
                    }
                }

            }
        }

        //// Check very near triangles
        //{
        //    const pointField& localPoints = surf.localPoints();
        //
        //    const boundBox bb(localPoints);
        //    scalar smallDim = 1e-6 * bb.mag();
        //    scalar smallDimSqr = Foam::sqr(smallDim);
        //
        //    const pointField& faceCentres = surf.faceCentres();
        //    forAll(faceCentres, faceI)
        //    {
        //        const point& fc = faceCentres[faceI];
        //        pointIndexHit hitInfo
        //        (
        //            tree.findNearest
        //            (
        //                fc,
        //                smallDimSqr,
        //                findSelfNearOp(tree, faceI)
        //            )
        //        );
        //
        //        if (hitInfo.hit() && intStreamPtr.valid())
        //        {
        //            intStreamPtr().write(hitInfo.hitPoint());
        //
        //            label nearFaceI = hitInfo.index();
        //            triPointRef nearTri(surf[nearFaceI].tri(surf.points()));
        //            triStreamPtr().write
        //            (
        //                surf[faceI].tri(surf.points()),
        //                false
        //            );
        //            triStreamPtr().write(nearTri, false);
        //            nInt++;
        //        }
        //    }
        //}


        if (nInt == 0)
        {
            Info<< "Surface is not self-intersecting" << endl;
        }
        else
        {
            Info<< "Surface is self-intersecting at " << nInt
                << " locations." << endl;

            if (intStreamPtr.valid())
            {
                Info<< "Writing intersection points to "
                    << intStreamPtr().name() << endl;
            }
        }
        Info<< endl;
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
