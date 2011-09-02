/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "AMIInterpolation.H"
#include "treeDataPrimitivePatch.H"
#include "indexedOctree.H"
#include "meshTools.H"
#include "mergePoints.H"

#include "vtkSurfaceWriter.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writeIntersectionOBJ
(
    const scalar area,
    const face& f1,
    const face& f2,
    const pointField& f1Points,
    const pointField& f2Points
) const
{
    static label count = 1;

    const pointField f1pts = f1.points(f1Points);
    const pointField f2pts = f2.points(f2Points);

    Info<< "Face intersection area (" << count <<  "):" << nl
        << "    f1 face = " << f1 << nl
        << "    f1 pts  = " << f1pts << nl
        << "    f2 face = " << f2 << nl
        << "    f2 pts  = " << f2pts << nl
        << "    area    = " << area
        << endl;

    OFstream os("areas" + name(count) + ".obj");

    forAll(f1pts, i)
    {
        meshTools::writeOBJ(os, f1pts[i]);
    }
    os<< "l";
    forAll(f1pts, i)
    {
        os<< " " << i + 1;
    }
    os<< " 1" << endl;


    forAll(f2pts, i)
    {
        meshTools::writeOBJ(os, f2pts[i]);
    }
    os<< "l";
    forAll(f2pts, i)
    {
        os<< " " << f1pts.size() + i + 1;
    }
    os<< " " << f1pts.size() + 1 << endl;

    count++;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::checkPatches
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
)
{
    const scalar maxBoundsError = 0.05;

    // check bounds of source and target
    boundBox bbSrc(srcPatch.points(), srcPatch.meshPoints());
    boundBox bbTgt(tgtPatch.points(), tgtPatch.meshPoints());

    bbTgt.inflate(maxBoundsError);

    if (!bbTgt.contains(bbSrc))
    {
        WarningIn
        (
            "AMIInterpolation<SourcePatch, TargetPatch>::checkPatches"
            "("
                "const primitivePatch&, "
                "const primitivePatch&"
            ")"
        )   << "Source and target patch bounding boxes are not similar" << nl
            << "    src span : " << bbSrc.span() << nl
            << "    tgt span : " << bbTgt.span() << nl
            << "    source: " << bbSrc << nl
            << "    target: " << bbTgt << endl;
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label
Foam::AMIInterpolation<SourcePatch, TargetPatch>::calcOverlappingProcs
(
    const List<treeBoundBoxList>& procBb,
    const treeBoundBox& bb,
    boolList& overlaps
)
{
    overlaps.setSize(procBb.size());
    overlaps = false;

    label nOverlaps = 0;

    forAll(procBb, procI)
    {
        const List<treeBoundBox>& bbs = procBb[procI];

        forAll(bbs, bbI)
        {
            if (bbs[bbI].overlaps(bb))
            {
                overlaps[procI] = true;
                nOverlaps++;
                break;
            }
        }
    }
    return nOverlaps;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::distributePatches
(
    const mapDistribute& map,
    const primitivePatch& pp,
    const globalIndex& gi,
    List<faceList>& faces,
    List<pointField>& points,
    List<labelList>& faceIDs
)
{
    PstreamBuffers pBufs(Pstream::nonBlocking);

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& sendElems = map.subMap()[domain];

        if (domain != Pstream::myProcNo() && sendElems.size())
        {
            labelList globalElems(sendElems.size());
            forAll(sendElems, i)
            {
                globalElems[i] = gi.toGlobal(sendElems[i]);
            }

            faceList subFaces(UIndirectList<face>(pp, sendElems));
            primitivePatch subPatch
            (
                SubList<face>(subFaces, subFaces.size()),
                pp.points()
            );

            if (debug)
            {
                Pout<< "distributePatches: to processor " << domain
                    << " sending faces " << subPatch.faceCentres() << endl;
            }

            UOPstream toDomain(domain, pBufs);
            toDomain
                << subPatch.localFaces() << subPatch.localPoints()
                << globalElems;
        }
    }

    // Start receiving
    pBufs.finishedSends();

    faces.setSize(Pstream::nProcs());
    points.setSize(Pstream::nProcs());
    faceIDs.setSize(Pstream::nProcs());

    {
        // Set up 'send' to myself
        const labelList& sendElems = map.subMap()[Pstream::myProcNo()];
        faceList subFaces(UIndirectList<face>(pp, sendElems));
        primitivePatch subPatch
        (
            SubList<face>(subFaces, subFaces.size()),
            pp.points()
        );

        // Receive
        if (debug)
        {
            Pout<< "distributePatches: to processor " << Pstream::myProcNo()
                << " sending faces " << subPatch.faceCentres() << endl;
        }

        faces[Pstream::myProcNo()] = subPatch.localFaces();
        points[Pstream::myProcNo()] = subPatch.localPoints();

        faceIDs[Pstream::myProcNo()].setSize(sendElems.size());
        forAll(sendElems, i)
        {
            faceIDs[Pstream::myProcNo()][i] = gi.toGlobal(sendElems[i]);
        }
    }

    // Consume
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
    {
        const labelList& recvElems = map.constructMap()[domain];

        if (domain != Pstream::myProcNo() && recvElems.size())
        {
            UIPstream str(domain, pBufs);

            str >> faces[domain]
                >> points[domain]
                >> faceIDs[domain]; 
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::
distributeAndMergePatches
(
    const mapDistribute& map,
    const primitivePatch& tgtPatch,
    const globalIndex& gi,
    faceList& tgtFaces,
    pointField& tgtPoints,
    labelList& tgtFaceIDs
)
{
    // Exchange per-processor data
    List<faceList> allFaces;
    List<pointField> allPoints;
    List<labelList> allTgtFaceIDs;
    distributePatches(map, tgtPatch, gi, allFaces, allPoints, allTgtFaceIDs);

    // Renumber and flatten
    label nFaces = 0;
    label nPoints = 0;
    forAll(allFaces, procI)
    {
        nFaces += allFaces[procI].size();
        nPoints += allPoints[procI].size();
    }

    tgtFaces.setSize(nFaces);
    tgtPoints.setSize(nPoints);
    tgtFaceIDs.setSize(nFaces);

    nFaces = 0;
    nPoints = 0;

    // My own data first
    {
        const labelList& faceIDs = allTgtFaceIDs[Pstream::myProcNo()];
        SubList<label>(tgtFaceIDs, faceIDs.size()).assign(faceIDs);

        const faceList& fcs = allFaces[Pstream::myProcNo()];
        forAll(fcs, i)
        {
            const face& f = fcs[i];
            face& newF = tgtFaces[nFaces++];
            newF.setSize(f.size());
            forAll(f, fp)
            {
                newF[fp] = f[fp] + nPoints;
            }
        }

        const pointField& pts = allPoints[Pstream::myProcNo()];
        forAll(pts, i)
        {
            tgtPoints[nPoints++] = pts[i];
        }
    }


    // Other proc data follows
    forAll(allFaces, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            const labelList& faceIDs = allTgtFaceIDs[procI];
            SubList<label>(tgtFaceIDs, faceIDs.size(), nFaces).assign(faceIDs);

            const faceList& fcs = allFaces[procI];
            forAll(fcs, i)
            {
                const face& f = fcs[i];
                face& newF = tgtFaces[nFaces++];
                newF.setSize(f.size());
                forAll(f, fp)
                {
                    newF[fp] = f[fp] + nPoints;
                }
            }

            const pointField& pts = allPoints[procI];
            forAll(pts, i)
            {
                tgtPoints[nPoints++] = pts[i];
            }
        }
    }

    // Merge
    labelList oldToNew;
    pointField newTgtPoints;
    bool hasMerged = mergePoints
    (
        tgtPoints,
        SMALL,
        false,
        oldToNew,
        newTgtPoints
    );

    if (hasMerged)
    {
        if (debug)
        {
            Pout<< "Merged from " << tgtPoints.size()
                << " down to " << newTgtPoints.size() << " points" << endl;
        }

        tgtPoints.transfer(newTgtPoints);
        forAll(tgtFaces, i)
        {
            inplaceRenumber(oldToNew, tgtFaces[i]);
        }
    }
}


template<class SourcePatch, class TargetPatch>
Foam::autoPtr<Foam::mapDistribute>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::calcProcMap
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
)
{
    // Get decomposition of patch
    List<treeBoundBoxList> procBb(Pstream::nProcs());

    if (srcPatch.size())
    {
        procBb[Pstream::myProcNo()] = treeBoundBoxList
        (
            1,  // For now single bounding box per proc
            treeBoundBox
            (
                srcPatch.points(),
                srcPatch.meshPoints()
            )
        );
    }
    else
    {
        procBb[Pstream::myProcNo()] = treeBoundBoxList();
    }

    // slightly increase size of bounding boxes to allow for cases where
    // bounding boxes are perfectly alligned
    forAll(procBb[Pstream::myProcNo()], bbI)
    {
        treeBoundBox& bb = procBb[Pstream::myProcNo()][bbI];
        bb.inflate(0.01);
    }

    Pstream::gatherList(procBb);
    Pstream::scatterList(procBb);

    // Determine which faces of tgtPatch overlaps srcPatch per proc
    const faceList& faces = tgtPatch.localFaces();
    const pointField& points = tgtPatch.localPoints();

    labelListList sendMap;

    {
        // Per processor indices into all segments to send
        List<DynamicList<label> > dynSendMap(Pstream::nProcs());

        // Work array - whether processor bb overlaps the face bounds
        boolList procBbOverlaps(Pstream::nProcs());

        forAll(faces, faceI)
        {
            if (faces[faceI].size())
            {
                treeBoundBox faceBb(points, faces[faceI]);

                // Find the processor this face overlaps
                calcOverlappingProcs(procBb, faceBb, procBbOverlaps);

                forAll(procBbOverlaps, procI)
                {
                    if (procBbOverlaps[procI])
                    {
                        dynSendMap[procI].append(faceI);
                    }
                }
            }
        }

        // Convert dynamicList to labelList
        sendMap.setSize(Pstream::nProcs());
        forAll(sendMap, procI)
        {
            sendMap[procI].transfer(dynSendMap[procI]);
        }
    }


    // Send over how many faces I need to receive
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendSizes[Pstream::myProcNo()][procI] = sendMap[procI].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);


    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    // My local segment first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label segmentI = constructMap[Pstream::myProcNo()].size();
    forAll(constructMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            // What I need to receive is what other processor is sending to me
            label nRecv = sendSizes[procI][Pstream::myProcNo()];
            constructMap[procI].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[procI][i] = segmentI++;
            }
        }
    }

    autoPtr<mapDistribute> mapPtr
    (
        new mapDistribute
        (
            segmentI,       // size after construction
            sendMap.xfer(),
            constructMap.xfer()
        )
    );

    return mapPtr;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::projectPointsToSurface
(
    const searchableSurface& surf,
    pointField& pts
) const
{
    if (debug)
    {
        Info<< "AMI: projecting points to surface" << endl;
    }

    List<pointIndexHit> nearInfo;

    surf.findNearest(pts, scalarField(pts.size(), GREAT), nearInfo);

    label nMiss = 0;
    forAll(nearInfo, i)
    {
        const pointIndexHit& pi = nearInfo[i];

        if (pi.hit())
        {
            pts[i] = pi.hitPoint();
        }
        else
        {
            pts[i] = pts[i];
            nMiss++;
        }
    }

    if (nMiss > 0)
    {
        FatalErrorIn
        (
            "void Foam::projectPointsToSurface"
            "("
                "const searchableSurface&, "
                "pointField&"
            ") const"
        )
            << "Error projecting points to surface: "
            << nMiss << " faces could not be determined"
            << abort(FatalError);
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::findTargetFace
(
    const label srcFaceI,
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
) const
{
    label targetFaceI = -1;

    treeBoundBox bb(tgtPatch.points());
    bb.inflate(0.01);

    typedef treeDataPrimitivePatch<face, SubList, const pointField&> treeType;

    indexedOctree<treeType> tree
    (
        treeType(false, tgtPatch),
        bb,                         // overall search domain
        8,                          // maxLevel
        10,                         // leaf size
        3.0                         // duplicity
    );

    const pointField& srcPts = srcPatch.points();
    const face& srcFace = srcPatch[srcFaceI];
    const point& srcPt = srcFace.centre(srcPts);
    const scalar srcFaceArea = srcFace.mag(srcPts);

//    pointIndexHit sample = tree.findNearest(srcPt, sqr(0.1*bb.mag()));
    pointIndexHit sample = tree.findNearest(srcPt, 10.0*srcFaceArea);


    if (debug)
    {
        Info<< "Source point = " << srcPt << ", Sample point = "
            << sample.hitPoint() << ", Sample index = " << sample.index()
            << endl;
    }

    if (sample.hit())
    {
        targetFaceI = sample.index();
    }

    return targetFaceI;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::appendNbrFaces
(
    const label faceI,
    const primitivePatch& patch,
    const DynamicList<label>& visitedFaces,
    DynamicList<label>& faceIDs
) const
{
//    const labelList& nbrFaces = patch.pointFaces()[faceI];
    const labelList& nbrFaces = patch.faceFaces()[faceI];

    // filter out faces already visited from src face neighbours
    forAll(nbrFaces, i)
    {
        label nbrFaceI = nbrFaces[i];
        bool valid = true;
        forAll(visitedFaces, j)
        {
            if (nbrFaceI == visitedFaces[j])
            {
                valid = false;
                break;
            }
        }

        if (valid)
        {
            forAll(faceIDs, j)
            {
                if (nbrFaceI == faceIDs[j])
                {
                    valid = false;
                    break;
                }
            }
        }

        if (valid)
        {
            faceIDs.append(nbrFaceI);
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::setNextFaces
(
    label& srcFaceI,
    label& tgtFaceI,
    const primitivePatch& srcPatch0,
    const primitivePatch& tgtPatch0,
    const boolList& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces
)
{
//    const labelList& srcNbrFaces = srcPatch0.pointFaces()[srcFaceI];
    const labelList& srcNbrFaces = srcPatch0.faceFaces()[srcFaceI];

    // set possible seeds for later use
    bool valuesSet = false;
    forAll(srcNbrFaces, i)
    {
        label faceS = srcNbrFaces[i];

        if (mapFlag[faceS] && seedFaces[faceS] == -1)
        {
            forAll(visitedFaces, j)
            {
                label faceT = visitedFaces[j];
                scalar area = interArea(faceS, faceT, srcPatch0, tgtPatch0);

                if (area > 0)
                {
                    // TODO - throwing area away - re-use in next iteration?

                    seedFaces[faceS] = faceT;

                    if (!valuesSet)
                    {
                        srcFaceI = faceS;
                        tgtFaceI = faceT;
                        valuesSet = true;
                    }
                }
            }
        }
    }

    // set next src and tgt faces if not set above
    if (valuesSet)
    {
        return;
    }
    else
    {
        // try to use existing seed
        bool foundNextSeed = false;
        for (label faceI = startSeedI_; faceI < mapFlag.size(); faceI++)
        {
            if (mapFlag[faceI])
            {
                if (!foundNextSeed)
                {
                    startSeedI_ = faceI;
                    foundNextSeed = true;
                }

                if (seedFaces[faceI] != -1)
                {
                    srcFaceI = faceI;
                    tgtFaceI = seedFaces[faceI];

                    return;
                }
            }
        }

        // perform new search to find match
        if (debug)
        {
            Info<< "Advancing front stalled: searching for new "
                << "target face" << endl;
        }

//        foundNextSeed = false;
        for (label faceI = startSeedI_; faceI < mapFlag.size(); faceI++)
        {
            if (mapFlag[faceI])
            {
                if (!foundNextSeed)
                {
                    startSeedI_ = faceI + 1;
                    foundNextSeed = true;
                }

                srcFaceI = faceI;
                tgtFaceI = findTargetFace(srcFaceI, srcPatch0, tgtPatch0);

                if (tgtFaceI >= 0)
                {
                    return;
                }
            }
        }

        FatalErrorIn
        (
            "void Foam::cyclicAMIPolyPatch::setNextFaces"
            "("
                "label&, "
                "label&, "
                "const primitivePatch&, "
                "const primitivePatch&, "
                "const boolList&, "
                "labelList&, "
                "const DynamicList<label>&"
            ") const"
        )  << "Unable to set source and target faces" << abort(FatalError);
    }
}


template<class SourcePatch, class TargetPatch>
Foam::scalar Foam::AMIInterpolation<SourcePatch, TargetPatch>::interArea
(
    const label srcFaceI,
    const label tgtFaceI,
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
) const
{
    const pointField& srcPoints = srcPatch.points();
    const pointField& tgtPoints = tgtPatch.points();

    const face& src = srcPatch[srcFaceI];
    const face& tgt = tgtPatch[tgtFaceI];

    // quick reject if either face has zero area
    if ((src.mag(srcPoints) < ROOTVSMALL) || (tgt.mag(tgtPoints) < ROOTVSMALL))
    {
        return 0.0;
    }

    // create intersection object
    faceAreaIntersect inter(srcPoints, tgtPoints);

    // crude resultant norm
    const vector n = 0.5*(tgt.normal(tgtPoints) - src.normal(srcPoints));

    scalar area = inter.calc(src, tgt, n, triMode_);

    if ((debug > 1) && (area > 0))
    {
        writeIntersectionOBJ(area, src, tgt, srcPoints, tgtPoints);
    }

    return area;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::calcAddressing
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    label srcFaceI,
    label tgtFaceI
)
{
    if (!srcPatch.size() || !tgtPatch.size())
    {
        if (debug)
        {
            Pout<< "AMI: Patches not on processor: Source faces = "
                << srcPatch.size() << ", target faces = " << tgtPatch.size()
                << endl;
        }

        return;
    }


    if (debug)
    {
        Info<< "AMI: calcAddressing" << endl;
        writePatch(srcPatch, "VTK", "source");
        writePatch(tgtPatch, "VTK", "target");
    }

    // temporary storage for addressing and weights
    List<DynamicList<label> > srcAddr(srcPatch.size());
    List<DynamicList<scalar> > srcWght(srcPatch.size());
    List<DynamicList<label> > tgtAddr(tgtPatch.size());
    List<DynamicList<scalar> > tgtWght(tgtPatch.size());


    // find initial face match using brute force/octree search
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ((srcFaceI == -1) || (tgtFaceI == -1))
    {
        srcFaceI = 0;
        tgtFaceI = 0;
        bool foundFace = false;
        forAll(srcPatch, faceI)
        {
            tgtFaceI = findTargetFace(faceI, srcPatch, tgtPatch);
            if (tgtFaceI >= 0)
            {
                srcFaceI = faceI;
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorIn
            (
                "void Foam::AMIInterpolation<SourcePatch, TargetPatch>::"
                "calcAddressing"
                "("
                    "const primitivePatch&, "
                    "const primitivePatch&, "
                    "label, "
                    "label"
                ")"
            )   << "Unable to find initial target face" << abort(FatalError);
        }
    }

    if (debug)
    {
        Info<< "AMI: initial target face = " << tgtFaceI << endl;
    }


    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    label nFacesRemaining = srcPatch.size();

    // list of tgt face neighbour faces
    DynamicList<label> nbrFaces(10);

    // list of faces currently visited for srcFaceI to avoid multiple hits
    DynamicList<label> visitedFaces(10);

    // list to keep track of tgt faces used to seed src faces
    labelList seedFaces(nFacesRemaining, -1);
    seedFaces[srcFaceI] = tgtFaceI;

    // list to keep track of whether src face can be mapped
    boolList mapFlag(nFacesRemaining, true);

    // reset starting seed
    startSeedI_ = 0;

    label nNonOverlap = 0;
    do
    {
        nbrFaces.clear();
        visitedFaces.clear();

        // append initial target face and neighbours
        nbrFaces.append(tgtFaceI);
        appendNbrFaces(tgtFaceI, tgtPatch, visitedFaces, nbrFaces);

        bool faceProcessed = false;

        do
        {
            // process new target face
            tgtFaceI = nbrFaces.remove();
            visitedFaces.append(tgtFaceI);
            scalar area = interArea(srcFaceI, tgtFaceI, srcPatch, tgtPatch);

            // store when intersection area > 0
            if (area > 0)
            {
                srcAddr[srcFaceI].append(tgtFaceI);
                srcWght[srcFaceI].append(area);

                tgtAddr[tgtFaceI].append(srcFaceI);
                tgtWght[tgtFaceI].append(area);

                appendNbrFaces(tgtFaceI, tgtPatch, visitedFaces, nbrFaces);

                faceProcessed = true;
            }

        } while (nbrFaces.size() > 0);

        mapFlag[srcFaceI] = false;

        nFacesRemaining--;

        if (!faceProcessed)
        {
            nNonOverlap++;
        }

        // choose new src face from current src face neighbour
        if (nFacesRemaining > 0)
        {
            setNextFaces
            (
                srcFaceI,
                tgtFaceI,
                srcPatch,
                tgtPatch,
                mapFlag,
                seedFaces,
                visitedFaces
            );
        }
    } while (nFacesRemaining > 0);

    if (nNonOverlap != 0)
    {
        Pout<< "AMI: " << nNonOverlap << " non-overlap faces identified"
            << endl;
    }

    // transfer data to persistent storage
    srcAddress_.setSize(srcPatch.size());
    srcWeights_.setSize(srcPatch.size());
    forAll(srcAddr, i)
    {
        srcAddress_[i].transfer(srcAddr[i]);
        srcWeights_[i].transfer(srcWght[i]);
    }

    tgtAddress_.setSize(tgtPatch.size());
    tgtWeights_.setSize(tgtPatch.size());
    forAll(tgtAddr, i)
    {
        tgtAddress_[i].transfer(tgtAddr[i]);
        tgtWeights_[i].transfer(tgtWght[i]);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::normaliseWeights
(
    const primitivePatch& patch,
    const word& patchName,
    const labelListList& addr,
    scalarListList& wght,
    const bool output
)
{
    scalarList wghtSum(wght.size(), 0.0);

    scalar minBound = VGREAT;
    scalar maxBound = -VGREAT;

    // Normalise the weights
    forAll(wght, faceI)
    {
        scalar s = sum(wght[faceI]);
        wghtSum[faceI] = s;

        scalar t = s/patch[faceI].mag(patch.points());
        if (t < minBound)
        {
            minBound = t;
        }

        if (t > maxBound)
        {
            maxBound = t;
        }

        forAll(addr[faceI], i)
        {
            wght[faceI][i] /= s;
        }
    }

    if (output)
    {
        Info<< "AMI: Patch " << patchName << " weights min/max = "
            << returnReduce(minBound, minOp<scalar>()) << ", "
            << returnReduce(maxBound, maxOp<scalar>()) << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode
)
:
    srcAddress_(),
    srcWeights_(),
    tgtAddress_(),
    tgtWeights_(),
    startSeedI_(0),
    triMode_(triMode),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    label srcSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtSize = returnReduce(tgtPatch.size(), sumOp<label>());

    Info<< "AMI: Creating addressing and weights between "
        << srcSize << " source faces and " << tgtSize << " target faces"
        << endl;

    if (surfPtr.valid())
    {
        // create new patches for source and target
        pointField srcPoints = srcPatch.points();
        primitivePatch srcPatch0
        (
            SubList<face>
            (
                srcPatch,
                srcPatch.size(),
                0
            ),
            srcPoints
        );

        if (debug)
        {
            OFstream os("amiSrcPoints.obj");
            forAll(srcPoints, i)
            {
                meshTools::writeOBJ(os, srcPoints[i]);
            }
        }

        pointField tgtPoints = tgtPatch.points();
        primitivePatch tgtPatch0
        (
            SubList<face>
            (
                tgtPatch,
                tgtPatch.size(),
                0
            ),
            tgtPoints
        );

        if (debug)
        {
            OFstream os("amiTgtPoints.obj");
            forAll(tgtPoints, i)
            {
                meshTools::writeOBJ(os, tgtPoints[i]);
            }
        }


        // map source and target patches onto projection surface
        projectPointsToSurface(surfPtr(), srcPoints);
        projectPointsToSurface(surfPtr(), tgtPoints);


        // calculate AMI interpolation
        update(srcPatch0, tgtPatch0);
    }
    else
    {
        update(srcPatch, tgtPatch);
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::~AMIInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::update
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
)
{
    static label patchI = 0;

    if (Pstream::parRun())
    {
        // convert local addressing to global addressing
        globalIndex globalSrcFaces(srcPatch.size());
        globalIndex globalTgtFaces(tgtPatch.size());

        // Create processor map of overlapping faces. This map gets
        // (possibly remote) faces from the tgtPatch such that they (together)
        // cover all of the srcPatch
        autoPtr<mapDistribute> mapPtr = calcProcMap(srcPatch, tgtPatch);
        const mapDistribute& map = mapPtr();

        // create new target patch that fully encompasses source patch

        // faces and points
        faceList newTgtFaces;
        pointField newTgtPoints;
        // original faces from tgtPatch (in globalIndexing since might be
        // remote)
        labelList tgtFaceIDs;
        distributeAndMergePatches
        (
            map,
            tgtPatch,
            globalTgtFaces,
            newTgtFaces,
            newTgtPoints,
            tgtFaceIDs
        );

        primitivePatch
            newTgtPatch
            (
                SubList<face>
                (
                    newTgtFaces,
                    newTgtFaces.size()
                ),
                newTgtPoints
            );

        checkPatches(srcPatch, newTgtPatch);


        // calculate AMI interpolation
        calcAddressing(srcPatch, newTgtPatch);

        // Now
        // ~~~
        //  srcAddress_ :   per srcPatch face a list of the newTgtPatch (not
        //                  tgtPatch) faces it overlaps
        //  tgtAddress_ :   per newTgtPatch (not tgtPatch) face a list of the
        //                  srcPatch faces it overlaps


        // Rework newTgtPatch indices into globalIndices of tgtPatch
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
        forAll(srcAddress_, i)
        {
            labelList& addressing = srcAddress_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = tgtFaceIDs[addressing[addrI]];
            }
        }

        forAll(tgtAddress_, i)
        {
            labelList& addressing = tgtAddress_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = globalSrcFaces.toGlobal(addressing[addrI]);
            }
        }

        // send data back to originating procs. Note that contributions
        // from different processors get added (ListPlusEqOp)

        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            map.subMap(),
            tgtAddress_,
            ListPlusEqOp<label>(),
            labelList()
        );

        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            map.subMap(),
            tgtWeights_,
            ListPlusEqOp<scalar>(),
            scalarList()
        );

        // weights normalisation
        normaliseWeights(srcPatch, "source", srcAddress_, srcWeights_, true);
        normaliseWeights(tgtPatch, "target", tgtAddress_, tgtWeights_, true);

        // cache maps and reset addresses
        List<Map<label> > cMap;
        srcMapPtr_.reset(new mapDistribute(globalSrcFaces, tgtAddress_, cMap));
        tgtMapPtr_.reset(new mapDistribute(globalTgtFaces, srcAddress_, cMap));


        if (debug)
        {
            writeWeights(srcWeights_, srcPatch, "VTK", "source");
            writeWeights(tgtWeights_, tgtPatch, "VTK", "target");
            writeFaceConnectivity(srcPatch, newTgtPatch, srcAddress_);
        }
    }
    else
    {
        checkPatches(srcPatch, tgtPatch);

        calcAddressing(srcPatch, tgtPatch);

        if (debug)
        {
            writeWeights(srcWeights_, srcPatch, "VTK", "source");
            writeWeights(tgtWeights_, tgtPatch, "VTK", "target");
        }

        normaliseWeights(srcPatch, "source", srcAddress_, srcWeights_, true);
        normaliseWeights(tgtPatch, "target", tgtAddress_, tgtWeights_, true);
    }

    patchI++;
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld
) const
{
    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorIn
        (
            "AMIInterpolation::interpolateToSource(const Field<Type>) const"
        )   << "Supplied field size is not equal to target patch size. "
            << "Target patch = " << tgtAddress_.size() << ", supplied field = "
            << fld.size() << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            srcAddress_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    if (Pstream::parRun())
    {
        const mapDistribute& map = tgtMapPtr_();

        Field<Type> work(fld);
        map.distribute(work);

        forAll(result, faceI)
        {
            const labelList& faces = srcAddress_[faceI];
            const scalarList& weights = srcWeights_[faceI];

            forAll(faces, i)
            {
                result[faceI] += work[faces[i]]*weights[i];
            }
        }
    }
    else
    {
        forAll(result, faceI)
        {
            const labelList& faces = srcAddress_[faceI];
            const scalarList& weights = srcWeights_[faceI];

            forAll(faces, i)
            {
                result[faceI] += fld[faces[i]]*weights[i];
            }
        }
    }

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type> >& tFld
) const
{
    return interpolateToSource(tFld());
}



template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld
) const
{
    if (fld.size() != srcAddress_.size())
    {
        FatalErrorIn
        (
            "AMIInterpolation::interpolateToSource(const Field<Type>&) const"
        )   << "Supplied field size is not equal to source patch size. "
            << "Source patch = " << srcAddress_.size() << ", supplied field = "
            << fld.size() << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            tgtAddress_.size(),
            pTraits<Type>::zero
        )
    );

    Field<Type>& result = tresult();

    if (Pstream::parRun())
    {
        const mapDistribute& map = srcMapPtr_();

        Field<Type> work(fld);
        map.distribute(work);

        forAll(result, faceI)
        {
            const labelList& faces = tgtAddress_[faceI];
            const scalarList& weights = tgtWeights_[faceI];

            forAll(faces, i)
            {
                result[faceI] += work[faces[i]]*weights[i];
            }
        }
    }
    else
    {
        forAll(result, faceI)
        {
            const labelList& faces = tgtAddress_[faceI];
            const scalarList& weights = tgtWeights_[faceI];

            forAll(faces, i)
            {
                result[faceI] += fld[faces[i]]*weights[i];
            }
        }
    }

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type> >& tFld
) const
{
    return interpolateToTarget(tFld());
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writeFaceConnectivity
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const labelListList& srcAddress
)
const
{
    OFstream os("faceConnectivity" + Foam::name(Pstream::myProcNo()) + ".obj");

    label ptI = 1;

    forAll(srcAddress, i)
    {
        const labelList& addr = srcAddress[i];
        const point& srcPt = srcPatch.faceCentres()[i];
        forAll(addr, j)
        {
            label tgtPtI = addr[j];
            const point& tgtPt = tgtPatch.faceCentres()[tgtPtI];

            meshTools::writeOBJ(os, srcPt);
            meshTools::writeOBJ(os, tgtPt);

            os  << "l " << ptI << " " << ptI + 1 << endl;

            ptI += 2;
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writeWeights
(
    const scalarListList& weights,
    const primitivePatch& patch,
    const word& folder,
    const word& prefix
)
const
{
    static label i = 0;

    scalarField wghtSum(weights.size(), 0.0);

    forAll(weights, faceI)
    {
        scalar s = sum(weights[faceI]);
        wghtSum[faceI] = s;
    }

    vtkSurfaceWriter writer;

    writer.write
    (
        folder,
        prefix
          + '_' + Foam::name(i) + "_proc" + Foam::name(Pstream::myProcNo()),
        patch.localPoints(),
        patch.localFaces(),
        "weights",
        wghtSum,
        false
    );

    i++;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writePatch
(
    const primitivePatch& patch,
    const word& folder,
    const word& prefix
)
const
{
    static label i = 0;

    vtkSurfaceWriter writer;
    writer.write
    (
        folder,
        prefix
          + '_' + Foam::name(i) + "_proc" + Foam::name(Pstream::myProcNo()),
        patch.localPoints(),
        patch.localFaces(),
        "AMIPatch",
        scalarField(patch.size(), Pstream::myProcNo()),
        false
    );

    i++;
}


// ************************************************************************* //
