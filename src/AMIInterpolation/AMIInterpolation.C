/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
#include "Random.H"
#include "treeDataPrimitivePatch.H"
#include "indexedOctree.H"
#include "primitivePatch.H"
#include "meshTools.H"

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

    // sanity checks
    boundBox bbSrc(srcPatch.points(), srcPatch.meshPoints());
    boundBox bbTgt(tgtPatch.points(), tgtPatch.meshPoints());
    boundBox bbSurf(bbTgt); // USE POINTS OF SURFACE!!!!!!!!!!!!!!!!!!


    // projection surface bounds - check against bounds of source and target

    bbSurf.inflate(maxBoundsError);

    if (!bbSurf.contains(bbSrc))
    {
        WarningIn
        (
            "AMIInterpolation<SourcePatch, TargetPatch>::checkPatches"
            "("
                "const primitivePatch&, "
                "const primitivePatch&"
            ")"
        )   << "Source patch is larger than, or misaligned with the "
            << "projection surface" << endl;
    }

    if (!bbSurf.contains(bbTgt))
    {
        WarningIn
        (
            "AMIInterpolation<SourcePatch, TargetPatch>::checkPatches"
            "("
                "const primitivePatch&, "
                "const primitivePatch&"
            ")"
        )   << "Target patch is larger than, or misaligned with the "
            << "projection surface" << endl;
    }


    // check bounds of source and target
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
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::projectPointsToSurface
(
    const searchableSurface& surf,
    pointField& pts
) const
{
    if (!projectPoints_)
    {
        return;
    }

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

    Random rndGen(123456);

    treeBoundBox bb(tgtPatch.points());

    typedef treeDataPrimitivePatch<face, SubList, const pointField&> treeType;

    indexedOctree<treeType> tree
    (
        treeType(false, tgtPatch),
        bb.extend(rndGen, 1E-4),                // overall search domain
        8,                                      // maxLevel
        10,                                     // leafsize
        3.0                                     // duplicity
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
    else
    {
        FatalErrorIn
        (
            "Foam::label Foam::cyclicAMIPolyPatch::findTargetFace"
            "("
                "const label, "
                "const primitivePatch&"
                "const primitivePatch&"
            ") const"
        )   << "Unable to find target face for source face" << srcFaceI
            << abort(FatalError);
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
                    startSeedI_ = faceI + 1;
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

        for (label faceI = startSeedI_; faceI < mapFlag.size(); faceI++)
        {
            if (mapFlag[faceI])
            {
                startSeedI_ = faceI + 1;

                srcFaceI = faceI;
                tgtFaceI = findTargetFace(srcFaceI, srcPatch0, tgtPatch0);

                return;
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
                "const labelList&, "
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
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::normaliseWeights
(
    const word& patchType,
    const primitivePatch& patch,
    const pointField& points,
    const List<DynamicList<label> >& addr,
    List<DynamicList<scalar> >& wght,
    const bool output
)
{
    scalarField wghtSum(patch.size(), 0.0);

    // normalise weights by face areas
    forAll(wght, faceI)
    {
        const DynamicList<label>& addressing = addr[faceI];
        forAll(addressing, addrI)
        {
            const scalar faceArea = patch[faceI].mag(points);
            wght[faceI][addrI] /= faceArea;
            wghtSum[faceI] += wght[faceI][addrI];
        }
    }

//    if (debug)
    {
        vtkSurfaceWriter writer;

        writer.write
        (
            "VTK",
            patchType,
            patch.localPoints(),
            patch.localFaces(),
            "weights",
            wghtSum,
            false
        );
    }

    if (output)
    {
        Info<< "AMI: Patch " << patchType << " weights min/max = "
            << min(wghtSum) << ", " << max(wghtSum) << endl;
    }


    forAll(wght, faceI)
    {
        const DynamicList<label>& addressing = addr[faceI];
        forAll(addressing, addrI)
        {
            wght[faceI][addrI] /= wghtSum[faceI];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const searchableSurface& surf,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool projectPoints
)
:
    srcAddress_(),
    srcWeights_(),
    tgtAddress_(),
    tgtWeights_(),
    startSeedI_(0),
    triMode_(triMode),
    projectPoints_(projectPoints)
{
    if (srcPatch.size() && tgtPatch.size())
    {
        checkPatches(srcPatch, tgtPatch);

        Info<< "AMI: Creating interpolation addressing and weights" << nl
            << "    Source faces = " << srcPatch.size() << nl
            << "    Target faces = " << tgtPatch.size() << endl;

        calcAddressing(srcPatch, tgtPatch, surf);
    }
}



// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::~AMIInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::calcAddressing
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const searchableSurface& surf
)
{
    if (debug)
    {
        Info<< "AMI: calcAddressing" << endl;
    }

    // temporary storage for addressing and weights
    List<DynamicList<label> > srcAddr(srcPatch.size());
    List<DynamicList<scalar> > srcWght(srcPatch.size());
    List<DynamicList<label> > tgtAddr(tgtPatch.size());
    List<DynamicList<scalar> > tgtWght(tgtPatch.size());


    // create new patches for source and target
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    projectPointsToSurface(surf, srcPoints);
    projectPointsToSurface(surf, tgtPoints);


    // find initial face match using brute force/octree search
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    label srcFaceI = 0;
    label tgtFaceI = findTargetFace(srcFaceI, srcPatch0, tgtPatch0);

    if (debug)
    {
        Info<< "AMI: initial target face = " << tgtFaceI << endl;
    }


    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    label nFacesRemaining = srcPatch0.size();

    // list of tgt face neighbour faces
    DynamicList<label> nbrFaces(10);

    // list of faces currently visited for srcFaceI to avoid multiple hits
    DynamicList<label> visitedFaces(10);

    // list to keep track of tgt faces used to seed src faces
    labelList seedFaces(nFacesRemaining, -1);
    seedFaces[srcFaceI] = tgtFaceI;

    // list to keep track of whether src face can be mapped
    boolList mapFlag(nFacesRemaining, true);

    label nFail = 0;

    do
    {
        nbrFaces.clear();
        visitedFaces.clear();

        // append initial target face and neighbours
        nbrFaces.append(tgtFaceI);
        appendNbrFaces(tgtFaceI, tgtPatch0, visitedFaces, nbrFaces);

        bool faceProcessed = false;

        do
        {
            // process new target face
            tgtFaceI = nbrFaces.remove();
            visitedFaces.append(tgtFaceI);
            scalar area = interArea(srcFaceI, tgtFaceI, srcPatch0, tgtPatch0);

            // store when intersection area > 0
            if (area > 0)
            {
                srcAddr[srcFaceI].append(tgtFaceI);
                srcWght[srcFaceI].append(area);

                tgtAddr[tgtFaceI].append(srcFaceI);
                tgtWght[tgtFaceI].append(area);

                appendNbrFaces(tgtFaceI, tgtPatch0, visitedFaces, nbrFaces);

                faceProcessed = true;
            }

        } while (nbrFaces.size() > 0);

        if (faceProcessed)
        {
            mapFlag[srcFaceI] = false;

            nFacesRemaining--;
        }
        else
        {
            nFail++;

            if (nFail > 1000)
            {
                Info<< "nFail = " << nFail << endl;
            }
        }

        // choose new src face from current src face neighbour
        if (nFacesRemaining > 0)
        {
            setNextFaces
            (
                srcFaceI,
                tgtFaceI,
                srcPatch0,
                tgtPatch0,
                mapFlag,
                seedFaces,
                visitedFaces
            );
        }
    } while (nFacesRemaining > 0);


    // weights normalisation
    normaliseWeights("source", srcPatch0, srcPoints, srcAddr, srcWght, true);
    normaliseWeights("target", tgtPatch0, tgtPoints, tgtAddr, tgtWght, true);

    // transfer data to persistent storage
    srcAddress_.setSize(srcAddr.size());
    srcWeights_.setSize(srcWght.size());
    forAll(srcAddr, i)
    {
        srcAddress_[i].transfer(srcAddr[i]);
        srcWeights_[i].transfer(srcWght[i]);
    }

    tgtAddress_.setSize(tgtAddr.size());
    tgtWeights_.setSize(tgtWght.size());
    forAll(tgtAddr, i)
    {
        tgtAddress_[i].transfer(tgtAddr[i]);
        tgtWeights_[i].transfer(tgtWght[i]);
    }

// ADD PARALLEL BITS
/*
    if (Pstream::parRun())
    {
        // Source fields - addesses and weights
        const label nSrc = srcAddess.size();
        const globalIndex globalSrcFaces(nSrc);

        List<Map<label> > compactMapSrcAddress;
        labelListList transSrcAddress;

        mapDistribute srcMapDist
        (
            globalSrcFaces,
            srcAddress_,
            globalIndexAndTransform(mesh), // mesh not available!!!!!!!!
            transSrcAddress,
            compactMapSrcAddress
        );


    }
*/
    // reset starting seed
    startSeedI_ = 0;
}


template<class SourcePatch, class TargetPatch>
bool Foam::AMIInterpolation<SourcePatch, TargetPatch>::movePoints()
{
    return true;
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
            "AMIInterpolation::interpolateToSource(const Field<Type>)"
        )   << "Supplied field size is not equal to target patch size. "
            << "Target patch = " << tgtAddress_.size() << ", supplied field: "
            << fld.size() << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>(srcAddress_.size(), pTraits<Type>::zero)
    );

    Field<Type>& result = tresult();

    forAll(result, faceI)
    {
        const labelList& faces = srcAddress_[faceI];
        const scalarList& weights = srcWeights_[faceI];

        forAll(faces, i)
        {
            result[faceI] += fld[faces[i]]*weights[i];
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
            "AMIInterpolation::interpolateToSource(const Field<Type>)"
        )   << "Supplied field size is not equal to source patch size. "
            << "Source patch = " << srcAddress_.size() << ", supplied field: "
            << fld.size() << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>(tgtAddress_.size(), pTraits<Type>::zero)
    );

    Field<Type>& result = tresult();

    forAll(result, faceI)
    {
        const labelList& faces = tgtAddress_[faceI];
        const scalarList& weights = tgtWeights_[faceI];

        forAll(faces, i)
        {
            result[faceI] += fld[faces[i]]*weights[i];
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
    const primitivePatch& tgtPatch
)
const
{
    OFstream os("faceConnectivity.obj");

    label ptI = 1;

    forAll(srcAddress_, i)
    {
        const labelList& addr = srcAddress_[i];
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


// ************************************************************************* //
