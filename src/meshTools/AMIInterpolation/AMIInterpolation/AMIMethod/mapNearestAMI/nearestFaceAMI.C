/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "nearestFaceAMI.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::autoPtr<Foam::mapDistribute>
Foam::nearestFaceAMI<SourcePatch, TargetPatch>::calcFaceMap
(
    const List<nearestAndDist>& localInfo,
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
) const
{
    // Generate the list of processor bounding boxes
    List<boundBox> procBbs(Pstream::nProcs());
    procBbs[Pstream::myProcNo()] =
        boundBox(srcPatch.points(), srcPatch.meshPoints(), true);
    Pstream::gatherList(procBbs);
    Pstream::scatterList(procBbs);

    // Identify which of my local tgt faces intersect with each processor bb
    // within the current match's search distance
    const pointField& tgtCcs = tgtPatch.faceCentres();
    List<DynamicList<label>> dynSendMap(Pstream::nProcs());

    forAll(localInfo, tgtFacei)
    {
        const scalar r2 = localInfo[tgtFacei].second();

        // Construct local bounding box to test against processor bb
        forAll(procBbs, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                if (procBbs[proci].overlaps(tgtCcs[tgtFacei], r2))
                {
                    dynSendMap[proci].append(tgtFacei);
                }
            }
        }
    }

    // Convert dynamicList to labelList
    labelListList sendMap(Pstream::nProcs());
    forAll(sendMap, proci)
    {
        dynSendMap[proci].shrink();
        sendMap[proci].transfer(dynSendMap[proci]);

        if (debug)
        {
            Pout<< "send map - to proc " << proci << " sending "
                << sendMap[proci].size() << " elements" << endl;
        }
    }

    return autoPtr<mapDistribute>::New(std::move(sendMap));
}


template<class SourcePatch, class TargetPatch>
Foam::autoPtr<Foam::mapDistribute>
Foam::nearestFaceAMI<SourcePatch, TargetPatch>::calcDistributed
(
    const SourcePatch& src,
    const TargetPatch& tgt,
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght
) const
{
    const auto tgtTreePtr = this->createTree(tgt);
    const auto& tgtTree = tgtTreePtr();

    // Create global indexing for each patch
    globalIndex globalTgtCells(src.size());

    // First pass: determine local matches

    // Identify local nearest matches
    pointField srcCcs(src.faceCentres());

    List<nearestAndDist> localInfo(src.size());
    forAll(srcCcs, srcCelli)
    {
        const point& srcCc = srcCcs[srcCelli];

        pointIndexHit& test = localInfo[srcCelli].first();
        test = tgtTree.findNearest(srcCc, GREAT);

        if (test.hit())
        {
            // With a search radius2 of GREAT all cells should receive a hit
            localInfo[srcCelli].second() = magSqr(srcCc - test.hitPoint());
            test.setIndex(globalTgtCells.toGlobal(test.index()));
        }
    }

    // Second pass: determine remote matches

    autoPtr<mapDistribute> mapPtr = calcFaceMap(localInfo, src, tgt);
    mapDistribute& map = mapPtr();

    List<nearestAndDist> remoteInfo(localInfo);
    map.distribute(remoteInfo);

    // Note: re-using srcCcs
    map.distribute(srcCcs);

    // Test remote target cells against local source cells
    nearestAndDist testInfo;
    pointIndexHit& test = testInfo.first();
    forAll(remoteInfo, i)
    {
        test = tgtTree.findNearest(srcCcs[i], remoteInfo[i].second());
        if (test.hit())
        {
            testInfo.first().setIndex
            (
                globalTgtCells.toGlobal(test.index())
            );
            testInfo.second() = magSqr(test.hitPoint() - srcCcs[i]);
            nearestEqOp()(remoteInfo[i], testInfo);
        }
    }

    // Send back to originating processor. Choose best if sent to multiple
    // processors. Note that afterwards all unused entries have the unique
    // value nearestZero (distance < 0). This is used later on to see if
    // the sample was sent to any processor.
    const nearestAndDist nearestZero(pointIndexHit(), -GREAT);
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(0),
        src.size(),
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        remoteInfo,
        nearestEqOp(),
        noOp(),             // no flipping
        nearestZero
    );

    // Third pass: combine local and remote info and filter out any
    // connections that are further away than threshold distance squared
    srcToTgtAddr.setSize(src.size());
    srcToTgtWght.setSize(src.size());
    forAll(srcToTgtAddr, srcCelli)
    {
        nearestEqOp()(localInfo[srcCelli], remoteInfo[srcCelli]);
        if (localInfo[srcCelli].second() < maxDistance2_)
        {
            const label tgtCelli = localInfo[srcCelli].first().index();
            srcToTgtAddr[srcCelli] = labelList(1, tgtCelli);
            srcToTgtWght[srcCelli] = scalarList(1, 1.0);
        }
    }

    List<Map<label>> cMap;
    return autoPtr<mapDistribute>::New(globalTgtCells, srcToTgtAddr, cMap);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::nearestFaceAMI<SourcePatch, TargetPatch>::nearestFaceAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch
)
:
    AMIMethod<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        triMode,
        reverseTarget,
        requireMatch
    ),
    maxDistance2_(GREAT)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
bool Foam::nearestFaceAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    pointListList& srcCentroids,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    scalarList& srcMagSf,
    scalarList& tgtMagSf,
    autoPtr<mapDistribute>& srcMapPtr,
    autoPtr<mapDistribute>& tgtMapPtr,
    label srcFacei,
    label tgtFacei
)
{
bool symmetric_ = true;

    if (this->distributed())
    {
        tgtMapPtr =
            calcDistributed
            (
                this->srcPatch0_,
                this->tgtPatch0_,
                srcAddress,
                srcWeights
            );

        if (symmetric_)
        {
            srcMapPtr =
                calcDistributed
                (
                    this->tgtPatch0_,
                    this->srcPatch0_,
                    tgtAddress,
                    tgtWeights
                );
            }
    }
    else
    {
        srcAddress.setSize(this->srcPatch0_.size());
        srcWeights.setSize(this->srcPatch0_.size());

        if (symmetric_)
        {
            tgtAddress.setSize(this->tgtPatch0_.size());
            tgtWeights.setSize(this->tgtPatch0_.size());
        }

        const pointField& srcCcs = this->srcPatch0_.faceCentres();
        const pointField& tgtCcs = this->tgtPatch0_.faceCentres();

        const auto tgtTreePtr = this->createTree(this->tgtPatch0_);
        const auto& tgtTree = tgtTreePtr();

        forAll(srcCcs, srcFacei)
        {
            const point& srcCc = srcCcs[srcFacei];
            const pointIndexHit hit = tgtTree.findNearest(srcCc, GREAT);

            if
            (
                hit.hit()
             && (magSqr(srcCc - tgtCcs[hit.index()]) < maxDistance2_)
            )
            {
                label tgtFacei = hit.index();
                srcAddress[srcFacei] = labelList(1, tgtFacei);
                srcWeights[srcFacei] = scalarList(1, 1.0);

                if (symmetric_)
                {
                    tgtAddress[tgtFacei] = labelList(1, srcFacei);
                    tgtWeights[tgtFacei] = scalarList(1, 1.0);
                }
            }
            else
            {
                if (debug)
                {
                    WarningInFunction
                        << "Unable to find target face for source face "
                        << srcFacei << endl;
                }
            }
        }

        if (symmetric_)
        {
            const auto srcTreePtr = this->createTree(this->srcPatch0_);
            const auto& srcTree = srcTreePtr();

            // Check that all source cells have connections and populate any
            // missing entries
            forAll(tgtWeights, tgtCelli)
            {
                if (tgtAddress[tgtCelli].empty())
                {
                    const point& tgtCc = tgtCcs[tgtCelli];
                    pointIndexHit hit = srcTree.findNearest(tgtCc, GREAT);

                    if
                    (
                        hit.hit()
                     && (magSqr(tgtCc - srcCcs[hit.index()]) < maxDistance2_)
                    )
                    {
                        tgtAddress[tgtCelli] = labelList(1, hit.index());
                        tgtWeights[tgtCelli] = scalarList(1, 1.0);
                    }
                }
                else
                {
                    if (debug)
                    {
                        WarningInFunction
                            << "Unable to find source face for target face "
                            << tgtCelli << endl;
                    }
                }
            }
        }
    }

    return true;
}


template<class SourcePatch, class TargetPatch>
void Foam::nearestFaceAMI<SourcePatch, TargetPatch>::normaliseWeights
(
    const bool verbose,
    AMIInterpolation<SourcePatch, TargetPatch>& inter
)
{
    // Do nothing - weights already 1-to-1 and normalised
}


// ************************************************************************* //
