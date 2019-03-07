/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "faceAreaWeightAMI.H"
#include "profiling.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::calcAddressing
(
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<point>>& srcCtr,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght,
    label srcFacei,
    label tgtFacei
)
{
    addProfiling(ami, "faceAreaWeightAMI::calcAddressing");

    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nFacesRemaining = srcAddr.size();

    // list of tgt face neighbour faces
    DynamicList<label> nbrFaces(10);

    // list of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label> visitedFaces(10);

    // list to keep track of tgt faces used to seed src faces
    labelList seedFaces(nFacesRemaining, -1);
    seedFaces[srcFacei] = tgtFacei;

    // list to keep track of whether src face can be mapped
    bitSet mapFlag(nFacesRemaining, true);

    // reset starting seed
    label startSeedi = 0;

    DynamicList<label> nonOverlapFaces;
    do
    {
        nbrFaces.clear();
        visitedFaces.clear();

        // Do advancing front starting from srcFacei,tgtFacei
        bool faceProcessed = processSourceFace
        (
            srcFacei,
            tgtFacei,

            nbrFaces,
            visitedFaces,

            srcAddr,
            srcWght,
            srcCtr,
            tgtAddr,
            tgtWght
        );

        mapFlag.unset(srcFacei);

        nFacesRemaining--;

        if (!faceProcessed)
        {
            nonOverlapFaces.append(srcFacei);
        }

        // choose new src face from current src face neighbour
        if (nFacesRemaining > 0)
        {
            setNextFaces
            (
                startSeedi,
                srcFacei,
                tgtFacei,
                mapFlag,
                seedFaces,
                visitedFaces
            );
        }
    } while (nFacesRemaining > 0);

    this->srcNonOverlap_.transfer(nonOverlapFaces);
}


template<class SourcePatch, class TargetPatch>
bool Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::processSourceFace
(
    const label srcFacei,
    const label tgtStartFacei,

    // list of tgt face neighbour faces
    DynamicList<label>& nbrFaces,
    // list of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label>& visitedFaces,

    // temporary storage for addressing, weights and centroid
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<point>>& srcCtr,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght
)
{
    addProfiling(ami, "faceAreaWeightAMI::processSourceFace");

    if (tgtStartFacei == -1)
    {
        return false;
    }

    // append initial target face and neighbours
    nbrFaces.append(tgtStartFacei);
    this->appendNbrFaces
    (
        tgtStartFacei,
        this->tgtPatch_,
        visitedFaces,
        nbrFaces
    );

    bool faceProcessed = false;

    label maxNeighbourFaces = nbrFaces.size();

    do
    {
        // process new target face
        label tgtFacei = nbrFaces.remove();
        visitedFaces.append(tgtFacei);

        scalar interArea = 0;
        vector interCentroid(Zero);
        calcInterArea(srcFacei, tgtFacei, interArea, interCentroid);

        // store when intersection fractional area > tolerance
        if
        (
            interArea/this->srcMagSf_[srcFacei]
          > faceAreaIntersect::tolerance()
        )
        {
            srcAddr[srcFacei].append(tgtFacei);
            srcWght[srcFacei].append(interArea);
            srcCtr[srcFacei].append(interCentroid);

            tgtAddr[tgtFacei].append(srcFacei);
            tgtWght[tgtFacei].append(interArea);

            this->appendNbrFaces
            (
                tgtFacei,
                this->tgtPatch_,
                visitedFaces,
                nbrFaces
            );

            faceProcessed = true;

            maxNeighbourFaces = max(maxNeighbourFaces, nbrFaces.size());
        }

    } while (nbrFaces.size() > 0);

    if (debug > 1)
    {
        DebugVar(maxNeighbourFaces);
    }

    return faceProcessed;
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::setNextFaces
(
    label& startSeedi,
    label& srcFacei,
    label& tgtFacei,
    const bitSet& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces,
    const bool errorOnNotFound
) const
{
    addProfiling(ami, "faceAreaWeightAMI::setNextFaces");

    const labelList& srcNbrFaces = this->srcPatch_.faceFaces()[srcFacei];

    // initialise tgtFacei
    tgtFacei = -1;

    // set possible seeds for later use
    bool valuesSet = false;
    for (label faceS: srcNbrFaces)
    {
        if (mapFlag.test(faceS) && seedFaces[faceS] == -1)
        {
            for (label faceT : visitedFaces)
            {
                const scalar threshold =
                    this->srcMagSf_[faceS]*faceAreaIntersect::tolerance();

                // Store when intersection fractional area > threshold
                if (overlaps(faceS, faceT, threshold))
                {
                    seedFaces[faceS] = faceT;

                    if (!valuesSet)
                    {
                        srcFacei = faceS;
                        tgtFacei = faceT;
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
        label facei = startSeedi;
        if (!mapFlag.test(startSeedi))
        {
            facei = mapFlag.find_next(facei);
        }
        const label startSeedi0 = facei;

        bool foundNextSeed = false;
        while (facei != -1)
        {
            if (!foundNextSeed)
            {
                startSeedi = facei;
                foundNextSeed = true;
            }

            if (seedFaces[facei] != -1)
            {
                srcFacei = facei;
                tgtFacei = seedFaces[facei];

                return;
            }

            facei = mapFlag.find_next(facei);
        }

        // perform new search to find match
        if (debug)
        {
            Pout<< "Advancing front stalled: searching for new "
                << "target face" << endl;
        }

        facei = startSeedi0;
        while (facei != -1)
        {
            srcFacei = facei;
            tgtFacei = this->findTargetFace(srcFacei, visitedFaces);

            if (tgtFacei >= 0)
            {
                return;
            }

            facei = mapFlag.find_next(facei);
        }

        if (errorOnNotFound)
        {
            FatalErrorInFunction
               << "Unable to set source and target faces" << abort(FatalError);
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::calcInterArea
(
    const label srcFacei,
    const label tgtFacei,
    scalar& area,
    vector& centroid
) const
{
    addProfiling(ami, "faceAreaWeightAMI::interArea");

    const pointField& srcPoints = this->srcPatch_.points();
    const pointField& tgtPoints = this->tgtPatch_.points();

    // references to candidate faces
    const face& src = this->srcPatch_[srcFacei];
    const face& tgt = this->tgtPatch_[tgtFacei];

    // quick reject if either face has zero area
    if
    (
        (this->srcMagSf_[srcFacei] < ROOTVSMALL)
     || (this->tgtMagSf_[tgtFacei] < ROOTVSMALL)
    )
    {
        return;
    }

    // create intersection object
    faceAreaIntersect inter
    (
        srcPoints,
        tgtPoints,
        this->srcTris_[srcFacei],
        this->tgtTris_[tgtFacei],
        this->reverseTarget_,
        AMIInterpolation<SourcePatch, TargetPatch>::cacheIntersections_
    );

    // crude resultant norm
    vector n(-this->srcPatch_.faceNormals()[srcFacei]);
    if (this->reverseTarget_)
    {
        n -= this->tgtPatch_.faceNormals()[tgtFacei];
    }
    else
    {
        n += this->tgtPatch_.faceNormals()[tgtFacei];
    }
    scalar magN = mag(n);

    if (magN > ROOTVSMALL)
    {
        inter.calc(src, tgt, n/magN, area, centroid);
    }
    else
    {
        WarningInFunction
            << "Invalid normal for source face " << srcFacei
            << " points " << UIndirectList<point>(srcPoints, src)
            << " target face " << tgtFacei
            << " points " << UIndirectList<point>(tgtPoints, tgt)
            << endl;
    }

    if
    (
        AMIInterpolation<SourcePatch, TargetPatch>::cacheIntersections_
     && debug
    )
    {
        static OBJstream tris("intersectionTris.obj");
        const auto& triPts = inter.triangles();
        for (const auto& tp : triPts)
        {
            tris.write(triPointRef(tp[0], tp[1], tp[2]), false);
        }
    }

    if ((debug > 1) && (area > 0))
    {
        this->writeIntersectionOBJ(area, src, tgt, srcPoints, tgtPoints);
    }
}


template<class SourcePatch, class TargetPatch>
bool Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::overlaps
(
    const label srcFacei,
    const label tgtFacei,
    const scalar threshold
) const
{
    const pointField& srcPoints = this->srcPatch_.points();
    const pointField& tgtPoints = this->tgtPatch_.points();

    // references to candidate faces
    const face& src = this->srcPatch_[srcFacei];
    const face& tgt = this->tgtPatch_[tgtFacei];

    // quick reject if either face has zero area
    if
    (
        (this->srcMagSf_[srcFacei] < ROOTVSMALL)
     || (this->tgtMagSf_[tgtFacei] < ROOTVSMALL)
    )
    {
        return false;
    }

    faceAreaIntersect inter
    (
        srcPoints,
        tgtPoints,
        this->srcTris_[srcFacei],
        this->tgtTris_[tgtFacei],
        this->reverseTarget_,
        AMIInterpolation<SourcePatch, TargetPatch>::cacheIntersections_
    );

    // crude resultant norm
    vector n(-this->srcPatch_.faceNormals()[srcFacei]);
    if (this->reverseTarget_)
    {
        n -= this->tgtPatch_.faceNormals()[tgtFacei];
    }
    else
    {
        n += this->tgtPatch_.faceNormals()[tgtFacei];
    }
    scalar magN = mag(n);

    if (magN > ROOTVSMALL)
    {
        return inter.overlaps(src, tgt, n/magN, threshold);
    }
    else
    {
        WarningInFunction
            << "Invalid normal for source face " << srcFacei
            << " points " << UIndirectList<point>(srcPoints, src)
            << " target face " << tgtFacei
            << " points " << UIndirectList<point>(tgtPoints, tgt)
            << endl;
    }

    return false;
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::
restartUncoveredSourceFace
(
    List<DynamicList<label>>& srcAddr,
    List<DynamicList<scalar>>& srcWght,
    List<DynamicList<point>>& srcCtr,
    List<DynamicList<label>>& tgtAddr,
    List<DynamicList<scalar>>& tgtWght
)
{
    addProfiling(ami, "faceAreaWeightAMI::restartUncoveredSourceFace");

    // Note: exclude faces in srcNonOverlap_ for ACMI?

    label nBelowMinWeight = 0;
    const scalar minWeight = 0.95;

    // list of tgt face neighbour faces
    DynamicList<label> nbrFaces(10);

    // list of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label> visitedFaces(10);

    forAll(srcWght, srcFacei)
    {
        const scalar s = sum(srcWght[srcFacei]);
        const scalar t = s/this->srcMagSf_[srcFacei];

        if (t < minWeight)
        {
            ++nBelowMinWeight;

            const face& f = this->srcPatch_[srcFacei];

            forAll(f, fpi)
            {
                const label tgtFacei =
                    this->findTargetFace(srcFacei, srcAddr[srcFacei], fpi);

                if (tgtFacei != -1)
                {
                    nbrFaces.clear();
                    visitedFaces = srcAddr[srcFacei];

                    (void)processSourceFace
                    (
                        srcFacei,
                        tgtFacei,

                        nbrFaces,
                        visitedFaces,

                        srcAddr,
                        srcWght,
                        srcCtr,
                        tgtAddr,
                        tgtWght
                    );
                }
            }
        }
    }

    if (debug && nBelowMinWeight)
    {
        WarningInFunction
            << "Restarted search on " << nBelowMinWeight
            << " faces since sum of weights < " << minWeight
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::faceAreaWeightAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch,
    const bool restartUncoveredSourceFace
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
    restartUncoveredSourceFace_(restartUncoveredSourceFace),
    srcTris_(),
    tgtTris_()
{
    this->triangulatePatch(srcPatch, srcTris_, this->srcMagSf_);
    this->triangulatePatch(tgtPatch, tgtTris_, this->tgtMagSf_);

    if (debug)
    {
        static label nAMI = 0;

        // Write out triangulated surfaces as OBJ files
        OBJstream srcTriObj("srcTris_" + Foam::name(nAMI) + ".obj");
        const pointField& srcPts = srcPatch.points();
        forAll(srcTris_, facei)
        {
            const DynamicList<face>& faces = srcTris_[facei];
            for (const face& f : faces)
            {
                srcTriObj.write
                (
                    triPointRef(srcPts[f[0]], srcPts[f[1]], srcPts[f[2]])
                );
            }
        }

        OBJstream tgtTriObj("tgtTris_" + Foam::name(nAMI) + ".obj");
        const pointField& tgtPts = tgtPatch.points();
        forAll(tgtTris_, facei)
        {
            const DynamicList<face>& faces = tgtTris_[facei];
            for (const face& f : faces)
            {
                tgtTriObj.write
                (
                    triPointRef(tgtPts[f[0]], tgtPts[f[1]], tgtPts[f[2]])
                );
            }
        }

        ++nAMI;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    pointListList& srcCentroids,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    label srcFacei,
    label tgtFacei
)
{
    addProfiling(ami, "faceAreaWeightAMI::calculate");

    bool ok =
        this->initialise
        (
            srcAddress,
            srcWeights,
            tgtAddress,
            tgtWeights,
            srcFacei,
            tgtFacei
        );

    srcCentroids.setSize(srcAddress.size());


    if (!ok)
    {
        return;
    }

    // temporary storage for addressing and weights
    List<DynamicList<label>> srcAddr(this->srcPatch_.size());
    List<DynamicList<scalar>> srcWght(srcAddr.size());
    List<DynamicList<point>> srcCtr(srcAddr.size());
    List<DynamicList<label>> tgtAddr(this->tgtPatch_.size());
    List<DynamicList<scalar>> tgtWght(tgtAddr.size());

    calcAddressing
    (
        srcAddr,
        srcWght,
        srcCtr,
        tgtAddr,
        tgtWght,
        srcFacei,
        tgtFacei
    );

    if (debug && !this->srcNonOverlap_.empty())
    {
        Pout<< "    AMI: " << this->srcNonOverlap_.size()
            << " non-overlap faces identified"
            << endl;
    }


    // Check for badly covered faces
    if (restartUncoveredSourceFace_)
    {
        restartUncoveredSourceFace
        (
            srcAddr,
            srcWght,
            srcCtr,
            tgtAddr,
            tgtWght
        );
    }

    // Transfer data to persistent storage
    forAll(srcAddr, i)
    {
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i].transfer(srcWght[i]);
        srcCentroids[i].transfer(srcCtr[i]);
    }
    forAll(tgtAddr, i)
    {
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i].transfer(tgtWght[i]);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::setMagSf
(
    const TargetPatch& tgtPatch,
    const mapDistribute& map,
    scalarList& srcMagSf,
    scalarList& tgtMagSf
) const
{
    srcMagSf = std::move(this->srcMagSf_);
    tgtMagSf = std::move(this->tgtMagSf_);
    map.reverseDistribute(tgtPatch.size(), tgtMagSf);
}


template<class SourcePatch, class TargetPatch>
void Foam::faceAreaWeightAMI<SourcePatch, TargetPatch>::normaliseWeights
(
    const bool verbose,
    AMIInterpolation<SourcePatch, TargetPatch>& inter
)
{
    inter.normaliseWeights(this->conformal(), verbose);
}


// ************************************************************************* //
