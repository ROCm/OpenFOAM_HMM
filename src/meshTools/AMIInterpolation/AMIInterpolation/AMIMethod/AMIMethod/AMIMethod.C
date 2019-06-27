/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "AMIMethod.H"
#include "meshTools.H"
#include "mapDistribute.H"
#include "unitConversion.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::checkPatches() const
{
    if (debug && (!srcPatch_.size() || !tgtPatch_.size()))
    {
        Pout<< "AMI: Patches not on processor: Source faces = "
            << srcPatch_.size() << ", target faces = " << tgtPatch_.size()
            << endl;
    }


    if (conformal())
    {
        const scalar maxBoundsError = 0.05;

        // check bounds of source and target
        boundBox bbSrc(srcPatch_.points(), srcPatch_.meshPoints(), true);
        boundBox bbTgt(tgtPatch_.points(), tgtPatch_.meshPoints(), true);

        boundBox bbTgtInf(bbTgt);
        bbTgtInf.inflate(maxBoundsError);

        if (!bbTgtInf.contains(bbSrc))
        {
            WarningInFunction
                << "Source and target patch bounding boxes are not similar"
                << nl
                << "    source box span     : " << bbSrc.span() << nl
                << "    target box span     : " << bbTgt.span() << nl
                << "    source box          : " << bbSrc << nl
                << "    target box          : " << bbTgt << nl
                << "    inflated target box : " << bbTgtInf << endl;
        }
    }
}


template<class SourcePatch, class TargetPatch>
bool Foam::AMIMethod<SourcePatch, TargetPatch>::initialise
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    label& srcFacei,
    label& tgtFacei
)
{
    checkPatches();

    // set initial sizes for weights and addressing - must be done even if
    // returns false below
    srcAddress.setSize(srcPatch_.size());
    srcWeights.setSize(srcPatch_.size());
    tgtAddress.setSize(tgtPatch_.size());
    tgtWeights.setSize(tgtPatch_.size());

    // check that patch sizes are valid
    if (!srcPatch_.size())
    {
        return false;
    }
    else if (!tgtPatch_.size())
    {
        WarningInFunction
            << srcPatch_.size() << " source faces but no target faces" << endl;

        return false;
    }

    // reset the octree
    resetTree();

    // find initial face match using brute force/octree search
    if ((srcFacei == -1) || (tgtFacei == -1))
    {
        srcFacei = 0;
        tgtFacei = 0;
        bool foundFace = false;
        forAll(srcPatch_, facei)
        {
            tgtFacei = findTargetFace(facei);
            if (tgtFacei >= 0)
            {
                srcFacei = facei;
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            if (requireMatch_)
            {
                FatalErrorInFunction
                    << "Unable to find initial target face"
                    << abort(FatalError);
            }

            return false;
        }
    }

    if (debug)
    {
        Pout<< "AMI: initial target face = " << tgtFacei << endl;
    }

    return true;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::writeIntersectionOBJ
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

    Pout<< "Face intersection area (" << count <<  "):" << nl
        << "    f1 face = " << f1 << nl
        << "    f1 pts  = " << f1pts << nl
        << "    f2 face = " << f2 << nl
        << "    f2 pts  = " << f2pts << nl
        << "    area    = " << area
        << endl;

    OFstream os("areas" + name(count) + ".obj");

    for (const point& pt : f1pts)
    {
        meshTools::writeOBJ(os, pt);
    }
    os<< "l";
    forAll(f1pts, i)
    {
        os<< " " << i + 1;
    }
    os<< " 1" << endl;

    for (const point& pt : f2pts)
    {
        meshTools::writeOBJ(os, pt);
    }
    os<< "l";
    const label n = f1pts.size();
    forAll(f2pts, i)
    {
        os<< " " << n + i + 1;
    }
    os<< " " << n + 1 << endl;

    ++count;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::resetTree()
{
    // Clear the old octree
    treePtr_.clear();

    treeBoundBox bb(tgtPatch_.points(), tgtPatch_.meshPoints());
    bb.inflate(0.01);

    if (!treePtr_.valid())
    {
        treePtr_.reset
        (
            new indexedOctree<treeType>
            (
                treeType
                (
                    false,
                    tgtPatch_,
                    indexedOctree<treeType>::perturbTol()
                ),
                bb,                         // overall search domain
                8,                          // maxLevel
                10,                         // leaf size
                3.0                         // duplicity
            )
        );
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIMethod<SourcePatch, TargetPatch>::findTargetFace
(
    const label srcFacei
) const
{
    label targetFacei = -1;

    const pointField& srcPts = srcPatch_.points();
    const face& srcFace = srcPatch_[srcFacei];
    const point srcPt = srcFace.centre(srcPts);
    const scalar srcFaceArea = srcMagSf_[srcFacei];

    pointIndexHit sample = treePtr_->findNearest(srcPt, 10.0*srcFaceArea);

    if (sample.hit())
    {
        targetFacei = sample.index();

        if (debug)
        {
            Pout<< "Source point = " << srcPt << ", Sample point = "
                << sample.hitPoint() << ", Sample index = " << sample.index()
                << endl;
        }
    }

    return targetFacei;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::appendNbrFaces
(
    const label facei,
    const TargetPatch& patch,
    const DynamicList<label>& visitedFaces,
    DynamicList<label>& faceIDs
) const
{
    static const scalar thetaCos = Foam::cos(degToRad(89.0));

    const labelList& nbrFaces = patch.faceFaces()[facei];

    // filter out faces already visited from face neighbours
    for (const label nbrFacei : nbrFaces)
    {
        bool valid = true;
        for (const label visitedFacei : visitedFaces)
        {
            if (nbrFacei == visitedFacei)
            {
                valid = false;
                break;
            }
        }

        if (valid)
        {
            for (const label testFacei : faceIDs)
            {
                if (nbrFacei == testFacei)
                {
                    valid = false;
                    break;
                }
            }
        }

        // prevent addition of face if it is not on the same plane-ish
        if (valid)
        {
            const vector& n1 = patch.faceNormals()[facei];
            const vector& n2 = patch.faceNormals()[nbrFacei];

            const scalar cosI = n1 & n2;

            if (cosI > thetaCos)
            {
                faceIDs.append(nbrFacei);
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class PatchType>
void Foam::AMIMethod<SourcePatch, TargetPatch>::triangulatePatch
(
    const PatchType& patch,
    List<DynamicList<face>>& tris,
    List<scalar>& magSf
) const
{
    const pointField& points = patch.points();
    tris.setSize(patch.size());
    magSf.setSize(patch.size());

    // Using methods that index into existing points
    forAll(patch, facei)
    {
        switch (triMode_)
        {
            case faceAreaIntersect::tmFan:
            {
                faceAreaIntersect::triangleFan(patch[facei], tris[facei]);
                break;
            }
            case faceAreaIntersect::tmMesh:
            {
                patch[facei].triangles(points, tris[facei]);
                break;
            }
        }

        const DynamicList<face>& triFaces = tris[facei];
        magSf[facei] = 0;
        for (const face& f : triFaces)
        {
            magSf[facei] +=
                triPointRef
                (
                    points[f[0]],
                    points[f[1]],
                    points[f[2]]
                ).mag();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIMethod<SourcePatch, TargetPatch>::AMIMethod
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget,
    const bool requireMatch
)
:
    srcPatch_(srcPatch),
    tgtPatch_(tgtPatch),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    srcMagSf_(srcPatch.size(), 1.0),
    tgtMagSf_(tgtPatch.size(), 1.0),
    srcNonOverlap_(),
    triMode_(triMode)
{
    // Note: setting srcMagSf and tgtMagSf to 1 by default for 1-to-1 methods
    // - others will need to overwrite as necessary
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
bool Foam::AMIMethod<SourcePatch, TargetPatch>::conformal() const
{
    return true;
}


// ************************************************************************* //
