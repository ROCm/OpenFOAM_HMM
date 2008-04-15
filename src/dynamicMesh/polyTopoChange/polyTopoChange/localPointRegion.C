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

\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "localPointRegion.H"
#include "polyMesh.H"
#include "indirectPrimitivePatch.H"
#include "mapPolyMesh.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(localPointRegion, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Set face to global region. Do its connected points.
void Foam::localPointRegion::setFaceRegion
(
    const polyMesh& mesh,
    const Map<label>& localPoints,
    const label faceI,

    labelListList& pointFaceRegion // in pointFaces addressing.
)
{
    label localFaceI = localFaces_[faceI];

    if (faceRegion_[localFaceI] == -1)
    {
        faceRegion_[localFaceI] = nRegions_;

        // And visit all points on face
        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            label pointI = f[fp];

            Map<label>::const_iterator iter = localPoints.find(pointI);

            if (iter != localPoints.end())
            {
                setPointRegion
                (
                    mesh,
                    localPoints,
                    pointI,     // global point index
                    iter(),     // local point index
                    faceI,

                    pointFaceRegion
                );
            }
        }
    }
    else if (faceRegion_[localFaceI] != nRegions_)
    {
        FatalErrorIn
        (
            "localPointRegion::setFaceRegion"
            "(const polyMesh&, const Map<label>&, const label, const label)"
        )   << "face:" << faceI << " fc:" << mesh.faceCentres()[faceI]
            << " has global region " << faceRegion_[localFaceI]
            << " but is also connected to region " << nRegions_
            << abort(FatalError);
    }
}


// Set faces using point to global region 
void Foam::localPointRegion::setPointRegion
(
    const polyMesh& mesh,
    const Map<label>& localPoints,

    const label pointI,
    const label localPointI,
    const label faceI,

    labelListList& pointFaceRegion
)
{
    const labelList& pFaces = mesh.pointFaces()[pointI];

    // faceI is set to global region. Find the corresponding local region.

    label facePos = findIndex(pFaces, faceI);

    if (facePos != -1)
    {
        // Face uses point. Check what local region (from pointFaceRegion)
        // corresponds to the global region and renumber all local regions.

        labelList& pRegions = pointFaceRegion[localPointI];

        if (pRegions.size() > 0 && pRegions[facePos] != -1)
        {
            label localRegionI = pRegions[facePos];

            // Renumber all pRegions with same local region
            forAll(pFaces, i)
            {
                if (pRegions[i] == localRegionI)
                {
                    // Mark point region (prevents extra layer of recursion)
                    pRegions[i] = -1;

                    // Set the corresponding face region
                    setFaceRegion
                    (
                        mesh,
                        localPoints,
                        pFaces[i],          // faceI

                        pointFaceRegion
                    );
                }
            }
        }
    }
}


// Gets local face regions (local to each point) and merges them into global
// point region (always possible since faces cannot cross a baffle.
void Foam::localPointRegion::mergeLocalRegions
(
    const polyMesh& mesh,
    const labelList& pointsToBeDuplicated,
    labelListList& pointFaceRegion
)
{
    // Build inverse addressing
    Map<label> localPoints(2*pointsToBeDuplicated.size());

    forAll(pointsToBeDuplicated, i)
    {
        localPoints.insert(pointsToBeDuplicated[i], i);
    }


    nRegions_ = 0;

    forAll(meshFaces_, localFaceI)
    {
        label meshFaceI = meshFaces_[localFaceI];

        if (faceRegion_[localFaceI] == -1)
        {
            setFaceRegion
            (
                mesh,
                localPoints,
                meshFaceI,

                pointFaceRegion // in pointFaces addressing.
            );

            nRegions_++;
        }
    }
}


// Change whole region to different region number
void Foam::localPointRegion::changeRegion
(
    const polyMesh& mesh,
    const label faceI,
    const label oldRegionI,
    const label newRegionI,
    labelHashSet& changedFaces
)
{
    Map<label>::const_iterator iter = localFaces_.find(faceI);

    if (iter != localFaces_.end())
    {
        if (faceRegion_[iter()] == oldRegionI)
        {
            faceRegion_[iter()] = newRegionI;
            changedFaces.insert(iter());

            // And visit all points on face
            const face& f = mesh.faces()[faceI];

            forAll(f, fp)
            {
                label pointI = f[fp];

                const labelList& pFaces = mesh.pointFaces()[pointI];

                forAll(pFaces, i)
                {
                    changeRegion
                    (
                        mesh,
                        pFaces[i],
                        oldRegionI,
                        newRegionI,
                        changedFaces
                    );
                }
            }
        }
    }
}


// Walk cell-face-cell whilst staying on a point. Used to determine
// whether point is a non-manifold point.
Foam::label Foam::localPointRegion::walkCellFaceCell
(
    const primitiveMesh& mesh,
    const label startCellI,
    const label startPointI,
    const label regionI,
    labelList& regionPerFace    // in order of pFaces.
)
{
    const cell& cFaces = mesh.cells()[startCellI];
    const labelList& pFaces = mesh.pointFaces()[startPointI];

    label nChanged = 0;

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];

        label index = findIndex(pFaces, faceI);

        if (index != -1)
        {
            if (regionPerFace[index] == -1)
            {
                regionPerFace[index] = regionI;
                nChanged++;

                if (mesh.isInternalFace(faceI))
                {
                    label otherCellI = mesh.faceOwner()[faceI];

                    if (otherCellI == startCellI)
                    {
                        otherCellI = mesh.faceNeighbour()[faceI];
                    }

                    nChanged += walkCellFaceCell
                    (
                        mesh,
                        otherCellI,
                        startPointI,
                        regionI,
                        regionPerFace
                    );
                }
            }
        }
    }
    return nChanged;
}


// Are two lists identical either in forward or in reverse order.
bool Foam::localPointRegion::isDuplicate
(
    const face& f0,
    const face& f1,
    const bool forward
)
{
    label fp1 = findIndex(f1, f0[0]);

    if (fp1 == -1)
    {
        return false;
    }

    forAll(f0, fp0)
    {
        if (f0[fp0] != f1[fp1])
        {
            return false;
        }

        if (forward)
        {
            fp1 = f1.fcIndex(fp1);
        }
        else
        {
            fp1 = f1.rcIndex(fp1);
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::localPointRegion::localPointRegion
(
    const polyMesh& mesh,
    const labelList& pointsToBeDuplicated
)
:
    nRegions_(0),
    localFaces_(2*pointsToBeDuplicated.size()),
    meshFaces_(0),
    faceRegion_(0),
    pointRegions_(pointsToBeDuplicated.size())
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    //if (debug)
    {
        //// Check that pointsToBeDuplicated are not on coupled boundaries.
        //labelHashSet dupPoints(pointsToBeDuplicated);
        //
        //forAll(patches, patchI)
        //{
        //    const polyPatch& pp = patches[patchI];
        //
        //    if (pp.coupled())
        //    {
        //        const labelList& meshPoints = pp.meshPoints();
        //
        //        forAll(meshPoints, i)
        //        {
        //            label pointI = meshPoints[i];
        //
        //            if (dupPoints.found(pointI))
        //            {
        //                FatalErrorIn
        //                (
        //                    "localPointRegion::localPointRegion"
        //                    "(const polyMesh&, const labelList&)"
        //                )   << "Point-to-be-duplicated " << pointI
        //                    << " coord:" << mesh.points()[pointI]
        //                    << " is on coupled patch " << pp.name() << '.' << nl
        //                    << "This is not allowed." << abort(FatalError);
        //            }
        //        }
        //    }
        //}

        // Check that pointsToBeDuplicated are non-manifold
        forAll(pointsToBeDuplicated, i)
        {
            label pointI = pointsToBeDuplicated[i];

            if (isSingleCellRegion(mesh, pointI))
            {
                FatalErrorIn
                (
                    "localPointRegion::localPointRegion"
                    "(const polyMesh&, const labelList&)"
                )   << "Point-to-be-duplicated " << pointI
                    << " coord:" << mesh.points()[pointI]
                    << " is on a single cell region."
                    << abort(FatalError);
            }
        }
    }


    // Precalculate all faces connected to pointsToBeDuplicated and assign
    // local face labels.

    label localFaceI = 0;

    forAll(pointsToBeDuplicated, i)
    {
        label pointI = pointsToBeDuplicated[i];

        const labelList& pFaces = mesh.pointFaces()[pointI];

        forAll(pFaces, j)
        {
            label faceI = pFaces[j];

            if (localFaces_.insert(faceI, localFaceI))
            {
                localFaceI++;
            }
        }
    }


    // Calculate reverse addressing (so from local face to mesh face)
    meshFaces_.setSize(localFaceI);

    forAllConstIter(Map<label>, localFaces_, iter)
    {
        meshFaces_[iter()] = iter.key();
    }


    // Calculate local point wise regions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Per points-to-be-duplicated
    // either size 0 or give per connected face the local region
    labelListList pointFaceRegion(pointsToBeDuplicated.size());
    {
        forAll(pointsToBeDuplicated, i)
        {
            label pointI = pointsToBeDuplicated[i];

            const labelList& pFaces = mesh.pointFaces()[pointI];

            labelList& pRegions = pointFaceRegion[i];
            pRegions.setSize(pFaces.size());
            pRegions = -1;

            // Walk cell face cell on the point
            label regionI = 0;

            forAll(pRegions, j)
            {
                label faceI = pFaces[j];

                label nChanged = walkCellFaceCell
                (
                    mesh,
                    mesh.faceOwner()[faceI],
                    pointI,
                    regionI,
                    pRegions
                );

                if (nChanged > 0)
                {
                    regionI++;
                }
            }

            if (regionI <= 1)
            {
                // Single region so save some storage.
                pRegions.setSize(0);
            }
        }
    }


    // Set the face regions from connected point regions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    faceRegion_.setSize(meshFaces_.size());
    faceRegion_ = -1;

    mergeLocalRegions(mesh, pointsToBeDuplicated, pointFaceRegion);

    if (debug)
    {
        Pout<< "localPointRegion :"
            << "pointsToBeDuplicated:" << pointsToBeDuplicated.size()
            << " point connected faces:" << meshFaces_.size()
            << " disconnected face regions:" << nRegions_
            << endl;

        forAll(faceRegion_, localFaceI)
        {
            label faceI = meshFaces_[localFaceI];

            Pout<< "    face:" << faceI
                <<  " fc:" << mesh.faceCentres()[faceI]
                << " region:" << faceRegion_[localFaceI]
                << endl;
        }
        Pout<< endl;
    }


    bool hasCoupledPatches = false;    
    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            hasCoupledPatches = true;
            break;
        }
    }
    reduce(hasCoupledPatches, orOp<bool>());

    if (hasCoupledPatches)
    {
        // Calculate global offsets
        labelList nLocalRegions(Pstream::nProcs());
        nLocalRegions[Pstream::myProcNo()] = nRegions_;
        Pstream::gatherList(nLocalRegions);
        Pstream::scatterList(nLocalRegions);

        label myOffset = 0;
        for (label procI = 0; procI < Pstream::myProcNo(); procI++)
        {
            myOffset += nLocalRegions[procI];
        }
        label nGlobalRegions = sum(nLocalRegions);

        if (debug)
        {
            Pout<< "nLocalRegions:" << nLocalRegions << endl
                << "My region offset:" << myOffset
                << " Global regions:" << nGlobalRegions
                << endl;
        }

        // Make all regions global
        forAll(faceRegion_, localFaceI)
        {
            if (faceRegion_[localFaceI] != -1)
            {
                faceRegion_[localFaceI] += myOffset;
            }
        }

        // Merge regions
        // ~~~~~~~~~~~~~
        // Does iteratively merge connected regions. Stops if no regions
        // merged.

        while (true)
        {
            labelList boundaryRegion
            (
                mesh.nFaces()-mesh.nInternalFaces(),
                labelMax
            );

            forAll(faceRegion_, localFaceI)
            {
                if (faceRegion_[localFaceI] != -1)
                {
                    label bFaceI = meshFaces_[localFaceI]-mesh.nInternalFaces();

                    if (bFaceI >= 0)
                    {
                        boundaryRegion[bFaceI] = faceRegion_[localFaceI];
                    }
                }
            }

            syncTools::syncBoundaryFaceList
            (
                mesh,
                boundaryRegion,
                minEqOp<label>(),
                false
            );

            // Now on coupled faces the boundaryRegion will be a valid region
            // or labelMax.


            // Count of merged regions so we know when to stop.
            label nMerged = 0;

            // To keep track of regions already changed. Makes sure we
            // don't go into loop. Scenario:
            // - at some face myRegion=2, bRegion=0
            // - all faces with myRegion=2 get changed
            // - some other face which had myRegion 2 but now has myRegion 0
            //   now sees bRegion=2 and starts changing back.

            labelHashSet changedLocalFaces(2*faceRegion_.size());

            forAll(faceRegion_, localFaceI)
            {
                if (!changedLocalFaces.found(localFaceI))
                {
                    label meshFaceI = meshFaces_[localFaceI];
                    label bFaceI = meshFaceI-mesh.nInternalFaces();

                    if (bFaceI >= 0)
                    {
                        label bRegion = boundaryRegion[bFaceI];
                        if (bRegion == labelMax)
                        {
                            bRegion = -1;
                        }
                        label myRegion = faceRegion_[localFaceI];

                        if (myRegion == -1 || bRegion  == -1)
                        {
                            // Both sides should have no region

                            if (myRegion != bRegion)
                            {
                                FatalErrorIn
                                (
                                    "localPointRegion::localPointRegion"
                                    "(const polyMesh&, const labelList&)"
                                )   << "Two coupled faces of which one has"
                                    << " an illegal region." << endl
                                    << "Face:" << meshFaceI
                                    << " has global region:" << myRegion
                                    << " and is coupled to a face"
                                    << " with global region " << bRegion
                                    << abort(FatalError);
                            }
                        }
                        else if (myRegion != bRegion)
                        {
                            // Both sides have valid and different region

                            if (debug)
                            {
                                Pout<< "Merging my region:"
                                    << myRegion << " and " << bRegion
                                    << " into " << bRegion << endl;
                            }
                            changeRegion
                            (
                                mesh,
                                meshFaceI,
                                myRegion,
                                bRegion,
                                changedLocalFaces
                            );

                            nMerged++;
                        }
                    }
                }
            }

            if (returnReduce(nMerged, sumOp<label>()) == 0)
            {
                break;
            }
        }


        // Compact region numbering
        // ~~~~~~~~~~~~~~~~~~~~~~~~


        // Determine (global) number of faces per region

        labelList nFaces(nGlobalRegions, 0);
        forAll(faceRegion_, localFaceI)
        {
            label region = faceRegion_[localFaceI];

            if (region != -1)
            {
                nFaces[region]++;
            }
        }
        Pstream::listCombineGather(nFaces, plusEqOp<label>());
        Pstream::listCombineScatter(nFaces);

        // Make compaction list.
        labelList oldToNew(nGlobalRegions, -1);
        nRegions_ = 0;
        forAll(nFaces, regionI)
        {
            if (nFaces[regionI] > 0)
            {
                oldToNew[regionI] = nRegions_++;
            }
        }
        inplaceRenumber(oldToNew, faceRegion_);

        if (debug)
        {
            Pout<< "Compacted from local regions:"
                << nLocalRegions[Pstream::myProcNo()]
                << " global regions:" << nGlobalRegions
                << " down to:" << nRegions_
                << endl;
        }
    }


    if (debug)
    {
        Pout<< "localPointRegion :"
            << "pointsToBeDuplicated:" << pointsToBeDuplicated.size()
            << " point connected faces:" << meshFaces_.size()
            << " disconnected face regions:" << nRegions_
            << endl;

        forAll(faceRegion_, localFaceI)
        {
            label faceI = meshFaces_[localFaceI];

            Pout<< "    face:" << faceI
                <<  " fc:" << mesh.faceCentres()[faceI]
                << " region:" << faceRegion_[localFaceI]
                << endl;
        }
        Pout<< endl;
    }



    // Collect regions per point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Done in one pass since usually only one or two regions.

    forAll(pointsToBeDuplicated, i)
    {
        label pointI = pointsToBeDuplicated[i];

        labelList& pRegions = pointRegions_[i];

        const labelList& pFaces = mesh.pointFaces()[pointI];

        forAll(pFaces, pFaceI)
        {
            label faceI = pFaces[pFaceI];

            label region = faceRegion_[localFaces_[faceI]];

            if (findIndex(pRegions, region) == -1)
            {
                label sz = pRegions.size();
                pRegions.setSize(sz+1);
                pRegions[sz] = region;
            }
        }

        if (debug)
        {
            Pout<< "    point:" << pointI << " coord:" << mesh.points()[pointI]
                << " regions:" << pRegions << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::localPointRegion::isSingleCellRegion
(
    const primitiveMesh& mesh,
    const label pointI
)
{
    const labelList& pFaces = mesh.pointFaces()[pointI];

    // Mark off visited faces.
    labelList regionPerFace(pFaces.size(), -1);

    // Check if all faces using pointI can be reached through face-cell-face
    // walk
    label nVisited = walkCellFaceCell
    (
        mesh,
        mesh.faceOwner()[pFaces[0]],
        pointI,
        1,              // arbitrary value to mark visited faces with
        regionPerFace   // status of face
    );

    return (nVisited == pFaces.size());
}


// Return a list (in allPatch indices) with either -1 or the face label
// of the face that uses the same vertices.
Foam::labelList Foam::localPointRegion::findDuplicateFaces
(
    const primitiveMesh& mesh,
    const labelList& boundaryFaces
)
{
    // Addressing engine for all boundary faces.
    indirectPrimitivePatch allPatch
    (
        IndirectList<face>(mesh.faces(), boundaryFaces),
        mesh.points()
    );

    labelList duplicateFace(allPatch.size(), -1);
    label nDuplicateFaces = 0;

    // Find all duplicate faces.
    forAll(allPatch, bFaceI)
    {
        const face& f = allPatch.localFaces()[bFaceI];

        // Get faces connected to f[0].
        // Check whether share all points with f.
        const labelList& pFaces = allPatch.pointFaces()[f[0]];

        forAll(pFaces, i)
        {
            label otherFaceI = pFaces[i];

            if (otherFaceI > bFaceI)
            {
                const face& otherF = allPatch.localFaces()[otherFaceI];

                if (isDuplicate(f, otherF, true))
                {
                    FatalErrorIn
                    (
                        "findDuplicateFaces(const primitiveMesh&"
                        ", const labelList&)"
                    )   << "Face:" << bFaceI + mesh.nInternalFaces()
                        << " has local points:" << f
                        << " which are in same order as face:"
                        << otherFaceI + mesh.nInternalFaces()
                        << " with local points:" << otherF
                        << abort(FatalError);
                }
                else if (isDuplicate(f, otherF, false))
                {
                    label meshFace0 = bFaceI + mesh.nInternalFaces();
                    label meshFace1 = otherFaceI + mesh.nInternalFaces();

                    if
                    (
                        duplicateFace[bFaceI] != -1
                     || duplicateFace[otherFaceI] != -1
                    )
                    {
                        FatalErrorIn
                        (
                            "findDuplicateFaces(const primitiveMesh&"
                            ", const labelList&)"
                        )   << "One of two duplicate faces already marked"
                            << " as duplicate." << nl
                            << "This means that three or more faces share"
                            << " the same points and this is illegal." << nl
                            << "Face:" << meshFace0
                            << " with local points:" << f
                            << " which are in same order as face:"
                            << meshFace1
                            << " with local points:" << otherF
                            << abort(FatalError);
                    }

                    duplicateFace[bFaceI] = otherFaceI;
                    duplicateFace[otherFaceI] = bFaceI;
                    nDuplicateFaces++;
                }
            }
        }
    }

    return duplicateFace;
}


void Foam::localPointRegion::updateMesh(const mapPolyMesh& map)
{
    Map<label> newLocalFaces(localFaces_.size());

    forAllConstIter(Map<label>, localFaces_, iter)
    {
        label newFaceI = map.reverseFaceMap()[iter.key()];

        if (newFaceI >= 0)
        {
            newLocalFaces.insert(newFaceI, iter());
            meshFaces_[iter()] = newFaceI;
        }
    }

    localFaces_.transfer(newLocalFaces);
}


// ************************************************************************* //
