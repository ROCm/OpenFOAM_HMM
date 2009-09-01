/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "referredCellList.H"
#include "interactionLists.H"
#include "polyBoundaryMeshEntries.H"
#include "PstreamCombineReduceOps.H"
#include "Time.H"
#include "globalMeshData.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "sendingReferralList.H"
#include "receivingReferralList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::referredCellList<ParticleType>::buildReferredCellList
(
    bool pointPointListBuild
)
{
    Info << nl << "Building list of referred interaction neighbours" << endl;

    const polyMesh& mesh(il_.mesh());

    DynamicList<referredCell<ParticleType> > referredInteractionList;

    // realCellsWithinRangeOfAnyReferringPatch
    DynamicList<label> rCellsWRRP;

    // realFacesWithinRangeOfAnyReferringPatch
    DynamicList<label> rFacesWRRP;

    // realEdgesWithinRangeOfAnyReferringPatch
    DynamicList<label> rEdgesWRRP;

    // realPointsWithinRangeOfAnyReferringPatch
    DynamicList<label> rPointsWRRP;

    labelListList processorPatchSegmentMapping
    (
        mesh.globalData().processorPatches().size()
    );

    List<vectorList> allNeighbourFaceCentres
    (
        mesh.globalData().processorPatches().size()
    );

    List<vectorList> allNeighbourFaceAreas
    (
        mesh.globalData().processorPatches().size()
    );

    label nUndecomposedPatches = 0;

    if (Pstream::parRun())
    {
        dictionary patchDictionary;

        DynamicList<word> patchNames;

        Time undecomposedTime
        (
            Time::controlDictName,
            mesh.time().rootPath(),
            mesh.time().caseName().path()
        );

        IOobject undecomposedBoundaryHeader
        (
            "boundary",
            undecomposedTime.constant(),
            polyMesh::meshSubDir,
            undecomposedTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (undecomposedBoundaryHeader.headerOk())
        {
            polyBoundaryMeshEntries undecomposedPatchEntries
            (
                undecomposedBoundaryHeader
            );

            forAll(undecomposedPatchEntries, patchi)
            {
                patchNames.append
                (
                    undecomposedPatchEntries[patchi].keyword()
                );

                patchDictionary.add
                (
                    undecomposedPatchEntries[patchi]
                );
            }
        }
        else
        {
            FatalErrorIn ("referredCellList.C")
                << nl << "unable to read undecomposed boundary file from "
                << "constant/polyMesh" << nl
                << abort(FatalError);
        }

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                mesh.time().constant(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        labelList procPatches(mesh.globalData().processorPatches());

        nUndecomposedPatches = patchNames.size();

        // processorPatchSegmentMapping works by mapping the original patch and
        // half that a face on a processor patch was on before decomposition.
        // This creates a patch segment for each half of each original (cyclic)
        // patch which can be assessed separately.  There are n =
        // patchNames.size() original patches, k = 0 to n-1.  The mapping is:
        // value = 0: originally an internal face value = k, was originally on
        // the on the patch k-1, 1st half value = -k, was originally on the on
        // the patch k-1, 2nd half

        forAll(procPatches,pP)
        {
            const processorPolyPatch& patch = refCast<const processorPolyPatch>
            (
                mesh.boundaryMesh()[procPatches[pP]]
            );

            labelList& procPatchSegMap = processorPatchSegmentMapping[pP];

            procPatchSegMap.setSize(patch.size());

            forAll (patch, pI)
            {
                label decomposedMeshFace = patch.start() + pI;

                label faceProcAdd = faceProcAddressing[decomposedMeshFace];

                label globalFace = abs(faceProcAdd)-1;

                label minStart = -1;

                forAll(patchNames, pN)
                {
                    if (patchDictionary.found(patchNames[pN]))
                    {
                        const dictionary& patchDict =
                        patchDictionary.subDict(patchNames[pN]);

                        label startFace
                        (
                            readLabel
                            (
                                patchDict.lookup("startFace")
                            )
                        );

                        label nFaces(readLabel(patchDict.lookup("nFaces")));

                        if (minStart < 0 || startFace < minStart)
                        {
                            minStart = startFace;
                        }

                        if
                        (
                            globalFace >= startFace
                         && globalFace < startFace + nFaces/2
                        )
                        {
                            procPatchSegMap[pI] = pN + 1;
                        }
                        else if
                        (
                            globalFace >= startFace + nFaces/2
                         && globalFace < startFace + nFaces
                        )
                        {
                            procPatchSegMap[pI] = -(pN + 1);
                        }
                    }
                }

                if (globalFace < minStart)
                {
                    procPatchSegMap[pI] = 0;
                }
            }
        }

        forAll(procPatches,pP)
        {
            const processorPolyPatch& patch = refCast<const processorPolyPatch>
            (
                mesh.boundaryMesh()[procPatches[pP]]
            );

            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    patch.neighbProcNo()
                );

                toNeighbProc << patch.faceCentres() << patch.faceAreas();
            }
        }

        forAll(procPatches,pP)
        {
            const processorPolyPatch& patch = refCast<const processorPolyPatch>
            (
                mesh.boundaryMesh()[procPatches[pP]]
            );

            vectorList& neighbFaceCentres = allNeighbourFaceCentres[pP];

            neighbFaceCentres.setSize(patch.size());

            vectorList& neighbFaceAreas = allNeighbourFaceAreas[pP];

            neighbFaceAreas.setSize(patch.size());

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    patch.neighbProcNo()
                );

                fromNeighbProc >> neighbFaceCentres >> neighbFaceAreas;
            }
        }

        // *************************************************************
        // Tests that all processor patch segments from different
        // original patches prior to decomposition produce the same
        // transform. Check before 1st iteration.
        // *************************************************************

        forAll(procPatches,pP)
        {
            const processorPolyPatch& patch = refCast<const processorPolyPatch>
            (
                mesh.boundaryMesh()[procPatches[pP]]
            );

            const vectorList& neighbFaceCentres = allNeighbourFaceCentres[pP];

            const vectorList& neighbFaceAreas = allNeighbourFaceAreas[pP];

            label nUP;

            for
            (
                nUP = -nUndecomposedPatches;
                nUP <= nUndecomposedPatches;
                nUP++
            )
            {
                DynamicList<vector> refOff;

                DynamicList<tensor> refTrans;

                forAll (patch, faceI)
                {
                    if (processorPatchSegmentMapping[pP][faceI] == nUP)
                    {
                        referredCell<ParticleType> testRefCell
                        (
                            mesh,
                            -1,
                            -1,
                            patch.faceCentres()[faceI],
                            neighbFaceCentres[faceI],
                            patch.faceNormals()[faceI],
                            neighbFaceAreas[faceI]
                           /(mag(neighbFaceAreas[faceI]) + VSMALL)
                        );

                        refOff.append(testRefCell.offset());

                        refTrans.append(testRefCell.rotation());
                    }
                }

                refOff.shrink();

                refTrans.shrink();

                if (refOff.size())
                {
                    if
                    (
                        sum(mag(refOff-refOff[0]))/refOff.size()
                            > interactionLists<ParticleType>::transTol
                     || sum(mag(refTrans-refTrans[0]))/refTrans.size()
                            > interactionLists<ParticleType>::transTol
                    )
                    {
                        FatalErrorIn ("referredCellList.C")
                            << nl << "Face pairs on patch "
                            << patch.name()
                            << ", segment " << patchNames[nUP]
                            << " do not give the same referring "
                            << " transformations to within tolerance of "
                            << interactionLists<ParticleType>::transTol << nl
                            << " Referring offsets:" << refOff << nl
                            << " Average sum of mag difference: "
                            << sum(mag(refOff-refOff[0]))/refOff.size() << nl
                            << " Referring transforms:" << refTrans << nl
                            << " Average sum of mag difference: "
                            << sum(mag(refTrans-refTrans[0]))/refTrans.size()
                            << nl << abort(FatalError);
                    }
                }
            }
        }
    }

    label cellsReferredThisIteration = 1;

    label iterationNo = 0;

    while (cellsReferredThisIteration > 0)
    {
        label refIntListStartSize = referredInteractionList.size();

        forAll (mesh.boundaryMesh(), patchI)
        {
            // Treat local cyclics on each processor before processor
            // boundaries.  Separate treatment allows the serial version to run
            // transparently.

            if (mesh.boundaryMesh()[patchI].type() == "cyclic")
            {
                const cyclicPolyPatch& patch = refCast<const cyclicPolyPatch>
                (
                    mesh.boundaryMesh()[patchI]
                );

                if (patch.size())
                {
                    if (iterationNo == 0)
                    {
                        // Tests that all combinations of face pairs produce the
                        // same transform.  Only check on 1st iteration.

                        label faceL;
                        // A face in the 1st half of the patch

                        label faceM;
                        // Face corresponding to faceL in the 2nd half of the
                        // patch. Assumes correct face ordering.

                        vectorList refOff(patch.size()/2);

                        List<tensor> refTrans(patch.size()/2);

                        for
                        (
                            faceL = 0, faceM = patch.size()/2;
                            faceL < patch.size()/2;
                            faceL++, faceM++
                        )
                        {
                            referredCell<ParticleType> testRefCell
                            (
                                mesh,
                                -1,
                                -1,
                                patch.faceCentres()[faceL],
                                patch.faceCentres()[faceM],
                                patch.faceNormals()[faceL],
                                patch.faceNormals()[faceM]
                            );

                            refOff[faceL] = testRefCell.offset();

                            refTrans[faceL] = testRefCell.rotation();
                        }

                        if
                        (
                            sum(mag(refOff - refOff[0]))/(patch.size()/2)
                                > interactionLists<ParticleType>::transTol
                         || sum(mag(refTrans - refTrans[0]))/(patch.size()/2)
                                > interactionLists<ParticleType>::transTol
                        )
                        {
                            FatalErrorIn ("referredCellList.C")
                                << nl << "Face pairs on patch "
                                << patch.name()
                                << " do not give the same referring "
                                << " transformations to within tolerance of "
                                << interactionLists<ParticleType>::transTol
                                << nl << " Referring offsets:" << refOff << nl
                                << " Average sum of mag difference: "
                                << sum(mag(refOff - refOff[0]))/refOff.size()
                                << nl
                                << " Referring transforms:" << refTrans << nl
                                << " Average sum of mag difference: "
                                << sum(mag(refTrans - refTrans[0]))
                                  /refTrans.size()
                                << nl << abort(FatalError);
                        }
                    }

                    // *********************************************************
                    // 1st half of face list - 1st side of boundary
                    // *********************************************************

                    label faceI;

                    DynamicList<label> meshFacesOnThisSegment;

                    for (faceI = 0; faceI < patch.size()/2; faceI++)
                    {
                        // unable to use the normal accessors for the polyPatch
                        // because points on separate halves need used
                        // separately

                        meshFacesOnThisSegment.append(faceI + patch.start());
                    }

                    meshFacesOnThisSegment.shrink();

                    DynamicList<label> meshEdgesOnThisSegment;

                    DynamicList<label> meshPointsOnThisSegment;

                    forAll(meshFacesOnThisSegment, mFOTS)
                    {
                        const label segFace = meshFacesOnThisSegment[mFOTS];

                        const labelList& faceEdges = mesh.faceEdges()[segFace];

                        forAll (faceEdges, fE)
                        {
                            const label faceEdge(faceEdges[fE]);

                            if
                            (
                                findIndex
                                (
                                    meshEdgesOnThisSegment,
                                    faceEdge
                                ) == -1
                            )
                            {
                                meshEdgesOnThisSegment.append(faceEdge);
                            }
                        }

                        const face& facePoints(mesh.faces()[segFace]);

                        forAll (facePoints, fP)
                        {
                            const label facePoint(facePoints[fP]);

                            if
                            (
                                findIndex
                                (
                                    meshPointsOnThisSegment,
                                    facePoint
                                )
                             ==
                                -1
                            )
                            {
                                meshPointsOnThisSegment.append(facePoint);
                            }
                        }
                    }

                    meshEdgesOnThisSegment.shrink();

                    meshPointsOnThisSegment.shrink();

                    if (iterationNo == 0)
                    {
                        // Assessing real cells in range is only required on
                        // the 1st iteration because they do not change from
                        // iteration to iteration.

                        labelList realCellsFoundInRange
                        (
                            il_.realCellsInRangeOfSegment
                            (
                                meshFacesOnThisSegment,
                                meshEdgesOnThisSegment,
                                meshPointsOnThisSegment
                            )
                        );

                        forAll(realCellsFoundInRange,cFIR)
                        {
                            const label realCell = realCellsFoundInRange[cFIR];

                            referredCell<ParticleType> cellToRefer
                            (
                                mesh,
                                Pstream::myProcNo(),
                                realCell,
                                patch.faceCentres()[0],
                                patch.faceCentres()[patch.size()/2],
                                patch.faceNormals()[0],
                                patch.faceNormals()[patch.size()/2]
                            );

                            // Test all existing referred and real cells to
                            // check duplicates are not being made or cells
                            // aren't being referred back onto themselves

                            bool addCellToRefer = true;

                            // Check if cellToRefer is an existing referred cell

                            forAll(referredInteractionList, rIL)
                            {
                                if
                                (
                                    cellToRefer.duplicate
                                    (
                                        referredInteractionList[rIL]
                                    )
                                )
                                {
                                    addCellToRefer = false;

                                    break;
                                }
                            }

                            // Check for cellToRefer being referred back
                            // ontop of a real cell

                            if
                            (
                                cellToRefer.duplicate
                                (
                                    Pstream::myProcNo(),
                                    mesh.nCells()
                                )
                            )
                            {
                                addCellToRefer = false;
                            }

                            if (addCellToRefer)
                            {
                                referredInteractionList.append(cellToRefer);
                            }

                            // add real cells found in range of cyclic patch
                            // to whole mesh list

                            if (findIndex (rCellsWRRP, realCell) == -1)
                            {
                                rCellsWRRP.append(realCell);
                            }
                        }
                    }

                    referredInteractionList.shrink();

                    labelList referredCellsFoundInRange
                    (
                        il_.referredCellsInRangeOfSegment
                        (
                            referredInteractionList,
                            meshFacesOnThisSegment,
                            meshEdgesOnThisSegment,
                            meshPointsOnThisSegment
                        )
                    );

                    forAll(referredCellsFoundInRange,cFIR)
                    {
                        referredCell<ParticleType>& existingRefCell =
                            referredInteractionList
                            [
                                referredCellsFoundInRange[cFIR]
                            ];

                        referredCell<ParticleType> cellToReRefer =
                            existingRefCell.reRefer
                            (
                                patch.faceCentres()[0],
                                patch.faceCentres()[patch.size()/2],
                                patch.faceNormals()[0],
                                patch.faceNormals()[patch.size()/2]
                            );

                        // Test all existing referred and real cells to check
                        // duplicates are not being made or cells aren't being
                        // referred back onto themselves

                        bool addCellToReRefer = true;

                        // Check if cellToRefer is an existing referred cell

                        forAll(referredInteractionList, rIL)
                        {
                            if
                            (
                                cellToReRefer.duplicate
                                (
                                    referredInteractionList[rIL]
                                )
                            )
                            {
                                addCellToReRefer = false;

                                break;
                            }
                        }

                        // Check for cellToRefer being referred back
                        // ontop of a real cell

                        if
                        (
                            cellToReRefer.duplicate
                            (
                                Pstream::myProcNo(),
                                mesh.nCells()
                            )
                        )
                        {
                            addCellToReRefer = false;
                        }

                        if (addCellToReRefer)
                        {
                            referredInteractionList.append(cellToReRefer);
                        }
                    }

                    // *********************************************************
                    // 2nd half of face list - 2nd side of boundary
                    // *********************************************************

                    meshFacesOnThisSegment.clear();

                    for (faceI = patch.size()/2; faceI < patch.size(); faceI++)
                    {
                        // unable to use the normal accessors for the polyPatch
                        // because points on separate halves need used
                        // separately

                        meshFacesOnThisSegment.append(faceI + patch.start());
                    }

                    meshFacesOnThisSegment.shrink();

                    meshEdgesOnThisSegment.clear();

                    meshPointsOnThisSegment.clear();

                    forAll(meshFacesOnThisSegment, mFOTS)
                    {
                        const label segFace = meshFacesOnThisSegment[mFOTS];

                        const labelList& faceEdges = mesh.faceEdges()[segFace];

                        forAll (faceEdges, fE)
                        {
                            const label faceEdge(faceEdges[fE]);

                            if
                            (
                                findIndex
                                (
                                    meshEdgesOnThisSegment,
                                    faceEdge
                                )
                             ==
                                -1
                            )
                            {
                                meshEdgesOnThisSegment.append(faceEdge);
                            }
                        }

                        const face& facePoints(mesh.faces()[segFace]);

                        forAll (facePoints, fP)
                        {
                            const label facePoint(facePoints[fP]);

                            if
                            (
                                findIndex
                                (
                                    meshPointsOnThisSegment,
                                    facePoint
                                )
                             ==
                                -1
                            )
                            {
                                meshPointsOnThisSegment.append(facePoint);
                            }
                        }
                    }

                    meshEdgesOnThisSegment.shrink();

                    meshPointsOnThisSegment.shrink();

                    if (iterationNo == 0)
                    {
                        // Assessing real cells in range is only required on
                        // the 1st iteration because they do not change from
                        // iteration to iteration.

                        labelList realCellsFoundInRange
                        (
                            il_.realCellsInRangeOfSegment
                            (
                                meshFacesOnThisSegment,
                                meshEdgesOnThisSegment,
                                meshPointsOnThisSegment
                            )
                        );

                        forAll(realCellsFoundInRange,cFIR)
                        {
                            const label realCell = realCellsFoundInRange[cFIR];

                            referredCell<ParticleType> cellToRefer
                            (
                                mesh,
                                Pstream::myProcNo(),
                                realCell,
                                patch.faceCentres()[patch.size()/2],
                                patch.faceCentres()[0],
                                patch.faceNormals()[patch.size()/2],
                                patch.faceNormals()[0]
                            );

                            // Test all existing referred and real cells to
                            // check duplicates are not being made or cells
                            // aren't being referred back onto themselves

                            bool addCellToRefer = true;

                            // Check if cellToRefer is an existing referred cell

                            forAll(referredInteractionList, rIL)
                            {
                                if
                                (
                                    cellToRefer.duplicate
                                    (
                                        referredInteractionList[rIL]
                                    )
                                )
                                {
                                    addCellToRefer = false;

                                    break;
                                }
                            }

                            // Check for cellToRefer being referred back
                            // ontop of a real cell

                            if
                            (
                                cellToRefer.duplicate
                                (
                                    Pstream::myProcNo(),
                                    mesh.nCells()
                                )
                            )
                            {
                                addCellToRefer = false;
                            }

                            if (addCellToRefer)
                            {
                                referredInteractionList.append(cellToRefer);
                            }

                            // add real cells found in range of cyclic patch
                            // to whole mesh list

                            if (findIndex (rCellsWRRP, realCell) == -1)
                            {
                                rCellsWRRP.append(realCell);
                            }
                        }
                    }

                    referredInteractionList.shrink();

                    referredCellsFoundInRange =
                        il_.referredCellsInRangeOfSegment
                        (
                            referredInteractionList,
                            meshFacesOnThisSegment,
                            meshEdgesOnThisSegment,
                            meshPointsOnThisSegment
                        );

                    forAll(referredCellsFoundInRange,cFIR)
                    {
                        referredCell<ParticleType>& existingRefCell =
                            referredInteractionList
                            [
                                referredCellsFoundInRange[cFIR]
                            ];

                        referredCell<ParticleType> cellToReRefer =
                            existingRefCell.reRefer
                            (
                                patch.faceCentres()[patch.size()/2],
                                patch.faceCentres()[0],
                                patch.faceNormals()[patch.size()/2],
                                patch.faceNormals()[0]
                            );

                        // Test all existing referred and real cells to check
                        // duplicates are not being made or cells aren't being
                        // referred back onto themselves

                        bool addCellToReRefer = true;

                        // Check if cellToRefer is an existing referred cell

                        forAll(referredInteractionList, rIL)
                        {
                            if
                            (
                                cellToReRefer.duplicate
                                (
                                    referredInteractionList[rIL]
                                )
                            )
                            {
                                addCellToReRefer = false;

                                break;
                            }
                        }

                        // Check for cellToRefer being referred back
                        // ontop of a real cell

                        if
                        (
                            cellToReRefer.duplicate
                            (
                                Pstream::myProcNo(),
                                mesh.nCells()
                            )
                        )
                        {
                            addCellToReRefer = false;
                        }

                        if (addCellToReRefer)
                        {
                            referredInteractionList.append(cellToReRefer);
                        }
                    }
                }
            }
        }

        if (Pstream::parRun())
        {
            labelList procPatches(mesh.globalData().processorPatches());

            forAll(procPatches,pP)
            {
                const processorPolyPatch& patch =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[procPatches[pP]]
                    );

                DynamicList<referredCell<ParticleType> > referredCellsToTransfer;

                const vectorList& neighbFaceCentres =
                    allNeighbourFaceCentres[pP];

                const vectorList& neighbFaceAreas = allNeighbourFaceAreas[pP];

                label nUP;

                for
                (
                    nUP = -nUndecomposedPatches;
                    nUP <= nUndecomposedPatches;
                    nUP++
                )
                {
                    // faceT is used to specify one face on this patch segment
                    // that will be used to calculate the transformation values.
                    // All faces are guaranteed to produce the same transform
                    // because of the checks carried out at the start of the
                    // function.  Setting to -1 until the 1st face on this
                    // segment is found.

                    label faceT = -1;

                    DynamicList<label> meshFacesOnThisSegment;

                    forAll (patch, faceI)
                    {
                        if (processorPatchSegmentMapping[pP][faceI] == nUP)
                        {
                            if (faceT == -1)
                            {
                                faceT = faceI;
                            }

                            meshFacesOnThisSegment.append
                            (
                                faceI + patch.start()
                            );
                        }
                    }

                    meshFacesOnThisSegment.shrink();

                    DynamicList<label> meshEdgesOnThisSegment;

                    DynamicList<label> meshPointsOnThisSegment;

                    forAll(meshFacesOnThisSegment, mFOTS)
                    {
                        const label segFace = meshFacesOnThisSegment[mFOTS];

                        const labelList& faceEdges = mesh.faceEdges()[segFace];

                        forAll (faceEdges, fE)
                        {
                            const label faceEdge(faceEdges[fE]);

                            if
                            (
                                findIndex
                                (
                                    meshEdgesOnThisSegment,
                                    faceEdge
                                )
                             ==
                                -1
                            )
                            {
                                meshEdgesOnThisSegment.append(faceEdge);
                            }
                        }

                        const face& facePoints(mesh.faces()[segFace]);

                        forAll (facePoints, fP)
                        {
                            const label facePoint(facePoints[fP]);

                            if
                            (
                                findIndex
                                (
                                    meshPointsOnThisSegment,
                                    facePoint
                                )
                             ==
                                -1
                            )
                            {
                                meshPointsOnThisSegment.append(facePoint);
                            }
                        }
                    }

                    meshEdgesOnThisSegment.shrink();

                    meshPointsOnThisSegment.shrink();

                    if (meshFacesOnThisSegment.size())
                    {
                        if (faceT == -1)
                        {
                            FatalErrorIn ("referredCellList.C")
                                << nl << "faceT == -1 encountered but "
                                << meshFacesOnThisSegment.size()
                                << " faces found on patch segment."
                                << abort(FatalError);
                        }

                        if (iterationNo == 0)
                        {
                            // Assessing real cells in range is only required on
                            // the 1st iteration because they do not change from
                            // iteration to iteration.

                            labelList realCellsFoundInRange
                            (
                                il_.realCellsInRangeOfSegment
                                (
                                    meshFacesOnThisSegment,
                                    meshEdgesOnThisSegment,
                                    meshPointsOnThisSegment
                                )
                            );

                            forAll(realCellsFoundInRange,cFIR)
                            {
                                const label realCell =
                                    realCellsFoundInRange[cFIR];

                                referredCell<ParticleType> cellToRefer
                                (
                                    mesh,
                                    Pstream::myProcNo(),
                                    realCell,
                                    patch.faceCentres()[faceT],
                                    neighbFaceCentres[faceT],
                                    patch.faceNormals()[faceT],
                                    neighbFaceAreas[faceT]
                                   /(mag(neighbFaceAreas[faceT]) + VSMALL)
                                );

                                referredCellsToTransfer.append(cellToRefer);

                                // add real cells found in range of processor
                                // patch to whole mesh list

                                if (findIndex (rCellsWRRP, realCell) == -1)
                                {
                                    rCellsWRRP.append(realCell);
                                }
                            }
                        }

                        referredInteractionList.shrink();

                        labelList referredCellsFoundInRange
                        (
                            il_.referredCellsInRangeOfSegment
                            (
                                referredInteractionList,
                                meshFacesOnThisSegment,
                                meshEdgesOnThisSegment,
                                meshPointsOnThisSegment
                            )
                        );

                        forAll(referredCellsFoundInRange,cFIR)
                        {
                            referredCell<ParticleType>& existingRefCell =
                                referredInteractionList
                                [
                                    referredCellsFoundInRange[cFIR]
                                ];

                            referredCell<ParticleType> cellToReRefer =
                                existingRefCell.reRefer
                                (
                                    patch.faceCentres()[faceT],
                                    neighbFaceCentres[faceT],
                                    patch.faceNormals()[faceT],
                                    neighbFaceAreas[faceT]
                                   /(mag(neighbFaceAreas[faceT]) + VSMALL)
                                );

                            referredCellsToTransfer.append(cellToReRefer);
                        }
                    }
                }

                referredCellsToTransfer.shrink();

                // Send these cells to the neighbouring processor.

                {
                    OPstream toNeighbProc
                    (
                        Pstream::blocking,
                        patch.neighbProcNo()
                    );

                    toNeighbProc << referredCellsToTransfer;
                }
            }

            forAll(procPatches,pP)
            {
                const processorPolyPatch& patch =
                refCast<const processorPolyPatch>
                (
                    mesh.boundaryMesh()[procPatches[pP]]
                );

                // Receive the cells from neighbour

                List<referredCell<ParticleType> >
                    referredCellsFromNeighbour(patch.size());

                {
                    IPstream fromNeighbProc
                    (
                        Pstream::blocking,
                        patch.neighbProcNo()
                    );

                    fromNeighbProc >> referredCellsFromNeighbour;
                }

                // Check to see if they are duplicates, if not append
                // them to the referredInteractionList

                forAll(referredCellsFromNeighbour,rCFN)
                {
                    referredCell<ParticleType>& cellToRefer =
                        referredCellsFromNeighbour[rCFN];

                    // Test all existing referred and real cells to check
                    // duplicates are not being made or cells aren't being
                    // referred back onto themselves

                    bool addCellToRefer = true;

                    // Check if cellToRefer is an existing referred cell

                    forAll(referredInteractionList, rIL)
                    {
                        if (cellToRefer.duplicate(referredInteractionList[rIL]))
                        {
                            addCellToRefer = false;

                            break;
                        }
                    }

                    // Check for cellToRefer being referred back ontop of a real
                    // cell

                    if
                    (
                        cellToRefer.duplicate
                        (
                            Pstream::myProcNo(),
                            mesh.nCells()
                        )
                    )
                    {
                        addCellToRefer = false;
                    }

                    if (addCellToRefer)
                    {
                        referredInteractionList.append(cellToRefer);
                    }
                }
            }
        }

        if (iterationNo == 0)
        {
            // record all real cells in range of any referring patch (cyclic or
            // processor) on the first iteration when the real cells are
            // evaluated.

            rCellsWRRP.shrink();

            // construct {faces, edges, points}WithinRangeOfAnyReferringPatch

            forAll(rCellsWRRP, rCWR)
            {
                const label realCell(rCellsWRRP[rCWR]);

                const labelList& rCFaces
                (
                    mesh.cells()[realCell]
                );

                forAll(rCFaces, rCF)
                {
                    const label f(rCFaces[rCF]);

                    if (findIndex(rFacesWRRP,f) == -1)
                    {
                        rFacesWRRP.append(f);
                    }
                }

                const labelList& rCEdges
                (
                    mesh.cellEdges()[realCell]
                );

                forAll(rCEdges, rCE)
                {
                    const label e(rCEdges[rCE]);

                    if (findIndex(rEdgesWRRP,e) == -1)
                    {
                        rEdgesWRRP.append(e);
                    }
                }

                const labelList& rCPoints
                (
                    mesh.cellPoints()[realCell]
                );

                forAll(rCPoints, rCP)
                {
                    const label p(rCPoints[rCP]);

                    if (findIndex(rPointsWRRP,p) == -1)
                    {
                        rPointsWRRP.append(p);
                    }
                }
            }

            rFacesWRRP.shrink();

            rEdgesWRRP.shrink();

            rPointsWRRP.shrink();
        }

        iterationNo++;

        cellsReferredThisIteration =
            referredInteractionList.size() - refIntListStartSize;

        reduce(cellsReferredThisIteration, sumOp<label>());

        Info<< tab << "Cells added this iteration: "
            << cellsReferredThisIteration << endl;
    }

    referredInteractionList.shrink();

    (*this).setSize
    (
        referredInteractionList.size()
    );

    forAll(referredInteractionList, rIL)
    {
        (*this)[rIL] = referredInteractionList[rIL];
    }

    Info<< nl << "Finding real cells in range of referred cells" << endl;

    forAll(*this, rC)
    {
        const polyMesh& mesh(il_.mesh());

        referredCell<ParticleType>& refCell = (*this)[rC];

        DynamicList<label> realCellsFoundInRange;

        const vectorList& refCellPoints = refCell.vertexPositions();

        forAll(rFacesWRRP, rCF)
        {
            const label f(rFacesWRRP[rCF]);

            if (il_.testPointFaceDistance(refCellPoints,f))
            {
                const label cellO(mesh.faceOwner()[f]);

                if (findIndex(realCellsFoundInRange, cellO) == -1)
                {
                    realCellsFoundInRange.append(cellO);
                }

                if (mesh.isInternalFace(f))
                {
                    // boundary faces will not have neighbour information

                    const label cellN(mesh.faceNeighbour()[f]);

                    if (findIndex(realCellsFoundInRange, cellN) == -1)
                    {
                        realCellsFoundInRange.append(cellN);
                    }
                }
            }
        }

        forAll(rPointsWRRP, rCP)
        {
            const label p(rPointsWRRP[rCP]);

            if (il_.testPointFaceDistance(p,refCell))
            {
                const labelList& pCells(mesh.pointCells()[p]);

                forAll(pCells, pC)
                {
                    const label cellI(pCells[pC]);

                    if (findIndex(realCellsFoundInRange, cellI) == -1)
                    {
                        realCellsFoundInRange.append(cellI);
                    }
                }
            }
        }


        const edgeList& refCellEdges = refCell.edges();

        forAll(rEdgesWRRP, rCE)
        {
            const label edgeIIndex(rEdgesWRRP[rCE]);

            const edge& eI(mesh.edges()[edgeIIndex]);

            forAll(refCellEdges, rCE)
            {
                const edge& eJ(refCellEdges[rCE]);

                if
                (
                    il_.testEdgeEdgeDistance
                    (
                        eI,
                        refCellPoints[eJ.start()],
                        refCellPoints[eJ.end()]
                    )
                )
                {
                    const labelList& eICells(mesh.edgeCells()[edgeIIndex]);

                    forAll(eICells, eIC)
                    {
                        const label cellI(eICells[eIC]);

                        if (findIndex(realCellsFoundInRange, cellI) == -1)
                        {
                            realCellsFoundInRange.append(cellI);
                        }
                    }
                }
            }
        }

//         scalar rCutMaxSqr = molCloud_.rCutMax()*molCloud_.rCutMax();
//
//         forAll (molCloud_.mesh().points(), pointIIndex)
//         {
//             const point& ptI
//             (
//                 molCloud_.mesh().points()[pointIIndex]
//             );
//
//             forAll(refCellPoints, rCP)
//             {
//                 if (magSqr(ptI - refCellPoints[rCP]) <= rCutMaxSqr)
//                 {
//                     const labelList& ptICells
//                     (
//                         molCloud_.mesh().pointCells()[pointIIndex]
//                     );
//
//                     forAll(ptICells, pIC)
//                     {
//                         const label cellI(ptICells[pIC]);
//
//                         if (findIndex(realCellsFoundInRange, cellI) == -1)
//                         {
//                             realCellsFoundInRange.append(cellI);
//                         }
//                     }
//                 }
//             }
//         }

        refCell.realCells() = realCellsFoundInRange.shrink();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::referredCellList<ParticleType>::referredCellList
(
    interactionLists<ParticleType>& il,
    bool pointPointListBuild
)
:
    List<referredCell<ParticleType> >(),
    il_(il),
    cloud_(il_.mesh(), "referredCellListDummyCloud", IDLList<ParticleType>())

{
    buildReferredCellList(pointPointListBuild);
}


template<class ParticleType>
Foam::referredCellList<ParticleType>::referredCellList
(
    interactionLists<ParticleType>& il
)
:
    List<referredCell<ParticleType> >(),
    il_(il),
    cloud_(il_.mesh(), IDLList<ParticleType>())
{
    Info<< "Read referredCellList from disk not implemented" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::referredCellList<ParticleType>::~referredCellList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::referredCellList<ParticleType>::referParticles
(
    const List<DynamicList<ParticleType*> >& cellOccupancy
)
{
    // Clear all existing referred particles

    forAll(*this, i)
    {
        (*this)[i].clear();
    }

    cloud_.clear();

    // Create referred particles for sending using cell occupancy and
    // cellSendingReferralLists

    forAll(il_.cellSendingReferralLists(), cSRL)
    {
        const sendingReferralList& sRL
        (
            il_.cellSendingReferralLists()[cSRL]
        );

        IDLList<ParticleType> particlesToReferOut;

        DynamicList<label> destinationReferredCell;

        forAll(sRL, sRLI)
        {
            List<ParticleType*> realParticles = cellOccupancy[sRL[sRLI]];

            forAll (realParticles, rM)
            {
                const ParticleType& particle = *realParticles[rM];

                particlesToReferOut.append(particle.clone().ptr());

                destinationReferredCell.append(sRLI);
            }
        }

        // Send lists of referred particles to other processors

        if (sRL.destinationProc() != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                OPstream toInteractingProc
                (
                    Pstream::blocking,
                    sRL.destinationProc()
                );

                toInteractingProc
                    << destinationReferredCell
                    << particlesToReferOut;
            }
        }
        else
        {
            // Refer particles directly for referred cells on the same
            // processor.

            const receivingReferralList& rRL
            (
                il_.cellReceivingReferralLists()[cSRL]
            );

            // forAll(rRL, rRLI)
            // {
            //     forAll(rRL[rRLI], rC)
            //     {
            //         // referredCell<ParticleType>& refCellToRefParticlesTo =
            //         //     (*this)[rRL[rRLI][rC]];

            //         // refCellToRefParticlesTo.referInParticles
            //         // (
            //         //     particlesToReferOut[rRLI]
            //         // );
            //     }
            // }

            forAllConstIter
            (
                typename IDLList<ParticleType>,
                particlesToReferOut,
                iter
            )
            {
                cloud_.addParticle(iter().clone().ptr());
            }

            label particleI = 0;

            forAllIter
            (
                typename IDLList<ParticleType>,
                particlesToReferOut,
                iter
            )
            {
                ParticleType& p = iter();

                labelList refCellsToReferTo =
                    rRL[destinationReferredCell[particleI]];

                forAll(refCellsToReferTo, refCellI)
                {
                    referredCell<ParticleType>& refCellToRefParticlesTo =
                        (*this)[refCellsToReferTo[refCellI]];

                    refCellToRefParticlesTo.referInParticle(p);
                }

                particleI++;
            }
        }
    }

    // Receive referred particle lists to and distribute to referredCells
    // according to cellReceivingReferralLists, referredCells deal with the
    // transformations themselves

    forAll(il_.cellReceivingReferralLists(), cRRL)
    {
        const receivingReferralList& rRL
        (
            il_.cellReceivingReferralLists()[cRRL]
        );

        IDLList<ParticleType> particlesToReferIn;

        List<label> destinationReferredCell;

        if (rRL.sourceProc() != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                IPstream fromInteractingProc
                (
                    Pstream::blocking,
                    rRL.sourceProc()
                );

                fromInteractingProc
                    >> destinationReferredCell;

                particlesToReferIn = IDLList<ParticleType>
                (
                    fromInteractingProc,
                    typename ParticleType::iNew(cloud_)
                );

                forAllConstIter
                (
                    typename IDLList<ParticleType>,
                    particlesToReferIn,
                    iter
                )
                {
                    cloud_.addParticle(iter().clone().ptr());
                }
            }

            // forAll(rRL, rRLI)
            // {
            //     forAll(rRL[rRLI], rC)
            //     {
            //         referredCell<ParticleType>& refCellToRefParticlesTo =
            //             (*this)[rRL[rRLI][rC]];

            //         refCellToRefParticlesTo.referInParticles
            //         (
            //             particlesToReferIn[rRLI]
            //         );
            //     }
            // }

            label particleI = 0;

            forAllIter
            (
                typename IDLList<ParticleType>,
                particlesToReferIn,
                iter
            )
            {
                ParticleType& p = iter();

                labelList refCellsToReferTo =
                    rRL[destinationReferredCell[particleI]];

                forAll(refCellsToReferTo, refCellI)
                {
                    referredCell<ParticleType>& refCellToRefParticlesTo =
                        (*this)[refCellsToReferTo[refCellI]];

                    refCellToRefParticlesTo.referInParticle(p);
                }

                particleI++;
            }
        }
    }

    if (il_.mesh().time().outputTime())
    {
        Particle<ParticleType>::writeFields(cloud_);
    }
}


// ************************************************************************* //
