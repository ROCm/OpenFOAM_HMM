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

#include "ReferredCellList.H"
#include "InteractionLists.H"
#include "polyBoundaryMeshEntries.H"
#include "PstreamCombineReduceOps.H"
#include "Time.H"
#include "globalMeshData.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::ReferredCellList<ParticleType>::buildReferredCellList
(
    bool pointPointListBuild
)
{
    Info<< "    Building list of referred interaction neighbours" << endl;

    const polyMesh& mesh(il_.mesh());

    DynamicList<ReferredCell<ParticleType> > referredInteractionList;

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
            FatalErrorIn ("ReferredCellList.C")
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
                        ReferredCell<ParticleType> testRefCell
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
                            > InteractionLists<ParticleType>::transTol
                     || sum(mag(refTrans-refTrans[0]))/refTrans.size()
                            > InteractionLists<ParticleType>::transTol
                    )
                    {
                        FatalErrorIn ("ReferredCellList.C")
                            << nl << "Face pairs on patch "
                            << patch.name()
                            << ", segment " << patchNames[nUP]
                            << " do not give the same referring "
                            << " transformations to within tolerance of "
                            << InteractionLists<ParticleType>::transTol << nl
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
                            ReferredCell<ParticleType> testRefCell
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
                                > InteractionLists<ParticleType>::transTol
                         || sum(mag(refTrans - refTrans[0]))/(patch.size()/2)
                                > InteractionLists<ParticleType>::transTol
                        )
                        {
                            FatalErrorIn ("ReferredCellList.C")
                                << nl << "Face pairs on patch "
                                << patch.name()
                                << " do not give the same referring "
                                << " transformations to within tolerance of "
                                << InteractionLists<ParticleType>::transTol
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

                            ReferredCell<ParticleType> cellToRefer
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

                    labelList ReferredCellsFoundInRange
                    (
                        il_.ReferredCellsInRangeOfSegment
                        (
                            referredInteractionList,
                            meshFacesOnThisSegment,
                            meshEdgesOnThisSegment,
                            meshPointsOnThisSegment
                        )
                    );

                    forAll(ReferredCellsFoundInRange,cFIR)
                    {
                        ReferredCell<ParticleType>& existingRefCell =
                            referredInteractionList
                            [
                                ReferredCellsFoundInRange[cFIR]
                            ];

                        ReferredCell<ParticleType> cellToReRefer =
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

                            ReferredCell<ParticleType> cellToRefer
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

                    ReferredCellsFoundInRange =
                        il_.ReferredCellsInRangeOfSegment
                        (
                            referredInteractionList,
                            meshFacesOnThisSegment,
                            meshEdgesOnThisSegment,
                            meshPointsOnThisSegment
                        );

                    forAll(ReferredCellsFoundInRange,cFIR)
                    {
                        ReferredCell<ParticleType>& existingRefCell =
                            referredInteractionList
                            [
                                ReferredCellsFoundInRange[cFIR]
                            ];

                        ReferredCell<ParticleType> cellToReRefer =
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

                DynamicList<ReferredCell<ParticleType> >
                    ReferredCellsToTransfer;

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
                            FatalErrorIn ("ReferredCellList.C")
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

                                ReferredCell<ParticleType> cellToRefer
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

                                ReferredCellsToTransfer.append(cellToRefer);

                                // add real cells found in range of processor
                                // patch to whole mesh list

                                if (findIndex (rCellsWRRP, realCell) == -1)
                                {
                                    rCellsWRRP.append(realCell);
                                }
                            }
                        }

                        referredInteractionList.shrink();

                        labelList ReferredCellsFoundInRange
                        (
                            il_.ReferredCellsInRangeOfSegment
                            (
                                referredInteractionList,
                                meshFacesOnThisSegment,
                                meshEdgesOnThisSegment,
                                meshPointsOnThisSegment
                            )
                        );

                        forAll(ReferredCellsFoundInRange,cFIR)
                        {
                            ReferredCell<ParticleType>& existingRefCell =
                                referredInteractionList
                                [
                                    ReferredCellsFoundInRange[cFIR]
                                ];

                            ReferredCell<ParticleType> cellToReRefer =
                                existingRefCell.reRefer
                                (
                                    patch.faceCentres()[faceT],
                                    neighbFaceCentres[faceT],
                                    patch.faceNormals()[faceT],
                                    neighbFaceAreas[faceT]
                                   /(mag(neighbFaceAreas[faceT]) + VSMALL)
                                );

                            ReferredCellsToTransfer.append(cellToReRefer);
                        }
                    }
                }

                ReferredCellsToTransfer.shrink();

                // Send these cells to the neighbouring processor.

                {
                    OPstream toNeighbProc
                    (
                        Pstream::blocking,
                        patch.neighbProcNo()
                    );

                    toNeighbProc << ReferredCellsToTransfer;
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

                List<ReferredCell<ParticleType> >
                    ReferredCellsFromNeighbour(patch.size());

                {
                    IPstream fromNeighbProc
                    (
                        Pstream::blocking,
                        patch.neighbProcNo()
                    );

                    fromNeighbProc >> ReferredCellsFromNeighbour;
                }

                // Check to see if they are duplicates, if not append
                // them to the referredInteractionList

                forAll(ReferredCellsFromNeighbour,rCFN)
                {
                    ReferredCell<ParticleType>& cellToRefer =
                        ReferredCellsFromNeighbour[rCFN];

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

        Info<< "        Cells added this iteration: "
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

    Info<< "    Finding real cells in range of referred cells" << endl;

    forAll(*this, rC)
    {
        const polyMesh& mesh(il_.mesh());

        ReferredCell<ParticleType>& refCell = (*this)[rC];

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

        refCell.realCells() = realCellsFoundInRange.shrink();
    }
}


template<class ParticleType>
void Foam::ReferredCellList<ParticleType>::buildCellReferralLists()
{
    Info<< "    Determining particle referring schedule" << endl;

    DynamicList<label> referralProcs;

    // Run through all ReferredCells to build list of interacting processors

    forAll(*this, refCellI)
    {
        const ReferredCell<ParticleType>& refCell((*this)[refCellI]);

        if (findIndex(referralProcs, refCell.sourceProc()) == -1)
        {
            referralProcs.append(refCell.sourceProc());
        }
    }

    List<DynamicList<label> > cellSendingReferralLists(referralProcs.size());

    List<DynamicList<DynamicList<label> > >
        cellReceivingReferralLists(referralProcs.size());

    // Run through all ReferredCells again building up send and receive info

    forAll(*this, refCellI)
    {
        const ReferredCell<ParticleType>& rC((*this)[refCellI]);

        label rPI = findIndex(referralProcs, rC.sourceProc());

        DynamicList<DynamicList<label> >& rRL(cellReceivingReferralLists[rPI]);

        DynamicList<label>& sRL(cellSendingReferralLists[rPI]);

        label existingSource = findIndex(sRL, rC.sourceCell());

        // Check to see if this source cell has already been allocated to
        // come to this processor.  If not, add the source cell to the sending
        // list and add the current referred cell to the receiving list.

        // It shouldn't be possible for the sending and receiving lists to be
        // different lengths, because their append operations happen at the
        // same time.

        if (existingSource == -1)
        {
            sRL.append(rC.sourceCell());

            rRL.append(DynamicList<label>(labelList(1, refCellI)));
        }
        else
        {
            rRL[existingSource].append(refCellI);
        }
    }

    // It is assumed that cell exchange is reciprocal, if proc A has cells to
    // send to proc B, then proc B must have some to send to proc A.

    cellReceivingReferralLists_.setSize(referralProcs.size());

    cellSendingReferralLists_.setSize(referralProcs.size());

    forAll(referralProcs, rPI)
    {
        DynamicList<DynamicList<label> >& rRL(cellReceivingReferralLists[rPI]);

        labelListList translLL(rRL.size());

        forAll(rRL, rRLI)
        {
            translLL[rRLI] = rRL[rRLI];
        }

        cellReceivingReferralLists_[rPI] = receivingReferralList
        (
            referralProcs[rPI],
            translLL
        );
    }

    // Send sendingReferralLists to each interacting processor.

    forAll(referralProcs, rPI)
    {
        DynamicList<label>& sRL(cellSendingReferralLists[rPI]);

        if (referralProcs[rPI] != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                OPstream toInteractingProc
                (
                    Pstream::blocking,
                    referralProcs[rPI]
                );

                toInteractingProc << sendingReferralList
                (
                    Pstream::myProcNo(),
                    sRL
                );
            }
        }
    }

    // Receive sendingReferralLists from each interacting processor.

    forAll(referralProcs, rPI)
    {
        if (referralProcs[rPI] != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                IPstream fromInteractingProc
                (
                    Pstream::blocking,
                    referralProcs[rPI]
                );

                fromInteractingProc >> cellSendingReferralLists_[rPI];
            }
        }
        else
        {
            cellSendingReferralLists_[rPI] = sendingReferralList
            (
                Pstream::myProcNo(),
                cellSendingReferralLists[rPI]
            );
        }
    }
}


template<class ParticleType>
void Foam::ReferredCellList<ParticleType>::storeParticles
(
    const receivingReferralList& rRL,
    const labelList& sourceReferredCell,
    IDLList<ParticleType>& particlesToReferIn
)
{
    label particleI = 0;

    forAllIter
    (
        typename IDLList<ParticleType>,
        particlesToReferIn,
        referInIter
    )
    {
        ParticleType& p = referInIter();

        labelList refCellsToReferTo =
            rRL[sourceReferredCell[particleI]];

        forAll(refCellsToReferTo, refCellI)
        {
            ReferredCell<ParticleType>& refCellToRefParticlesTo =
                (*this)[refCellsToReferTo[refCellI]];

            refCellToRefParticlesTo.referInParticle
            (
                p.clone().ptr()
            );
        }

        particleI++;
    }

    particlesToReferIn.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ReferredCellList<ParticleType>::ReferredCellList
(
    InteractionLists<ParticleType>& il,
    bool pointPointListBuild
)
:
    List<ReferredCell<ParticleType> >(),
    il_(il),
    cloud_(il_.mesh(), "referredParticleCloud", IDLList<ParticleType>()),
    cellSendingReferralLists_(),
    cellReceivingReferralLists_()
{
    buildReferredCellList(pointPointListBuild);

    buildCellReferralLists();

    writeReferredCells();
}


template<class ParticleType>
Foam::ReferredCellList<ParticleType>::ReferredCellList
(
    InteractionLists<ParticleType>& il
)
:
    List<ReferredCell<ParticleType> >(),
    il_(il),
    cloud_(il_.mesh(), IDLList<ParticleType>())
{
    Info<< "    Read ReferredCellList from disk not implemented" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::ReferredCellList<ParticleType>::~ReferredCellList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::ReferredCellList<ParticleType>::referParticles
(
    const List<DynamicList<ParticleType*> >& cellOccupancy
)
{
    Info<< "    Refer particles" << endl;

    // Clear all existing referred particles

    forAll(*this, i)
    {
        (*this)[i].clear();
    }

    // Create referred particles for sending using cell occupancy and
    // cellSendingReferralLists

    forAll(cellSendingReferralLists_, cSRL)
    {
        const sendingReferralList& sRL
        (
            cellSendingReferralLists_[cSRL]
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

            // Remove particles after transfer
            particlesToReferOut.clear();
        }
        else
        {
            // Refer particles directly for referred cells on the same
            // processor.

            const receivingReferralList& rRL
            (
                cellReceivingReferralLists_[cSRL]
            );

            storeParticles(rRL, destinationReferredCell, particlesToReferOut);
        }

    }

    // Receive referred particle lists to and distribute to ReferredCells
    // according to cellReceivingReferralLists, ReferredCells deal with the
    // transformations themselves

    forAll(cellReceivingReferralLists_, cRRL)
    {
        const receivingReferralList& rRL
        (
            cellReceivingReferralLists_[cRRL]
        );

        IDLList<ParticleType> particlesToReferIn;

        labelList sourceReferredCell;

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
                    >> sourceReferredCell;

                particlesToReferIn = IDLList<ParticleType>
                (
                    fromInteractingProc,
                    typename ParticleType::iNew(cloud_)
                );
            }

            storeParticles(rRL, sourceReferredCell, particlesToReferIn);
        }
    }

    bool writeCloud = false;

    if (il_.mesh().time().outputTime() && writeCloud)
    {
        cloud_.clear();

        forAll(*this, refCellI)
        {
            ReferredCell<ParticleType>& refCell = (*this)[refCellI];

            forAllIter
            (
                typename IDLList<ParticleType>,
                refCell,
                iter
            )
            {
                cloud_.addParticle(iter().clone().ptr());
            }
        }

        Particle<ParticleType>::writeFields(cloud_);

        cloud_.clear();
    }
}


template<class ParticleType>
void Foam::ReferredCellList<ParticleType>::writeReferredCells() const
{
    fileName fName = il_.mesh().time().path()/"referredCells.obj";

    Info<< "    Writing " << fName.name() << endl;

    OFstream file(fName);

    label vertexOffset = 1;

    forAll(*this, refCellI)
    {
        const ReferredCell<ParticleType>& refCell = (*this)[refCellI];

        const vectorList& refCellPts = refCell.vertexPositions();

        const labelListList& refCellFaces = refCell.faces();

        forAll(refCellPts, ptI)
        {
            file<< "v "
                << refCellPts[ptI].x() << " "
                << refCellPts[ptI].y() << " "
                << refCellPts[ptI].z()
                << nl;
        }

        forAll(refCellFaces, faceI)
        {
            file<< "f ";

            forAll(refCellFaces[faceI], fPtI)
            {
                file<< " " << refCellFaces[faceI][fPtI] + vertexOffset;
            }

            file<< nl;
        }

        vertexOffset += refCellPts.size();
    }

    file.flush();
}

// ************************************************************************* //
