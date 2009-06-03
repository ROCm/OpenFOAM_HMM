/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "extendedStencil.H"
#include "globalIndex.H"
#include "syncTools.H"
#include "SortableList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistribute> Foam::extendedStencil::calcDistributeMap
(
    const globalIndex& globalNumbering,
    labelListList& faceStencil
)
{
    const label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();


    // Convert stencil to schedule
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // We now know what information we need from other processors. This needs
    // to be converted into what information I need to send as well
    // (mapDistribute)


    // 1. Construct per processor compact addressing of the global cells
    //    needed. The ones from the local processor are not included since
    //    these are always all needed.
    List<Map<label> > globalToProc(Pstream::nProcs());
    {
        const labelList& procPatchMap = mesh_.globalData().procPatchMap();
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Presize with (as estimate) size of patch to neighbour.
        forAll(procPatchMap, procI)
        {
            if (procPatchMap[procI] != -1)
            {
                globalToProc[procI].resize
                (
                    patches[procPatchMap[procI]].size()
                );
            }
        }

        // Collect all (non-local) globalcells/faces needed.
        forAll(faceStencil, faceI)
        {
            const labelList& stencilCells = faceStencil[faceI];

            forAll(stencilCells, i)
            {
                label globalCellI = stencilCells[i];
                label procI = globalNumbering.whichProcID(stencilCells[i]);

                if (procI != Pstream::myProcNo())
                {
                    label nCompact = globalToProc[procI].size();
                    globalToProc[procI].insert(globalCellI, nCompact);
                }
            }
        }
        // Sort global cells needed (not really necessary)
        forAll(globalToProc, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                Map<label>& globalMap = globalToProc[procI];

                SortableList<label> sorted(globalMap.toc().xfer());

                forAll(sorted, i)
                {
                    Map<label>::iterator iter = globalMap.find(sorted[i]);
                    iter() = i;
                }
            }
        }


        //forAll(globalToProc, procI)
        //{
        //    Pout<< "From processor:" << procI << " want cells/faces:" << endl;
        //    forAllConstIter(Map<label>, globalToProc[procI], iter)
        //    {
        //        Pout<< "    global:" << iter.key()
        //            << " local:" << globalNumbering.toLocal(procI, iter.key())
        //            << endl;
        //    }
        //    Pout<< endl;
        //}
    }


    // 2. The overall compact addressing is
    // - myProcNo first
    // - all other processors consecutively

    labelList compactStart(Pstream::nProcs());
    compactStart[Pstream::myProcNo()] = 0;
    label nCompact = mesh_.nCells()+nBnd;
    forAll(compactStart, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            compactStart[procI] = nCompact;
            nCompact += globalToProc[procI].size();
        }
    }


    // 3. Find out what to receive/send in compact addressing.
    labelListList recvCompact(Pstream::nProcs());
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            labelList wantedGlobals(globalToProc[procI].size());
            recvCompact[procI].setSize(globalToProc[procI].size());

            label i = 0;
            forAllConstIter(Map<label>, globalToProc[procI], iter)
            {
                wantedGlobals[i] = iter.key();
                recvCompact[procI][i] = compactStart[procI]+iter();
                i++;
            }

            // Send the global cell numbers I need from procI
            OPstream str(Pstream::blocking, procI);
            str << wantedGlobals;
        }
        else
        {
            recvCompact[procI] =
                compactStart[procI]
              + identity(mesh_.nCells()+nBnd);
        }
    }
    labelListList sendCompact(Pstream::nProcs());
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            // See what neighbour wants to receive (= what I need to send)

            IPstream str(Pstream::blocking, procI);
            labelList globalCells(str);

            labelList& procCompact = sendCompact[procI];
            procCompact.setSize(globalCells.size());

            // Convert from globalCells (all on my processor!) into compact
            // addressing
            forAll(globalCells, i)
            {
                label cellI = globalNumbering.toLocal(globalCells[i]);
                procCompact[i] = compactStart[Pstream::myProcNo()]+cellI;
            }
        }
        else
        {
            sendCompact[procI] = recvCompact[procI];
        }
    }

    // Convert stencil to compact numbering
    forAll(faceStencil, faceI)
    {
        labelList& stencilCells = faceStencil[faceI];

        forAll(stencilCells, i)
        {
            label globalCellI = stencilCells[i];
            label procI = globalNumbering.whichProcID(globalCellI);
            if (procI != Pstream::myProcNo())
            {
                label localCompact = globalToProc[procI][globalCellI];
                stencilCells[i] = compactStart[procI]+localCompact;
            }
            else
            {
                label localCompact = globalNumbering.toLocal(globalCellI);
                stencilCells[i] = compactStart[procI]+localCompact;
            }

        }
    }

    // Constuct map for distribution of compact data.
    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            nCompact,
            sendCompact,
            recvCompact,
            true            // reuse send/recv maps.
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedStencil::extendedStencil(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// ************************************************************************* //
