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

#include "domainDecomposition.H"
#include "decompositionMethod.H"
#include "cpuTime.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void domainDecomposition::distributeCells()
{
    Info<< "\nCalculating distribution of cells" << endl;

    cpuTime decompositionTime;


    // See if any faces need to have owner and neighbour on same processor
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelHashSet sameProcFaces;

    if (decompositionDict_.found("preservePatches"))
    {
        wordList pNames(decompositionDict_.lookup("preservePatches"));

        Info<< "Keeping owner and neighbour of faces in patches " << pNames
            << " on same processor" << endl;

        const polyBoundaryMesh& patches = boundaryMesh();

        forAll(pNames, i)
        {
            label patchI = patches.findPatchID(pNames[i]);

            if (patchI == -1)
            {
                FatalErrorIn("domainDecomposition::distributeCells()")
                    << "Unknown preservePatch " << pNames[i]
                    << endl << "Valid patches are " << patches.names()
                    << exit(FatalError);
            }

            const polyPatch& pp = patches[patchI];

            forAll(pp, i)
            {
                sameProcFaces.insert(pp.start() + i);
            }
        }
    }
    if (decompositionDict_.found("preserveFaceZones"))
    {
        wordList zNames(decompositionDict_.lookup("preserveFaceZones"));

        Info<< "Keeping owner and neighbour of faces in zones " << zNames
            << " on same processor" << endl;

        const faceZoneMesh& fZones = faceZones();

        forAll(zNames, i)
        {
            label zoneI = fZones.findZoneID(zNames[i]);

            if (zoneI == -1)
            {
                FatalErrorIn("domainDecomposition::distributeCells()")
                    << "Unknown preserveFaceZone " << zNames[i]
                    << endl << "Valid faceZones are " << fZones.names()
                    << exit(FatalError);
            }

            const faceZone& fz = fZones[zoneI];

            forAll(fz, i)
            {
                sameProcFaces.insert(fz[i]);
            }
        }
    }

    if (sameProcFaces.size() > 0)
    {
        Info<< "Selected " << sameProcFaces.size()
            << " faces whose owner and neighbour cell should be kept on the"
            << " same processor" << endl;
    }



    // Construct decomposition method and either do decomposition on
    // cell centres or on agglomeration


    autoPtr<decompositionMethod> decomposePtr = decompositionMethod::New
    (
        decompositionDict_,
        *this
    );

    if (sameProcFaces.size() == 0)
    {
        cellToProc_ = decomposePtr().decompose(cellCentres());
    }
    else
    {
        

        // Work the faces whose neighbours need to be kept together into an
        // agglomeration.

        // Per cell the region/agglomeration it is in
        labelList cellToRegion(nCells(), -1);

        // Current region
        label regionI = 0;

        labelHashSet freeRegions;

        forAllConstIter(labelHashSet, sameProcFaces, iter)
        {
            label patchI = boundaryMesh().whichPatch(iter.key());

            label own = faceOwner()[iter.key()];
            label nei = -1;

            if (patchI == -1)
            {
                nei = faceNeighbour()[iter.key()];
            }
            else if (isA<cyclicPolyPatch>(boundaryMesh()[patchI]))
            {
                const cyclicPolyPatch& pp =
                    refCast<const cyclicPolyPatch>(boundaryMesh()[patchI]);

                nei = faceOwner()[pp.transformGlobalFace(iter.key())];
            }

            if (nei != -1)
            {
                label ownRegion = cellToRegion[own];
                label neiRegion = cellToRegion[nei];

                if (ownRegion == -1 && neiRegion == -1)
                {
                    // Allocate new agglomeration
                    cellToRegion[own] = regionI;
                    cellToRegion[nei] = regionI;
                    regionI++;
                }
                else if (ownRegion != -1)
                {
                    // Owner already part of agglomeration. Add nei to it.
                    cellToRegion[nei] = ownRegion;
                }
                else if (neiRegion != -1)
                {
                    // nei already part of agglomeration. Add own to it.
                    cellToRegion[own] = neiRegion;
                }
                else if (ownRegion < neiRegion)
                {
                    // Renumber neiRegion
                    forAll(cellToRegion, cellI)
                    {
                        if (cellToRegion[cellI] == neiRegion)
                        {
                            cellToRegion[cellI] = ownRegion;
                        }
                    }
                    freeRegions.insert(neiRegion);
                }
                else if (ownRegion > neiRegion)
                {
                    // Renumber ownRegion
                    forAll(cellToRegion, cellI)
                    {
                        if (cellToRegion[cellI] == ownRegion)
                        {
                            cellToRegion[cellI] = neiRegion;
                        }
                    }
                    freeRegions.insert(ownRegion);
                }
            }
        }


        // Do all other cells
        forAll(cellToRegion, cellI)
        {
            if (cellToRegion[cellI] == -1)
            {
                cellToRegion[cellI] = regionI++;
            }
        }


        // Compact out freeRegions
        // ~~~~~~~~~~~~~~~~~~~~~~~

        {
            labelList compactRegion(regionI, -1);

            regionI = 0;

            forAll(compactRegion, i)
            {
                if (!freeRegions.found(compactRegion[i]))
                {
                    compactRegion[i] = regionI++;
                }
            }

            inplaceRenumber(compactRegion, cellToRegion);
        }


        // Determine region cell centres
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // This just takes the first cell in the region. Otherwise the problem
        // is with cyclics - if we'd average the region centre might be
        // somewhere in the middle of the domain which might not be anywhere
        // near any of the cells.

        const point greatPoint(GREAT, GREAT, GREAT);

        pointField regionCentres(regionI, greatPoint);

        forAll(cellToRegion, cellI)
        {
            label regionI = cellToRegion[cellI];

            if (regionCentres[regionI] == greatPoint)
            {
                regionCentres[regionI] = cellCentres()[cellI];
            }
        }


        // Do decomposition on agglomeration
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        cellToProc_ = decomposePtr().decompose(cellToRegion, regionCentres);
    }

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;
}


// ************************************************************************* //
