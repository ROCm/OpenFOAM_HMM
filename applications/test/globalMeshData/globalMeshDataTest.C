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

Application
    globalMeshDataTest

Description
    Test global point communication

\*---------------------------------------------------------------------------*/

#include "globalMeshData.H"
#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "mapDistribute.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    const globalMeshData& globalData = mesh.globalData();
    const indirectPrimitivePatch& coupledPatch = globalData.coupledPatch();


    // Test:print shared points
    {
        const labelListList& globalPointSlaves =
            globalData.globalPointSlaves();
        const mapDistribute& globalPointSlavesMap =
            globalData.globalPointSlavesMap();

        pointField coords(globalPointSlavesMap.constructSize());
        SubList<point>(coords, coupledPatch.nPoints()).assign
        (
            coupledPatch.localPoints()
        );

        // Exchange data
        globalPointSlavesMap.distribute(coords);

        // Print
        forAll(globalPointSlaves, pointI)
        {
            const labelList& slavePoints = globalPointSlaves[pointI];

            if (slavePoints.size() > 0)
            {
                Pout<< "Master point:" << pointI
                    << " coord:" << coords[pointI]
                    << " connected to slave points:" << endl;

                forAll(slavePoints, i)
                {
                    Pout<< "    " << coords[slavePoints[i]] << endl;
                }
            }
        }
    }



    // Test: point to faces addressing
    {
        const labelListList& globalPointBoundaryFaces =
            globalData.globalPointBoundaryFaces();
        const mapDistribute& globalPointBoundaryFacesMap =
            globalData.globalPointBoundaryFacesMap();

        label nBnd = mesh.nFaces()-mesh.nInternalFaces();

        pointField fc(globalPointBoundaryFacesMap.constructSize());
        SubList<point>(fc, nBnd).assign
        (
            primitivePatch
            (
                SubList<face>
                (
                    mesh.faces(),
                    nBnd,
                    mesh.nInternalFaces()
                ),
                mesh.points()
            ).faceCentres()
        );

        // Exchange data
        globalPointBoundaryFacesMap.distribute(fc);

        // Print
        forAll(globalPointBoundaryFaces, pointI)
        {
            const labelList& bFaces = globalPointBoundaryFaces[pointI];

            Pout<< "Point:" << pointI
                << " at:" << coupledPatch.localPoints()[pointI]
                << " connected to faces:" << endl;

            forAll(bFaces, i)
            {
                Pout<< "    " << fc[bFaces[i]] << endl;
            }
        }
    }





    // Test:point to cells addressing
    {
        const labelList& boundaryCells = globalData.boundaryCells();
        const labelListList& globalPointBoundaryCells =
            globalData.globalPointBoundaryCells();
        const mapDistribute& globalPointBoundaryCellsMap =
            globalData.globalPointBoundaryCellsMap();

        pointField cc(globalPointBoundaryCellsMap.constructSize());
        forAll(boundaryCells, i)
        {
            cc[i] = mesh.cellCentres()[boundaryCells[i]];
        }

        // Exchange data
        globalPointBoundaryCellsMap.distribute(cc);

        // Print
        forAll(globalPointBoundaryCells, pointI)
        {
            const labelList& bCells = globalPointBoundaryCells[pointI];

            Pout<< "Point:" << pointI
                << " at:" << coupledPatch.localPoints()[pointI]
                << " connected to cells:" << endl;

            forAll(bCells, i)
            {
                Pout<< "    " << cc[bCells[i]] << endl;
            }
        }
    }



    // Test:print shared edges
    {
        const labelListList& globalEdgeSlaves =
            globalData.globalEdgeSlaves();
        const mapDistribute& globalEdgeSlavesMap =
            globalData.globalEdgeSlavesMap();

        // Test: distribute edge centres
        pointField ec(globalEdgeSlavesMap.constructSize());
        forAll(coupledPatch.edges(), edgeI)
        {
            ec[edgeI] = coupledPatch.edges()[edgeI].centre
            (
                coupledPatch.localPoints()
            );
        }

        // Exchange data
        globalEdgeSlavesMap.distribute(ec);

        // Print
        forAll(globalEdgeSlaves, edgeI)
        {
            const labelList& slaveEdges = globalEdgeSlaves[edgeI];

            if (slaveEdges.size() > 0)
            {
                Pout<< "Master edge:" << edgeI
                    << " centre:" << ec[edgeI]
                    << " connected to slave edges:" << endl;

                forAll(slaveEdges, i)
                {
                    Pout<< "    " << ec[slaveEdges[i]] << endl;
                }
            }
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
