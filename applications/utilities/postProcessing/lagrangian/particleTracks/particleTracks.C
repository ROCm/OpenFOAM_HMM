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
    particleTracks

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "passiveParticle.H"
#include "SortableList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject origProcHeader
        (
            "origProc",
            runTime.timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::MUST_READ
        );
        IOobject idHeader
        (
            "id",
            runTime.timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::MUST_READ
        );
        if (idHeader.headerOk() && origProcHeader.headerOk())
        {
            IOField<label> origProc(origProcHeader);
            IOField<label> id(idHeader);
            forAll(id, i)
            {
                maxIds[origProc[i]] = max(maxIds[origProc[i]], id[i]);
            }
        }
    }
    Pstream::listCombineGather(maxIds, maxOp<label>());
    Pstream::listCombineScatter(maxIds);
    labelList numIds = maxIds + 1;

    Info<< "numIds = " << numIds << endl;

    // calc starting ids for particles on each processor
    List<label> startIds(numIds.size(), 0);
    for (label i = 0; i < numIds.size()-1; i++)
    {
        startIds[i+1] += startIds[i] + numIds[i];
    }
    label nParticle = startIds[startIds.size()-1] + numIds[startIds.size()-1];

    // number of tracks to generate
    label nTracks = nParticle/sampleFrequency;

    // storage for all particle tracks
    List<DynamicList<vector> > allTracks(nTracks);

    Info<< "\nGenerating " << nTracks << " particle tracks" << nl << endl;
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject positionsHeader
        (
            "positions",
            runTime.timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        IOobject origProcHeader
        (
            "origProc",
            runTime.timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        IOobject idHeader
        (
            "id",
            runTime.timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if
        (
            positionsHeader.headerOk()
         && origProcHeader.headerOk()
         && idHeader.headerOk()
        )
        {
            Info<< "    Reading particle positions" << endl;
            Cloud<passiveParticle> myCloud(mesh, cloudName, false);

            pointField positions(myCloud.size(), vector::zero);
            label i = 0;
            forAllConstIter(Cloud<passiveParticle>, myCloud, iter)
            {
                positions[i++] = iter().position();
            }
            IOField<label> id(idHeader);
            IOField<label> origProc(origProcHeader);

            List<pointField> allPositions(Pstream::nProcs());
            allPositions[Pstream::myProcNo()] = positions;
            Pstream::gatherList(allPositions);

            List<labelList> allIds(Pstream::nProcs());
            allIds[Pstream::myProcNo()] = id;
            Pstream::gatherList(allIds);

            List<labelList> allOrigProcs(Pstream::nProcs());
            allOrigProcs[Pstream::myProcNo()] = origProc;
            Pstream::gatherList(allOrigProcs);

            if (Pstream::master())
            {
                forAll(allPositions, procI)
                {
                    forAll(allPositions[procI], i)
                    {
                        label globalId =
                            startIds[allOrigProcs[procI][i]]
                        + allIds[procI][i];

                        if (globalId % sampleFrequency == 0)
                        {
                            label trackId = globalId/sampleFrequency;
                            if (allTracks[trackId].size() < maxPositions)
                            {
                                allTracks[trackId].append
                                (
                                    allPositions[procI][i]
                                );
                            }
                        }
                    }
                }
            }
            else
            {
                Info<< "    No particles read" << nl << endl;
            }
        }
    }

    if (Pstream::master())
    {
        Info<< "Writing particle tracks" << nl << endl;

        OFstream vtkTracks("particleTracks.vtk");

        // Total number of points in tracks + 1 per track
        label nPoints = 0;
        forAll(allTracks, trackI)
        {
            nPoints += allTracks[trackI].size();
        }

        vtkTracks
            << "# vtk DataFile Version 2.0\n"
            << "particleTracks" << nl
            << "ASCII\n"
            << "DATASET POLYDATA\n"
            << "POINTS " << nPoints << " float\n";

        // Write track points to file
        forAll(allTracks, trackI)
        {
            forAll(allTracks[trackI], i)
            {
                const vector& pt = allTracks[trackI][i];
                vtkTracks << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
            }
        }

        // write track (line) connectivity to file
        vtkTracks << "LINES " << nTracks << ' ' << nPoints+nTracks << nl;

        // Write ids of track points to file
        label globalPtI = 0;
        forAll(allTracks, trackI)
        {
            vtkTracks << allTracks[trackI].size();

            forAll(allTracks[trackI], i)
            {
                vtkTracks << ' ' << globalPtI;
                globalPtI++;
            }

            vtkTracks << nl;
        }

        Info<< "end" << endl;
    }

    return 0;
}


// ************************************************************************* //
