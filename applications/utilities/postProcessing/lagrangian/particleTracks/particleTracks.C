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

Foam::label checkHeader(IOobject& io)
{
    if (io.headerOk())
    {
        return 0;
    }
    else
    {
        Info<< "    no " << io.name() << " field" << endl;
        return 1;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"
#   include "createFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
    Parallel mode
      - nProcs known
      - scan through times to determine max id for each proc

      - scan through times to populate for each particle (indexed to procId):
          - id
          - time
      - construct single flat list from all procs of ids and times:
          - starting id = sum of ids of lower procs
      - sort particles according to time
      - join the dots...

    Serial mode
      - don't know if case was run in serial mode, or is a reconstructed
        parallel run
      - scan times to determine max origProcId and max id for each origProc

      - scan through times to populate for each particle (indexed to procId):
          - id
          - time
      - construct single flat list from all procs of ids and times:
          - starting id = sum of ids of lower procs
      - sort particles according to time
      - join the dots...


*/

    DynamicList<label> maxIds;
    labelHashSet validIds;
    if (Pstream::parRun())
    {
        maxIds(Pstream::myProcNo()) = -1;
        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);
            Info<< "Time = " << runTime.timeName() << endl;

            IOobject idHeader
            (
                "id",
                runTime.timeName(),
                cloud::prefix/cloudName,
                mesh,
                IOobject::MUST_READ
            );
            if (idHeader.headerOk())
            {
                IOField<label> id(idHeader);
                maxIds[Pstream::myProcNo()]
                    = max(maxIds[Pstream::myProcNo()], max(id));
            }
        }
        Pstream::listCombineGather(maxIds, maxOp<label>());
//        reduce(maxIds, maxOp<DynamicList<label> >());
    }
    else
    {
        Info<< "Assuming case is a reconstructed parallel run" << nl;
        Info<< "Scanning time folders for tracking information" << nl << endl;
        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);
            Info<< "Time = " << runTime.timeName() << endl;

            label proceed = 1;
            IOobject origProcHeader
            (
                "origProc",
                runTime.timeName(),
                cloud::prefix/cloudName,
                mesh,
                IOobject::MUST_READ
            );
            proceed -= checkHeader(origProcHeader);

            IOobject idHeader
            (
                "id",
                runTime.timeName(),
                cloud::prefix/cloudName,
                mesh,
                IOobject::MUST_READ
            );
            proceed -= checkHeader(idHeader);

            if (proceed == 1)
            {
                IOField<label> origProc(origProcHeader);
                IOField<label> id(idHeader);
                forAll(origProc, i)
                {
                    validIds.insert(origProc[i]);
                    if (maxIds.size() > origProc[i])
                    {
                        // proc group already exists
                        maxIds[origProc[i]] = max(maxIds[origProc[i]], id[i]);
                    }
                    else
                    {
                        maxIds(origProc[i]) = id[i];
                    }
                }
            }
        }
    }

    maxIds.shrink();

    // determine number of unique particles
    label nParticle = 0;
    forAll(maxIds, procI)
    {
        if (validIds.found(procI))
        {
            nParticle += maxIds[procI] + 1;
        }
    }

    // calc starting ids for particles on each processor
    List<label> startIds(validIds.size(), 0);
    for (label i = 0; i < validIds.size(); i++)
    {
        if (validIds.found(i))
        {
            startIds[i+1] += startIds[i] + maxIds[i];
        }
    }

    // number of tracks to generate
    label nTracks = nParticle/sampleFrequency;

    // storage for all particle tracks
    List<DynamicList<vector> > allTracks(nTracks);

    // storage for all particle times
    List<DynamicList<scalar> > allTimes(nTracks);

    Info<< "\nGenerating " << nTracks << " particle tracks" << nl << endl;
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        label proceed = 1;
        IOobject positionsHeader
        (
            "positions",
            runTime.timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::MUST_READ
        );
        proceed -= checkHeader(positionsHeader);

        IOobject idHeader
        (
            "id",
            runTime.timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::MUST_READ
        );
        proceed -= checkHeader(idHeader);

        if (proceed == 1)
        {
            Info<< "    Reading particle positions" << endl;
            Cloud<passiveParticle> positions(mesh, cloudName, false);

            if (positions.size() > 0)
            {
                IOField<label> id(idHeader);
                Info<< "    Reading particle ids" << nl << endl;

                label i = 0;
                forAllConstIter(Cloud<passiveParticle>, positions, iter)
                {
                    label globalId = startIds[Pstream::myProcNo()] + id[i++];
                    if (globalId % sampleFrequency == 0)
                    {
                        label trackId = globalId/sampleFrequency;
                        if (allTracks[trackId].size() < maxPositions)
                        {
                            allTracks[trackId].append(iter().position());
                            allTimes[trackId].append(timeDirs[timeI].value());
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

    if (Pstream::master)
    {
        Info<< "Writing particle tracks" << nl << endl;

        forAll(allTracks, trackI)
        {
            Pstream::listCombineGather(allTracks[trackI], plusEqOp<vector>());
            Pstream::listCombineGather(allTimes[trackI], plusEqOp<scalar>());
        }

        OFstream vtkTracks("particleTracks.vtk");

        // Total number of points in tracks + 1 per track
        label nPoints = 0;
        forAll(allTracks, trackI)
        {
            nPoints += allTracks[trackI].size();
        }
        reduce(nPoints, sumOp<label>());

        vtkTracks
            << "# vtk DataFile Version 2.0\n"
            << "particleTracks" << nl
            << "ASCII\n"
            << "DATASET POLYDATA\n"
            << "POINTS " << nPoints << " float\n";

        // Write track points to file
        forAll(allTracks, trackI)
        {
            SortableList<scalar> sortedTimes(allTimes[trackI]);

            const DynamicList<vector>& trackPts = allTracks[trackI];
            forAll(trackPts, i)
            {
                label ptI = sortedTimes.indices()[i];
                const vector& pt = trackPts[ptI];
                vtkTracks << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
            }
        }

        // write track (line) connectivity to file
        vtkTracks << "LINES " << nTracks << ' ' << nPoints+nTracks << nl;

        // Write track points to file
        label globalPtI = 0;
        forAll(allTracks, trackI)
        {
            const DynamicList<vector>& trackPts = allTracks[trackI];

            vtkTracks << trackPts.size();

            forAll(trackPts, i)
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
