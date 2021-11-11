/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Application
    particleTracks

Group
    grpPostProcessingUtilities

Description
    Generate particle tracks for cases that were computed using a
    tracked-parcel-type cloud.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "passiveParticleCloud.H"
#include "writer.H"
#include "ListOps.H"

#define createTrack(field, trackValues)                                        \
createTrackField                                                               \
(                                                                              \
    field,                                                                     \
    sampleFrequency,                                                           \
    maxPositions,                                                              \
    startIds,                                                                  \
    allOrigProcs,                                                              \
    allOrigIds,                                                                \
    trackValues                                                                \
);

#define setFields(fields, fieldNames)                                          \
setTrackFields                                                                 \
(                                                                              \
    obr,                                                                       \
    fields,                                                                    \
    fieldNames,                                                                \
    nTracks,                                                                   \
    startIds,                                                                  \
    allOrigProcs,                                                              \
    allOrigIds,                                                                \
    maxPositions,                                                              \
    sampleFrequency                                                            \
);

#define writeFields(fields, fieldNames, tracks, times, dirs)                   \
writeTrackFields                                                               \
(                                                                              \
    fields,                                                                    \
    fieldNames,                                                                \
    tracks,                                                                    \
    times,                                                                     \
    dirs,                                                                      \
    setFormat,                                                                 \
    formatOptions,                                                             \
    cloudName                                                                  \
);

using namespace Foam;

template<class Type>
void createTrackField
(
    const Field<Type>& values,
    const label sampleFrequency,
    const label maxPositions,
    const labelList& startIds,
    const List<labelField>& allOrigProcs,
    const List<labelField>& allOrigIds,
    List<DynamicList<Type>>& trackValues
)
{
    List<Field<Type>> procField(Pstream::nProcs());
    procField[Pstream::myProcNo()] = values;
    Pstream::gatherList(procField);

    if (!Pstream::master())
    {
        return;
    }

    const label nTracks = trackValues.size();

    forAll(procField, proci)
    {
        forAll(procField[proci], i)
        {
            const label globalId =
                startIds[allOrigProcs[proci][i]] + allOrigIds[proci][i];

            if (globalId % sampleFrequency == 0)
            {
                const label trackId = globalId/sampleFrequency;

                if
                (
                    trackId < nTracks
                 && trackValues[trackId].size() < maxPositions
                )
                {
                    trackValues[trackId].append(procField[proci][i]);
                }
            }
        }
    }
}


template<class Type>
void writeTrackFields
(
    List<List<DynamicList<Type>>>& fieldValues,
    const wordList& fieldNames,
    const PtrList<coordSet>& tracks,
    const List<scalarField>& times,
    const List<vectorField>& dirs,
    const word& setFormat,
    const dictionary& formatOptions,
    const word& cloudName
)
{
    if (fieldValues.empty())
    {
        return;
    }

    auto writerPtr = writer<Type>::New(setFormat, formatOptions);

    const fileName outFile(writerPtr().getFileName(tracks[0], wordList(0)));

    const fileName outPath
    (
        functionObject::outputPrefix/cloud::prefix/cloudName/"particleTracks"
    );
    mkDir(outPath);

    OFstream os(outPath/(pTraits<Type>::typeName & "tracks." + outFile.ext()));

    Info<< "Writing " << pTraits<Type>::typeName << " particle tracks in "
        << setFormat << " format to " << os.name() << endl;


    List<List<Field<Type>>> fields(fieldValues.size());
    forAll(fields, fieldi)
    {
        fields[fieldi].setSize(fieldValues[fieldi].size());
        forAll(fields[fieldi], tracki)
        {
            fields[fieldi][tracki].transfer(fieldValues[fieldi][tracki]);
        }
    }

    writerPtr().write(true, times, tracks, fieldNames, fields, os);
}


template<class Type>
Foam::label setTrackFields
(
    const objectRegistry& obr,
    List<List<DynamicList<Type>>>& fields,
    List<word>& fieldNames,
    const label nTracks,
    const labelList& startIds,
    const List<labelField>& allOrigProcs,
    const List<labelField>& allOrigIds,
    const label maxPositions,
    const label sampleFrequency
)
{
    const auto availableFieldPtrs = obr.lookupClass<IOField<Type>>();

    fieldNames = availableFieldPtrs.toc();

    if (Pstream::parRun())
    {
        Pstream::combineGather(fieldNames, ListOps::uniqueEqOp<word>());
        Pstream::combineScatter(fieldNames);
        Foam::sort(fieldNames);
    }

    const label nFields = fieldNames.size();

    if (fields.empty())
    {
        fields.setSize(nFields);
        fieldNames.setSize(nFields);
        forAll(fields, i)
        {
            fields[i].setSize(nTracks);
        }
    }

    forAll(fieldNames, fieldi)
    {
        const word& fieldName = fieldNames[fieldi];

        const auto* fldPtr = obr.cfindObject<IOField<Type>>(fieldName);

        createTrack
        (
            fldPtr ? static_cast<const Field<Type>>(*fldPtr) : Field<Type>(),
            fields[fieldi]
        );
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate a file of particle tracks for cases that were"
        " computed using a tracked-parcel-type cloud"
    );

    timeSelector::addOptions();
    #include "addRegionOption.H"

    #include "setRootCase.H"

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const fileName vtkPath(runTime.rootPath()/runTime.globalCaseName()/"VTK");
    mkDir(vtkPath);

    Info<< "Scanning times to determine track data for cloud " << cloudName
        << nl << endl;

    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "    Reading particle positions" << endl;
        passiveParticleCloud myCloud(mesh, cloudName);

        Info<< "    Read " << returnReduce(myCloud.size(), sumOp<label>())
            << " particles" << endl;

        for (const passiveParticle& p : myCloud)
        {
            const label origId = p.origId();
            const label origProc = p.origProc();

            // Handle case where we are processing particle data generated in
            // parallel using more cores than when running this application.
            if (origProc >= maxIds.size())
            {
                maxIds.setSize(origProc+1, -1);
            }

            maxIds[origProc] = max(maxIds[origProc], origId);
        }
    }

    label maxNProcs = returnReduce(maxIds.size(), maxOp<label>());

    Info<< "Detected particles originating from " << maxNProcs
        << " processors." << nl << endl;

    maxIds.setSize(maxNProcs, -1);

    Pstream::listCombineGather(maxIds, maxEqOp<label>());
    Pstream::listCombineScatter(maxIds);

    labelList numIds = maxIds + 1;

    Info<< nl << "Particle statistics:" << endl;
    forAll(maxIds, proci)
    {
        Info<< "    Found " << numIds[proci] << " particles originating"
            << " from processor " << proci << endl;
    }
    Info<< nl << endl;


    // Calculate starting ids for particles on each processor
    labelList startIds(numIds.size(), Zero);
    for (label i = 0; i < numIds.size()-1; ++i)
    {
        startIds[i+1] += startIds[i] + numIds[i];
    }
    label nParticle = startIds.last() + numIds[startIds.size()-1];


    // Number of tracks to generate
    const label nTracks =
        maxTracks > 0
      ? min(nParticle/sampleFrequency, maxTracks)
      : nParticle/sampleFrequency;

    // Storage for all particle tracks
    List<DynamicList<vector>> allTracks(nTracks);
    List<DynamicList<scalar>> allTrackTimes(nTracks);

    // Lists of field, tracki, trackValues
    //List<List<DynamicList<label>>> labelFields;
    List<List<DynamicList<scalar>>> scalarFields;
    List<List<DynamicList<vector>>> vectorFields;
    List<List<DynamicList<sphericalTensor>>> sphTensorFields;
    List<List<DynamicList<symmTensor>>> symTensorFields;
    List<List<DynamicList<tensor>>> tensorFields;
    //List<word> labelFieldNames;
    List<word> scalarFieldNames;
    List<word> vectorFieldNames;
    List<word> sphTensorFieldNames;
    List<word> symTensorFieldNames;
    List<word> tensorFieldNames;

    Info<< "\nGenerating " << nTracks << " particle tracks for cloud "
        << cloudName << nl << endl;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;

        // Read particles. Will be size 0 if no particles.
        Info<< "    Reading particle positions" << endl;
        passiveParticleCloud myCloud(mesh, cloudName);

        pointField localPositions(myCloud.size(), Zero);
        scalarField localTimes(myCloud.size(), Zero);

        List<labelField> allOrigIds(Pstream::nProcs());
        List<labelField> allOrigProcs(Pstream::nProcs());

        // Collect the track data on all processors that have positions
        allOrigIds[Pstream::myProcNo()].setSize(myCloud.size(), Zero);
        allOrigProcs[Pstream::myProcNo()].setSize(myCloud.size(), Zero);

        label i = 0;
        for (const passiveParticle& p : myCloud)
        {
            allOrigIds[Pstream::myProcNo()][i] = p.origId();
            allOrigProcs[Pstream::myProcNo()][i] = p.origProc();
            localPositions[i] = p.position();
            localTimes[i] = runTime.value();
            ++i;
        }

        // Collect the track data on the master processor
        Pstream::gatherList(allOrigIds);
        Pstream::gatherList(allOrigProcs);

        objectRegistry obr
        (
            IOobject
            (
                "cloudFields",
                runTime.timeName(),
                runTime
            )
        );

        myCloud.readFromFiles(obr, fieldNames);

        // Create track positions and track time fields
        // (not registered as IOFields)
        // Note: createTrack is a local #define to reduce arg count...
        createTrack(localPositions, allTracks);
        createTrack(localTimes, allTrackTimes);

        // Create the track fields
        // Note: setFields is a local #define to reduce arg count...
        //setFields(labelFields, labelFieldNames);
        setFields(scalarFields, scalarFieldNames);
        setFields(vectorFields, vectorFieldNames);
        setFields(sphTensorFields, sphTensorFieldNames);
        setFields(symTensorFields, symTensorFieldNames);
        setFields(tensorFields, tensorFieldNames);
    }


    if (Pstream::master())
    {
        PtrList<coordSet> tracks(allTracks.size());
        List<scalarField> times(tracks.size());
        forAll(tracks, tracki)
        {
            tracks.set
            (
                tracki,
                new coordSet("track" + Foam::name(tracki), "distance")
            );
            tracks[tracki].transfer(allTracks[tracki]);
            times[tracki].transfer(allTrackTimes[tracki]);
        }

        Info<< nl;

        const label Uid = vectorFieldNames.find(UName);
        List<vectorField> dirs(nTracks);

        if (Uid != -1)
        {
            const auto& UTracks = vectorFields[Uid];
            forAll(UTracks, tracki)
            {
                const auto& U = UTracks[tracki];
                dirs[tracki] = U/(mag(U) + ROOTVSMALL);
            }
        }

        // Write track fields
        // Note: writeFields is a local #define to reduce arg count...
        //writeFields(allLabelFields, labelFieldNames, tracks);
        writeFields(scalarFields, scalarFieldNames, tracks, times, dirs);
        writeFields(vectorFields, vectorFieldNames, tracks, times, dirs);
        writeFields(sphTensorFields, sphTensorFieldNames, tracks, times, dirs);
        writeFields(symTensorFields, symTensorFieldNames, tracks, times, dirs);
        writeFields(tensorFields, tensorFieldNames, tracks, times, dirs);
    }

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
