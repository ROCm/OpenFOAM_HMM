/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "surfaceNoise.H"
#include "surfaceReader.H"
#include "surfaceWriter.H"
#include "noiseFFT.H"
#include "graph.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace noiseModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceNoise, 0);
addToRunTimeSelectionTable(noiseModel, surfaceNoise, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void surfaceNoise::initialise(const dictionary& dict)
{
    // All reading performed on the master processor only
    if (Pstream::master())
    {
        // Create the surface reader
        const word readerType(dict.lookup("reader"));
        readerPtr_.reset(surfaceReader::New(readerType, inputFileName_).ptr());

        // Find the index of the pressure data
        const word pName(dict.lookupOrDefault<word>("pName", "p"));
        const List<word> fieldNames(readerPtr_->fieldNames(0));
        pIndex_ = findIndex(fieldNames, pName);
        if (pIndex_ == -1)
        {
            FatalErrorInFunction
                << "Unable to find pressure field name " << pName
                << " in list of available fields: " << fieldNames
                << exit(FatalError);
        }

        // Create the surface writer
        // - Could be done later, but since this utility can process a lot of
        //   data we can ensure that the user-input is correct prior to doing
        //   the heavy lifting
        const word writerType(dict.lookup("writer"));
        dictionary optDict
        (
            dict.subOrEmptyDict("writeOptions").subOrEmptyDict(writerType)
        );
        writerPtr_.reset(surfaceWriter::New(writerType, optDict).ptr());

        // Set the time range
        const instantList allTimes = readerPtr_->times();
        startTimeIndex_ = findStartTimeIndex(allTimes, startTime_);

        // Determine the windowing
        label nAvailableTimes = allTimes.size() - startTimeIndex_;
        label nRequiredTimes = windowModelPtr_->validate(nAvailableTimes);

        // Restrict times
        times_.setSize(nRequiredTimes);
        forAll(times_, timeI)
        {
            times_[timeI] = allTimes[timeI + startTimeIndex_].value();
        }
        deltaT_ = checkUniformTimeStep(times_);

        const meshedSurface& surf = readerPtr_->geometry();
        nFace_ = surf.size();
    }

    Pstream::scatter(pIndex_);
    Pstream::scatter(times_);
    Pstream::scatter(startTimeIndex_);
    Pstream::scatter(deltaT_);
    Pstream::scatter(nFace_);
}


void surfaceNoise::readSurfaceData
(
    const labelList& procFaceOffset,
    List<scalarField>& pData
)
{
    // Data is stored as pressure on surface at a given time.  Now we need to
    // pivot the data so that we have the complete pressure time history per
    // surface face.  In serial mode, this results in all pressure data being
    // loaded into memory (!)

    Info << "Reading pressure data" << endl;

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::nonBlocking);

        // Procedure:
        // 1. Master processor reads pressure data for all faces for all times
        // 2. Master sends each processor a subset of faces
        // 3. Local processor reconstructs the full time history for the subset
        //    of faces
        // Note: reading all data on master to avoid potential NFS problems...

        const label myProcI = Pstream::myProcNo();
        const label nLocalFace =
            procFaceOffset[myProcI + 1] - procFaceOffset[myProcI];

        // Complete pressure time history data for subset of faces
        pData.setSize(nLocalFace);
        const label nTimes = times_.size();
        forAll(pData, faceI)
        {
            pData[faceI].setSize(nTimes);
        }

        // Read and send data
        forAll(times_, i)
        {
            pBufs.clear();

            if (Pstream::master())
            {
                label timeI = i + startTimeIndex_;

                Info<< "    time: " << times_[i] << endl;

                // Read pressure at all faces for time timeI
                scalarField p(readerPtr_->field(timeI, pIndex_, scalar(0)));

                // Send subset of faces to each processor
                for (label procI = 0; procI < Pstream::nProcs(); procI++)
                {
                    label face0 = procFaceOffset[procI];
                    label nLocalFace =
                        procFaceOffset[procI + 1] - procFaceOffset[procI];

                    UOPstream toProc(procI, pBufs);
                    toProc << SubList<scalar>(p, nLocalFace, face0);
                }
            }

            pBufs.finishedSends();

            // Receive data from the master
            UIPstream fromProc(0, pBufs);

            scalarList pSlice(fromProc);

            forAll(pSlice, faceI)
            {
                pData[faceI][i] = pSlice[faceI];
            }
        }
    }
    else
    {
        const label nLocalFace = procFaceOffset[0];

        pData.setSize(nLocalFace);
        forAll(times_, timeI)
        {
            forAll(pData, faceI)
            {
                pData[faceI].setSize(times_.size());
            }
        }

        forAll(times_, i)
        {
            label timeI = i + startTimeIndex_;

            Info<< "    time: " << times_[i] << endl;
            const scalarField p(readerPtr_->field(timeI, pIndex_, scalar(0)));

            forAll(p, faceI)
            {
                pData[faceI][i] = p[faceI];
            }
        }
    }

    Info<< "Read "
        << returnReduce(pData.size(), sumOp<label>())
        << " pressure traces with "
        << times_.size()
        << " time values" << nl << endl;
}


void surfaceNoise::writeSurfaceData
(
    const word& fName,
    const word& groupName,
    const word& title,
    const scalar freq,
    const scalarField& data,
    const labelList& procFaceOffset
) const
{
    Info<< "    processing " << title << " for frequency " << freq << endl;

    fileName outDir
    (
        fileName("postProcessing")/"noise"/groupName/Foam::name(freq)
    );

    if (Pstream::parRun())
    {
        // Collect the surface data so that we can output the surfaces

        PstreamBuffers pBufs(Pstream::nonBlocking);

        if (!Pstream::master())
        {
            UOPstream toProc(0, pBufs);
            toProc << data;
        }

        pBufs.finishedSends();

        if (Pstream::master())
        {
            const meshedSurface& surf = readerPtr_->geometry();

            scalarField allData(surf.size());

            forAll(data, faceI)
            {
                // Master procFaceOffset is zero...
                allData[faceI] = data[faceI];
            }

            for (label procI = 1; procI < Pstream::nProcs(); procI++)
            {
                UIPstream fromProc(procI, pBufs);
                scalarList dataSlice(fromProc);
                forAll(dataSlice, i)
                {
                    label faceI = procFaceOffset[procI] + i;
                    allData[faceI] = dataSlice[i];
                }
            }

            fileName outFileName = writerPtr_->write
            (
                outDir,
                fName,
                surf.points(),
                surf.faces(),
                title,
                allData,
                false
            );
        }
    }
    else
    {
        const meshedSurface& surf = readerPtr_->geometry();

        writerPtr_->write
        (
            outDir,
            fName,
            surf.points(),
            surf.faces(),
            title,
            data,
            false
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceNoise::surfaceNoise(const dictionary& dict)
:
    noiseModel(dict),
    inputFileName_(dict.lookup("inputFile")),
    pIndex_(0),
    times_(),
    deltaT_(0),
    startTimeIndex_(0),
    nFace_(0),
    fftWriteInterval_(dict.lookupOrDefault("fftWriteInterval", 1))
{
    initialise(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

surfaceNoise::~surfaceNoise()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void surfaceNoise::calculate()
{
    // Container for pressure time history data per face
    List<scalarField> pData;

    // Processor procFaceOffsets
    labelList procFaceOffset;
    if (Pstream::parRun())
    {
        const label nProcs = Pstream::nProcs();
        const label nFacePerProc = floor(nFace_/nProcs) + 1;

        procFaceOffset.setSize(nProcs + 1, 0);
        for (label i = 1; i < procFaceOffset.size(); i++)
        {
            procFaceOffset[i] = min(i*nFacePerProc, nFace_);
        }
    }
    else
    {
        procFaceOffset.setSize(1, nFace_);
    }

    // Read pressure data from file
    readSurfaceData(procFaceOffset, pData);

    // Process the pressure data, and store results as surface values per
    // frequency so that it can be output using the surface writer

    Info<< "Creating noise FFTs" << endl;

    // Storage for FFT data
    const label nLocalFace = pData.size();
    const scalarField freq1(noiseFFT::frequencies(nSamples_, deltaT_));
    const label nFFT = freq1.size()/fftWriteInterval_;
    List<scalarField> surfPf(nFFT);
    List<scalarField> surfLf(nFFT);
    List<scalarField> surfPSD(nFFT);
    forAll(surfPf, freqI)
    {
        surfPf[freqI].setSize(nLocalFace);
        surfLf[freqI].setSize(nLocalFace);
        surfPSD[freqI].setSize(nLocalFace);
    }

    // Storage for 1/3 octave data
    labelList octave13BandIDs;
    noiseFFT::octaveFrequenciesIDs(freq1, fLower_, fUpper_, 3, octave13BandIDs);

    List<scalarField> surfPdelta(octave13BandIDs.size() - 1);
    List<scalarField> surfLdelta(octave13BandIDs.size() - 1);
    forAll(surfPdelta, freqI)
    {
        surfPdelta[freqI].setSize(nLocalFace);
        surfLdelta[freqI].setSize(nLocalFace);
    }

    forAll(pData, faceI)
    {
        const scalarField& p = pData[faceI];

        noiseFFT nfft(deltaT_, p);

        nfft -= pRef_;
        graph Pf(nfft.RMSmeanPf(windowModelPtr_()));
        graph Lf(nfft.Lf(Pf));
        graph PSDf(nfft.PSDf(windowModelPtr_()));
        graph PSD(nfft.PSD(PSDf));

        // Store the frequency results in slot for face of surface
        forAll(surfPf, i)
        {
            label freqI = (i + 1)*fftWriteInterval_ - 1;
            surfPf[i][faceI] = Pf.y()[freqI];
            surfLf[i][faceI] = Lf.y()[freqI];
            surfPSD[i][faceI] = PSD.y()[freqI];
        }

        graph Pdelta(nfft.Pdelta(Pf, octave13BandIDs));
        graph Ldelta(nfft.Ldelta(Lf, octave13BandIDs));

        // Store the 1/3 octave results in slot for face of surface
        forAll(surfPdelta, freqI)
        {
            surfPdelta[freqI][faceI] = Pdelta.y()[freqI];
            surfLdelta[freqI][faceI] = Ldelta.y()[freqI];
        }

        // Free the storage for p
//        p.clear();
    }

    Info<< "Writing fft surface data" << endl;

    forAll(surfPf, i)
    {
        label freqI = i*fftWriteInterval_;
        const word& fName = inputFileName_.name(true);
        const word gName = "fft";
        writeSurfaceData
        (
            fName,
            gName,
            "Pf",
            freq1[freqI],
            surfPf[i],
            procFaceOffset
        );
        writeSurfaceData
        (
            fName,
            gName,
            "Lf",
            freq1[freqI],
            surfLf[i],
            procFaceOffset
        );
        writeSurfaceData
        (
            fName,
            gName,
            "PSD",
            freq1[freqI],
            surfPSD[i],
            procFaceOffset
        );
    }

    Info<< "Writing one-third octave surface data" << endl;

    forAll(surfPdelta, i)
    {
        const word& fName = inputFileName_.name(true);
        const word gName = "oneThirdOctave";
        writeSurfaceData
        (
            fName,
            gName,
            "Pdelta",
            octave13BandIDs[i],
            surfPdelta[i],
            procFaceOffset
        );
        writeSurfaceData
        (
            fName,
            gName,
            "Ldelta",
            octave13BandIDs[i],
            surfLdelta[i],
            procFaceOffset
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// ************************************************************************* //



