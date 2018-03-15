/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenCFD Ltd.
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

void surfaceNoise::initialise(const fileName& fName)
{
    Info<< "Reading data file " << fName << endl;

    label nAvailableTimes = 0;

    // All reading performed on the master processor only
    if (Pstream::master())
    {
        // Create the surface reader
        readerPtr_ = surfaceReader::New(readerType_, fName);

        // Find the index of the pressure data
        const List<word> fieldNames(readerPtr_->fieldNames(0));
        pIndex_ = fieldNames.find(pName_);
        if (pIndex_ == -1)
        {
            FatalErrorInFunction
                << "Unable to find pressure field name " << pName_
                << " in list of available fields: " << fieldNames
                << exit(FatalError);
        }

        // Create the surface writer
        // - Could be done later, but since this utility can process a lot of
        //   data we can ensure that the user-input is correct prior to doing
        //   the heavy lifting

        // Set the time range
        const instantList allTimes = readerPtr_->times();
        startTimeIndex_ = findStartTimeIndex(allTimes, startTime_);

        // Determine the windowing
        nAvailableTimes = allTimes.size() - startTimeIndex_;
    }

    Pstream::scatter(pIndex_);
    Pstream::scatter(startTimeIndex_);
    Pstream::scatter(nAvailableTimes);


    // Note: all processors should call the windowing validate function
    label nRequiredTimes = windowModelPtr_->validate(nAvailableTimes);

    if (Pstream::master())
    {
        // Restrict times
        const instantList allTimes = readerPtr_->times();

        times_.setSize(nRequiredTimes);
        forAll(times_, timeI)
        {
            times_[timeI] = allTimes[timeI + startTimeIndex_].value();
        }
        deltaT_ = checkUniformTimeStep(times_);

        // Read the surface geometry
        const meshedSurface& surf = readerPtr_->geometry();
        nFace_ = surf.size();
    }

    Pstream::scatter(times_);
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
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

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
        for (scalarField& pf : pData)
        {
            pf.setSize(nTimes);
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

                // Apply conversions
                p *= rhoRef_;

                // Send subset of faces to each processor
                for (label procI = 0; procI < Pstream::nProcs(); ++procI)
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
        for (scalarField& pf : pData)
        {
            pf.setSize(times_.size());
        }

        forAll(times_, i)
        {
            label timeI = i + startTimeIndex_;

            Info<< "    time: " << times_[i] << endl;
            const scalarField p(readerPtr_->field(timeI, pIndex_, scalar(0)));

            forAll(p, faceI)
            {
                pData[faceI][i] = p[faceI]*rhoRef_;
            }
        }
    }

    Info<< "Read "
        << returnReduce(pData.size(), sumOp<label>())
        << " pressure traces with "
        << times_.size()
        << " time values" << nl << endl;
}


Foam::scalar surfaceNoise::writeSurfaceData
(
    const fileName& outDirBase,
    const word& fName,
    const word& title,
    const scalar freq,
    const scalarField& data,
    const labelList& procFaceOffset,
    const bool writeSurface
) const
{
    Info<< "    processing " << title << " for frequency " << freq << endl;

    const fileName outDir(outDirBase/Foam::name(freq));

    if (Pstream::parRun())
    {
        // Collect the surface data so that we can output the surfaces

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (!Pstream::master())
        {
            UOPstream toProc(0, pBufs);
            toProc << data;
        }

        pBufs.finishedSends();

        scalar areaAverage = 0;
        if (Pstream::master())
        {
            const meshedSurface& surf = readerPtr_->geometry();

            scalarField allData(surf.size());

            forAll(data, faceI)
            {
                // Master procFaceOffset is zero...
                allData[faceI] = data[faceI];
            }

            for (label procI = 1; procI < Pstream::nProcs(); ++procI)
            {
                UIPstream fromProc(procI, pBufs);
                scalarList dataSlice(fromProc);
                forAll(dataSlice, i)
                {
                    label faceI = procFaceOffset[procI] + i;
                    allData[faceI] = dataSlice[i];
                }
            }

            // Could also have meshedSurface implement meshedSurf
            if (writeSurface)
            {
                fileName outFileName = writerPtr_->write
                (
                    outDir,
                    fName,
                    meshedSurfRef
                    (
                        surf.points(),
                        surf.surfFaces()
                    ),
                    title,
                    allData,
                    false
                );
            }

            // TO BE VERIFIED: area-averaged values
            // areaAverage = sum(allData*surf.magSf())/sum(surf.magSf());
            areaAverage = sum(allData)/allData.size();
        }
        Pstream::scatter(areaAverage);

        return areaAverage;
    }
    else
    {
        const meshedSurface& surf = readerPtr_->geometry();

        // Could also have meshedSurface implement meshedSurf
        if (writeSurface)
        {
            writerPtr_->write
            (
                outDir,
                fName,
                meshedSurfRef
                (
                    surf.points(),
                    surf.surfFaces()
                ),
                title,
                data,
                false
            );
        }

        // TO BE VERIFIED: area-averaged values
        // return sum(data*surf.magSf())/sum(surf.magSf());
        return sum(data)/data.size();
    }
}


Foam::scalar surfaceNoise::surfaceAverage
(
    const scalarField& data,
    const labelList& procFaceOffset
) const
{
    if (Pstream::parRun())
    {
        // Collect the surface data so that we can output the surfaces

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (!Pstream::master())
        {
            UOPstream toProc(0, pBufs);
            toProc << data;
        }

        pBufs.finishedSends();

        scalar areaAverage = 0;
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

            // TO BE VERIFIED: area-averaged values
            // areaAverage = sum(allData*surf.magSf())/sum(surf.magSf());
            areaAverage = sum(allData)/allData.size();
        }
        Pstream::scatter(areaAverage);

        return areaAverage;
    }
    else
    {

        // TO BE VERIFIED: area-averaged values
        // const meshedSurface& surf = readerPtr_->geometry();
        // return sum(data*surf.magSf())/sum(surf.magSf());
        return sum(data)/data.size();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceNoise::surfaceNoise(const dictionary& dict, const bool readFields)
:
    noiseModel(dict, false),
    inputFileNames_(),
    pName_("p"),
    pIndex_(0),
    times_(),
    deltaT_(0),
    startTimeIndex_(0),
    nFace_(0),
    fftWriteInterval_(1),
    readerType_(word::null),
    readerPtr_(nullptr),
    writerPtr_(nullptr)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

surfaceNoise::~surfaceNoise()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool surfaceNoise::read(const dictionary& dict)
{
    if (noiseModel::read(dict))
    {
        if (dict.found("file"))
        {
            inputFileNames_.setSize(1);
            dict.lookup("file") >> inputFileNames_[0];
        }
        else
        {
            dict.lookup("files") >> inputFileNames_;
        }

        dict.readIfPresent("fftWriteInterval", fftWriteInterval_);
        dict.readIfPresent("p", pName_);

        dict.lookup("reader") >> readerType_;

        word writerType(dict.lookup("writer"));
        dictionary optDict
        (
            dict.subOrEmptyDict("writeOptions").subOrEmptyDict(writerType)
        );

        writerPtr_ = surfaceWriter::New(writerType, optDict);

        return true;
    }

    return false;
}


void surfaceNoise::calculate()
{
    forAll(inputFileNames_, filei)
    {
        fileName fName = inputFileNames_[filei];
        fName.expand();

        if (!fName.isAbsolute())
        {
            fName = "$FOAM_CASE"/fName;
        }

        initialise(fName.expand());

        // Container for pressure time history data per face
        List<scalarField> pData;

        // Processor procFaceOffsets
        labelList procFaceOffset;
        if (Pstream::parRun())
        {
            const label nProcs = Pstream::nProcs();
            const label nFacePerProc = floor(nFace_/nProcs) + 1;

            procFaceOffset.setSize(nProcs + 1, 0);
            for (label i = 1; i < procFaceOffset.size(); ++i)
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

        const scalarField freq1(noiseFFT::frequencies(nSamples_, deltaT_));

        // Reset desired frequency range if outside actual frequency range
        fLower_ = min(fLower_, max(freq1));
        fUpper_ = min(fUpper_, max(freq1));

        // Storage for FFT data
        const label nLocalFace = pData.size();
        const label nFFT = ceil(freq1.size()/scalar(fftWriteInterval_));

        List<scalarField> surfPrmsf(nFFT);
        List<scalarField> surfPSDf(nFFT);
        forAll(surfPrmsf, freqI)
        {
            surfPrmsf[freqI].setSize(nLocalFace);
            surfPSDf[freqI].setSize(nLocalFace);
        }

        // Storage for 1/3 octave data
        labelList octave13BandIDs;
        scalarField octave13FreqCentre;
        noiseFFT::octaveBandInfo
        (
            freq1,
            fLower_,
            fUpper_,
            3,
            octave13BandIDs,
            octave13FreqCentre
        );

        label bandSize = 0;
        if (octave13BandIDs.empty())
        {
            WarningInFunction
                << "Octave band calculation failed (zero sized). "
                << "Please check your input data"
                << endl;
        }
        else
        {
            bandSize = octave13BandIDs.size() - 1;
        }

        List<scalarField> surfPSD13f(bandSize);
        List<scalarField> surfPrms13f2(bandSize);
        forAll(surfPSD13f, freqI)
        {
            surfPSD13f[freqI].setSize(nLocalFace);
            surfPrms13f2[freqI].setSize(nLocalFace);
        }

        const windowModel& win = windowModelPtr_();

        {
            noiseFFT nfft(deltaT_, win.nSamples());

            forAll(pData, faceI)
            {
                // Note: passes storage from pData to noiseFFT
                nfft.setData(pData[faceI]);

                // Generate the FFT-based data
                graph Prmsf(nfft.RMSmeanPf(win));
                graph PSDf(nfft.PSDf(win));

                // Store the frequency results in slot for face of surface
                forAll(surfPrmsf, i)
                {
                    label freqI = i*fftWriteInterval_;
                    surfPrmsf[i][faceI] = Prmsf.y()[freqI];
                    surfPSDf[i][faceI] = PSDf.y()[freqI];
                }

                // PSD [Pa^2/Hz]
                graph PSD13f(nfft.octaves(PSDf, octave13BandIDs, false));

                // Integrated PSD = P(rms)^2 [Pa^2]
                graph Prms13f2(nfft.octaves(PSDf, octave13BandIDs, true));

                // Store the 1/3 octave results in slot for face of surface
                forAll(surfPSD13f, freqI)
                {
                    surfPSD13f[freqI][faceI] = PSD13f.y()[freqI];
                    surfPrms13f2[freqI][faceI] = Prms13f2.y()[freqI];
                }
            }
        }

        const word fNameBase = fName.nameLessExt();

        // Output directory for graphs
        fileName outDirBase(baseFileDir(filei)/fNameBase);

        const scalar deltaf = 1.0/(deltaT_*win.nSamples());
        Info<< "Writing fft surface data";
        if (fftWriteInterval_ == 1)
        {
            Info<< endl;
        }
        else
        {
            Info<< " at every " << fftWriteInterval_ << " frequency points"
                << endl;
        }

        {
            fileName outDir(outDirBase/"fft");

            // Determine frequency range of interest
            // Note: frequencies have fixed interval, and are in the range
            //       0 to fftWriteInterval_*(n-1)*deltaf
            label f0 = ceil(fLower_/deltaf/scalar(fftWriteInterval_));
            label f1 = floor(fUpper_/deltaf/scalar(fftWriteInterval_));
            label nFreq = f1 - f0;

            scalarField PrmsfAve(nFreq, 0);
            scalarField PSDfAve(nFreq, 0);
            scalarField fOut(nFreq, 0);

            if (nFreq == 0)
            {
                WarningInFunction
                    << "No surface data available using a fftWriteInterval of "
                    << fftWriteInterval_ << endl;
            }
            else
            {
                forAll(fOut, i)
                {
                    label freqI = (i + f0)*fftWriteInterval_;
                    fOut[i] = freq1[freqI];


                    PrmsfAve[i] = writeSurfaceData
                    (
                        outDir,
                        fNameBase,
                        "Prmsf",
                        freq1[freqI],
                        surfPrmsf[i + f0],
                        procFaceOffset,
                        writePrmsf_
                    );

                    PSDfAve[i] = writeSurfaceData
                    (
                        outDir,
                        fNameBase,
                        "PSDf",
                        freq1[freqI],
                        surfPSDf[i + f0],
                        procFaceOffset,
                        writePSDf_
                    );
                    writeSurfaceData
                    (
                        outDir,
                        fNameBase,
                        "PSD",
                        freq1[freqI],
                        noiseFFT::PSD(surfPSDf[i + f0]),
                        procFaceOffset,
                        writePSD_
                    );
                    writeSurfaceData
                    (
                        outDir,
                        fNameBase,
                        "SPL",
                        freq1[freqI],
                        noiseFFT::SPL(surfPSDf[i + f0]*deltaf),
                        procFaceOffset,
                        writeSPL_
                    );
                }
            }

            graph Prmsfg
            (
                "Average Prms(f)",
                "f [Hz]",
                "P(f) [Pa]",
                fOut,
                PrmsfAve
            );
            Prmsfg.write(outDir, graph::wordify(Prmsfg.title()), graphFormat_);

            graph PSDfg
            (
                "Average PSD_f(f)",
                "f [Hz]",
                "PSD(f) [PaPa_Hz]",
                fOut,
                PSDfAve
            );
            PSDfg.write(outDir, graph::wordify(PSDfg.title()), graphFormat_);

            graph PSDg
            (
                "Average PSD_dB_Hz(f)",
                "f [Hz]",
                "PSD(f) [dB_Hz]",
                fOut,
                noiseFFT::PSD(PSDfAve)
            );
            PSDg.write(outDir, graph::wordify(PSDg.title()), graphFormat_);

            graph SPLg
            (
                "Average SPL_dB(f)",
                "f [Hz]",
                "SPL(f) [dB]",
                fOut,
                noiseFFT::SPL(PSDfAve*deltaf)
            );
            SPLg.write(outDir, graph::wordify(SPLg.title()), graphFormat_);
        }


        Info<< "Writing one-third octave surface data" << endl;
        {
            fileName outDir(outDirBase/"oneThirdOctave");

            scalarField PSDfAve(surfPSD13f.size(), 0);
            scalarField Prms13f2Ave(surfPSD13f.size(), 0);

            forAll(surfPSD13f, i)
            {
                PSDfAve[i] = writeSurfaceData
                (
                    outDir,
                    fNameBase,
                    "PSD13f",
                    octave13FreqCentre[i],
                    surfPSD13f[i],
                    procFaceOffset,
                    writeOctaves_
                );
                writeSurfaceData
                (
                    outDir,
                    fNameBase,
                    "PSD13",
                    octave13FreqCentre[i],
                    noiseFFT::PSD(surfPSD13f[i]),
                    procFaceOffset,
                    writeOctaves_
                );
                writeSurfaceData
                (
                    outDir,
                    fNameBase,
                    "SPL13",
                    octave13FreqCentre[i],
                    noiseFFT::SPL(surfPrms13f2[i]),
                    procFaceOffset,
                    writeOctaves_
                );

                Prms13f2Ave[i] =
                    surfaceAverage(surfPrms13f2[i], procFaceOffset);
            }

            graph PSD13g
            (
                "Average PSD13_dB_Hz(fm)",
                "fm [Hz]",
                "PSD(fm) [dB_Hz]",
                octave13FreqCentre,
                noiseFFT::PSD(PSDfAve)
            );
            PSD13g.write(outDir, graph::wordify(PSD13g.title()), graphFormat_);

            graph SPL13g
            (
                "Average SPL13_dB(fm)",
                "fm [Hz]",
                "SPL(fm) [dB]",
                octave13FreqCentre,
                noiseFFT::SPL(Prms13f2Ave)
            );
            SPL13g.write(outDir, graph::wordify(SPL13g.title()), graphFormat_);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// ************************************************************************* //
