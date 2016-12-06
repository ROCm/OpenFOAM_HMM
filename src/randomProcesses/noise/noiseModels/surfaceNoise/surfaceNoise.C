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
    label nAvailableTimes = 0;

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

                // Apply conversions
                p *= rhoRef_;

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

            // TODO: Move faceAreas to demand-driven function in MeshedSurface
            // scalarField faceAreas(surf.faces().size());
            // forAll(faceAreas, i)
            // {
            //     faceAreas[i] = surf.faces()[i].mag(surf.points());
            // }
            //
            // areaAverage = sum(allData*faceAreas)/sum(faceAreas);
            areaAverage = sum(allData)/allData.size();
        }
        Pstream::scatter(areaAverage);

        return areaAverage;
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

        // TODO: Move faceAreas to demand-driven function in MeshedSurface
        // scalarField faceAreas(surf.faces().size());
        // forAll(faceAreas, i)
        // {
        //     faceAreas[i] = surf.faces()[i].mag(surf.points());
        // }
        //
        // return sum(data*faceAreas)/sum(faceAreas);
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

        PstreamBuffers pBufs(Pstream::nonBlocking);

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

            // TODO: Move faceAreas to demand-driven function in MeshedSurface
            scalarField faceAreas(surf.faces().size());
            forAll(faceAreas, i)
            {
                faceAreas[i] = surf.faces()[i].mag(surf.points());
            }

//            areaAverage = sum(allData*faceAreas)/sum(faceAreas);
            areaAverage = sum(allData)/allData.size();
        }
        Pstream::scatter(areaAverage);

        return areaAverage;
    }
    else
    {
        const meshedSurface& surf = readerPtr_->geometry();

        // TODO: Move faceAreas to demand-driven function in MeshedSurface
        scalarField faceAreas(surf.faces().size());
        forAll(faceAreas, i)
        {
            faceAreas[i] = surf.faces()[i].mag(surf.points());
        }

//        return sum(data*faceAreas)/sum(faceAreas);
        return sum(data)/data.size();
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
            << "Ocatve band calculation failed (zero sized). "
            << "please check your input data"
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

    forAll(pData, faceI)
    {
        const scalarField& p = pData[faceI];

        noiseFFT nfft(deltaT_, p);
        graph Prmsf(nfft.RMSmeanPf(win));
        graph PSDf(nfft.PSDf(win));

        // Store the frequency results in slot for face of surface
        forAll(surfPrmsf, i)
        {
            label freqI = (i + 1)*fftWriteInterval_ - 1;
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

        // Free the storage for p
//        p.clear();
    }

    // Output directory for graphs
    fileName outDir(fileName("postProcessing")/"noise"/typeName);

    const scalar deltaf = 1.0/(deltaT_*win.nSamples());
    Info<< "Writing fft surface data" << endl;
    {
        scalarField PrmsfAve(surfPrmsf.size(), 0);
        scalarField PSDfAve(surfPrmsf.size(), 0);
        scalarField fOut(surfPrmsf.size(), 0);

        forAll(surfPrmsf, i)
        {
            label freqI = (i + 1)*fftWriteInterval_ - 1;
            fOut[i] = freq1[freqI];
            const word& fName = inputFileName_.name(true);
            const word gName = "fft";
            PrmsfAve[i] = writeSurfaceData
            (
                fName,
                gName,
                "Prmsf",
                freq1[freqI],
                surfPrmsf[i],
                procFaceOffset
            );

            PSDfAve[i] = writeSurfaceData
            (
                fName,
                gName,
                "PSDf",
                freq1[freqI],
                surfPSDf[i],
                procFaceOffset
            );
            writeSurfaceData
            (
                fName,
                gName,
                "PSD",
                freq1[freqI],
                noiseFFT::PSD(surfPSDf[i]),
                procFaceOffset
            );
            writeSurfaceData
            (
                fName,
                gName,
                "SPL",
                freq1[freqI],
                noiseFFT::SPL(surfPSDf[i]*deltaf),
                procFaceOffset
            );
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
        scalarField PSDfAve(surfPSD13f.size(), 0);
        scalarField Prms13f2Ave(surfPSD13f.size(), 0);

        forAll(surfPSD13f, i)
        {
            const word& fName = inputFileName_.name(true);
            const word gName = "oneThirdOctave";
            PSDfAve[i] = writeSurfaceData
            (
                fName,
                gName,
                "PSD13f",
                octave13FreqCentre[i],
                surfPSD13f[i],
                procFaceOffset
            );
            writeSurfaceData
            (
                fName,
                gName,
                "PSD13",
                octave13FreqCentre[i],
                noiseFFT::PSD(surfPSD13f[i]),
                procFaceOffset
            );
            writeSurfaceData
            (
                fName,
                gName,
                "SPL13",
                octave13FreqCentre[i],
                noiseFFT::SPL(surfPrms13f2[i]),
                procFaceOffset
            );

            Prms13f2Ave[i] = surfaceAverage(surfPrms13f2[i], procFaceOffset);
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// ************************************************************************* //
