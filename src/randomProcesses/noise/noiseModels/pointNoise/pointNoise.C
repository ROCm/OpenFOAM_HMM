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

#include "pointNoise.H"
#include "noiseFFT.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace noiseModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointNoise, 0);
addToRunTimeSelectionTable(noiseModel, pointNoise, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void pointNoise::filterTimeData
(
    const scalarField& t0,
    const scalarField& p0,
    scalarField& t,
    scalarField& p
) const
{
    DynamicList<scalar> tf(t0.size());
    DynamicList<scalar> pf(t0.size());

    forAll(t0, timeI)
    {
        if (t0[timeI] >= startTime_)
        {
            tf.append(t0[timeI]);
            pf.append(p0[timeI]);
        }
    }

    t.transfer(tf);
    p.transfer(pf);
}


void pointNoise::processData
(
    const label dataseti,
    const Function1Types::CSV<scalar>& data
)
{
    Info<< "Reading data file " << data.fName() << endl;

    const word fNameBase = data.fName().nameLessExt();

    // Time and pressure history data
    scalarField t, p;
    filterTimeData(data.x(), data.y(), t, p);
    p *= rhoRef_;

    Info<< "    read " << t.size() << " values" << nl << endl;

    Info<< "Creating noise FFT" << endl;

    const scalar deltaT = checkUniformTimeStep(t);

    if (!validateBounds(p))
    {
        Info<< "No noise data generated" << endl;
        return;
    }


    // Determine the windowing
    windowModelPtr_->validate(t.size());
    const windowModel& win = windowModelPtr_();
    const scalar deltaf = 1.0/(deltaT*win.nSamples());
    fileName outDir(baseFileDir(dataseti)/fNameBase);

    // Create the fft
    noiseFFT nfft(deltaT, win.nSamples());
    nfft.setData(p);


    // Narrow band data
    // ----------------

    // RMS pressure [Pa]
    graph Prmsf(nfft.RMSmeanPf(win));
    if (customBounds_)
    {
        Prmsf.setXRange(fLower_, fUpper_);
    }
    if (writePrmsf_)
    {
        Info<< "    Creating graph for " << Prmsf.title() << endl;
        Prmsf.write(outDir, graph::wordify(Prmsf.title()), graphFormat_);
    }

    // PSD [Pa^2/Hz]
    graph PSDf(nfft.PSDf(win));
    if (customBounds_)
    {
        PSDf.setXRange(fLower_, fUpper_);
    }
    if (writePSDf_)
    {
        Info<< "    Creating graph for " << PSDf.title() << endl;
        PSDf.write(outDir, graph::wordify(PSDf.title()), graphFormat_);
    }

    // PSD [dB/Hz]
    graph PSDg
    (
        "PSD_dB_Hz(f)",
        "f [Hz]",
        "PSD(f) [dB_Hz]",
        Prmsf.x(),
        noiseFFT::PSD(PSDf.y())
    );

    if (writePSD_)
    {
        Info<< "    Creating graph for " << PSDg.title() << endl;
        PSDg.write(outDir, graph::wordify(PSDg.title()), graphFormat_);
    }

    // SPL [dB]
    graph SPLg
    (
        "SPL_dB(f)",
        "f [Hz]",
        "SPL(f) [dB]",
        Prmsf.x(),
        noiseFFT::SPL(PSDf.y()*deltaf)
    );

    if (writeSPL_)
    {
        Info<< "    Creating graph for " << SPLg.title() << endl;
        SPLg.write(outDir, graph::wordify(SPLg.title()), graphFormat_);
    }

    if (writeOctaves_)
    {
        labelList octave13BandIDs;
        scalarField octave13FreqCentre;
        noiseFFT::octaveBandInfo
        (
            Prmsf.x(),
            fLower_,
            fUpper_,
            3,
            octave13BandIDs,
            octave13FreqCentre
        );


        // 1/3 octave data
        // ---------------

        // PSD [Pa^2/Hz]
        graph PSD13f(nfft.octaves(PSDf, octave13BandIDs, false));

        // Integrated PSD = P(rms)^2 [Pa^2]
        graph Prms13f2(nfft.octaves(PSDf, octave13BandIDs, true));

        graph PSD13g
        (
            "PSD13_dB_Hz(fm)",
            "fm [Hz]",
            "PSD(fm) [dB_Hz]",
            octave13FreqCentre,
            noiseFFT::PSD(PSD13f.y())
        );
        Info<< "    Creating graph for " << PSD13g.title() << endl;
        PSD13g.write(outDir, graph::wordify(PSD13g.title()), graphFormat_);

        graph SPL13g
        (
            "SPL13_dB(fm)",
            "fm [Hz]",
            "SPL(fm) [dB]",
            octave13FreqCentre,
            noiseFFT::SPL(Prms13f2.y())
        );
        Info<< "    Creating graph for " << SPL13g.title() << endl;
        SPL13g.write(outDir, graph::wordify(SPL13g.title()), graphFormat_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pointNoise::pointNoise(const dictionary& dict, const bool readFields)
:
    noiseModel(dict, false)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pointNoise::~pointNoise()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointNoise::calculate()
{
    // Point data only handled by master
    if (!Pstream::master())
    {
        return;
    }


    forAll(inputFileNames_, filei)
    {
        fileName fName = inputFileNames_[filei];
        fName.expand();
        if (!fName.isAbsolute())
        {
            fName = "$FOAM_CASE"/fName;
            fName.expand();
        }
        Function1Types::CSV<scalar> data("pressure", dict_, fName);
        processData(filei, data);
    }
}


bool pointNoise::read(const dictionary& dict)
{
    if (noiseModel::read(dict))
    {
        if (!dict.readIfPresent("files", inputFileNames_))
        {
            inputFileNames_.setSize(1);

            // Note: lookup uses same keyword as used by the CSV constructor
            dict.lookup("file") >> inputFileNames_[0];
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// ************************************************************************* //
