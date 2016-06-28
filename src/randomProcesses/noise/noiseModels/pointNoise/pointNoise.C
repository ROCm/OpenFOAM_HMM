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
    const CSV<scalar>& pData,
    scalarField& t,
    scalarField& p
)
{
    const scalarField t0(pData.x());
    const scalarField p0(pData.y());

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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointNoise::calculate()
{
    // Point data only handled by master
    if (!Pstream::master())
    {
        return;
    }

    Info<< "Reading data file" << endl;

    CSV<scalar> pData("pressure", dict_, "Data");
    fileName baseFileName(pData.fName().lessExt());

    // Time and pressure history data
    scalarField t, p;
    filterTimeData(pData, t, p);
    Info<< "    read " << t.size() << " values" << nl << endl;

    Info<< "Creating noise FFT" << endl;
    const scalar deltaT = checkUniformTimeStep(t);

    // Determine the windowing
    windowModelPtr_->validate(t.size());

    // Create the fft
    noiseFFT nfft(deltaT, p);

    nfft -= pRef_;

    graph Pf(nfft.RMSmeanPf(windowModelPtr_()));
    Info<< "    Creating graph for " << Pf.title() << endl;
    Pf.write(baseFileName + graph::wordify(Pf.title()), graphFormat_);

    graph Lf(nfft.Lf(Pf));
    Info<< "    Creating graph for " << Lf.title() << endl;
    Lf.write(baseFileName + graph::wordify(Lf.title()), graphFormat_);

    graph PSDf(nfft.PSDf(windowModelPtr_()));
    Info<< "    Creating graph for " << PSDf.title() << endl;
    PSDf.write(baseFileName + graph::wordify(PSDf.title()), graphFormat_);

    graph PSD(nfft.PSD(PSDf));
    Info<< "    Creating graph for " << PSD.title() << endl;
    PSD.write(baseFileName + graph::wordify(PSD.title()), graphFormat_);

    labelList octave13BandIDs;
    noiseFFT::octaveFrequenciesIDs(Pf.x(), fLower_, fUpper_, 3, octave13BandIDs);

    graph Ldelta(nfft.Ldelta(Lf, octave13BandIDs));
    Info<< "    Creating graph for " << Ldelta.title() << endl;
    Ldelta.write(baseFileName + graph::wordify(Ldelta.title()), graphFormat_);

    graph Pdelta(nfft.Pdelta(Pf, octave13BandIDs));
    Info<< "    Creating graph for " << Pdelta.title() << endl;
    Pdelta.write(baseFileName + graph::wordify(Pdelta.title()), graphFormat_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pointNoise::pointNoise(const dictionary& dict)
:
    noiseModel(dict),
    graphFormat_(dict.lookupOrDefault<word>("graphFormat", "raw"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pointNoise::~pointNoise()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// ************************************************************************* //
