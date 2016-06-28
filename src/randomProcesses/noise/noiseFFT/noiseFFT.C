/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "noiseFFT.H"
#include "IFstream.H"
#include "DynamicList.H"
#include "fft.H"
#include "SubField.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * //

Foam::scalar Foam::noiseFFT::p0 = 2e-5;


Foam::tmp<Foam::scalarField> Foam::noiseFFT::frequencies
(
    const label N,
    const scalar deltaT
)
{
    tmp<scalarField> tf(new scalarField(N/2, 0));
    scalarField& f = tf();

    scalar deltaf = 1.0/(N*deltaT);
    forAll(f, i)
    {
        f[i] = i*deltaf;
    }

    return tf;
}


void Foam::noiseFFT::octaveFrequenciesIDs
(
    const scalarField& f,
    const scalar fLower,
    const scalar fUpper,
    const scalar octave,
    labelList& freqBandIDs
)
{
    // Set the indices of to the lower frequency bands for the input frequency
    // range. Ensure that the centre frequency passes though 1000 Hz

    // Low frequency bound given by:
    //     fLow = f0*(2^(bandI/octave/2))
    // Centre frequency given by:
    //     fCentre = f0*(2^(bandI/octave))

    scalar f0 = 1000;
    scalar minFrequency = max(fLower, min(f));

    // Lower frequency band limit
    label band0Low = ceil(2*octave*log(minFrequency/f0)/log(2.0));

    // Centre frequency band limit
    //label band0Centre = ceil(octave*log(fLower/f0)/log(2.0));

    scalar fLowerBand = f0*pow(2, band0Low/octave/2);
    scalar fRatio = pow(2, 1.0/octave);

    bool complete = false;
    DynamicList<label> bandIDs(f.size());
    forAll(f, i)
    {
        while (f[i] >= fLowerBand)
        {
            bandIDs.append(i);
            fLowerBand *= fRatio;

            if (fLowerBand > fUpper)
            {
                complete = true;
                break;
            }
        }

        if (complete) break;
    }

    freqBandIDs.transfer(bandIDs);
}


Foam::tmp<Foam::scalarField> Foam::noiseFFT::octaveFrequencies
(
    const scalarField& f,
    const scalar fLower,
    const scalar fUpper,
    const scalar octave
)
{
    labelList freqBandIDs;
    octaveFrequenciesIDs(f, fLower, fUpper, octave, freqBandIDs);
    tmp<scalarField> tf(new scalarField(f, freqBandIDs));
    return tf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noiseFFT::noiseFFT
(
    const scalar deltaT,
    const scalarField& pressure
)
:
    scalarField(pressure),
    deltaT_(deltaT)
{}


Foam::noiseFFT::noiseFFT(const fileName& pFileName, const label skip)
:
    scalarField(),
    deltaT_(0.0)
{
    // Construct pressure data file
    IFstream pFile(pFileName);

    // Check pFile stream is OK
    if (!pFile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << pFileName
            << exit(FatalError);
    }

    if (skip)
    {
        scalar dummyt, dummyp;

        for (label i = 0; i < skip; i++)
        {
            pFile >> dummyt;

            if (!pFile.good() || pFile.eof())
            {
                FatalErrorInFunction
                    << "Number of points in file " << pFileName
                    << " is less than the number to be skipped = " << skip
                    << exit(FatalError);
            }

            pFile >> dummyp;
        }
    }

    scalar t = 0, T0 = 0, T1 = 0;
    DynamicList<scalar> pData(100000);
    label i = 0;

    while (!(pFile >> t).eof())
    {
        if (i == 0)
        {
            T0 = t;
        }

        T1 = t;
        pFile >> pData(i++);
    }

    // Note: assumes fixed time step
    deltaT_ = (T1 - T0)/pData.size();

    this->transfer(pData);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::graph Foam::noiseFFT::pt() const
{
    scalarField t(size());
    forAll(t, i)
    {
        t[i] = i*deltaT_;
    }

    return graph
    (
        "p(t)",
        "t [s]",
        "p(t) [Pa]",
        t,
        *this
    );
}


Foam::tmp<Foam::scalarField> Foam::noiseFFT::Pf
(
    const tmp<scalarField>& tpn
) const
{
    tmp<scalarField> tPn2
    (
        mag
        (
            fft::reverseTransform
            (
                ReComplexField(tpn),
                labelList(1, tpn().size())
            )
        )
    );

    tpn.clear();

    tmp<scalarField> tPn
    (
        new scalarField
        (
            scalarField::subField(tPn2(), tPn2().size()/2)
        )
    );
    scalarField& Pn = tPn();

    Pn *= 2.0/sqrt(scalar(tPn2().size()));
    Pn[0] /= 2.0;

    return tPn;
}


Foam::graph Foam::noiseFFT::meanPf(const windowModel& window) const
{
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    scalarField meanPf(N/2, 0.0);

    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        meanPf += Pf(window.apply<scalar>(*this, windowI));
    }

    meanPf /= scalar(nWindow);

    scalar deltaf = 1.0/(N*deltaT_);
    scalarField f(meanPf.size());
    forAll(f, i)
    {
        f[i] = i*deltaf;
    }

    return graph
    (
        "P(f)",
        "f [Hz]",
        "P(f) [Pa]",
        f,
        meanPf
    );
}


Foam::graph Foam::noiseFFT::RMSmeanPf(const windowModel& window) const
{
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    scalarField RMSMeanPf(N/2, 0.0);

    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        RMSMeanPf += sqr(Pf(window.apply<scalar>(*this, windowI)));
    }

    RMSMeanPf = sqrt(RMSMeanPf/scalar(nWindow));

    scalar deltaf = 1.0/(N*deltaT_);
    scalarField f(RMSMeanPf.size());
    forAll(f, i)
    {
        f[i] = i*deltaf;
    }

    return graph
    (
        "P(f)",
        "f [Hz]",
        "P(f) [Pa]",
        f,
        RMSMeanPf
    );
}


Foam::graph Foam::noiseFFT::PSDf(const windowModel& window) const
{
    const label N = window.nSamples();
    const label nWindow = window.nWindow();

    scalarField psd(N/2, 0.0);

    for (label windowI = 0; windowI < nWindow; ++windowI)
    {
        psd += 0.5*sqr(Pf(window.apply<scalar>(*this, windowI)));
    }

    scalar deltaf = 1.0/(N*deltaT_);

    psd /= nWindow*deltaf;

    scalarField f(psd.size());

    forAll(f, i)
    {
        f[i] = i*deltaf;
    }

    return graph
    (
        "PSD(f)",
        "f [Hz]",
        "PSD(f) [PaPa_Hz]",
        f,
        psd
    );
}


Foam::graph Foam::noiseFFT::PSD(const graph& gPSD) const
{
    return graph
    (
        "PSD(dB)",
        "f [Hz]",
        "PSD_dB(f) [dB]",
        gPSD.x(),
        10*log10(gPSD.y()/sqr(p0))
    );
}


Foam::graph Foam::noiseFFT::Lf(const graph& gPf) const
{
    return graph
    (
        "L(f)",
        "f [Hz]",
        "L(f) [dB]",
        gPf.x(),
        20*log10(gPf.y()/p0)
    );
}


Foam::graph Foam::noiseFFT::Ldelta
(
    const graph& gLf,
    const labelList& freqBandIDs
) const
{
    if (freqBandIDs.size() < 2)
    {
        WarningInFunction
            << "Octave frequency bands are not defined "
            << "- skipping Ldelta calculation"
            << endl;

        return graph
        (
            "Ldelta",
            "fm [Hz]",
            "Ldelta [dB]",
            scalarField(),
            scalarField()
        );
    }

    const scalarField& f = gLf.x();
    const scalarField& Lf = gLf.y();

    scalarField ldelta(freqBandIDs.size() - 1, 0.0);
    scalarField fm(freqBandIDs.size() -1, 0.0);

    for (label bandI = 0; bandI < freqBandIDs.size() - 1; bandI++)
    {
        label f0 = freqBandIDs[bandI];
        label f1 = freqBandIDs[bandI+1];
        fm[bandI] = f[f0];

        if (f0 == f1) continue;

        for (label freqI = f0; freqI < f1; freqI++)
        {
            ldelta[bandI] += pow(10, Lf[freqI]/10.0);
        }
        ldelta[bandI] = 10*log10(ldelta[bandI]);
    }

    return graph
    (
        "Ldelta",
        "fm [Hz]",
        "Ldelta [dB]",
        fm,
        ldelta
    );
}


Foam::graph Foam::noiseFFT::Pdelta
(
    const graph& gPf,
    const labelList& freqBandIDs
) const
{
    if (freqBandIDs.size() < 2)
    {
        WarningInFunction
            << "Octave frequency bands are not defined "
            << "- skipping Pdelta calculation"
            << endl;

        return graph
        (
            "Pdelta",
            "fm [Hz]",
            "Pdelta [dB]",
            scalarField(),
            scalarField()
        );
    }

    const scalarField& f = gPf.x();
    const scalarField& Pf = gPf.y();

    scalarField pdelta(freqBandIDs.size() - 1, 0.0);
    scalarField fm(pdelta.size());

    for (label bandI = 0; bandI < freqBandIDs.size() - 1; bandI++)
    {
        label f0 = freqBandIDs[bandI];
        label f1 = freqBandIDs[bandI+1];
        fm[bandI] = f[f0];

        if (f0 == f1) continue;

        for (label freqI = f0; freqI < f1; freqI++)
        {
            pdelta[bandI] += sqr(Pf[freqI]);
        }
        pdelta[bandI] = sqrt((2.0/3.0)*pdelta[bandI]);
    }

    return graph
    (
        "Pdelta",
        "fm [Hz]",
        "Pdelta [dB]",
        fm,
        pdelta
    );
}


Foam::scalar Foam::noiseFFT::Lsum(const graph& gLf) const
{
    const scalarField& Lf = gLf.y();

    scalar lsum = 0.0;

    forAll(Lf, i)
    {
        lsum += pow(10, Lf[i]/10.0);
    }

    lsum = 10*log10(lsum);

    return lsum;
}


Foam::scalar Foam::noiseFFT::dbToPa(const scalar db) const
{
    return p0*pow(10.0, db/20.0);
}


Foam::tmp<Foam::scalarField> Foam::noiseFFT::dbToPa
(
    const tmp<scalarField>& db
) const
{
    return p0*pow(10.0, db/20.0);
}


// ************************************************************************* //
