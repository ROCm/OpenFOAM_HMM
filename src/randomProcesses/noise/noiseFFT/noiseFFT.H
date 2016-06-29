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

Class
    Foam::noiseFFT

Description
    Performs FFT of pressure field to generate noise data.

    General functionality:
    - Pf: fft of the pressure data
    - meanPf: multi-window mean fft
    - RMSmeanPf: multi-window RMS mean fft
    - PSDf: multi-window power spectral density (PSD) in frequency domain
    - PSD: power spectral density in dB/Hz
    - SPL: sound pressure level in dB

    Octave-based data:
    - PSD spectrum
    - SPL spectrum

SourceFiles
    noiseFFT.C

SeeAlso
    windowModel.H

\*---------------------------------------------------------------------------*/

#ifndef noiseFFT_H
#define noiseFFT_H

#include "scalarField.H"
#include "graph.H"
#include "windowModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class noiseFFT Declaration
\*---------------------------------------------------------------------------*/

class noiseFFT
:
    public scalarField
{
    // Private data

        //- Time spacing of the raw data
        scalar deltaT_;

    struct octaveBandInfo
    {
        label octave;

        // IDs of bin boundaries in pressure data
        labelList binIDs;

        // Centre frequencies for each bin
        scalarField centreFreq;
    };


public:

    //- Reference pressure
    static scalar p0;


    // Constructors

        //- Construct from pressure field
        noiseFFT
        (
            const scalar deltat,
            const scalarField& pressure
        );

        //- Construct from Istream
        noiseFFT(Istream&);

        //- Construct from pressure field file name
        noiseFFT(const fileName& pFileName, const label skip = 0);


    // Member Functions

        //- Return the FFT frequencies
        static tmp<scalarField> frequencies
        (
            const label N,
            const scalar deltaT
        );

        //- Return the PSD [dB/Hz]
        //  Input PSD in [Pa^2/Hz]
        static tmp<scalarField> PSD(const scalarField& PSDf);

        //- Return the SPL [dB]
        //  Input P(rms)^2 in [Pa^2]
        static tmp<scalarField> SPL(const scalarField& Prms2);

        //- Return a list of the frequency indices wrt f field that
        //  correspond to the bands limits for a given octave
        static void octaveBandInfo
        (
            const scalarField& f,
            const scalar fLower,
            const scalar fUpper,
            const scalar octave,
            labelList& fBandIDs,
            scalarField& fCentre
        );

        //- Return the graph of pressure as a function of time
        graph pt() const;

        //- Return the fft of the given pressure data
        tmp<scalarField> Pf(const tmp<scalarField>& pn) const;

        //- Return the multi-window mean fft of the complete pressure data [Pa]
        graph meanPf(const windowModel& window) const;

        //- Return the multi-window RMS mean fft of the complete pressure
        //  data [Pa]
        graph RMSmeanPf(const windowModel& window) const;

        //- Return the multi-window PSD (power spectral density) of the complete
        //  pressure data [Pa^2/Hz]
        graph PSDf(const windowModel& window) const;

        //- Return the PSD [dB/Hz]
        //  Takes PSD in [Pa^2/Hz]
        graph PSD(const graph& gPSDf) const;

        //- Generate octave data
        graph octaves
        (
            const graph& g,
            const labelList& freqBandIDs,
            bool integrate
        ) const;

        //- Convert the db into Pa
        scalar dbToPa(const scalar db) const;

        //- Convert the db-field into Pa
        tmp<scalarField> dbToPa(const tmp<scalarField>& db) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
