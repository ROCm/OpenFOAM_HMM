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
    - PSD: multi-window power spectral density (PSD) in dB
    - Lf: narrow-band pressure-fluctuation level (PFL) in frequency domain

    Octave-based data:
    - Ldelta: PFL spectrum
    - Pdelta: pressure spectrum

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

        //- Return a list of the frequency indices wrt f field that
        //  correspond to the bands limits for a given octave
        static void octaveFrequenciesIDs
        (
            const scalarField& f,
            const scalar fLower,
            const scalar fUpper,
            const scalar octave,
            labelList& freqBandIDs
        );

        //- Return the 1/octave octave frequency bounds
        static tmp<scalarField> octaveFrequencies
        (
            const scalarField& f,
            const scalar fLower,
            const scalar fUpper,
            const scalar octave
        );

        //- Return the graph of p(t)
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

        //- Return the PSD [dB]
        graph PSD(const graph& gPSD) const;

        //- Return the narrow-band PFL (pressure-fluctuation level)
        //  spectrum [dB]
        graph Lf(const graph& gPf) const;

        //- Return the octave-band PFL spectrum starting at octave
        //  frequencies given by the supplied frequency bands [dB]
        graph Ldelta
        (
            const graph& gLf,
            const labelList& freqBandIDs
        ) const;

        //- Return the octave-band pressure spectrum at octave
        //  frequencies given by the supplied frequency bands [dB]
        graph Pdelta
        (
            const graph& gLf, 
            const labelList& freqBandIDs
        ) const;

        //- Return the total PFL as the sum of Lf over all frequencies
        scalar Lsum(const graph& gLf) const;

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
