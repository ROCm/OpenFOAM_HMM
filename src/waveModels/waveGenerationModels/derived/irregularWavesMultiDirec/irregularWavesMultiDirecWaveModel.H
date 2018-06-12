/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2015 IH-Cantabria
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
    Foam::waveModels::irregularWavesMultiDirec

Description
    Multi-directional irregular wave model

\*---------------------------------------------------------------------------*/

#ifndef waveModels_irregularWavesMultiDirec_H
#define waveModels_irregularWavesMultiDirec_H

#include "irregularWaveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                            Class irregularWavesMultiDirec Declaration
\*---------------------------------------------------------------------------*/

class irregularWavesMultiDirec
:
    public irregularWaveModel
{
private:

    // Proteced Data

        //- Wave periods for irregularWavesMultiDirec case (seconds)
        List< List<scalar> > irregWavePeriods_;

        //- Wave heights for irregularWavesMultiDirec case (meters)
	List< List<scalar> > irregWaveHeights_;

        //- Wave lengths for irregularWavesMultiDirec case (meters)
	List< List<scalar> > irregWaveLengths_;

        //- Wave phases for irregularWavesMultiDirec case (radians)
	List< List<scalar> > irregWavePhases_;

        //- Direction of propagation for irregularWavesMultiDirec case (degrees, from X axis)
	List< List<scalar> > irregWaveDirs_;


    // Private Member Functions

        //- First Order Wave height (same as StokesI)
        virtual scalar eta
        (
            const scalar H,
            const scalar Kx,
            const scalar x,
            const scalar Ky,
            const scalar y,
            const scalar omega,
            const scalar t,
            const scalar phase
        ) const;


protected:

    // Protected Member Functions

        //- Return the wavelength
        virtual scalar waveLength
	(
	     const scalar h, 
	     const scalar T
	) const;

        //- Wave velocity
        virtual vector U
        (
            const scalar d,
            const scalar x,
            const scalar y,
            const scalar t,
            const scalar z
        ) const;

        //-
        virtual vector uMultiDirec
        (
            const scalar irregH,
            const scalar irregWaveOmega,
            const scalar phaseTot,
            const scalar irregWaveKs,
            const scalar z,
            const scalar h,
	    const scalar irregDir
        ) const;

        //- Set the water level
        virtual void setLevel
        (
            const scalar t,
            const scalar tCoeff,
            scalarField& level
        ) const;

        //- Calculate the wave model velocity
        virtual void setVelocity
        (
            const scalar t,
            const scalar tCoeff,
            const scalarField& level
        );

public:

    //- Runtime type information
    TypeName("irregularWavesMultiDirectional");

    //- Constructor
    irregularWavesMultiDirec
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );

    //- Destructor
    virtual ~irregularWavesMultiDirec();


    // Public Member Functions

        //- Read from dictionary
        virtual bool readDict(const dictionary& overrideDict);

        //- Info
        virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
