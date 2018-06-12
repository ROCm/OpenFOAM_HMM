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

\*---------------------------------------------------------------------------*/

#include "irregularWavesMultiDirecWaveModel.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(irregularWavesMultiDirec, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        irregularWavesMultiDirec,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// First order wave height
Foam::scalar Foam::waveModels::irregularWavesMultiDirec::eta
(
    const scalar H,
    const scalar Kx,
    const scalar x,
    const scalar Ky,
    const scalar y,
    const scalar omega,
    const scalar t,
    const scalar phase
) const
{
    scalar phaseTot = Kx*x + Ky*y - omega*t + phase;
    return H*0.5*cos(phaseTot);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

//Calculate waveLength
Foam::scalar Foam::waveModels::irregularWavesMultiDirec::waveLength
(
    const scalar h,
    const scalar T
) const
{
    scalar L0 = mag(g_)*T*T/(2.0*mathematical::pi);
    scalar L = L0;

    for (int iii=1; iii<=100; iii++)
    {
        L = L0*tanh(2.0*mathematical::pi*h/L);
    }

    return L;
}


Foam::vector Foam::waveModels::irregularWavesMultiDirec::U
(
    const scalar h,
    const scalar x,
    const scalar y,
    const scalar t,
    const scalar z
) const
{
    scalar u = 0.0;
    scalar v = 0.0;
    scalar w = 0.0;

    scalar phaseTot = 0.0;
    scalar waveKs_ = 0.0;
    scalar waveOmegas_ = 0.0;

    int ii;
    int jj;
    scalar COLUMNAS = 0;
    scalar FILAS = irregWaveHeights_.size();

    for (ii=0; ii<FILAS; ++ii)
    {		
	COLUMNAS= irregWaveHeights_[ii].size();	

	for (jj=0; jj<COLUMNAS; ++jj)
    	{
	     waveKs_ = mathematical::twoPi/irregWaveLengths_[ii][jj];
	     waveOmegas_ = mathematical::twoPi/irregWavePeriods_[ii][jj];

	     phaseTot = waveKs_*x*cos(irregWaveDirs_[ii][jj]) + waveKs_*y*sin(irregWaveDirs_[ii][jj]) - waveOmegas_*t + irregWavePhases_[ii][jj];

	     const vector Uf = uMultiDirec
	            (
                    irregWaveHeights_[ii][jj],
		    waveOmegas_,
		    phaseTot,
		    waveKs_,
		    z,
		    h,
		    irregWaveDirs_[ii][jj]
              );
	      u = u + Uf[0];
	      v = v + Uf[1];
	      w = w + Uf[2];
	}	    
    }

    return vector(u, v, w);
}


void Foam::waveModels::irregularWavesMultiDirec::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    scalar eta = 0.0;

    scalar COLUMNAS = 0;
    scalar FILAS = 0;

    scalar waveKs_ = 0.0;
    scalar waveOmegas_ = 0.0;

    int ii;
    int jj;

    forAll(level, paddlei)
    {
        eta = 0.0;
	FILAS= irregWaveHeights_.size();

        for (ii=0; ii<FILAS; ++ii)
        {	
	    COLUMNAS= irregWaveHeights_[ii].size();

	    for (jj=0; jj<COLUMNAS; ++jj)
    	    {
		waveKs_ = mathematical::twoPi/irregWaveLengths_[ii][jj];
		waveOmegas_ = mathematical::twoPi/irregWavePeriods_[ii][jj];

                eta +=
                    this->eta
	                (
                        irregWaveHeights_[ii][jj],
			waveKs_*cos(irregWaveDirs_[ii][jj]),
		        xPaddle_[paddlei],
			waveKs_*sin(irregWaveDirs_[ii][jj]),
		        yPaddle_[paddlei],
			waveOmegas_,
                        t,
                        irregWavePhases_[ii][jj]
                );
	     }
	}

	level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}

Foam::vector Foam::waveModels::irregularWavesMultiDirec::uMultiDirec
(
    const scalar irregH,
    const scalar irregWaveOmega,
    const scalar pha,
    const scalar irregWaveKs,
    const scalar zz,
    const scalar hh,
    const scalar irregDir
) const
{

    scalar u =
             irregH*0.5*irregWaveOmega*
             cos(pha)*
             (cosh(irregWaveKs*zz)/
             sinh(irregWaveKs*hh))*
	     cos(irregDir);

    scalar v =     
             irregH*0.5*irregWaveOmega*
             cos(pha)*
             (cosh(irregWaveKs*zz)/
             sinh(irregWaveKs*hh))*
	     sin(irregDir);

    scalar w =
             irregH*0.5*irregWaveOmega*
             sin(pha)*
             (sinh(irregWaveKs*zz)/
             sinh(irregWaveKs*hh));

    return vector(u, v, w);
}

void Foam::waveModels::irregularWavesMultiDirec::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{
    forAll(U_, facei)
    {
        // Fraction of geometry represented by paddle - to be set
        scalar fraction = 1;

        // Height - to be set
        scalar z = 0;

        setPaddlePropeties(level, facei, fraction, z);

        if (fraction > 0)
        {
            const label paddlei = faceToPaddle_[facei];

            const vector Uf = U
	        (
                waterDepthRef_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                t,
                z
            );

            U_[facei] = fraction*Uf*tCoeff;

        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::irregularWavesMultiDirec::irregularWavesMultiDirec
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    irregularWaveModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::irregularWavesMultiDirec::~irregularWavesMultiDirec()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::irregularWavesMultiDirec::readDict(const dictionary& overrideDict)
{
    int ii;
    int jj;

    if (irregularWaveModel::readDict(overrideDict))
    {

	lookup("irregWavePeriods") >> irregWavePeriods_;
	lookup("irregWaveHeights") >> irregWaveHeights_;
	lookup("irregWavePhases") >> irregWavePhases_;
	lookup("irregWaveDirs") >> irregWaveDirs_;

	 irregWaveLengths_ = irregWaveHeights_;
         scalar COLUMNAS = 0;
         scalar FILAS = irregWaveHeights_.size();

         for (ii=0; ii<FILAS; ++ii)
         {		
	      COLUMNAS= irregWaveHeights_[ii].size();	

	      for (jj=0; jj<COLUMNAS; ++jj)
    	      {
		  irregWaveLengths_[ii][jj] = waveLength (waterDepthRef_, irregWavePeriods_[ii][jj]);
		  irregWaveDirs_[ii][jj]  = irregWaveDirs_[ii][jj]  * (mathematical::pi/180);
	      }
         }

        return true;

    }

    return false;
}

void Foam::waveModels::irregularWavesMultiDirec::info(Ostream& os) const
{
    irregularWaveModel::info(os);

    os  << "    WavePeriods (s) coefficients : " << irregWavePeriods_ << nl
        << "    WaveHeights (m) coefficients : " << irregWaveHeights_ << nl
        << "    WavePhases (rad) coefficients : " << irregWavePhases_ << nl
        << "    WaveLengths (m) coefficients : " << irregWaveLengths_ << nl
        << "    WaveDirections (rad) coefficients : " << irregWaveDirs_ << nl;
}

// ************************************************************************* //
