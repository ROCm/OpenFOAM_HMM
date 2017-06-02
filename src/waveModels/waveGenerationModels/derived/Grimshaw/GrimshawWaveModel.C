/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "GrimshawWaveModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Grimshaw, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        Grimshaw,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Grimshaw::alfa
(
    const scalar H,
    const scalar h
) const
{
    scalar eps = H/h;
    scalar alfa = sqrt(0.75*eps)*(1.0 - (5.0/8.0)*eps + (71.0/128.0)*eps*eps);

    return alfa;
}

//- Wave height
Foam::scalar Foam::waveModels::Grimshaw::eta
(
    const scalar H,
    const scalar h,
    const scalar x,
    const scalar y,
    const scalar theta,
    const scalar t,
    const scalar X0
) const
{
    scalar eps = H/h;
    scalar C = sqrt(mag(g_)*h)*sqrt(1.0+eps-(1.0/20.0)*eps*eps-(3.0/70.0)*eps*eps*eps);

    scalar ts = 3.5*h/sqrt(H/h);
    scalar xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);

    scalar alfa = this->alfa(H, h);

    scalar s = (1.0)/(cosh(alfa*(xa/h)));
    scalar q = tanh(alfa*(xa/h));

    return h*(eps*s*s - 0.75*eps*eps*s*s*q*q + eps*eps*eps*((5.0/8.0)*s*s*q*q - (101.0/80.0)*s*s*s*s*q*q));
}

//- Wave velocity
Foam::vector Foam::waveModels::Grimshaw::U
(
    const scalar H,
    const scalar h,
    const scalar x,
    const scalar y,
    const scalar theta,
    const scalar t,
    const scalar X0,
    const scalar z
) const
{
    scalar eps = H/h;
    scalar C = sqrt(mag(g_)*h)*sqrt(1.0+eps-(1.0/20.0)*eps*eps-(3.0/70.0)*eps*eps*eps);
    scalar ts = 3.5*h/sqrt(H/h);
    scalar xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    scalar alfa = this->alfa(H, h);

    scalar s = (1.0)/(cosh(alfa*(xa/h)));

    scalar outa = eps*s*s - eps*eps*(-(1.0/4.0)*s*s + s*s*s*s + ((z/h)*(z/h))*((3.0/2.0)*s*s - (9.0/4.0)*s*s*s*s));
    scalar outb = (19.0/40.0)*s*s + (1.0/5.0)*s*s*s*s - (6.0/5.0)*s*s*s*s*s*s;
    scalar outc = ((z/h)*(z/h)) * ( -(3.0/2.0)*s*s - (15.0/4.0)*s*s*s*s + (15.0/2.0)*s*s*s*s*s*s);
    scalar outd = ((z/h)*(z/h)*(z/h)*(z/h)) * ((-3.0/8.0)*s*s + (45.0/16.0)*s*s*s*s - (45.0/16.0)*s*s*s*s*s*s);

    scalar u = sqrt(mag(g_)*h)*(outa - eps*eps*eps*(outb+outc+outd));

    outa = eps*s*s - eps*eps*((3.0/8.0)*s*s + 2.0*s*s*s*s + ((z/h)*(z/h))*(0.5*s*s - (3.0/2.0)*s*s*s*s));
    outb = (49.0/640.0)*s*s - (17.0/20.0)*s*s*s*s - (18.0/5.0)*s*s*s*s*s*s;
    outc = ((z/h)*(z/h)) * ((-13.0/16.0)*s*s -(25.0/16.0)*s*s*s*s + (15.0/2.0)*s*s*s*s*s*s);
    outd = ((z/h)*(z/h)*(z/h)*(z/h)) * ((-3.0/40.0)*s*s -(9.0/8.0)*s*s*s*s - (27.0/16.0)*s*s*s*s*s*s);

    scalar w = sqrt(mag(g_)*h)*(outa - eps*eps*eps*(outb+outc+outd));

    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}

void Foam::waveModels::Grimshaw::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    forAll(level, paddlei)
    {
        const scalar eta =
            this->eta
            (
                waveHeight_,
                waterDepthRef_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveAngle_,
                t,
                x0_
            );

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Grimshaw::Grimshaw
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    solitaryWaveModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::Grimshaw::~Grimshaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::Grimshaw::read(const dictionary& overrideDict)
{
    if (solitaryWaveModel::read(overrideDict))
    {
        return true;
    }

    return false;
}

void Foam::waveModels::Grimshaw::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level,
    const scalar tg
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

	    if ( (tg<0) || (t >= tg) )
	    {
		    const label paddlei = faceToPaddle_[facei];

		    const vector Uf = U
			(
		        waveHeight_,
		        waterDepthRef_,
		        xPaddle_[paddlei],
		        yPaddle_[paddlei],
		        waveAngle_,
		        t,
		        x0_,
		        z
		    );

                U_[facei] = fraction*Uf*tCoeff + fraction*UCurrent_;

	    }
	    else if ( tg>=t )
	    {
	        U_[facei] = fraction*UCurrent_;	
	    }
        }
    }
}

void Foam::waveModels::Grimshaw::setVelocityAbsorption
(
    const scalarField& calculatedLevel,
    const scalarField& activeLevel
) 
{

    forAll(U_, facei)
    {
	const label paddlei = faceToPaddle_[facei];

	scalar activeLevelMBL=activeLevel[paddlei];

	scalar zMin = zMin_[facei];

//------ not needed anymore in new release
	if (fabs(zMinGb_)>1.0e-3)
    	{
	      zMin = zMin - zMinGb_;
    	}
//------
	
	if (zMin < activeLevelMBL)

	{
	    scalar UCorr =
	    	 (calculatedLevel[paddlei] - activeLevel[paddlei])
	    	 *sqrt(mag(g_)/activeLevel[paddlei]);
					
	    U_[facei].x() += UCorr;
	}
	else
	{
	    U_[facei].x() = 0;
	}
    }
}

void Foam::waveModels::Grimshaw::info(Ostream& os) const
{
    solitaryWaveModel::info(os);
}


// ************************************************************************* //
