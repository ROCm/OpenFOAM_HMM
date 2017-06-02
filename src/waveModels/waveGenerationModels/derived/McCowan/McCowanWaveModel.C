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

#include "McCowanWaveModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(McCowan, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        McCowan,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::McCowan::eta
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
    vector vec = this->mn(H,h);
    scalar mm = vec[0];
    scalar nn = vec[1];

    scalar C = sqrt(((mag(g_)*h)/mm)*tan(mm));
    scalar ts = 3.5*h/sqrt(H/h);
    scalar Xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);

    scalar xin = 0.5*H;
    scalar etas = newtonRapsonF2(xin,H,h,Xa,mm,nn);
    return etas;
}

Foam::vector Foam::waveModels::McCowan::mn
(
    const scalar H,
    const scalar h
) const
{
    //m
    scalar xin = 1.0;
    scalar m = newtonRapsonF1(xin,H,h);

    //n
    scalar c1=sin(m+(1.0+((2.0*H)/(3.0*h))));
    scalar n = (2.0/3.0)*pow(c1,2);

    scalar aux=n;

    return vector(m, n, aux);
}

Foam::scalar Foam::waveModels::McCowan::newtonRapsonF1
(
    const scalar x0,
    const scalar H,
    const scalar h
) const
{
    scalar N=10000; 
    scalar eps=1.e-5; 
    scalar maxval = 10000.0;     

    scalar xn=0;
    scalar x=0;
    scalar c1=0;
    scalar c2=0;
    scalar c3=0;
    scalar fx=0;
    scalar fprime=0;
    scalar fxn=0;
    scalar fxx=0;

    //define value for divergence
    scalar xx=x0; 
    while (N>0) 
    {
       //f
       c1=sin(xx+(1.0+((2.0*H)/(3.0*h))));
       fx = (2.0/3.0)*pow(c1,2) - (xx*H)/(h*tan(0.5*xx*(1.0+(H/h))));

       //fprime
       c2=1/tan(0.5*xx*(h/H + 1.0));
       c3=1/sin(0.5*xx*(H/h + 1.0));
       fprime=(4.0/3.0)*sin((2.0*H)/(3.0*h) + xx + 1.0)*cos((2.0*H)/(3.0*h) + xx + 1.0)-(h*c2)/H - (0.5*h*xx*(H/h + 1.0)*pow(c3,2))/h;
       xn = xx-fx/fprime;

       c1=sin(xn+(1.0+((2.0*H)/(3.0*h))));
       fxn = (2.0/3.0)*pow(c1,2) - (xn*H)/(h*tan(0.5*xn*(1.0+(H/h)))); 
       if (fabs(fxn)<eps)
       {
           x=xn; 
           return x; 
       }
   
       c1=sin(xx+(1.0+((2.0*H)/(3.0*h))));
       fxx = (2.0/3.0)*pow(c1,2) - (xx*H)/(h*tan(0.5*xx*(1.0+(H/h))));
       if (fabs(fxx)>maxval)
       {
          FatalIOErrorInFunction(*this)
            << "fxx > maxval !!!"                    
            << exit(FatalIOError); 
       }

       N = N - 1; 
       xx = xn; 
    }

    return x;
}

Foam::scalar Foam::waveModels::McCowan::newtonRapsonF2
(
    const scalar x0,
    const scalar H,
    const scalar h,
    const scalar xa,
    const scalar m,
    const scalar n
) const
{
    scalar N=10000; 
    scalar eps=1.e-5; 
    scalar maxval = 10000.0;     

    scalar xn=0;
    scalar x=0;
    scalar c2=0;
    scalar c3=0;
    scalar fx=0;
    scalar fprime=0;
    scalar fxn=0;
    scalar fxx=0;

    //define value for divergence
    scalar xx=x0; 
    while (N>0) 
    {
       //f
       fx = xx-(h*(n/m)*((sin(m*(1.0+(xx/h))))/(cos(m*(1.0+(xx/h)))+cosh(m*(xa/h)))));

       //fprime
       c2=sin((m*(h + x))/h);
       c3=cosh((xa*m)/h) + cos((m*(h + x))/h);
       fprime =  1 - (n*cos((m*(h + x))/h))/(cosh((xa*m)/h) + cos((m*(h + x))/h)) - (n*pow(c2,2))/(pow(c3,2));

       xn = xx-fx/fprime;

       fxn = xn-(h*(n/m)*((sin(m*(1.0+(xn/h))))/(cos(m*(1.0+(xn/h)))+cosh(m*(xa/h)))));
       if (fabs(fxn)<eps)
       {
           x=xn;
           return x; 
       }
   
       fxx = xx-(h*(n/m)*((sin(m*(1.0+(xx/h))))/(cos(m*(1.0+(xx/h)))+cosh(m*(xa/h)))));
       if (fabs(fxx)>maxval)
       {
          FatalIOErrorInFunction(*this)
            << "fxx > maxval !!!"                    
            << exit(FatalIOError); 
       }

       N = N - 1; 
       xx = xn; 
    }

    return x;
}

Foam::vector Foam::waveModels::McCowan::U
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
    vector vec = this->mn(H,h);
    scalar mm = vec[0];
    scalar nn = vec[1];

    scalar C = sqrt(((mag(g_)*h)/mm)*tan(mm));
    scalar ts = 3.5*h/sqrt(H/h);
    scalar Xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    
    scalar outa = C*nn*(1.0+cos(mm*(z/h))*cosh(mm*(Xa/h)));
    scalar outb = (cos(mm*(z/h))+cosh(mm*(Xa/h))) * (cos(mm*(z/h))+cosh(mm*(Xa/h)));

    scalar u = outa/outb;

    outa = C*nn*(sin(mm*(z/h))*sinh(mm*(Xa/h)));
    outb = (cos(mm*(z/h))+cosh(mm*(Xa/h))) * (cos(mm*(z/h))+cosh(mm*(Xa/h)));

    scalar w = outa/outb;

    scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


void Foam::waveModels::McCowan::setLevel
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

Foam::waveModels::McCowan::McCowan
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

Foam::waveModels::McCowan::~McCowan()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::McCowan::read(const dictionary& overrideDict)
{
    if (solitaryWaveModel::read(overrideDict))
    {
        return true;
    }

    return false;
}

void Foam::waveModels::McCowan::setVelocity
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

void Foam::waveModels::McCowan::setVelocityAbsorption
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

void Foam::waveModels::McCowan::setCurrent
(
    const scalarField& levelMBO
)
{
    //No needed for generation
}

void Foam::waveModels::McCowan::info(Ostream& os) const
{
    solitaryWaveModel::info(os);
}


// ************************************************************************* //
