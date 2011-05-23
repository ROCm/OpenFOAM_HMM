/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "PilchErdman.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::PilchErdman<CloudType>::PilchErdman
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict,owner, typeName),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    B1_(readScalar(coeffsDict_.lookup("B1"))),
    B2_(readScalar(coeffsDict_.lookup("B2")))
{}

    /*
        These are the default values for this model...
        static const scalar B1     = 0.375;
        static const scalar B2     = 0.236;
    */


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::PilchErdman<CloudType>::~PilchErdman()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PilchErdman<CloudType>::update
(
    const scalar& dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
    scalar& y,
    scalar& yDot,
    const scalar& d0,
    const scalar& rho,
    const scalar& mu,
    const scalar& sigma,
    const vector& U,
    const scalar& rhoc,
    const scalar& muc,
    const vector& Urel,
    const scalar& Urmag,
    const scalar& tMom,
    const scalar& averageParcelMass,
    scalar& dChild,
    scalar& massChild,
    Random& rndGen
) const
{

    scalar semiMass = nParticle*pow(d, 3);
    scalar We = 0.5*rhoc*pow(Urmag, 2)*d/sigma;
    //scalar nuc = muc/rhoc;
    //scalar Re = Urmag*d/nuc;
    scalar Oh = mu/pow(rho*d*sigma, 0.5);

    scalar Wec = 6.0*(1.0 + 1.077*pow(Oh, 1.6));

    if (We > Wec)
    {

        // We > 1335, wave crest stripping
        scalar taubBar = 5.5;

        if (We > 175.0)
        {
            // sheet stripping
            taubBar = 0.766*pow(2.0*We - 12.0, 0.25);
        }
        else if (We > 22.0)
        {
            // Bag-and-stamen breakup
            taubBar = 14.1*pow(2.0*We - 12.0, -0.25);
        }
        else if (We > 9.0)
        {
            // Bag breakup
            taubBar = 2.45*pow(2.0*We - 12.0, 0.25);
        }
        else if (We > 6.0)
        {
            // Vibrational breakup
            taubBar = 6.0*pow(2.0*We - 12.0, -0.25);
        }

        scalar rho12 = pow(rhoc/rho, 0.5);

        scalar Vd = Urmag*rho12*(B1_*taubBar * B2_*taubBar*taubBar);
        scalar Vd1 = pow(1.0 - Vd/Urmag, 2.0);
        Vd1 = max(Vd1, SMALL);
        scalar Ds = 2.0*Wec*sigma*Vd1/(Vd1*rhoc*pow(Urmag, 2.0));
        scalar A = Urmag*rho12/d;

        scalar taub = taubBar/A;

        scalar frac = dt/taub;

        // update the droplet diameter according to the rate eq. (implicitly)
        d = (d + frac*Ds)/(1.0 + frac);

        // correct the number of particles to conserve mass
        nParticle = semiMass/pow(d, 3);
    }

    return false;

}


// ************************************************************************* //
