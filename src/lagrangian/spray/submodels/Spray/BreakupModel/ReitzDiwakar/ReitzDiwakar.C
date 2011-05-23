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

#include "ReitzDiwakar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ReitzDiwakar<CloudType>::ReitzDiwakar
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict,owner, typeName),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    Cbag_(readScalar(coeffsDict_.lookup("Cbag"))),
    Cb_(readScalar(coeffsDict_.lookup("Cb"))),
    Cstrip_(readScalar(coeffsDict_.lookup("Cstrip"))),
    Cs_(readScalar(coeffsDict_.lookup("Cs")))
{}

    /*
        These are the default values for this model...
        static const scalar Cbag   = 6.0;
        static const scalar Cb     = 0.785;
        static const scalar Cstrip = 0.5;
        static const scalar Cs     = 10.0;
    */


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ReitzDiwakar<CloudType>::~ReitzDiwakar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ReitzDiwakar<CloudType>::update
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

    scalar d1 = d;
    scalar nuc = muc/rhoc;
    scalar We = 0.5*rhoc*pow(Urmag, 2)*d/sigma;
    scalar Re = Urmag*d/nuc;

    scalar sqRey = sqrt(Re);

    if (We > Cbag_)
    {
        if (We > Cstrip_*sqRey)
        {
            scalar dStrip = pow(2.0*Cstrip_*sigma, 2.0)/
            (
                rhoc*pow(Urmag, 3.0)*muc
            );

            scalar tauStrip = Cs_*d*sqrt(rho/rhoc)/Urmag;
            scalar fraction = dt/tauStrip;

            // new droplet diameter, implicit calculation
            d = (fraction*dStrip + d)/(1.0 + fraction);
        }
        else
        {
            scalar dBag =
                2.0 * Cbag_ * sigma
              / (
                    rhoc * pow(Urmag, 2.0)
                );

            scalar tauBag =
                Cb_ * d
                * sqrt
                  (
                      rho * d / sigma
                  );

            scalar fraction = dt/tauBag;

            // new droplet diameter, implicit calculation
            d = (fraction*dBag + d)/(1.0 + fraction);
        }

        // preserve the total mass/volume, by increasing the number of particles in parcels due to breakup
        nParticle *= pow(d1/d, 3.0);

    }

    return false;

}


// ************************************************************************* //
