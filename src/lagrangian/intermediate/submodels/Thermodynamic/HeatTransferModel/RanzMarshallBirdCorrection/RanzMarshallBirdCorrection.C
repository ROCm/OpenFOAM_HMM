/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "RanzMarshallBirdCorrection.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::RanzMarshallBirdCorrection<CloudType>::RanzMarshallBirdCorrection
(
    const dictionary& dict,
    CloudType& cloud
)
:
    RanzMarshall<CloudType>(dict, cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::RanzMarshallBirdCorrection<CloudType>::~RanzMarshallBirdCorrection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::RanzMarshallBirdCorrection<CloudType>::htc
(
    const scalar dp,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW
) const
{
    scalar htc = RanzMarshall<CloudType>::htc(dp, Re, Pr, kappa, NCpW);

    // Bird correction
    if (mag(htc) > ROOTVSMALL && mag(NCpW) > ROOTVSMALL)
    {
        const scalar phit = min(NCpW/htc, 50);
        scalar fBird = 1.0;
        if (phit > 0.001)
        {
            fBird = phit/(exp(phit) - 1.0);
        }
        htc *= fBird;
    }

    return htc;
}


// ************************************************************************* //
