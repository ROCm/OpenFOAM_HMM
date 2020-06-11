/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "WeberNumberReacting.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WeberNumberReacting<CloudType>::WeberNumberReacting
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName)
{}


template<class CloudType>
Foam::WeberNumberReacting<CloudType>::WeberNumberReacting
(
    const WeberNumberReacting<CloudType>& we
)
:
    CloudFunctionObject<CloudType>(we)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::WeberNumberReacting<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    const auto& c = this->owner();

    if (!c.template foundObject<IOField<scalar>>("We"))
    {
        IOField<scalar>* WePtr =
            new IOField<scalar>
            (
                IOobject
                (
                    "We",
                    c.time().timeName(),
                    c,
                    IOobject::NO_READ
                )
            );

        WePtr->store();
    }

    auto& We = c.template lookupObjectRef<IOField<scalar>>("We");
    We.setSize(c.size());

    const auto& thermo = c.db().template lookupObject<SLGThermo>("SLGThermo");
    const auto& liquids = thermo.liquids();

    const auto& UInterp = td.UInterp();
    const auto& pInterp = td.pInterp();
    const auto& rhoInterp = td.rhoInterp();

    label parceli = 0;
    forAllConstIters(c, parcelIter)
    {
        const parcelType& p = parcelIter();

        const auto& coords = p.coordinates();
        const auto& tetIs = p.currentTetIndices();

        const vector Uc(UInterp.interpolate(coords, tetIs));

        const scalar pc =
            max
            (
                pInterp.interpolate(coords, tetIs),
                c.constProps().pMin()
            );

        const scalar rhoc(rhoInterp.interpolate(coords, tetIs));
        const scalarField X(liquids.X(p.YLiquid()));
        const scalar sigma = liquids.sigma(pc, p.T(), X);

        We[parceli++] = rhoc*magSqr(p.U() - Uc)*p.d()/sigma;
    }


    if (c.size() && c.time().writeTime())
    {
        We.write();
    }
}


// ************************************************************************* //
