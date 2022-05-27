/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "ParticleDose.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleDose<CloudType>::ParticleDose
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    GName_(this->coeffDict().getWord("GName"))
{}


template<class CloudType>
Foam::ParticleDose<CloudType>::ParticleDose
(
    const ParticleDose<CloudType>& re
)
:
    CloudFunctionObject<CloudType>(re),
    GName_(re.GName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleDose<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    auto& c = this->owner();

    if (!c.template foundObject<IOField<scalar>>("D"))
    {
        auto* DPtr =
            new IOField<scalar>
            (
                IOobject
                (
                    "D",
                    c.time().timeName(),
                    c,
                    IOobject::NO_READ
                )
            );

        DPtr->store();
    }

    auto& D = c.template lookupObjectRef<IOField<scalar>>("D");

    D.resize(c.size(), Zero);

    const fvMesh& mesh = this->owner().mesh();

    const auto& G = mesh.lookupObject<volScalarField>(GName_);

    label parceli = 0;
    forAllConstIters(c, parcelIter)
    {
        const parcelType& p = parcelIter();

        D[parceli] += G[p.cell()]*mesh.time().deltaTValue();
        parceli++;
    }

    if (c.size() && c.time().writeTime())
    {
        D.write();
    }
}


// ************************************************************************* //
