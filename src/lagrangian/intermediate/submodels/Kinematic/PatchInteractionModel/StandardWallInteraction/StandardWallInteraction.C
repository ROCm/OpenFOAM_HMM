/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

#include "StandardWallInteraction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StandardWallInteraction<CloudType>::StandardWallInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    interactionType_
    (
        this->wordToInteractionType(this->coeffDict().lookup("type"))
    ),
    e_(0.0),
    mu_(0.0)
{
    switch (interactionType_)
    {
        case PatchInteractionModel<CloudType>::itOther:
        {
            const word interactionTypeName(this->coeffDict().lookup("type"));

            FatalErrorIn
            (
                "StandardWallInteraction<CloudType>::StandardWallInteraction"
                "("
                    "const dictionary&, "
                    "CloudType& cloud"
                ")"
            )   << "Unknown interaction result type "
                << interactionTypeName
                << ". Valid selections are:" << this->interactionTypeNames_
                << endl << exit(FatalError);

            break;
        }
        case PatchInteractionModel<CloudType>::itRebound:
        {
            e_ = this->coeffDict().lookupOrDefault("e", 1.0);
            mu_ = this->coeffDict().lookupOrDefault("mu", 0.0);
            break;
        }
        default:
        {
            // do nothing
        }
    }
}


template<class CloudType>
Foam::StandardWallInteraction<CloudType>::StandardWallInteraction
(
    const StandardWallInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    interactionType_(pim.interactionType_),
    e_(pim.e_),
    mu_(pim.mu_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StandardWallInteraction<CloudType>::~StandardWallInteraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::StandardWallInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle,
    const scalar trackFraction,
    const tetIndices& tetIs
) const
{
    vector& U = p.U();

    bool& active = p.active();

    if (isA<wallPolyPatch>(pp))
    {
        switch (interactionType_)
        {
            case PatchInteractionModel<CloudType>::itEscape:
            {
                keepParticle = false;
                active = false;
                U = vector::zero;
                break;
            }
            case PatchInteractionModel<CloudType>::itStick:
            {
                keepParticle = true;
                active = false;
                U = vector::zero;
                break;
            }
            case PatchInteractionModel<CloudType>::itRebound:
            {
                keepParticle = true;
                active = true;

                vector nw;
                vector Up;

                this->patchData(p, pp, trackFraction, tetIs, nw, Up);

                // Calculate motion relative to patch velocity
                U -= Up;

                scalar Un = U & nw;
                vector Ut = U - Un*nw;

                if (Un > 0)
                {
                    U -= (1.0 + e_)*Un*nw;
                }

                U -= mu_*Ut;

                // Return velocity to global space
                U += Up;

                break;
            }
            default:
            {
                FatalErrorIn
                (
                    "bool StandardWallInteraction<CloudType>::correct"
                    "("
                        "const polyPatch&, "
                        "const label, "
                        "bool&, "
                        "vector&"
                    ") const"
                )   << "Unknown interaction type "
                    << this->interactionTypeToWord(interactionType_)
                    << "(" << interactionType_ << ")" << endl
                    << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
