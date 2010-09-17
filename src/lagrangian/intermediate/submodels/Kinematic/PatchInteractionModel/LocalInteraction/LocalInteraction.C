/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "LocalInteraction.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template <class CloudType>
Foam::label Foam::LocalInteraction<CloudType>::applyToPatch
(
    const label globalPatchI
) const
{
    forAll(patchIds_, patchI)
    {
        if (patchIds_[patchI] == globalPatchI)
        {
            return patchI;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    patchProperties_(this->coeffDict().lookup("patches")),
    patchIds_(patchProperties_.size())
{
    const polyMesh& mesh = cloud.mesh();
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // check that user patches are valid region patches
    forAll(patchProperties_, patchI)
    {
        const word& patchName = patchProperties_[patchI].patchName();
        patchIds_[patchI] = bMesh.findPatchID(patchName);
        if (patchIds_[patchI] < 0)
        {
            FatalErrorIn("LocalInteraction(const dictionary&, CloudType&)")
                << "Patch " << patchName << " not found. Available patches "
                << "are: " << bMesh.names() << nl << exit(FatalError);
        }
    }

    // check that all walls are specified
    DynamicList<word> badWalls;
    forAll(bMesh, patchI)
    {
        if
        (
            isA<wallPolyPatch>(bMesh[patchI])
         && applyToPatch(bMesh[patchI].index()) < 0
        )
        {
            badWalls.append(bMesh[patchI].name());
        }
    }

    if (badWalls.size() > 0)
    {
        FatalErrorIn("LocalInteraction(const dictionary&, CloudType&)")
            << "All wall patches must be specified when employing local patch "
            << "interaction. Please specify data for patches:" << nl
            << badWalls << nl << exit(FatalError);
    }

    // check that interactions are valid/specified
    forAll(patchProperties_, patchI)
    {
        const word& interactionTypeName =
            patchProperties_[patchI].interactionTypeName();
        const typename PatchInteractionModel<CloudType>::interactionType& it =
            this->wordToInteractionType(interactionTypeName);

        if (it == PatchInteractionModel<CloudType>::itOther)
        {
            const word& patchName = patchProperties_[patchI].patchName();
            FatalErrorIn("LocalInteraction(const dictionary&, CloudType&)")
                << "Unknown patch interaction type "
                << interactionTypeName << " for patch " << patchName
                << ". Valid selections are:"
                << this->PatchInteractionModel<CloudType>::interactionTypeNames_
                << nl << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LocalInteraction<CloudType>::~LocalInteraction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::LocalInteraction<CloudType>::active() const
{
    return true;
}


template <class CloudType>
bool Foam::LocalInteraction<CloudType>::correct
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

    label patchI = applyToPatch(pp.index());

    if (patchI >= 0)
    {
        typename PatchInteractionModel<CloudType>::interactionType it =
            this->wordToInteractionType
            (
                patchProperties_[patchI].interactionTypeName()
            );

        switch (it)
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
                    U -= (1.0 + patchProperties_[patchI].e())*Un*nw;
                }

                U -= patchProperties_[patchI].mu()*Ut;

                // Return velocity to global space
                U += Up;

                break;
            }
            default:
            {
                FatalErrorIn
                (
                    "bool LocalInteraction<CloudType>::correct"
                    "("
                        "const polyPatch&, "
                        "const label, "
                        "bool&, "
                        "vector&"
                    ") const"
                )   << "Unknown interaction type "
                    << patchProperties_[patchI].interactionTypeName()
                    << "(" << it << ") for patch "
                    << patchProperties_[patchI].patchName()
                    << ". Valid selections are:" << this->interactionTypeNames_
                    << endl << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
