/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    patchData_(cloud.mesh(), this->coeffDict()),
    nEscape_(patchData_.size()),
    massEscape_(patchData_.size()),
    nStick_(patchData_.size()),
    massStick_(patchData_.size()),
    writeFields_(this->coeffDict().lookupOrDefault("writeFields", false)),
    outputByInjectorId_(this->coeffDict().lookupOrDefault("outputByInjectorId", false)),
    injIdToIndex_(cloud.injectors().size()),
    massEscapePtr_(nullptr),
    massStickPtr_(nullptr)
{
    if (writeFields_)
    {
        word massEscapeName(this->owner().name() + ":massEscape");
        word massStickName(this->owner().name() + ":massStick");
        Info<< "    Interaction fields will be written to " << massEscapeName
            << " and " << massStickName << endl;

        (void)massEscape();
        (void)massStick();
    }
    else
    {
        Info<< "    Interaction fields will not be written" << endl;
    }

    // check that interactions are valid/specified
    forAll(patchData_, patchi)
    {
        const word& interactionTypeName =
            patchData_[patchi].interactionTypeName();
        const typename PatchInteractionModel<CloudType>::interactionType& it =
            this->wordToInteractionType(interactionTypeName);

        if (it == PatchInteractionModel<CloudType>::itOther)
        {
            const word& patchName = patchData_[patchi].patchName();
            FatalErrorInFunction
                << "Unknown patch interaction type "
                << interactionTypeName << " for patch " << patchName
                << ". Valid selections are:"
                << this->PatchInteractionModel<CloudType>::interactionTypeNames_
                << nl << exit(FatalError);
        }

        label nInjectors(1);
        if (outputByInjectorId_)
        {
            nInjectors = cloud.injectors().size();
            for (label i=0; i<nInjectors; i++)
            {
                injIdToIndex_.insert(cloud.injectors()[i].injectorID(), i);
            }
        }

        nEscape_[patchi].setSize(nInjectors, 0);
        massEscape_[patchi].setSize(nInjectors, 0.0);
        nStick_[patchi].setSize(nInjectors, 0);
        massStick_[patchi].setSize(nInjectors, 0.0);
    }
}


template<class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const LocalInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    patchData_(pim.patchData_),
    nEscape_(pim.nEscape_),
    massEscape_(pim.massEscape_),
    nStick_(pim.nStick_),
    massStick_(pim.massStick_),
    writeFields_(pim.writeFields_),
    outputByInjectorId_(pim.outputByInjectorId_),
    injIdToIndex_(pim.injIdToIndex_),
    massEscapePtr_(nullptr),
    massStickPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LocalInteraction<CloudType>::~LocalInteraction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::volScalarField& Foam::LocalInteraction<CloudType>::massEscape()
{
    if (!massEscapePtr_.valid())
    {
        const fvMesh& mesh = this->owner().mesh();

        massEscapePtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + ":massEscape",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
    }

    return massEscapePtr_();
}


template<class CloudType>
Foam::volScalarField& Foam::LocalInteraction<CloudType>::massStick()
{
    if (!massStickPtr_.valid())
    {
        const fvMesh& mesh = this->owner().mesh();

        massStickPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + ":massStick",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
    }

    return massStickPtr_();
}


template<class CloudType>
bool Foam::LocalInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    label patchi = patchData_.applyToPatch(pp.index());

    if (patchi >= 0)
    {
        vector& U = p.U();

        typename PatchInteractionModel<CloudType>::interactionType it =
            this->wordToInteractionType
            (
                patchData_[patchi].interactionTypeName()
            );

        switch (it)
        {
            case PatchInteractionModel<CloudType>::itNone:
            {
                return false;
            }
            case PatchInteractionModel<CloudType>::itEscape:
            {
                scalar dm = p.mass()*p.nParticle();

                keepParticle = false;
                p.active(false);
                U = Zero;
                if (outputByInjectorId_)
                {
                    nEscape_[patchi][injIdToIndex_[p.typeId()]]++;
                    massEscape_[patchi][injIdToIndex_[p.typeId()]] += dm;
                }
                else
                {
                    nEscape_[patchi][0]++;
                    massEscape_[patchi][0] += dm;
                }
                if (writeFields_)
                {
                    label pI = pp.index();
                    label fI = pp.whichFace(p.face());
                    massEscape().boundaryFieldRef()[pI][fI] += dm;
                }
                break;
            }
            case PatchInteractionModel<CloudType>::itStick:
            {
                scalar dm = p.mass()*p.nParticle();

                keepParticle = true;
                p.active(false);
                U = Zero;
                if (outputByInjectorId_)
                {
                    nStick_[patchi][injIdToIndex_[p.typeId()]]++;
                    massStick_[patchi][injIdToIndex_[p.typeId()]] += dm;
                }
                else
                {
                    nStick_[patchi][0]++;
                    massStick_[patchi][0] += dm;
                }
                if (writeFields_)
                {
                    label pI = pp.index();
                    label fI = pp.whichFace(p.face());
                    massStick().boundaryFieldRef()[pI][fI] += dm;
                }
                break;
            }
            case PatchInteractionModel<CloudType>::itRebound:
            {
                keepParticle = true;
                p.active(true);

                vector nw;
                vector Up;

                this->owner().patchData(p, pp, nw, Up);

                // Calculate motion relative to patch velocity
                U -= Up;

                scalar Un = U & nw;
                vector Ut = U - Un*nw;

                if (Un > 0)
                {
                    U -= (1.0 + patchData_[patchi].e())*Un*nw;
                }

                U -= patchData_[patchi].mu()*Ut;

                // Return velocity to global space
                U += Up;

                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown interaction type "
                    << patchData_[patchi].interactionTypeName()
                    << "(" << it << ") for patch "
                    << patchData_[patchi].patchName()
                    << ". Valid selections are:" << this->interactionTypeNames_
                    << endl << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


template<class CloudType>
void Foam::LocalInteraction<CloudType>::info(Ostream& os)
{
    PatchInteractionModel<CloudType>::info(os);

    // retrieve any stored data
    labelListList npe0(nEscape_);
    this->getModelProperty("nEscape", npe0);

    scalarListList mpe0(massEscape_);
    this->getModelProperty("massEscape", mpe0);

    labelListList nps0(nStick_);
    this->getModelProperty("nStick", nps0);

    scalarListList mps0(massStick_);
    this->getModelProperty("massStick", mps0);

    // accumulate current data
    labelListList npe(nEscape_);
    forAll(npe, i)
    {
        Pstream::listCombineGather(npe[i], plusEqOp<label>());
        npe[i] = npe[i] + npe0[i];
    }

    scalarListList mpe(massEscape_);
    forAll(mpe, i)
    {
        Pstream::listCombineGather(mpe[i], plusEqOp<scalar>());
        mpe[i] = mpe[i] + mpe0[i];
    }

    labelListList nps(nStick_);
    forAll(nps, i)
    {
        Pstream::listCombineGather(nps[i], plusEqOp<label>());
        nps[i] = nps[i] + nps0[i];
    }

    scalarListList mps(massStick_);
    forAll(nps, i)
    {
        Pstream::listCombineGather(mps[i], plusEqOp<scalar>());
        mps[i] = mps[i] + mps0[i];
    }


    if (outputByInjectorId_)
    {
        forAll(patchData_, i)
        {
            forAll (mpe[i], injId)
            {
                os  << "    Parcel fate: patch " <<  patchData_[i].patchName()
                    << " (number, mass)" << nl
                    << "      - escape  (injector " << injIdToIndex_.toc()[injId]
                    << " )  = " << npe[i][injId]
                    << ", " << mpe[i][injId] << nl
                    << "      - stick   (injector " << injIdToIndex_.toc()[injId]
                    << " )  = " << nps[i][injId]
                    << ", " << mps[i][injId] << nl;
            }
        }
    }
    else
    {
        forAll(patchData_, i)
        {
            os  << "    Parcel fate: patch " <<  patchData_[i].patchName()
                << " (number, mass)" << nl
                << "      - escape                      = " << npe[i][0]
                << ", " << mpe[i][0] << nl
                << "      - stick                       = " << nps[i][0]
                << ", " << mps[i][0] << nl;
        }
    }

    if (this->writeTime())
    {
        this->setModelProperty("nEscape", npe);
        nEscape_ = Zero;

        this->setModelProperty("massEscape", mpe);
        massEscape_ = Zero;

        this->setModelProperty("nStick", nps);
        nStick_ = Zero;

        this->setModelProperty("massStick", mps);
        massStick_ = Zero;
    }
}


// ************************************************************************* //
