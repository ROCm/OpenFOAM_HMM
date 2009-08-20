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

#include "DeterministicPairForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::DeterministicPairForce<CloudType>::buildCellOccupancy()
{
    Info<< "Build cell occupancy" << endl;

    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(typename CloudType, this->owner(), iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }

}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::DeterministicPairForce<CloudType>::evaluatePair
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    Info<< "PARTICLE FORCES HARD CODED AS TEST" << endl;
    pA.f() += -20.0*pB.mass()*pB.U();
    pB.f() += -20.0*pA.mass()*pA.U();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DeterministicPairForce<CloudType>::DeterministicPairForce
(
    const dictionary& dict,
    CloudType& owner
)
:
    CollisionModel<CloudType>(dict, owner, typeName),
    cellOccupancy_(owner.mesh().nCells()),
    il_(owner.mesh(), 1e-6, false)
{
    Info<< "SEARCH DISTANCE SQR HARD CODED" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DeterministicPairForce<CloudType>::~DeterministicPairForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::DeterministicPairForce<CloudType>::active() const
{
    return true;
}


template<class CloudType>
void Foam::DeterministicPairForce<CloudType>::collide()
{
    Info<< "Calculating collisions" << endl;

    forAllIter(typename CloudType, this->owner(), iter)
    {
        typename CloudType::parcelType& p = iter();
        p.f() = vector::zero;
    }

    buildCellOccupancy();

    const directInteractionList& dil(il_.dil());

    typename CloudType::parcelType* pA_ptr = NULL;
    typename CloudType::parcelType* pB_ptr = NULL;

    forAll(dil, d)
    {
        // Loop over all Parcels in cell A (a)
        forAll(cellOccupancy_[d], a)
        {
            pA_ptr = cellOccupancy_[d][a];

            forAll(dil[d], interactingCells)
            {
                List<typename CloudType::parcelType*> cellBParcels =
                    cellOccupancy_[dil[d][interactingCells]];

                // Loop over all Parcels in cell B (b)
                forAll(cellBParcels, b)
                {
                    pB_ptr = cellBParcels[b];

                    evaluatePair(*pA_ptr, *pB_ptr);
                }
            }

            // Loop over the other Parcels in cell A (aO)
            forAll(cellOccupancy_[d], aO)
            {
                pB_ptr = cellOccupancy_[d][aO];

                // Do not double-evaluate, compare pointers, arbitrary
                // order
                if (pB_ptr > pA_ptr)
                {
                    evaluatePair(*pA_ptr, *pB_ptr);
                }
            }
        }
    }

    Info<< "ADD COLLISIONS WITH WALLS HERE, DOES NOT NEED TO BE A TRACKING"
        << "OPERATION.  CALCULATE DISTANCE TO SURFACES OF WALL TYPE AND APPLY"
        << "WALL FORCE MODEL" << endl;
}


// ************************************************************************* //
